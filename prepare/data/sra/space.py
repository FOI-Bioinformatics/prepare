"""Free-space guard for downloads: start an accession only when the disk has room for it.

A download writes to three places: the FASTQ output folder (a project's ``fastq/`` or the
store's ``tmp/``, where the uncompressed files are built before gzip), ``fasterq-dump``'s
temporary folder, and, with ``prefetch``, the ``.sra`` cache. Before each accession starts, the
guard estimates what it needs on each filesystem and checks it against the free space left
once the downloads already running are counted, so parallel workers cannot each see the same
free space and together fill the disk. An accession that would fit once the downloads in
progress release their reservations waits for them; one that would not fit even then is refused
on its own (``insufficient-space: ...``) and the rest of the run goes on.

The estimate for an accession with a known run size (NCBI's ``.sra`` size, from the project
registry) is ``FASTQ_EXPANSION`` plus ``GZIP_EXPANSION`` times that size for the output folder
(compression writes each ``.gz`` before it removes the uncompressed file, so both are there
together), ``FASTQ_EXPANSION`` times for the temporary folder (an assumption: equal to the
output, not measured), plus the ``.sra`` size for the cache when ``prefetch`` is used. Needs of
locations that share a filesystem are added. An accession without a known size needs
``floor_bytes`` (``--min-free-gb``) on each filesystem instead. A filesystem whose free space
cannot be read is assumed to have room, as ``store_adopt`` does: refusing on an unreadable
``disk_usage`` would be worse than trying.

The free space ``disk_usage`` reports already falls as running downloads write, while their full
reservations are still subtracted, so a download in progress is partly counted twice; the guard
errs toward starting the next accession a little later.
"""

import logging
import os
import shutil
import threading
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple, Union

from prepare.data.sra.cleanup import PROJECT_STAGING_FOLDER

if TYPE_CHECKING:  # pragma: no cover - import cycle: prepare.store imports this package
    from prepare.store.layout import StorePaths

logger = logging.getLogger(__name__)

# Uncompressed FASTQ bytes per byte of .sra archive. Measured on seven runs of the crispatus
# test store (Illumina paired-end, 2x150 bp): 6.98, 7.28, 7.33, 7.34, 7.63, 7.69 and 7.73;
# rounded up to the next whole number.
FASTQ_EXPANSION = 8
# gzip-compressed FASTQ bytes per byte of .sra archive, on the same seven runs: 1.49 to 1.66;
# rounded up. Counted on the output folder, where the compressed files are written.
GZIP_EXPANSION = 2

GB = 1024**3

# The locations a download writes to, as keys of the ``locations`` mapping.
OUTPUT, TEMP, CACHE = "output", "temp", "cache"

# The start of a refusal for an accession that would not fit even with no other download running.
# Deliberately not a disk-full message: it fails that accession alone and does not stop the run.
INSUFFICIENT_SPACE_PREFIX = "insufficient-space: not enough free space"

# Seconds between checks of the free space while waiting for running downloads; ``release`` also
# wakes a waiter at once. Read at call time so tests can shorten it.
WAIT_POLL_SECONDS = 5.0


def _existing(path: Path) -> Path:
    """``path`` itself, or its nearest parent that exists (a folder a download will create)."""
    path = Path(path).absolute()
    while not path.exists() and path != path.parent:
        path = path.parent
    return path


def _device(path: Path) -> Optional[int]:
    """The filesystem (``st_dev``) ``path`` lies on, or None when it cannot be read."""
    try:
        return os.stat(_existing(path)).st_dev
    except OSError:
        return None


def mount_point(path: Path) -> Path:
    """The top folder of the filesystem ``path`` lies on (used to name it in messages)."""
    current = _existing(path)
    try:
        device = os.stat(current).st_dev
        while current != current.parent and os.stat(current.parent).st_dev == device:
            current = current.parent
    except OSError:
        pass
    return current


def _size(value: Any) -> Optional[int]:
    """A run size as a positive whole number of bytes, or None (registry values may be strings)."""
    try:
        size = int(value)
    except (TypeError, ValueError):
        return None
    return size if size > 0 else None


def _gb(value: int) -> str:
    """Bytes as gigabytes with one decimal."""
    return f"{value / GB:.1f}"


def download_locations(
    fastq_path: Union[str, Path],
    temp_folder: Optional[Union[str, Path]] = None,
    sra_cache: Optional[Union[str, Path]] = None,
    store: Optional["StorePaths"] = None,
) -> Dict[str, Path]:
    """Where a download writes: the output, temporary and ``.sra`` cache folders.

    Mirrors the defaults of the project and store downloaders: with a store the files are built
    in the store's ``tmp`` folder, which also holds the cache and temporary folder unless they
    are given; without one they are built under ``fastq_path``, the cache defaults to
    ``<fastq>/.sra-cache`` and the temporary folder to the staging folder
    ``<fastq>/.prepare-tmp``, which holds each accession's ``<ACC>_fqtmp`` scratch.
    """
    if store is not None:
        output = Path(store.tmp)
        return {
            OUTPUT: output,
            TEMP: Path(temp_folder) if temp_folder else output,
            CACHE: Path(sra_cache) if sra_cache else output / ".sra-cache",
        }
    fastq = Path(fastq_path)
    return {
        OUTPUT: fastq,
        TEMP: Path(temp_folder) if temp_folder else fastq / PROJECT_STAGING_FOLDER,
        CACHE: Path(sra_cache) if sra_cache else fastq / ".sra-cache",
    }


class SpaceGuard:
    """Reserve disk space for each accession before its download starts; see the module docstring.

    ``reserve`` and ``release`` are called from the worker threads; the reservations of
    downloads in progress are kept per filesystem under one lock. ``floor_bytes`` of zero
    turns the guard off.
    """

    def __init__(
        self,
        locations: Mapping[str, Path],
        floor_bytes: int,
        run_sizes: Mapping[str, Any],
        use_prefetch: bool,
        exempt: Iterable[str] = (),
    ) -> None:
        """Group ``locations`` by filesystem; ``exempt`` accessions (already in the store) need nothing."""
        self.floor_bytes = max(0, int(floor_bytes))
        self.use_prefetch = use_prefetch
        self._sizes = {acc: size for acc, value in (run_sizes or {}).items() if (size := _size(value)) is not None}
        self._exempt = frozenset(exempt)
        self._cond = threading.Condition()
        self._reserved: Dict[int, int] = {}
        self._held: Dict[str, Dict[int, int]] = {}
        # One representative folder and the roles it serves, per filesystem.
        self._paths: Dict[int, Path] = {}
        self._roles: Dict[int, List[str]] = {}
        for role, path in locations.items():
            device = _device(Path(path))
            if device is None:
                logger.debug("Free-space guard: cannot read the filesystem of %s; not checked", path)
                continue
            self._paths.setdefault(device, _existing(Path(path)))
            self._roles.setdefault(device, []).append(role)

    @property
    def enabled(self) -> bool:
        """False when the guard was turned off with a floor of zero."""
        return self.floor_bytes > 0

    def needs(self, accession: str) -> Dict[int, int]:
        """Bytes ``accession`` is expected to need, per filesystem."""
        if accession in self._exempt:
            return {}
        size = self._sizes.get(accession)
        if size is None:
            return {device: self.floor_bytes for device in self._paths}
        per_role = {
            OUTPUT: (FASTQ_EXPANSION + GZIP_EXPANSION) * size,
            TEMP: FASTQ_EXPANSION * size,
            CACHE: size if self.use_prefetch else 0,
        }
        needs = {device: sum(per_role.get(role, 0) for role in roles) for device, roles in self._roles.items()}
        return {device: need for device, need in needs.items() if need > 0}

    def _free(self, device: int) -> Optional[int]:
        """Free bytes on ``device``, or None when they cannot be read."""
        path = self._paths[device]
        try:
            return shutil.disk_usage(path).free
        except OSError as e:
            logger.warning("Could not check free space on %s: %s", path, e)
            return None

    def _shortfall(self, needs: Mapping[int, int]) -> Optional[Tuple[int, int, int]]:
        """The first filesystem without room for ``needs``: (device, free bytes, bytes needed), or None.

        Called with the lock held; in-flight reservations are subtracted from the free space.
        """
        for device, need in needs.items():
            free = self._free(device)
            if free is not None and free - self._reserved.get(device, 0) < need:
                return device, free, need
        return None

    def reserve(self, accession: str, should_stop: Optional[Callable[[], bool]] = None) -> Optional[str]:
        """Reserve the space ``accession`` needs; None once it is reserved, else why it was not.

        When the space is short only because downloads in progress hold reservations, this waits
        for them to release theirs (``release`` wakes it; the free space is read again each time)
        instead of refusing. It refuses at once, with a message starting ``insufficient-space:``,
        when the accession would not fit even if every reservation were returned: that fails
        this accession alone and leaves the rest of the run going. A disk that actually fills up
        during a download is reported by the tool itself and handled as disk-full by the loops.
        ``should_stop`` is checked while waiting; once it returns True the wait ends with
        ``"interrupted"`` and nothing is reserved. Reserving an accession that is still held
        raises ``RuntimeError``: it would be counted twice and only one release could return it.
        """
        if not self.enabled:
            return None
        needs = self.needs(accession)
        waiting = False
        with self._cond:
            if accession in self._held:
                raise RuntimeError(f"{accession} is already reserved")
            while (short := self._shortfall(needs)) is not None:
                device, free, need = short
                reserved = self._reserved.get(device, 0)
                if need > free + reserved:
                    return (
                        f"{INSUFFICIENT_SPACE_PREFIX} on {mount_point(self._paths[device])}: "
                        f"{_gb(free)} GB free, about {_gb(need)} GB needed"
                    )
                if should_stop is not None and should_stop():
                    return "interrupted"
                if not waiting:
                    waiting = True
                    logger.info(
                        "%s: waiting for running downloads to free space on %s (about %s GB needed, %s GB free, "
                        "%s GB reserved by downloads in progress)",
                        accession,
                        mount_point(self._paths[device]),
                        _gb(need),
                        _gb(free),
                        _gb(reserved),
                    )
                self._cond.wait(WAIT_POLL_SECONDS)
            for device, need in needs.items():
                self._reserved[device] = self._reserved.get(device, 0) + need
            self._held[accession] = needs
        return None

    def release(self, accession: str) -> None:
        """Return the space reserved for ``accession`` (its download ended, however it ended)."""
        with self._cond:
            for device, need in self._held.pop(accession, {}).items():
                self._reserved[device] = max(0, self._reserved.get(device, 0) - need)
            self._cond.notify_all()

    def preflight(self, accessions: Iterable[str]) -> List[str]:
        """Warnings for filesystems that cannot hold every download of ``accessions`` at once.

        Only accessions with a known run size are added up; the run still starts, and each
        download is checked again by ``reserve`` before it begins.
        """
        if not self.enabled:
            return []
        totals: Dict[int, int] = {}
        counted = 0
        for accession in accessions:
            if accession in self._exempt or accession not in self._sizes:
                continue
            counted += 1
            for device, need in self.needs(accession).items():
                totals[device] = totals.get(device, 0) + need
        warnings = []
        for device, total in totals.items():
            free = self._free(device)
            if free is not None and free < total:
                warnings.append(
                    f"The {counted} downloads with a known size together need up to about {_gb(total)} GB on "
                    f"{mount_point(self._paths[device])}, which has {_gb(free)} GB free; "
                    "downloads wait for earlier ones to finish, and one that cannot fit is not started"
                )
        return warnings
