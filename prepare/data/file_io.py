"""
File I/O utilities for PREPARE.

This module provides abstract file operations to handle different file types and formats.
"""

import logging
import os
import secrets
import shutil
import socket
import stat
import time
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any, Iterator, List, Optional, Union

import pandas as pd

from prepare.core.exceptions import DataAccessError

logger = logging.getLogger(__name__)

# os.replace can fail with PermissionError on Windows while another process (a virus scanner, an
# indexer, a reader opened without FILE_SHARE_DELETE) holds the target open; the rename is tried
# this many times, spread over REPLACE_RETRY_SECONDS, before the error is raised.
REPLACE_ATTEMPTS = 5
REPLACE_RETRY_SECONDS = 0.5

# pandas infers compression from the file name, which a temporary name hides; this maps the real
# target's suffix to the method pandas would have inferred for it.
_COMPRESSION_BY_SUFFIX = {".gz": "gzip", ".bz2": "bz2", ".zip": "zip", ".xz": "xz", ".zst": "zstd"}


def ensure_directory(path: Union[str, Path]) -> Path:
    """
    Ensure a directory exists, creating it if necessary.

    Args:
        path: Path to the directory

    Returns:
        Path object for the directory

    Raises:
        DataAccessError: If the directory cannot be created
    """
    try:
        dir_path = Path(path)
        dir_path.mkdir(parents=True, exist_ok=True)
        return dir_path
    except Exception as e:
        raise DataAccessError(f"Failed to create directory {path}: {e}")


def is_hidden_name(name: str) -> bool:
    """True for a dotfile name, including the ``._<name>`` AppleDouble files macOS writes next to
    every file on a volume without native extended attributes (ExFAT, SMB, some NAS shares)."""
    return name.startswith(".")


def visible_files(directory: Union[str, Path], *patterns: str, dirs: bool = False) -> List[Path]:
    """Entries directly in ``directory`` matching any of ``patterns``, hidden names removed.

    Returns files (or directories when ``dirs`` is True), sorted by name, without duplicates.
    A directory that does not exist yields an empty list. Every folder listing in prepare
    goes through here so that ``._*`` and ``.DS_Store`` never count as data.
    """
    base = Path(directory)
    if not base.is_dir():
        return []
    found = set()
    for pattern in patterns or ("*",):
        for candidate in base.glob(pattern):
            if is_hidden_name(candidate.name):
                continue
            try:
                keep = candidate.is_dir() if dirs else candidate.is_file()
            except OSError:
                continue
            if keep:
                found.add(candidate)
    return sorted(found)


def list_files(directory: Union[str, Path], pattern: str = "*", include_hidden: bool = False) -> List[Path]:
    """List the files in ``directory`` matching ``pattern``; hidden names are dropped unless asked for."""
    try:
        directory = Path(directory)
        if include_hidden:
            return list(directory.glob(pattern))
        return visible_files(directory, pattern)
    except Exception as e:
        logger.warning(f"Error listing files in {directory} with pattern {pattern}: {e}")
        return []


def copy_file(source: Union[str, Path], destination: Union[str, Path]) -> Path:
    """
    Copy a file from source to destination, publishing it with one rename (``atomic_path``).

    Args:
        source: Source file path
        destination: Destination file path

    Returns:
        Path to the destination file

    Raises:
        DataAccessError: If the file cannot be copied
    """
    try:
        source_path = Path(source)
        dest_path = Path(destination)

        if dest_path.is_dir():
            dest_path = dest_path / source_path.name
        with atomic_path(dest_path) as tmp:
            shutil.copy2(source_path, tmp)
        return dest_path
    except Exception as e:
        raise DataAccessError(f"Failed to copy {source} to {destination}: {e}")


def read_csv(file_path: Union[str, Path], **kwargs) -> pd.DataFrame:
    """
    Read a CSV file to a pandas DataFrame.

    Args:
        file_path: Path to the CSV file
        **kwargs: Additional arguments to pass to pd.read_csv

    Returns:
        Pandas DataFrame

    Raises:
        DataAccessError: If the file cannot be read
    """
    try:
        return pd.read_csv(file_path, **kwargs)
    except Exception as e:
        raise DataAccessError(f"Failed to read CSV file {file_path}: {e}")


def _short_hostname() -> str:
    """The first label of this host's name, keeping only characters that are safe in a file name."""
    label = socket.gethostname().split(".")[0]
    cleaned = "".join(c for c in label if c.isascii() and (c.isalnum() or c in "-_"))
    return cleaned or "host"


def unique_temp_path(target: Union[str, Path]) -> Path:
    """A temporary name next to ``target`` that no other process on any host will choose.

    The name is ``.<name>.<host>.<pid>.<token>.tmp`` in the target's folder, so the final rename
    stays on one filesystem. The leading dot keeps the leftover of an interrupted write out of
    every folder listing (see ``visible_files``); host, pid and a random token keep two processes
    on different machines that share one NFS folder from writing to the same temporary file.
    """
    path = Path(target)
    return path.parent / f".{path.name}.{_short_hostname()}.{os.getpid()}.{secrets.token_hex(4)}.tmp"


def _replace_with_retry(source: Path, destination: Path) -> None:
    """``os.replace``, tried again on ``PermissionError`` (see REPLACE_ATTEMPTS)."""
    for attempt in range(REPLACE_ATTEMPTS):
        try:
            os.replace(source, destination)
            return
        except PermissionError:
            if attempt == REPLACE_ATTEMPTS - 1:
                raise
            time.sleep(REPLACE_RETRY_SECONDS / REPLACE_ATTEMPTS)


def _fsync_path(path: Path, directory: bool = False) -> None:
    """Flush ``path`` to stable storage; a folder that cannot be flushed (Windows, some shares) is skipped."""
    try:
        fd = os.open(path, os.O_RDONLY if directory else os.O_RDWR)
    except OSError:
        if directory:
            return
        raise
    try:
        os.fsync(fd)
    except OSError:
        if not directory:
            raise
    finally:
        os.close(fd)


@contextmanager
def atomic_path(target: Union[str, Path], fsync: bool = False) -> Iterator[Path]:
    """Yield a temporary path to write; on a clean exit it replaces ``target`` in one rename.

    Readers see either the previous content or the complete new content, never a partial file.
    The temporary file is created empty (``O_EXCL``, mode 0o666 less the umask, so the group
    permissions of a shared folder apply) and takes the mode of an existing target just before
    the rename, so a read-only target is replaced too. A symlinked
    target is resolved and the file it points to is replaced, so the link itself stays. When the
    body raises, the rename fails or the process is interrupted (``KeyboardInterrupt``), the
    temporary file is removed and the target is left as it was. With ``fsync`` the file and its
    folder are flushed to disk, for state files that should survive a power loss.
    """
    final = Path(target)
    if final.is_symlink():
        final = final.resolve()
    final.parent.mkdir(parents=True, exist_ok=True)
    tmp = unique_temp_path(final)
    os.close(os.open(tmp, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666))
    replaced = False
    try:
        yield tmp
        if fsync:
            _fsync_path(tmp)
        # The target's mode is copied only now, so a read-only target (0o444) can still be
        # replaced: the body and the fsync write to a temporary file this process can open.
        try:
            os.chmod(tmp, stat.S_IMODE(final.stat().st_mode))
        except FileNotFoundError:
            pass
        _replace_with_retry(tmp, final)
        replaced = True
        if fsync:
            _fsync_path(final.parent, directory=True)
    finally:
        if not replaced:
            try:
                tmp.unlink(missing_ok=True)
            except OSError as e:
                logger.warning(f"Could not remove temporary file {tmp}: {e}")


def write_text_atomic(path: Union[str, Path], text: str, encoding: str = "utf-8", fsync: bool = False) -> Path:
    """Write ``text`` to ``path`` through ``atomic_path`` and return ``path``."""
    with atomic_path(path, fsync=fsync) as tmp:
        with open(tmp, "w", encoding=encoding) as handle:
            handle.write(text)
    return Path(path)


@contextmanager
def open_atomic(
    path: Union[str, Path],
    mode: str = "w",
    encoding: Optional[str] = None,
    newline: Optional[str] = None,
    fsync: bool = False,
) -> Iterator[IO[Any]]:
    """Open a handle for writing ``path`` in ``mode`` ("w" or "wb"), published by ``atomic_path`` on close.

    For writers that stream (``csv.writer``, ``json.dump``, a gzip stream): the target keeps its
    old content until the ``with`` block finishes without an error.
    """
    if mode not in ("w", "wb"):
        raise ValueError(f"open_atomic writes a whole file; mode must be 'w' or 'wb', not {mode!r}")
    if mode == "w" and encoding is None:
        encoding = "utf-8"
    with atomic_path(path, fsync=fsync) as tmp:
        with open(tmp, mode, encoding=encoding, newline=newline) as handle:
            yield handle


def write_bytes_atomic(path: Union[str, Path], data: bytes, fsync: bool = False) -> Path:
    """Write ``data`` to ``path`` through ``atomic_path`` and return ``path``."""
    with atomic_path(path, fsync=fsync) as tmp:
        with open(tmp, "wb") as handle:
            handle.write(data)
    return Path(path)


def _compression_for(path: Path) -> Any:
    """The pandas ``compression`` argument implied by the suffix of ``path``, or None."""
    method = _COMPRESSION_BY_SUFFIX.get(path.suffix.lower())
    if method == "zip":
        # Without a name the archive member would be named after the temporary file.
        return {"method": "zip", "archive_name": path.stem}
    return method


def write_csv(df: pd.DataFrame, file_path: Union[str, Path], **kwargs) -> None:
    """
    Write a pandas DataFrame to a CSV file atomically (see ``atomic_path``).

    Compression follows the suffix of the real target (``.gz``, ``.bz2``, ``.zip``, ``.xz``,
    ``.zst``) as it would for a direct ``to_csv``, unless ``compression`` is passed explicitly.

    Args:
        df: Pandas DataFrame to write
        file_path: Path to the output CSV file
        **kwargs: Additional arguments to pass to df.to_csv

    Raises:
        DataAccessError: If the file cannot be written
    """
    output_path = Path(file_path)
    if "compression" not in kwargs:
        kwargs["compression"] = _compression_for(output_path)
    try:
        with atomic_path(output_path) as tmp:
            df.to_csv(tmp, **kwargs)
    except Exception as e:
        raise DataAccessError(f"Failed to write CSV file {file_path}: {e}")
