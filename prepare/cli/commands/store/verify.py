"""
`store_verify`: check a dataset's files against its recorded size, md5 or NCBI spot count.
"""

import argparse
import logging
import re
from pathlib import Path
from typing import Any, Collection, Dict, List, Optional

from prepare.cli.base import BaseCommand, refuse_bad_accessions
from prepare.cli.commands.store._shared import _no_store_hint
from prepare.core.constants import SRA_ACCESSION_PATTERN
from prepare.core.exceptions import DataAccessError
from prepare.data import registry_blocks as rb
from prepare.data import run_log
from prepare.data.file_io import visible_files
from prepare.data.registry import load_registry, project_root
from prepare.data.sra import count_fastq_reads, fastq_files, orphan_fastq, primary_fastq, verify_download
from prepare.store.catalog import catalog_write
from prepare.store.layout import StorePaths, sidecar_path, sra_dir, store_paths
from prepare.store.locks import LockHeld, dataset_lock
from prepare.store.resolve import resolve_store_root
from prepare.store.sidecar import (
    Sidecar,
    build_sidecar,
    md5_file,
    ncbi_from_metadata_xml,
    read_sidecar,
    write_sidecar,
)

logger = logging.getLogger(__name__)


class StoreVerifyCommand(BaseCommand):
    """Command to verify store datasets against their sidecars, optionally by md5 or spot count.

    A check never takes any lock (it only reads). ``--fix-state``'s write-back does, one
    accession at a time and without waiting: see ``_fix_state_locked``.
    """

    @property
    def name(self) -> str:
        """Return the command name."""
        return "store_verify"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return "Verify store datasets against their sidecars (bytes, optionally md5 and NCBI spot counts)"

    @property
    def group(self) -> str:
        """Return the pipeline-step group this command is listed under."""
        return "Store"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's arguments."""
        parser.add_argument("accessions", nargs="*", help="Accessions to verify (default: every dataset in the store)")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument(
            "--registry",
            default=None,
            help="Path to the project registry file (defaults to the nearest prepare_registry.json)",
        )
        parser.add_argument("--md5", action="store_true", help="Recompute and compare each file's md5 checksum")
        parser.add_argument(
            "--spots", action="store_true", help="Compare the read count against NCBI's recorded spot count"
        )
        parser.add_argument(
            "--fix-state",
            action="store_true",
            help="Rewrite the sidecar state and catalogue entry when a check finds a mismatch",
        )
        parser.add_argument(
            "--rescan",
            action="store_true",
            help=(
                "Rebuild each dataset's file list from the files on disk before checking "
                "(use after files were added or removed by hand)"
            ),
        )
        parser.add_argument("--json", action="store_true", help="Emit the report as JSON")

    @staticmethod
    def _accessions_to_check(args: argparse.Namespace, paths: StorePaths) -> List[str]:
        """The accessions named on the command line, else every accession folder in the store.

        A folder in ``sra/`` not named like an SRA accession is not a dataset: it is logged as
        skipped, as ``store_reindex`` does, rather than stopping the check of the others.
        """
        if args.accessions:
            return list(args.accessions)
        accessions = []
        for folder in visible_files(paths.sra, dirs=True):
            if re.fullmatch(SRA_ACCESSION_PATTERN, folder.name):
                accessions.append(folder.name)
            else:
                logger.warning("%s in %s is not named like an SRA accession; skipped", folder.name, paths.sra)
        return accessions

    @staticmethod
    def _check_bytes_and_md5(accession: str, store_dir: Path, sidecar: Sidecar, check_md5: bool):
        """Returns ``(bytes_ok, md5_ok, detail)``; ``detail`` names the first file and reason
        for any bytes or md5 mismatch (or missing/unreadable file), else None."""
        bytes_ok = True
        md5_ok: Optional[bool] = True if check_md5 else None
        detail: Optional[str] = None
        for entry in sidecar.files:
            name = str(entry.get("name"))
            file_path = store_dir / name
            if not file_path.is_file():
                bytes_ok = False
                if check_md5:
                    md5_ok = False
                detail = detail or f"{name}: missing"
                continue
            try:
                actual_bytes = file_path.stat().st_size
                if actual_bytes != entry.get("bytes"):
                    bytes_ok = False
                    detail = detail or f"{name}: size mismatch ({actual_bytes} vs {entry.get('bytes')} bytes)"
                if check_md5:
                    actual_md5 = md5_file(file_path)
                    if actual_md5 != entry.get("md5"):
                        md5_ok = False
                        detail = detail or f"{name}: md5 mismatch"
            except OSError as e:
                logger.warning("%s: could not read %s: %s", accession, file_path, e)
                bytes_ok = False
                if check_md5:
                    md5_ok = False
                detail = detail or f"{name}: {e}"
        return bytes_ok, md5_ok, detail

    @staticmethod
    def _missing_result(accession: str, detail: Optional[str] = None) -> Dict[str, Any]:
        """The result shape for an accession this check cannot find anything usable for:
        no sidecar at all, or (with ``--rescan``) a store folder with no FASTQ files left to
        describe. ``sidecar`` is left None so ``_fix_state`` leaves whatever is on disk alone.
        ``state_before`` mirrors ``state`` here and is never rewritten by a later fix, so a
        caller can always tell what this accession looked like before ``--fix-state`` ran."""
        return {
            "accession": accession,
            "state": "missing",
            "state_before": "missing",
            "bytes_ok": False,
            "md5_ok": None,
            "verdict": "missing",
            "sidecar": None,
            "spots_verdict": None,
            "spots_ratio": None,
            "mismatch_detail": detail,
            "rescanned": False,
            "ncbi_found": None,
        }

    def _find_ncbi_spots(self, accession: str, metadata_dirs: List[Path]) -> Optional[Dict[str, Any]]:
        """The ``ncbi`` block read from the first ``<accession>_metadata.xml`` in
        ``metadata_dirs`` that both exists and yields a spot count, else None."""
        for metadata_dir in metadata_dirs:
            xml_path = metadata_dir / f"{accession}_metadata.xml"
            if not xml_path.is_file():
                continue
            found = ncbi_from_metadata_xml(xml_path)
            if found.get("spots"):
                self.logger.info("%s: spot count read from %s", accession, xml_path)
                return found
        return None

    def _verify_one(
        self,
        accession: str,
        paths: StorePaths,
        check_md5: bool,
        check_spots: bool,
        metadata_dirs: Optional[List[Path]] = None,
        rescan: bool = False,
    ) -> Dict[str, Any]:
        store_dir = sra_dir(paths, accession)
        sc_path = sidecar_path(paths, accession)
        sidecar = read_sidecar(sc_path)
        if sidecar is None:
            return self._missing_result(accession)

        rescanned = False
        if rescan:
            if not fastq_files(store_dir):
                # The folder this sidecar describes has no FASTQ files left (deleted, moved,
                # or never populated): rebuilding from it would silently overwrite a real
                # file list with an empty one, so report the dataset missing instead and
                # leave whatever is on disk untouched.
                return self._missing_result(accession, f"{store_dir}: no FASTQ files found to rescan")
            rebuilt = build_sidecar(
                accession,
                store_dir,
                sidecar.ncbi,
                sidecar.tool_version,
                sidecar.compression,
                tool=sidecar.tool,
                downloaded=sidecar.downloaded,
            )
            rebuilt.stats, rebuilt.stats_computed = sidecar.stats, sidecar.stats_computed
            sidecar = rebuilt
            rescanned = True

        bytes_ok, md5_ok, detail = self._check_bytes_and_md5(accession, store_dir, sidecar, check_md5)

        spots_verdict = None
        spots_ratio = None
        ncbi_found = None
        if check_spots:
            expected_spots = sidecar.ncbi.get("spots")
            if not expected_spots:
                ncbi_found = self._find_ncbi_spots(accession, metadata_dirs or [])
                if ncbi_found is not None:
                    expected_spots = ncbi_found.get("spots")
            try:
                verify = verify_download(accession, store_dir, expected_spots)
                spots_verdict = verify["verdict"]
                spots_ratio = verify["ratio"]
            except (EOFError, OSError) as e:
                self.logger.warning("%s: could not verify read counts: %s", accession, e)
                bytes_ok = False
                spots_verdict = "corrupt"
                detail = detail or f"read count verification failed: {e}"

        if not bytes_ok or (check_md5 and md5_ok is False):
            verdict = "corrupt"
        elif check_spots and spots_verdict == "truncated":
            verdict = "truncated"
        elif check_spots and spots_verdict == "unverified":
            verdict = "unverified"
        else:
            verdict = "ok"

        return {
            "accession": accession,
            "state": sidecar.state,
            "state_before": sidecar.state,
            "bytes_ok": bytes_ok,
            "md5_ok": md5_ok,
            "verdict": verdict,
            "sidecar": sidecar,
            "spots_verdict": spots_verdict,
            "spots_ratio": spots_ratio,
            "mismatch_detail": detail,
            "rescanned": rescanned,
            "ncbi_found": ncbi_found,
        }

    def _mark_failed(self, result: Dict[str, Any], paths: StorePaths, sidecar: Sidecar, error: str) -> None:
        """Rewrite ``sidecar`` to ``state="failed"`` with ``error``, unless it already is.

        The result is reported as ``corrupt`` either way, so the printed table and the exit
        status agree with the state the sidecar is left in. Sets ``result["fix"]`` to
        ``"updated"`` when the sidecar was actually rewritten, else ``"unchanged"`` (it was
        already ``"failed"`` with this same error).
        """
        state_before = sidecar.state
        result["verdict"] = "corrupt"
        result["state"] = "failed"
        if sidecar.state == "failed" and sidecar.error == error:
            result["fix"] = {
                "action": "unchanged",
                "state_before": state_before,
                "state_after": "failed",
                "error": error,
            }
            return
        sidecar.state = "failed"
        sidecar.error = error
        write_sidecar(sidecar_path(paths, result["accession"]), sidecar)
        with catalog_write(paths) as catalog:
            catalog.upsert_dataset(sidecar)
        result["state"] = "failed"
        result["fix"] = {"action": "updated", "state_before": state_before, "state_after": "failed", "error": error}

    def _read_through_error(self, store_dir: Path, sidecar: Sidecar, skip: Collection[str] = ()) -> Optional[str]:
        """Open and decompress every file ``sidecar.files`` records, except those named in
        ``skip`` (already read by a spot comparison), returning the first one's error
        (name-prefixed) when it cannot be read, else None.

        A file that matches its recorded size (and md5, if checked) can still be a truncated
        or corrupt gzip stream, since neither check opens it, and a spot comparison reads only
        the mate 1 (or single-end) file and the unpaired-read file. ``--fix-state`` must not
        promote a dataset to ``"complete"`` before every recorded file has been read through.
        """
        for entry in sidecar.files:
            name = str(entry.get("name"))
            if name in skip:
                continue
            try:
                count_fastq_reads(store_dir / name)
            except (EOFError, OSError) as e:
                self.logger.warning("%s: could not read %s: %s", sidecar.accession, name, e)
                return f"{name}: {e}"
        return None

    def _compare_known_spots(self, result: Dict[str, Any], paths: StorePaths, sidecar: Sidecar) -> Optional[str]:
        """Compare the reads on disk against the spot count ``sidecar.ncbi`` already records,
        for a run without ``--spots``; the outcome is written into ``result`` as a ``--spots``
        check would have written it. Returns a read error (name-prefixed where possible) when
        a file cannot be read, else None."""
        accession = result["accession"]
        store_dir = sra_dir(paths, accession)
        try:
            verify = verify_download(accession, store_dir, sidecar.ncbi.get("spots"))
        except (EOFError, OSError) as e:
            self.logger.warning("%s: could not verify read counts: %s", accession, e)
            return self._read_through_error(store_dir, sidecar) or f"read count verification failed: {e}"
        result["spots_verdict"] = verify["verdict"]
        result["spots_ratio"] = verify["ratio"]
        if verify["verdict"] == "truncated":
            result["verdict"] = "truncated"
        return None

    def _unread_file_error(self, result: Dict[str, Any], paths: StorePaths, sidecar: Sidecar) -> Optional[str]:
        """Read through every recorded file that no spot comparison has read in this run."""
        store_dir = sra_dir(paths, result["accession"])
        skip: set = set()
        if result.get("spots_verdict") is not None:
            skip = {p.name for p in (primary_fastq(store_dir), orphan_fastq(store_dir)) if p is not None}
        return self._read_through_error(store_dir, sidecar, skip)

    def _fix_state(self, result: Dict[str, Any], paths: StorePaths) -> None:
        """Rewrite the sidecar's state (and completeness) to match what this check found.

        A bytes or md5 mismatch always wins: the file itself is wrong, so the dataset is
        ``"failed"`` with an error naming the first mismatch, regardless of what the spots
        check says (a corrupt file can still happen to contain the right number of reads). A
        sidecar recording no files at all is never promoted either, since there is nothing to
        have verified. A spot count found in a metadata XML (``result["ncbi_found"]``, from
        ``_verify_one``'s fallback search) is first written into the sidecar's ``ncbi`` block
        so later runs read it straight from there. When ``--spots`` was not requested but the
        sidecar records a spot count, the reads are compared against it here, exactly as
        ``--spots`` would have. Every recorded file the spot comparison did not read (mate 2,
        or every file when no comparison ran) is then read through; a read failure is recorded
        as the new error, the same as a bytes/md5 mismatch. Once all of that checks out, a
        spots verdict of ``complete``/``truncated`` decides ``complete``/``partial``; with no
        spot count anywhere, the dataset is promoted to ``complete`` with an ``unverified``
        completeness, clearing any stale error. Every return path sets ``result["fix"]``
        (``"updated"`` or ``"unchanged"``; ``_mark_failed`` sets its own on the failure paths).
        """
        sidecar = result.get("sidecar")
        if sidecar is None:
            result["fix"] = {"action": "unchanged", "state_before": "missing", "state_after": "missing", "error": None}
            return

        if not result.get("bytes_ok") or result.get("md5_ok") is False:
            self._mark_failed(result, paths, sidecar, result.get("mismatch_detail") or "verify: bytes or md5 mismatch")
            return

        if not sidecar.files:
            self._mark_failed(result, paths, sidecar, "no files recorded")
            return

        if result.get("ncbi_found") and not sidecar.ncbi.get("spots"):
            sidecar.ncbi = dict(result["ncbi_found"])
        if result.get("spots_verdict") == "corrupt":
            result["fix"] = {
                "action": "unchanged",
                "state_before": sidecar.state,
                "state_after": sidecar.state,
                "error": None,
            }
            return
        read_error = None
        if result.get("spots_verdict") is None and sidecar.ncbi.get("spots"):
            read_error = self._compare_known_spots(result, paths, sidecar)
        read_error = read_error or self._unread_file_error(result, paths, sidecar)
        if read_error is not None:
            self._mark_failed(result, paths, sidecar, read_error)
            return

        state_before = sidecar.state
        spots_verdict = result.get("spots_verdict")
        state_for_verdict = {"complete": "complete", "truncated": "partial"}
        if spots_verdict in state_for_verdict:
            new_state = state_for_verdict[spots_verdict]
            completeness = {"method": "spots", "ratio": result.get("spots_ratio"), "verdict": spots_verdict}
        else:
            new_state = "complete"
            completeness = {"method": "unverified", "ratio": None, "verdict": "unverified"}
        changed = (
            result.get("rescanned")
            or new_state != sidecar.state
            or sidecar.completeness != completeness
            or sidecar.error
        )
        if not changed:
            result["fix"] = {
                "action": "unchanged",
                "state_before": state_before,
                "state_after": state_before,
                "error": None,
            }
            return
        sidecar.state = new_state
        sidecar.completeness = completeness
        sidecar.error = None
        write_sidecar(sidecar_path(paths, result["accession"]), sidecar)
        with catalog_write(paths) as catalog:
            catalog.upsert_dataset(sidecar)
        result["state"] = new_state
        result["fix"] = {"action": "updated", "state_before": state_before, "state_after": new_state, "error": None}

    def _fix_state_locked(self, result: Dict[str, Any], paths: StorePaths) -> None:
        """Apply ``_fix_state`` only under ``accession``'s lock, and only if nothing else has
        touched its sidecar since ``_verify_one`` read it.

        Held without waiting (``blocking=False``): another run already working on this
        accession, or one that grabs the lock first, means the write-back would race a write
        of its own, so it is skipped rather than raced. The lock alone is not enough, since a
        write could have already happened and finished before this call even reaches the
        lock; the sidecar is re-read under it and compared against what ``_verify_one`` saw.
        ``downloaded`` is compared for every run; ``files`` only when this run did not pass
        ``--rescan``, since a rescan's whole point is to make ``result["sidecar"].files``
        differ from what is still on disk in ``downloaded``'s sidecar file. When the write-back
        is skipped here, ``result["fix"]`` records why (``"skipped-changed"`` or
        ``"skipped-in-use"``) instead of ``_fix_state``'s own ``"updated"``/``"unchanged"``.
        """
        accession = result["accession"]
        sidecar_before = result.get("sidecar")
        if sidecar_before is None:
            self._fix_state(result, paths)
            return
        try:
            with dataset_lock(paths, accession, blocking=False):
                current = read_sidecar(sidecar_path(paths, accession))
                if current is None or current.downloaded != sidecar_before.downloaded:
                    changed = True
                else:
                    changed = not result.get("rescanned") and current.files != sidecar_before.files
                if changed:
                    self.logger.warning("%s: sidecar changed since the check; not updated", accession)
                    result["fix"] = {
                        "action": "skipped-changed",
                        "state_before": sidecar_before.state,
                        "state_after": current.state if current is not None else None,
                        "error": None,
                    }
                    return
                self._fix_state(result, paths)
        except LockHeld:
            self.logger.warning("%s: in use, not updated", accession)
            result["fix"] = {
                "action": "skipped-in-use",
                "state_before": sidecar_before.state,
                "state_after": sidecar_before.state,
                "error": None,
            }

    def _print_table(self, results: List[Dict[str, Any]]) -> None:
        self.emit(f"{'accession':<15s} {'state':<10s} {'bytes_ok':<9s} {'md5_ok':<7s} verdict")
        for r in results:
            md5_col = "-" if r["md5_ok"] is None else str(r["md5_ok"])
            self.emit(f"{r['accession']:<15s} {r['state']:<10s} {str(r['bytes_ok']):<9s} {md5_col:<7s} {r['verdict']}")

    def _print_fixes(self, results: List[Dict[str, Any]]) -> None:
        """Print one line per dataset ``--fix-state`` actually rewrote (``fix.action ==
        "updated"``); a dataset left unchanged, or skipped because it was locked or had
        changed since the check, prints nothing here."""
        for r in results:
            fix = r.get("fix")
            if fix and fix["action"] == "updated":
                self.emit(f"{r['accession']}: fixed ({fix['state_before']} -> {fix['state_after']})")

    @staticmethod
    def _json_row(result: Dict[str, Any]) -> Dict[str, Any]:
        """One dataset's check result in JSON-safe form: every field but the ``Sidecar``
        object and the raw ``ncbi_found`` lookup, neither of which a script reading
        ``--json`` output needs or can decode."""
        return {key: value for key, value in result.items() if key not in ("sidecar", "ncbi_found")}

    def _json_document(self, root: Path, args: argparse.Namespace, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """The one document ``--json`` prints: the checks that ran, a tally of verdicts, every
        dataset's row, and the datasets ``--fix-state`` actually rewrote. ``counts`` and
        ``fixed`` are meant to be read back out of this document by a caller summarizing the
        run (a project's run log, say) without re-deriving them from ``datasets``."""
        counts: Dict[str, int] = {}
        for r in results:
            counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
        fixed = [
            {
                "accession": r["accession"],
                "state_before": r["fix"]["state_before"],
                "state_after": r["fix"]["state_after"],
                "error": r["fix"]["error"],
            }
            for r in results
            if r.get("fix", {}).get("action") == "updated"
        ]
        return {
            "root": str(root),
            "checks": {
                "md5": bool(args.md5),
                "spots": bool(args.spots),
                "rescan": bool(args.rescan),
                "fix_state": bool(args.fix_state),
            },
            "counts": counts,
            "datasets": [self._json_row(r) for r in results],
            "fixed": fixed,
        }

    def records_run(self, args: argparse.Namespace) -> bool:
        """Every run is added to the project's run log (one with no project registry writes nothing)."""
        return True

    @staticmethod
    def _note_run(args: argparse.Namespace, document: Dict[str, Any]) -> None:
        """The verdict tally and the number fixed for the run-log summary; one row per dataset for its detail."""
        datasets = document["datasets"]
        run_log.note_run(
            args, summary={"datasets": len(datasets), "counts": document["counts"], "fixed": len(document["fixed"])}
        )
        rows = {
            row["accession"]: {
                "verdict": row["verdict"],
                "state_before": row.get("state_before"),
                "state_after": row.get("state"),
                "fix": (row.get("fix") or {}).get("action"),
            }
            for row in datasets
        }
        run_log.note_rows(args, rows, "verify")

    def execute(self, args: argparse.Namespace) -> int:
        """Run the command; return the exit code."""
        if (rc := refuse_bad_accessions(self.logger, args.accessions)) is not None:
            return rc
        try:
            registry = load_registry(args.registry)
            root = resolve_store_root(args.data_root, rb.store_block(registry).root)
        except DataAccessError as e:
            return self.fail(e, self.name)

        if root is None:
            _no_store_hint(args.json)
            return 1

        try:
            paths = store_paths(root)
            metadata_dirs = [paths.metadata]
            if registry.path is not None and registry.path.is_file():
                metadata_dirs.append(project_root(registry) / "metadata")
            accessions = self._accessions_to_check(args, paths)
            results = []
            for accession in accessions:
                result = self._verify_one(
                    accession, paths, args.md5, args.spots, metadata_dirs=metadata_dirs, rescan=args.rescan
                )
                if args.fix_state:
                    self._fix_state_locked(result, paths)
                results.append(result)
        except DataAccessError as e:
            return self.fail(e, self.name)

        document = self._json_document(root, args, results)
        self._note_run(args, document)
        if args.json:
            self.emit_json(document)
        else:
            self._print_table(results)
            if args.fix_state:
                self._print_fixes(results)
        failed = any(r["verdict"] in ("corrupt", "truncated", "missing") for r in results)
        return 1 if failed else 0
