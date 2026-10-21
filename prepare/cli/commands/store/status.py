"""
`store_status`: dataset and byte counts, and stale projects, for the resolved store.
"""

import argparse
from pathlib import Path
from typing import Any, Dict, List

from prepare.cli.base import BaseCommand
from prepare.cli.commands.store._shared import _stale_project_row, _required_store
from prepare.core.exceptions import DataAccessError
from prepare.store.catalog import Catalog
from prepare.store.layout import read_marker, store_paths
from prepare.store.usage import stale_projects


class StoreStatusCommand(BaseCommand):
    """Command to report on the shared data store: dataset counts, bytes, and stale projects."""

    @property
    def name(self) -> str:
        """Return the command name."""
        return "store_status"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return "Report dataset counts, bytes and stale projects for the shared data store"

    @property
    def group(self) -> str:
        """Return the pipeline-step group this command is listed under."""
        return "Store"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's arguments."""
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument(
            "--registry",
            default=None,
            help="Path to the project registry file (defaults to the nearest prepare_registry.json)",
        )
        parser.add_argument("--json", action="store_true", help="Emit the report as JSON")
        parser.add_argument(
            "--details",
            action="store_true",
            help="Also list every dataset in the store (formerly --verbose; -v now sets the log level)",
        )
        # Before the global -v/--verbose existed, --verbose meant --details. -v is now the log level.
        parser.add_argument("--verbose", dest="details_alias", action="store_true", help=argparse.SUPPRESS)

    @staticmethod
    def _dataset_counts_and_bytes(catalog: Catalog) -> Any:
        """Dataset counts by state and total bytes, without the placeholder rows.

        A usage row recorded before its dataset was catalogued inserts a ``state="unknown"``
        row that stands for no files at all; counting it as a dataset would overstate what
        the store holds.
        """
        rows = catalog.conn.execute(
            "SELECT state, COUNT(*) AS n, COALESCE(SUM(bytes_total), 0) AS bytes FROM datasets "
            "WHERE state IS NOT 'unknown' GROUP BY state"
        ).fetchall()
        counts = {row["state"]: row["n"] for row in rows}
        bytes_total = sum(row["bytes"] for row in rows)
        return counts, bytes_total

    @staticmethod
    def _datasets_list(catalog: Catalog) -> List[Dict[str, Any]]:
        rows = catalog.conn.execute(
            "SELECT accession, state, bytes_total FROM datasets WHERE state IS NOT 'unknown' ORDER BY accession"
        ).fetchall()
        result = []
        for row in rows:
            project_count = len(catalog.projects_for(row["accession"]))
            result.append(
                {
                    "accession": row["accession"],
                    "state": row["state"],
                    "bytes": row["bytes_total"] or 0,
                    "projects": project_count,
                }
            )
        return result

    @staticmethod
    def _stale_project_rows(catalog: Catalog) -> List[Dict[str, Any]]:
        return [_stale_project_row(row) for row in stale_projects(catalog)]

    def _wants_details(self, args: argparse.Namespace) -> bool:
        """True for ``--details``; the old ``--verbose`` spelling still works and logs a deprecation line."""
        if getattr(args, "details_alias", False):
            if not getattr(args, "_alias_warned", False):
                args._alias_warned = True
                self.logger.warning("store_status --verbose is deprecated and will be removed; use --details")
            return True
        return bool(args.details)

    def _build_report(self, root: Path, args: argparse.Namespace) -> Dict[str, Any]:
        paths = store_paths(root)
        marker = read_marker(root) or {}
        with Catalog(paths) as catalog:
            counts, bytes_total = self._dataset_counts_and_bytes(catalog)
            project_count = catalog.conn.execute("SELECT COUNT(*) AS n FROM projects").fetchone()["n"]
            stale = self._stale_project_rows(catalog)
            report: Dict[str, Any] = {
                "root": str(root),
                "id": marker.get("id"),
                "datasets": counts,
                "bytes_total": bytes_total,
                "projects": project_count,
                "stale_projects": stale,
            }
            if self._wants_details(args):
                report["datasets_list"] = self._datasets_list(catalog)
        return report

    def _print_report(self, report: Dict[str, Any], verbose: bool) -> None:
        self.emit("Store")
        self.emit("=====")
        self.emit(f"  Root         : {report['root']}")
        self.emit(f"  Id           : {report['id']}")
        self.emit(f"  Bytes total  : {report['bytes_total']}")
        self.emit(f"  Projects     : {report['projects']}")
        for state, count in sorted(report["datasets"].items()):
            self.emit(f"  {state:<12s}: {count}")
        if report["stale_projects"]:
            self.emit("  Stale projects:")
            for entry in report["stale_projects"]:
                self.emit(f"    {entry['name']} on {entry['hostname']} ({entry['reason']}: {entry['registry']})")
        if verbose:
            self.emit("\nDatasets")
            self.emit("========")
            for entry in report.get("datasets_list", []):
                self.emit(
                    f"  {entry['accession']:<15s} {entry['state']:<10s} "
                    f"{entry['bytes']:>12} bytes  {entry['projects']} project(s)"
                )

    def execute(self, args: argparse.Namespace) -> int:
        """Run the command; return the exit code."""
        ctx = _required_store(args, self)
        if isinstance(ctx, int):
            return ctx
        root = ctx.root

        try:
            report = self._build_report(root, args)
        except DataAccessError as e:
            return self.fail(e, self.name)

        if args.json:
            self.emit_json(report)
        else:
            self._print_report(report, self._wants_details(args))
        return 0
