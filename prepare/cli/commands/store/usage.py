"""
`store_usage`: which projects and genomes used an accession, or which datasets are unused.
"""

import argparse
from typing import Any, Dict, List, Tuple

from prepare.cli.base import BaseCommand
from prepare.cli.commands.store._shared import _required_store
from prepare.core.exceptions import DataAccessError
from prepare.store.catalog import Catalog
from prepare.store.layout import store_paths


class StoreUsageCommand(BaseCommand):
    """Command to report catalogue usage: by accession, project, organism, or store-wide."""

    _COLUMNS: Dict[str, List[Tuple[str, str]]] = {
        "accession": [
            ("project_name", "project"),
            ("project_id", "id"),
            ("genome_id", "genome_id"),
            ("stage", "stage"),
            ("first_used", "first_used"),
            ("last_used", "last_used"),
        ],
        "project": [
            ("accession", "accession"),
            ("genome_id", "genome_id"),
            ("stage", "stage"),
            ("last_used", "last_used"),
        ],
        "organism": [("accession", "accession"), ("project_name", "project"), ("stage", "stage")],
        "unused": [("accession", "accession"), ("state", "state"), ("bytes", "bytes")],
        "bytes-by-organism": [("genome_id", "genome_id"), ("datasets", "datasets"), ("bytes", "bytes")],
    }

    @property
    def name(self) -> str:
        """Return the command name."""
        return "store_usage"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return "Report catalogue usage by accession, project, organism, or store-wide"

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
        selector = parser.add_mutually_exclusive_group(required=True)
        selector.add_argument("--accession", default=None, help="Every project that has used this accession")
        selector.add_argument("--project", default=None, help="Every dataset this project (by name or id) has used")
        selector.add_argument("--organism", default=None, help="Every dataset used for this target genome id")
        selector.add_argument(
            "--unused", action="store_true", help="Datasets in the store with no recorded usage at all"
        )
        selector.add_argument(
            "--bytes-by-organism", action="store_true", help="Total bytes and dataset counts, grouped by genome id"
        )
        parser.add_argument("--json", action="store_true", help="Emit the report as JSON")

    # --------------------------------------------------------------- queries

    @staticmethod
    def _resolve_project_id(catalog: Catalog, value: str) -> Tuple[str, List[str]]:
        """Resolve ``value`` (a project id or name) to the id to query.

        Returns ``(project_id, [])`` when ``value`` is itself a known project id, or matches
        exactly one project's name. Returns ``(value, [])`` unchanged when it matches no
        project at all (the caller's query then simply returns no rows). Returns
        ``(None-ish, ambiguous_ids)`` when ``value`` matches more than one project's name; the
        caller must treat a non-empty second element as an error.
        """
        row = catalog.conn.execute("SELECT project_id FROM projects WHERE project_id = ?", (value,)).fetchone()
        if row is not None:
            return row["project_id"], []
        rows = catalog.conn.execute(
            "SELECT project_id FROM projects WHERE name = ? ORDER BY project_id", (value,)
        ).fetchall()
        ids = [r["project_id"] for r in rows]
        if len(ids) == 1:
            return ids[0], []
        if len(ids) > 1:
            return value, ids
        return value, []

    @staticmethod
    def _rows_for_accession(catalog: Catalog, accession: str) -> List[Dict[str, Any]]:
        rows = catalog.conn.execute(
            """
            SELECT p.name AS project_name, p.project_id AS project_id, u.genome_id AS genome_id,
                   u.stage AS stage, u.first_used AS first_used, u.last_used AS last_used
            FROM usage u
            JOIN projects p ON p.project_id = u.project_id
            WHERE u.accession = ?
            ORDER BY p.project_id, u.genome_id, u.stage
            """,
            (accession,),
        ).fetchall()
        return [dict(row) for row in rows]

    @staticmethod
    def _rows_for_project(catalog: Catalog, project_id: str) -> List[Dict[str, Any]]:
        rows = catalog.conn.execute(
            """
            SELECT DISTINCT u.accession AS accession, u.genome_id AS genome_id, u.stage AS stage,
                   u.last_used AS last_used
            FROM usage u
            WHERE u.project_id = ?
            ORDER BY u.accession, u.genome_id, u.stage
            """,
            (project_id,),
        ).fetchall()
        return [dict(row) for row in rows]

    @staticmethod
    def _rows_for_organism(catalog: Catalog, genome_id: str) -> List[Dict[str, Any]]:
        rows = catalog.conn.execute(
            """
            SELECT DISTINCT u.accession AS accession, p.name AS project_name, u.stage AS stage
            FROM usage u
            JOIN projects p ON p.project_id = u.project_id
            WHERE u.genome_id = ?
            ORDER BY u.accession, p.name, u.stage
            """,
            (genome_id,),
        ).fetchall()
        return [dict(row) for row in rows]

    @staticmethod
    def _rows_unused(catalog: Catalog) -> List[Dict[str, Any]]:
        rows = catalog.conn.execute("""
            SELECT d.accession AS accession, d.state AS state, COALESCE(d.bytes_total, 0) AS bytes
            FROM datasets d
            LEFT JOIN usage u ON u.accession = d.accession
            WHERE u.accession IS NULL
            ORDER BY d.accession
            """).fetchall()
        return [dict(row) for row in rows]

    @staticmethod
    def _rows_bytes_by_organism(catalog: Catalog) -> List[Dict[str, Any]]:
        rows = catalog.conn.execute("""
            SELECT pair.genome_id AS genome_id, COUNT(DISTINCT pair.accession) AS datasets,
                   COALESCE(SUM(d.bytes_total), 0) AS bytes
            FROM (SELECT DISTINCT genome_id, accession FROM usage) pair
            JOIN datasets d ON d.accession = pair.accession
            GROUP BY pair.genome_id
            ORDER BY pair.genome_id
            """).fetchall()
        return [dict(row) for row in rows]

    # ---------------------------------------------------------------- print

    def _print_rows(self, selector: str, rows: List[Dict[str, Any]]) -> None:
        columns = self._COLUMNS[selector]
        self.emit("  ".join(f"{label:<15s}" for _, label in columns))
        for row in rows:
            self.emit("  ".join(f"{str(row.get(key, '')):<15s}" for key, _ in columns))

    # --------------------------------------------------------------- execute

    def execute(self, args: argparse.Namespace) -> int:
        """Run the command; return the exit code."""
        ctx = _required_store(args, self)
        if isinstance(ctx, int):
            return ctx
        root = ctx.root

        try:
            paths = store_paths(root)
            with Catalog(paths) as catalog:
                if args.accession:
                    selector = "accession"
                    rows = self._rows_for_accession(catalog, args.accession)
                elif args.project:
                    selector = "project"
                    project_id, ambiguous = self._resolve_project_id(catalog, args.project)
                    if ambiguous:
                        self.logger.error(
                            "Project name %r is ambiguous: %s", args.project, ", ".join(sorted(ambiguous))
                        )
                        return 1
                    rows = self._rows_for_project(catalog, project_id)
                elif args.organism:
                    selector = "organism"
                    rows = self._rows_for_organism(catalog, args.organism)
                elif args.unused:
                    selector = "unused"
                    rows = self._rows_unused(catalog)
                else:
                    selector = "bytes-by-organism"
                    rows = self._rows_bytes_by_organism(catalog)
        except DataAccessError as e:
            return self.fail(e, self.name)

        if args.json:
            self.emit_json({"selector": selector, "rows": rows})
        else:
            self._print_rows(selector, rows)
        return 0
