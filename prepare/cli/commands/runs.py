"""``prepare runs``: list the project's run log, show one run, diff two runs, trace one accession.

The run log itself is ``prepare.data.run_log``; the comparisons are in
``prepare.processing.run_diff``. This command only reads the log and records nothing itself.
It exits 0 normally, 1 when the project has no run log and 1 (through ``fail``) for a run
selector that matches no run or more than one.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from prepare.cli.base import BaseCommand, emit_error_json
from prepare.core.exceptions import DataAccessError, PrepareError, ValidationError
from prepare.data import run_log
from prepare.data.run_log import RunRecord
from prepare.processing.run_diff import (
    accession_history,
    detail_kept,
    detail_rows,
    diff_details,
    diff_summaries,
)

DEFAULT_LIMIT = 20
SUMMARY_WIDTH = 60
NOT_KEPT_NOTE = (
    f"a run keeps one only when its command noted one, and only the last {run_log.DETAILS_KEPT_PER_COMMAND} "
    "per command are kept"
)


def _limit(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a whole number, got {text!r}")
    if value < 0:
        raise argparse.ArgumentTypeError(f"expected a whole number, 0 or more, got {value}")
    return value


def _text(value: Any) -> str:
    """One value for a text table: JSON for structures, ``-`` for None."""
    if value is None:
        return "-"
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return str(value)


def _summary_text(summary: Dict[str, Any]) -> str:
    text = " ".join(f"{key}={_text(value)}" for key, value in summary.items())
    return text if len(text) <= SUMMARY_WIDTH else text[: SUMMARY_WIDTH - 3] + "..."


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> List[str]:
    """Left-aligned columns separated by two spaces; the last column is not padded."""
    widths = [max(len(row[i]) for row in [header, *rows]) for i in range(len(header) - 1)]
    lines = []
    for row in [header, *rows]:
        cells = [cell.ljust(width) for cell, width in zip(row, widths)] + [row[-1]]
        lines.append("  ".join(cells).rstrip())
    return lines


class RunsCommand(BaseCommand):
    """List the project's run log, show one run, compare two runs or follow one accession."""

    @property
    def name(self) -> str:
        """The command name."""
        return "runs"

    @property
    def help(self) -> str:
        """One line for the command list."""
        return "List recorded runs of this project, show or compare runs, follow one accession"

    @property
    def group(self) -> str:
        """Listed under Environment in the main help."""
        return "Environment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the selection flags, the three views (mutually exclusive), ``--registry`` and ``--json``."""
        parser.add_argument("--command", dest="run_command", metavar="NAME", help="Only runs of this command")
        parser.add_argument(
            "--limit",
            type=_limit,
            default=DEFAULT_LIMIT,
            metavar="N",
            help="Number of most recent runs to list; 0 lists all",
        )
        views = parser.add_mutually_exclusive_group()
        views.add_argument(
            "--show",
            metavar="RUN",
            help="Show one run: a run ID, a unique prefix of one, 'latest' or 'previous'",
        )
        views.add_argument(
            "--diff",
            nargs=2,
            metavar=("RUN_A", "RUN_B"),
            help="Compare the summaries and details of two runs (selected as for --show)",
        )
        views.add_argument("--accession", metavar="ACC", help="Follow one accession through the recorded runs")
        parser.add_argument(
            "--registry",
            default=None,
            help="Registry file (default: prepare_registry.json found upwards from here)",
        )
        parser.add_argument("--json", action="store_true", help="Write the result as one JSON document")

    def execute(self, args: argparse.Namespace) -> int:
        """Run the selected view; 0 normally, 1 without a run log or for a selector that does not resolve."""
        project = run_log.project_for(args.registry)
        if project is None or not (run_log.runs_dir(project) / run_log.RUNS_FILE).is_file():
            where = "no project registry found" if project is None else f"nothing under {run_log.runs_dir(project)}"
            return self._failed(args, DataAccessError(f"This project has no run log ({where})"))
        try:
            records = run_log.read_runs(project, args.run_command)
            if args.show:
                self._show(args, project, run_log.resolve_run(records, args.show))
            elif args.diff:
                first, second = (run_log.resolve_run(records, selector) for selector in args.diff)
                self._diff(args, project, first, second)
            elif args.accession is not None:
                self._accession(args, project, records)
            else:
                self._list(args, project, records)
        except (ValidationError, DataAccessError) as e:
            return self._failed(args, e)
        return 0

    def _failed(self, args: argparse.Namespace, error: PrepareError) -> int:
        """Report ``error`` through ``fail``, and as a JSON error document with ``--json``."""
        if args.json:
            emit_error_json(str(error))
        return self.fail(error, self.name)

    @staticmethod
    def _run_dict(project: Path, record: RunRecord) -> Dict[str, Any]:
        return {**record.to_dict(), "detail_kept": detail_kept(project, record)}

    # --- list ----------------------------------------------------------------------------------

    def _list(self, args: argparse.Namespace, project: Path, records: List[RunRecord]) -> None:
        shown = list(reversed(records))
        if args.limit:
            shown = shown[: args.limit]
        if args.json:
            self.emit_json(
                {
                    "project": str(project),
                    "command": args.run_command,
                    "total": len(records),
                    "runs": [self._run_dict(project, record) for record in shown],
                }
            )
            return
        if not records:
            self.emit(f"No runs of {args.run_command} are recorded" if args.run_command else "No runs are recorded")
            return
        rows = [
            [
                record.run_id,
                record.command,
                record.started,
                f"{record.seconds:.1f}",
                str(record.exit_code),
                "kept" if detail_kept(project, record) else "-",
                _summary_text(record.summary),
            ]
            for record in shown
        ]
        for line in _table(["RUN", "COMMAND", "STARTED (UTC)", "SECONDS", "EXIT", "DETAIL", "SUMMARY"], rows):
            self.emit(line)
        self.emit(f"{len(shown)} of {len(records)} run(s) shown, newest first")

    # --- show ----------------------------------------------------------------------------------

    def _show(self, args: argparse.Namespace, project: Path, record: RunRecord) -> None:
        kept = detail_kept(project, record)
        detail = run_log.read_detail(project, record) if kept else None
        if args.json:
            self.emit_json({"run": self._run_dict(project, record), "detail_kept": kept, "detail": detail})
            return
        for label, value in (
            ("Run", record.run_id),
            ("Command", record.command),
            ("Started", record.started),
            ("Finished", record.finished),
            ("Seconds", f"{record.seconds:.3f}"),
            ("Exit code", str(record.exit_code)),
            ("Version", record.version),
            ("Host", f"{record.host} (PID {record.pid})"),
            ("Arguments", " ".join(record.argv)),
        ):
            self.emit(f"{label + ':':<11}{value}")
        self.emit("Summary:" if record.summary else "Summary:  none")
        for key, value in record.summary.items():
            self.emit(f"  {key}: {_text(value)}")
        if detail is None:
            self.emit(f"Detail: not kept ({NOT_KEPT_NOTE})")
            return
        self.emit("Detail: kept")
        rows = detail_rows(detail)
        if not rows:
            for line in json.dumps(detail, indent=2, sort_keys=True).splitlines():
                self.emit(f"  {line}")
            return
        for key in sorted(rows):
            fields = " ".join(f"{field}={_text(value)}" for field, value in rows[key].items())
            self.emit(f"  {key}  {fields}")

    # --- diff ----------------------------------------------------------------------------------

    def _diff(self, args: argparse.Namespace, project: Path, first: RunRecord, second: RunRecord) -> None:
        if first.command != second.command:
            self.logger.warning("Comparing runs of different commands: %s and %s", first.command, second.command)
        kept = {"run_a": detail_kept(project, first), "run_b": detail_kept(project, second)}
        summary = diff_summaries(first.summary, second.summary)
        detail: Optional[Dict[str, Any]] = None
        if all(kept.values()):
            detail = diff_details(run_log.read_detail(project, first), run_log.read_detail(project, second))
        if args.json:
            self.emit_json(
                {
                    "run_a": self._run_dict(project, first),
                    "run_b": self._run_dict(project, second),
                    "summary": summary,
                    "detail_kept": kept,
                    "detail": detail,
                }
            )
            return
        self.emit(f"Run A: {first.run_id}  ({first.command}, started {first.started}, exit {first.exit_code})")
        self.emit(f"Run B: {second.run_id}  ({second.command}, started {second.started}, exit {second.exit_code})")
        self.emit("Summary:")
        if summary:
            rows = [
                [
                    row["key"],
                    _text(row["before"]),
                    _text(row["after"]),
                    _text(row["delta"]),
                    "yes" if row["changed"] else "no",
                ]
                for row in summary
            ]
            for line in _table(["  KEY", "RUN A", "RUN B", "DELTA", "CHANGED"], [["  " + r[0], *r[1:]] for r in rows]):
                self.emit(line)
        else:
            self.emit("  neither run recorded a summary")
        for name, record in (("run_a", first), ("run_b", second)):
            if not kept[name]:
                self.emit(f"Detail: not kept for run {record.run_id} ({NOT_KEPT_NOTE}); only summaries compared")
        if detail is not None:
            self._emit_detail_diff(detail)

    def _emit_detail_diff(self, detail: Dict[str, Any]) -> None:
        self.emit("Detail:")
        self.emit(f"  only in run A: {', '.join(detail['removed']) or 'none'}")
        self.emit(f"  only in run B: {', '.join(detail['added']) or 'none'}")
        self.emit(f"  identical in both: {detail['unchanged']}")
        if not detail["changed"]:
            self.emit("  changed: none")
            return
        rows = [
            ["  " + key, field, _text(values[0]), _text(values[1])]
            for key, fields in detail["changed"].items()
            for field, values in fields.items()
        ]
        for line in _table(["  CHANGED", "FIELD", "RUN A", "RUN B"], rows):
            self.emit(line)

    # --- accession -----------------------------------------------------------------------------

    def _accession(self, args: argparse.Namespace, project: Path, records: List[RunRecord]) -> None:
        accession = args.accession.strip()
        if not accession:
            raise ValidationError("--accession needs an accession")
        history = accession_history(project, records, accession)
        listed = {entry["run_id"] for entry in history}
        unsearched = sum(1 for r in records if r.run_id not in listed and not detail_kept(project, r))
        if args.json:
            self.emit_json({"accession": accession, "runs": history, "runs_without_detail": unsearched})
            return
        if not history:
            self.emit(f"No recorded run with a kept detail holds values for {accession}")
        for entry in history:
            self.emit(f"{entry['run_id']}  {entry['command']}  started {entry['started']}  exit {entry['exit_code']}")
            if entry["values"] is None:
                self.emit(f"  detail not kept ({NOT_KEPT_NOTE})")
                continue
            for key, row in entry["values"].items():
                fields = " ".join(f"{field}={_text(value)}" for field, value in row.items())
                self.emit(f"  {key}  {fields}")
        if unsearched:
            self.emit(f"{unsearched} other run(s) have no kept detail and could not be searched")
