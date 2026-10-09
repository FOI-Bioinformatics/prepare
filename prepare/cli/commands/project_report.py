"""``prepare project_report``: one report of the whole project, from its registry and run log.

Writes ``project_report.md`` and ``project_report.json`` to ``--output-dir`` every time, and
``project_report.html`` when the ``interactive`` extra (plotly, jinja2) is installed. The report
itself is built by ``prepare.processing.project_report``; this command reads no FASTQ file. It
exits 1 without a project registry, and 3 with ``--html always`` when the extra is missing; in that
case nothing is written.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from prepare.cli.base import BaseCommand, emit_error_json
from prepare.core.exceptions import ConfigurationError, DataAccessError, PrepareError
from prepare.data import run_log
from prepare.data.file_io import write_text_atomic
from prepare.data.registry import load_registry, record_export, registry_path, registry_transaction
from prepare.processing.project_report import DEFAULT_MAX_ROWS, SECTIONS, build_project_report
from prepare.processing.project_report_markdown import render_markdown

EXPORT_NAME = "project_report"
MARKDOWN_FILE = "project_report.md"
JSON_FILE = "project_report.json"
HTML_FILE = "project_report.html"
HTML_CHOICES = ("auto", "always", "never")


def _max_rows(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a whole number, got {text!r}")
    if value < 0:
        raise argparse.ArgumentTypeError(f"expected a whole number, 0 or more, got {value}")
    return value


class ProjectReportCommand(BaseCommand):
    """Write the project report (Markdown, JSON and, with the interactive extra, HTML)."""

    @property
    def name(self) -> str:
        """The command name."""
        return "project_report"

    @property
    def help(self) -> str:
        """One line for the command list."""
        return "Write one report of the project from its registry: Markdown, JSON and HTML"

    @property
    def group(self) -> str:
        """Listed under Environment in the main help, next to ``doctor`` and ``runs``."""
        return "Environment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the output folder, HTML choice, row limit, environment and recording flags."""
        parser.add_argument("--output-dir", default="project_report", help="Folder for the report files")
        parser.add_argument(
            "--html",
            choices=HTML_CHOICES,
            default="auto",
            help="auto: write HTML when plotly and jinja2 are installed; always: exit 3 without them; never: no HTML",
        )
        parser.add_argument(
            "--max-rows",
            type=_max_rows,
            default=DEFAULT_MAX_ROWS,
            metavar="N",
            help="Rows kept per table and list; 0 keeps all (results_table writes every extraction row)",
        )
        parser.add_argument(
            "--no-environment",
            dest="no_environment",
            action="store_true",
            help="Leave out the environment checks of doctor",
        )
        parser.add_argument(
            "--no-record",
            dest="no_record",
            action="store_true",
            help="Write the report but do not record the export in the registry or the run in the run log",
        )
        parser.add_argument(
            "--registry",
            default=None,
            help="Registry file (default: prepare_registry.json found upwards from here)",
        )
        parser.add_argument("--json", action="store_true", help="Write the result as one JSON document")

    def execute(self, args: argparse.Namespace) -> int:
        """Build and write the report; 0 on success, 1 without a registry, 3 for a missing extra."""
        registry_file = registry_path(args.registry)
        if not registry_file.is_file():
            missing = DataAccessError(
                f"No project registry at {registry_file}; create one with 'prepare status --init' "
                "in the project folder, or name it with --registry"
            )
            return self._failed(args, missing)
        try:
            html = self._html_available(args.html)
            registry = load_registry(registry_file)
            report = build_project_report(registry, max_rows=args.max_rows, include_environment=not args.no_environment)
            texts = {MARKDOWN_FILE: render_markdown(report), JSON_FILE: json.dumps(report, indent=2) + "\n"}
            if html:
                from prepare.visualization.project_report import render_html

                texts[HTML_FILE] = render_html(report)
            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            for name, text in texts.items():
                write_text_atomic(output_dir / name, text)
            recorded = self._record(args, registry_file, output_dir, list(texts), report)
        except (PrepareError, OSError) as e:
            return self._failed(args, e)
        self._report_written(args, output_dir, list(texts), recorded)
        return 0

    def records_run(self, args: argparse.Namespace) -> bool:
        """Every run is added to the project's run log, unless ``--no-record``."""
        return not getattr(args, "no_record", False)

    def _html_available(self, choice: str) -> bool:
        """Whether to write HTML: ``never`` no; ``always`` yes or a ConfigurationError; ``auto`` when installed."""
        if choice == "never":
            return False
        from prepare.visualization.project_report import require_html

        try:
            require_html()
        except ConfigurationError as e:
            if choice == "always":
                raise
            self.logger.info("HTML report not written: %s", e)
            return False
        return True

    def _record(
        self, args: argparse.Namespace, registry_file: Path, output_dir: Path, files: List[str], report: Dict[str, Any]
    ) -> bool:
        """Record the export in the registry unless ``--no-record``; True when it was recorded."""
        if args.no_record:
            self.logger.info("Not recording this export in the registry (--no-record)")
            return False
        summary = {
            "files": files,
            "html": HTML_FILE in files,
            "datasets": report["project"]["datasets"],
            "extraction_rows": report["extractions"]["rows_total"],
            "failed_downloads": report["failures"]["rows_total"],
            "environment": report["environment"].get("status"),
            "max_rows": args.max_rows,
        }
        run_log.note_run(args, summary=summary)
        with registry_transaction(registry_file) as locked:
            record_export(locked, EXPORT_NAME, output_dir / MARKDOWN_FILE, summary)
        return True

    def _report_written(self, args: argparse.Namespace, output_dir: Path, files: List[str], recorded: bool) -> None:
        paths = {name: str(output_dir / name) for name in files}
        if args.json:
            html: Optional[str] = paths.get(HTML_FILE)
            self.emit_json(
                {
                    "output_dir": str(output_dir),
                    "files": {"markdown": paths[MARKDOWN_FILE], "json": paths[JSON_FILE], "html": html},
                    "recorded": recorded,
                    "sections": list(SECTIONS),
                }
            )
            return
        for path in paths.values():
            self.emit(f"Wrote {path}")

    def _failed(self, args: argparse.Namespace, error: BaseException) -> int:
        """Report ``error`` through ``fail``, and as a JSON error document with ``--json``."""
        if args.json:
            emit_error_json(str(error))
        return self.fail(error, "Error writing the project report")
