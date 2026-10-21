"""Tests for the log file, the quiet and verbose flags, host and PID in the file, and tracebacks.

``setup_logging`` installs a console handler (stderr) and, with a log file, a file handler that
appends. The file always receives INFO and the traceback of a failure; the console shows the
traceback only at DEBUG. The logging flags parse both before and after the command name.
"""

import logging
import os
import socket
import sys
from unittest.mock import patch

import pytest

from prepare.cli.main import create_parser, main
from prepare.utils.logging import ConsoleFormatter, setup_logging


def _tagged_handlers():
    return [h for h in logging.getLogger().handlers if getattr(h, "_prepare", False)]


@pytest.fixture(autouse=True)
def restore_root_logger():
    """Remove (and close) the handlers setup_logging installed and restore the root level."""
    root = logging.getLogger()
    level = root.level
    yield
    for handler in _tagged_handlers():
        root.removeHandler(handler)
        handler.close()
    root.setLevel(level)


def _raise_and_log(logger: logging.Logger) -> None:
    try:
        raise RuntimeError("the inner failure")
    except RuntimeError as e:
        logger.error("Step failed: %s", e, exc_info=True)


# --- setup_logging --------------------------------------------------------------


def test_setup_twice_leaves_one_console_handler():
    setup_logging(level=logging.INFO)
    setup_logging(level=logging.INFO)
    assert len(_tagged_handlers()) == 1


def test_setup_twice_with_a_file_leaves_one_handler_of_each_kind(tmp_path):
    setup_logging(level=logging.INFO, log_file=str(tmp_path / "run.log"))
    setup_logging(level=logging.INFO, log_file=str(tmp_path / "run.log"))
    kinds = sorted(type(h).__name__ for h in _tagged_handlers())
    assert kinds == ["FileHandler", "StreamHandler"]


def test_a_foreign_handler_survives():
    root = logging.getLogger()
    foreign = logging.NullHandler()
    root.addHandler(foreign)
    try:
        setup_logging(level=logging.INFO, log_file=None)
        assert foreign in root.handlers
    finally:
        root.removeHandler(foreign)


def test_the_file_is_appended_across_two_runs_and_its_folder_created(tmp_path):
    path = tmp_path / "logs" / "nested" / "run.log"
    setup_logging(level=logging.INFO, log_file=str(path))
    logging.getLogger("prepare.test").info("first run")
    setup_logging(level=logging.INFO, log_file=str(path))
    logging.getLogger("prepare.test").info("second run")
    lines = path.read_text(encoding="utf-8").splitlines()
    assert [line.endswith("first run") for line in lines] == [True, False]
    assert lines[1].endswith("second run")


def test_the_file_line_carries_host_and_pid(tmp_path):
    path = tmp_path / "run.log"
    setup_logging(level=logging.INFO, log_file=str(path))
    logging.getLogger("prepare.test").info("hello")
    line = path.read_text(encoding="utf-8").splitlines()[-1]
    assert f" {socket.gethostname()}[{os.getpid()}] INFO prepare.test: hello" in line


def test_the_console_shows_host_and_pid_only_when_asked(tmp_path, capsys):
    setup_logging(level=logging.INFO)
    logging.getLogger("prepare.test").info("plain")
    assert f"[{os.getpid()}]" not in capsys.readouterr().err
    setup_logging(level=logging.INFO, show_host=True)
    logging.getLogger("prepare.test").info("with host")
    assert f"{socket.gethostname()}[{os.getpid()}]" in capsys.readouterr().err


def test_a_quiet_console_still_leaves_info_in_the_file(tmp_path, capsys):
    path = tmp_path / "run.log"
    setup_logging(level=logging.WARNING, log_file=str(path))
    logging.getLogger("prepare.test").info("progress note")
    logging.getLogger("prepare.test").warning("a warning")
    err = capsys.readouterr().err
    assert "progress note" not in err and "a warning" in err
    text = path.read_text(encoding="utf-8")
    assert "progress note" in text and "a warning" in text


def test_the_traceback_goes_to_the_file_and_not_the_console(tmp_path, capsys):
    path = tmp_path / "run.log"
    setup_logging(level=logging.INFO, log_file=str(path))
    _raise_and_log(logging.getLogger("prepare.test"))
    err = capsys.readouterr().err
    assert "Step failed: the inner failure" in err
    assert "Traceback" not in err
    text = path.read_text(encoding="utf-8")
    assert "Traceback (most recent call last)" in text and "RuntimeError: the inner failure" in text


def test_the_console_shows_the_traceback_at_debug(capsys):
    setup_logging(level=logging.DEBUG, console_traceback=True)
    _raise_and_log(logging.getLogger("prepare.test"))
    assert "Traceback (most recent call last)" in capsys.readouterr().err


def test_the_console_formatter_leaves_the_record_as_it_found_it():
    try:
        raise ValueError("kept")
    except ValueError:
        record = logging.LogRecord("x", logging.ERROR, __file__, 1, "msg", None, exc_info=sys.exc_info())
    record.exc_text = "cached text from another handler"
    assert ConsoleFormatter("%(message)s").format(record) == "msg"
    assert record.exc_info is not None and record.exc_text == "cached text from another handler"


def test_a_log_file_that_cannot_be_opened_is_a_configuration_error(tmp_path):
    blocker = tmp_path / "a_file"
    blocker.write_text("x")
    assert main(["--log-file", str(blocker / "run.log"), "store_status", "--data-root", str(tmp_path)]) == 3


# --- the flags ------------------------------------------------------------------


@pytest.mark.parametrize(
    "argv",
    [
        ["--quiet", "--log-file", "run.log", "--progress-every", "7", "store_status"],
        ["store_status", "--quiet", "--log-file", "run.log", "--progress-every", "7"],
    ],
)
def test_both_flag_placements_parse(argv):
    args = create_parser().parse_args(argv)
    assert args.log_quiet is True and args.log_file == "run.log" and args.progress_every == 7


def test_a_subcommand_flag_left_out_keeps_the_main_parser_value():
    args = create_parser().parse_args(["--log-level", "ERROR", "store_status"])
    assert args.log_level == "ERROR" and args.log_quiet is False and args.log_file is None


def test_store_status_keeps_the_hidden_verbose_alias_and_takes_quiet():
    args = create_parser().parse_args(["-v", "store_status", "--verbose"])
    assert args.details_alias is True and args.log_verbose is True
    args = create_parser().parse_args(["store_status", "--quiet"])
    assert args.details_alias is False and args.log_verbose is False and args.log_quiet is True


@pytest.mark.parametrize(
    "argv", [["-q", "-v", "store_status"], ["store_status", "-q", "-v"], ["-q", "download_sra", "-v"]]
)
def test_quiet_and_verbose_together_are_a_usage_error(argv):
    with pytest.raises(SystemExit) as raised:
        main(argv)
    assert raised.value.code == 2


def test_a_negative_progress_interval_is_a_usage_error():
    with pytest.raises(SystemExit) as raised:
        create_parser().parse_args(["--progress-every", "-1", "store_status"])
    assert raised.value.code == 2


@pytest.mark.parametrize(
    "flags, level",
    [(["-q"], logging.WARNING), (["-v"], logging.DEBUG), ([], None)],
)
def test_quiet_and_verbose_override_the_log_level(tmp_path, flags, level):
    with patch("prepare.cli.main.setup_logging") as setup:
        main([*flags, "store_status", "--data-root", str(tmp_path / "none")])
    assert setup.call_args.kwargs["level"] == (level if level is not None else logging.INFO)


def test_main_passes_the_log_file_and_host_setting(tmp_path, monkeypatch):
    monkeypatch.setenv("PREPARE_LOG_HOST", "yes")
    path = tmp_path / "run.log"
    with patch("prepare.cli.main.setup_logging") as setup:
        main(["--log-file", str(path), "store_status", "--data-root", str(tmp_path / "none")])
    assert setup.call_args.kwargs["log_file"] == str(path)
    assert setup.call_args.kwargs["show_host"] is True
    assert setup.call_args.kwargs["console_traceback"] is False


def test_main_logs_version_argv_host_pid_and_slurm_ids_at_debug(tmp_path, monkeypatch, caplog):
    monkeypatch.setenv("SLURM_JOB_ID", "4242")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "7")
    with patch("prepare.cli.main.setup_logging"), caplog.at_level(logging.DEBUG):
        main(["-v", "store_status", "--data-root", str(tmp_path / "none")])
    text = caplog.text
    assert "PREPARE v" in text and "store_status" in text
    assert socket.gethostname() in text and str(os.getpid()) in text
    assert "SLURM_JOB_ID=4242" in text and "SLURM_ARRAY_TASK_ID=7" in text


def test_a_failing_command_writes_its_traceback_to_the_file_only(tmp_path, capsys):
    path = tmp_path / "run.log"
    code = main(["--quiet", "--log-file", str(path), "store_status", "--data-root", str(tmp_path / "none")])
    assert code != 0
    err = capsys.readouterr().err
    assert "Traceback" not in err
    assert "Use --log-level DEBUG" not in err
    text = path.read_text(encoding="utf-8")
    assert "Traceback (most recent call last)" in text


def test_the_logged_arguments_hide_the_api_key(caplog):
    from prepare.cli.main import _log_run_header, _masked_argv

    assert _masked_argv(["--api-key", "SECRET1", "--api-key=SECRET2", "x"]) == [
        "--api-key",
        "***",
        "--api-key=***",
        "x",
    ]
    with caplog.at_level(logging.DEBUG):
        _log_run_header(["download_metadata", "--api-key", "SECRET1", "--email", "a@b.c"])
    assert "SECRET1" not in caplog.text and "--api-key ***" in caplog.text


@pytest.mark.parametrize("flag", [["--api-k", "SECRETKEY1"], ["--api=SECRETKEY1"], ["--api-key", "SECRETKEY1"]])
def test_an_abbreviated_api_key_flag_is_masked_too(tmp_path, caplog, flag):
    from prepare.cli.main import main

    with patch("prepare.cli.main.setup_logging"), caplog.at_level(logging.DEBUG):
        main(["-v", "download_metadata", *flag, "--email", "a@b.c", "--matches-folder", str(tmp_path / "none")])
    assert "PREPARE v" in caplog.text
    assert "SECRETKEY1" not in caplog.text


def test_masking_without_a_namespace_covers_abbreviations():
    from prepare.cli.main import _masked_argv

    assert _masked_argv(["--api-k", "K1", "--api=K2", "--email", "a@b.c"]) == [
        "--api-k",
        "***",
        "--api=***",
        "--email",
        "a@b.c",
    ]


def test_store_status_accepts_the_global_v_flag(tmp_path):
    root = tmp_path / "store"
    with patch("prepare.cli.main.setup_logging") as setup:
        code = main(["store_status", "-v", "--data-root", str(root), "--registry", str(tmp_path / "r.json")])
    assert code in (0, 1)
    assert setup.call_args.kwargs["level"] == logging.DEBUG


def test_store_status_details_flag_and_hidden_verbose_alias_parse():
    parser = create_parser()
    assert parser.parse_args(["store_status", "--details"]).details is True
    alias = parser.parse_args(["store_status", "--verbose"])
    assert alias.details_alias is True
    assert alias.log_verbose is not True


@pytest.mark.parametrize("command", ["download_sra", "store_status", "find_by_taxonomy"])
def test_command_help_does_not_print_default_false(command, capsys):
    with pytest.raises(SystemExit):
        main([command, "--help"])
    assert "(default: False)" not in capsys.readouterr().out


def test_verify_downloads_help_keeps_one_default_note(capsys):
    with pytest.raises(SystemExit):
        main(["download_sra", "--help"])
    text = " ".join(capsys.readouterr().out.split())
    note = text.split("against NCBI spot counts", 1)[1].split(" --", 1)[0]
    assert note.count("(default") == 1 and "(default: on)" in note
    assert "(default: True)" not in text


def test_log_level_is_case_insensitive():
    assert create_parser().parse_args(["--log-level", "debug", "store_status"]).log_level == "DEBUG"


@pytest.mark.parametrize("flag", ["-v", "-q"])
@pytest.mark.parametrize("placement", ["before", "after"])
def test_log_level_with_quiet_or_verbose_is_a_usage_error(flag, placement):
    argv = ["--log-level", "ERROR", flag, "store_status"]
    if placement == "after":
        argv = ["store_status", "--log-level", "ERROR", flag]
    with pytest.raises(SystemExit) as raised:
        main(argv)
    assert raised.value.code == 2
