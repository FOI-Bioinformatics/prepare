"""``prepare doctor``: environment checks with fake tools on PATH, faked disk space and no network.

HOME, XDG_CONFIG_HOME and PREPARE_DATA are isolated as in the store tests, ``PATH`` holds only
``tests/helpers_tools.fake_tool`` scripts, and ``shutil.disk_usage`` and ``requests.get`` are
patched, so no real tool, store, config file or network is touched.
"""

import json
from collections import namedtuple
from unittest.mock import Mock, patch

import pytest
import requests

from helpers_tools import fake_tool
from prepare.cli.main import main
from prepare.core import settings
from prepare.processing import doctor_report
from prepare.processing.doctor_report import FAIL, OK, WARN, Check, overall_status, tools_needed_for
from prepare.store.layout import init_store

Usage = namedtuple("Usage", "total used free")
GB = 1024**3
PLENTY = Usage(1000 * GB, 100 * GB, 900 * GB)

VERSIONS = {
    "fasterq-dump": "fasterq-dump : 3.1.1",
    "prefetch": "prefetch : 3.1.1",
    "minimap2": "2.28-r1209",
    "samtools": "samtools 1.21",
    "megahit": "MEGAHIT v1.2.9",
    "pigz": "pigz 2.8",
    "datasets": "datasets version: 16.27.0",
    "seqkit": "seqkit v2.8.2",
}


@pytest.fixture(autouse=True)
def _isolated(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    monkeypatch.delenv("PREPARE_DATA", raising=False)
    # conftest sets this to 0 for every test; doctor's free-space checks need the default of 10.
    monkeypatch.delenv("PREPARE_MIN_FREE_GB", raising=False)
    for name in ("SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "PREPARE_TEMP_FOLDER"):
        monkeypatch.delenv(name, raising=False)
    project = tmp_path / "project"
    project.mkdir()
    monkeypatch.chdir(project)
    settings.reset_for_tests()
    with patch("prepare.processing.doctor_report.shutil.disk_usage", return_value=PLENTY):
        yield
    settings.reset_for_tests()


def _tools_on_path(tmp_path, monkeypatch, leave_out=()):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    for name, text in VERSIONS.items():
        if name not in leave_out:
            if name == "pigz":
                fake_tool(tmp_path, name, "", stderr=text)
            else:
                fake_tool(tmp_path, name, text)
    monkeypatch.setenv("PATH", str(bin_dir))


def _run_json(capsys, *extra):
    rc = main(["doctor", "--json", *extra])
    out = capsys.readouterr().out
    return rc, json.loads(out)


def _by_name(report):
    return {check["name"]: check for check in report["checks"]}


def test_every_tool_present_gives_exit_0_and_one_json_document(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    with patch("prepare.processing.doctor_report.requests.get") as get:
        rc, report = _run_json(capsys)
    assert rc == 0
    assert report["status"] == OK
    checks = _by_name(report)
    for name, text in VERSIONS.items():
        tool = checks[f"tool {name}"]
        assert tool["status"] == OK, tool
        assert tool["data"]["version"], tool
        assert tool["data"]["path"].endswith(name)
    assert checks["tool minimap2"]["data"]["min_version"] == "2.17"
    assert checks["python"]["data"]["prepare"]
    assert checks["config"]["status"] == OK
    assert checks["config"]["data"]["settings"]["min_free_gb"] == {"value": 10.0, "source": "default"}
    assert checks["store"]["status"] == OK and checks["store"]["data"]["root"] is None
    assert checks["registry"]["status"] == OK and checks["registry"]["data"]["exists"] is False
    assert checks["resources"]["data"]["cpus_available"] >= 1
    # No network check without --network, and nothing was requested.
    assert not [name for name in checks if name.startswith("network")]
    get.assert_not_called()


def test_a_missing_tool_is_a_warning_without_for(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch, leave_out=("minimap2",))
    rc, report = _run_json(capsys)
    assert rc == 0
    assert report["status"] == WARN
    minimap2 = _by_name(report)["tool minimap2"]
    assert minimap2["status"] == WARN
    assert "not found on PATH" in minimap2["detail"] and "conda install" in minimap2["detail"]


def test_for_a_command_that_needs_a_missing_tool_fails_with_exit_3(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch, leave_out=("minimap2", "megahit"))
    rc, report = _run_json(capsys, "--for", "extract_target_reads")
    assert rc == 3
    assert report["status"] == FAIL
    checks = _by_name(report)
    assert checks["tool minimap2"]["status"] == FAIL
    # megahit is used only with --assemble, so it stays a warning.
    assert checks["tool megahit"]["status"] == WARN
    # A tool the command does not use stays a warning when missing, and fine when present.
    assert checks["tool datasets"]["status"] == OK


def test_a_tool_below_its_floor_fails_even_without_for(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    fake_tool(tmp_path, "minimap2", "2.16-r922")
    rc, report = _run_json(capsys)
    assert rc == 3
    minimap2 = _by_name(report)["tool minimap2"]
    assert minimap2["status"] == FAIL and "2.17" in minimap2["detail"]


def test_for_an_unknown_command_is_a_usage_error(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    assert main(["doctor", "--for", "no_such_command"]) == 2


def test_for_rejects_a_hidden_former_command_name():
    assert main(["doctor", "--for", "sra_stats"]) == 2


def test_low_disk_space_warns(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    with patch("prepare.processing.doctor_report.shutil.disk_usage", return_value=Usage(100 * GB, 98 * GB, 2 * GB)):
        rc, report = _run_json(capsys)
    assert rc == 0
    checks = _by_name(report)
    for name in ("free space project", "free space temp", "free space sra-cache"):
        assert checks[name]["status"] == WARN, checks[name]
        assert "below" in checks[name]["detail"] and "10" in checks[name]["detail"]
    assert checks["free space project"]["data"]["free_bytes"] == 2 * GB


def test_min_free_gb_zero_turns_the_free_space_warning_off(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    monkeypatch.setenv("PREPARE_MIN_FREE_GB", "0")
    with patch("prepare.processing.doctor_report.shutil.disk_usage", return_value=Usage(100 * GB, 98 * GB, 2 * GB)):
        rc, report = _run_json(capsys)
    assert rc == 0
    assert _by_name(report)["free space project"]["status"] == OK


def test_malformed_config_is_a_failed_check_and_exit_3_through_main(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    config = tmp_path / "xdg" / "prepare" / "config.toml"
    config.parent.mkdir(parents=True)
    config.write_text("[runtime\nnot toml")
    rc, report = _run_json(capsys)
    assert rc == 3
    checks = _by_name(report)
    assert checks["config"]["status"] == FAIL
    assert "not valid TOML" in checks["config"]["detail"]
    # The rest of the report still ran, and the one problem is one failed check, not also the store's.
    assert checks["tool minimap2"]["status"] == OK
    assert checks["store"]["status"] == WARN and "config file does not parse" in checks["store"]["detail"]
    assert [check["name"] for check in report["checks"] if check["status"] == FAIL] == ["config"]


def test_malformed_config_still_stops_any_other_command(tmp_path, monkeypatch):
    config = tmp_path / "xdg" / "prepare" / "config.toml"
    config.parent.mkdir(parents=True)
    config.write_text("[runtime\nnot toml")
    assert main(["store_status"]) == 3


def test_a_bad_environment_value_is_reported_by_the_config_check(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    monkeypatch.setenv("PREPARE_TIMEOUT", "soon")
    rc, report = _run_json(capsys)
    assert rc == 3
    config = _by_name(report)["config"]
    assert config["status"] == FAIL and "PREPARE_TIMEOUT" in config["detail"]


def test_unknown_config_key_warns(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    config = tmp_path / "xdg" / "prepare" / "config.toml"
    config.parent.mkdir(parents=True)
    config.write_text("[runtime]\nmin_free_gb = 5\nno_such_key = 1\n")
    rc, report = _run_json(capsys)
    assert rc == 0
    checks = _by_name(report)
    assert checks["config"]["status"] == WARN and "no_such_key" in checks["config"]["detail"]
    assert checks["config"]["data"]["settings"]["min_free_gb"]["value"] == 5.0


def test_store_with_marker_is_checked_for_writing_and_space(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    root = tmp_path / "store"
    init_store(root)
    rc, report = _run_json(capsys, "--data-root", str(root))
    assert rc == 0
    store = _by_name(report)["store"]
    assert store["status"] == OK, store
    assert store["data"]["root"] == str(root.resolve())
    assert store["data"]["writable"] is True
    assert store["data"]["free_bytes"] == PLENTY.free
    assert list(root.iterdir())  # the probe file was removed; only the store's own entries remain
    assert not [p for p in root.iterdir() if p.name.startswith(".doctor")]


def test_store_without_marker_fails(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    (tmp_path / "not-a-store").mkdir()
    rc, report = _run_json(capsys, "--data-root", str(tmp_path / "not-a-store"))
    assert rc == 3
    store = _by_name(report)["store"]
    assert store["status"] == FAIL and "marker" in store["detail"]


def test_a_registry_that_does_not_parse_fails(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    (tmp_path / "project" / "prepare_registry.json").write_text("{not json")
    rc, report = _run_json(capsys)
    assert rc == 3
    registry = _by_name(report)["registry"]
    assert registry["status"] == FAIL and "not valid JSON" in registry["detail"]


def test_the_nearest_registry_above_the_project_folder_is_found(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    (tmp_path / "project" / "prepare_registry.json").write_text('{"version": 2, "datasets": {"SRR1": {}}}')
    sub = tmp_path / "project" / "fastq"
    sub.mkdir()
    rc, report = _run_json(capsys, "--project", str(sub))
    assert rc == 0
    registry = _by_name(report)["registry"]
    assert registry["data"]["exists"] is True and registry["data"]["datasets"] == 1


def test_slurm_variables_are_reported(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    monkeypatch.setenv("SLURM_JOB_ID", "4242")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    rc, report = _run_json(capsys)
    resources = _by_name(report)["resources"]
    assert resources["data"]["slurm"] == {"SLURM_JOB_ID": "4242", "SLURM_CPUS_PER_TASK": "8"}
    assert "SLURM job 4242" in resources["detail"]


def test_network_checks_run_only_with_network(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    with patch("prepare.processing.doctor_report.requests.get", return_value=Mock(status_code=200)) as get:
        rc, report = _run_json(capsys, "--network")
    assert rc == 0
    checks = _by_name(report)
    assert checks["network ncbi"]["status"] == OK
    assert checks["network branchwater"]["status"] == OK
    assert all(call.kwargs["timeout"] == 10 for call in get.call_args_list)
    assert len(get.call_args_list) == 2


def test_an_unreachable_service_fails_the_network_check(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    with patch(
        "prepare.processing.doctor_report.requests.get", side_effect=requests.ConnectionError("no route to host")
    ):
        rc, report = _run_json(capsys, "--network")
    assert rc == 3
    assert _by_name(report)["network ncbi"]["status"] == FAIL
    assert "no route to host" in _by_name(report)["network ncbi"]["detail"]


def test_text_output_lists_each_check_and_a_summary(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch, leave_out=("pigz",))
    rc = main(["doctor"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "[ok]" in out and "[warn]" in out
    assert "tool minimap2" in out and "2.28" in out
    assert out.strip().splitlines()[-1].startswith("Result:")


def test_tools_needed_for_lists_only_required_tools():
    assert tools_needed_for("extract_target_reads") == {"minimap2", "samtools"}
    assert tools_needed_for("download_sra") == {"fasterq-dump"}
    assert tools_needed_for("genome_download") == {"datasets"}
    assert tools_needed_for("status") == set()


def test_overall_status_is_the_worst_one():
    assert overall_status([Check("a", OK, "")]) == OK
    assert overall_status([Check("a", OK, ""), Check("b", WARN, "")]) == WARN
    assert overall_status([Check("a", WARN, ""), Check("b", FAIL, "")]) == FAIL
    assert overall_status([]) == OK


def test_check_serialises_to_a_plain_dict():
    assert Check("x", OK, "fine", {"n": 1}).to_dict() == {"name": "x", "status": OK, "detail": "fine", "data": {"n": 1}}


def test_run_checks_is_usable_without_the_cli(tmp_path, monkeypatch):
    _tools_on_path(tmp_path, monkeypatch)
    checks = doctor_report.run_checks(project=tmp_path / "project")
    assert {check.name for check in checks} >= {"python", "config", "store", "registry", "resources"}


def test_for_fails_a_required_tool_that_cannot_run(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    loader = "samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file"
    fake_tool(tmp_path, "samtools", "", stderr=loader, rc=1)
    rc, report = _run_json(capsys, "--for", "extract_target_reads")
    assert rc == 3
    samtools = _by_name(report)["tool samtools"]
    assert samtools["status"] == FAIL and "not runnable" in samtools["detail"]


def test_a_tool_that_cannot_run_is_a_warning_when_not_required(tmp_path, monkeypatch, capsys):
    _tools_on_path(tmp_path, monkeypatch)
    fake_tool(tmp_path, "samtools", "", stderr="broken", rc=1)
    rc, report = _run_json(capsys)
    assert _by_name(report)["tool samtools"]["status"] == WARN


def test_the_temp_location_is_where_a_download_puts_its_scratch(tmp_path, monkeypatch):
    project = tmp_path / "project"
    locations = dict(doctor_report.space_locations(project, None))
    assert locations["temp"] == project / "fastq" / ".prepare-tmp"
    assert locations["sra-cache"] == project / "fastq" / ".sra-cache"

    store_root = tmp_path / "store"
    with_store = dict(doctor_report.space_locations(project, str(store_root)))
    assert with_store["temp"] == store_root / "tmp"

    monkeypatch.setenv("PREPARE_TEMP_FOLDER", str(tmp_path / "local"))
    settings.reset_for_tests()
    assert dict(doctor_report.space_locations(project, None))["temp"] == tmp_path / "local"
