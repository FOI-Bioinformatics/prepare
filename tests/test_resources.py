"""Tests for CPU and memory detection (prepare.utils.resources) and the megahit memory flag."""

import argparse
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.cli.commands import extraction_assembly as assembly_cmd
from prepare.cli.commands.read_extraction import ExtractTargetReadsCommand
from prepare.core.exceptions import ConfigurationError
from prepare.data import assembly as assembly_mod
from prepare.utils import resources
from helpers_extraction import _fake_tools

GIB = 1024**3


@pytest.fixture(autouse=True)
def no_slurm(monkeypatch):
    """Start every test without SLURM job variables, whatever the host sets."""
    for name in ("SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU"):
        monkeypatch.delenv(name, raising=False)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _cgroup_v2(tmp_path: Path, group: str, limits: dict) -> tuple:
    """A fake /proc and /sys tree for cgroup v2: ``limits`` maps a group path to its memory.max text."""
    proc = tmp_path / "proc"
    sys_root = tmp_path / "sys"
    _write(proc / "self" / "cgroup", f"0::{group}\n")
    (sys_root / "fs" / "cgroup").mkdir(parents=True, exist_ok=True)
    for path, text in limits.items():
        _write(sys_root / "fs" / "cgroup" / path.lstrip("/") / "memory.max", text + "\n")
    return proc, sys_root


class TestAvailableCpus:
    def test_affinity_mask_wins(self, monkeypatch):
        monkeypatch.setattr(resources.os, "sched_getaffinity", lambda pid: {0, 1, 2}, raising=False)
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        assert resources.available_cpus() == 3

    def test_slurm_when_no_affinity(self, monkeypatch):
        monkeypatch.delattr(resources.os, "sched_getaffinity", raising=False)
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "6")
        monkeypatch.setattr(resources.os, "cpu_count", lambda: 64)
        assert resources.available_cpus() == 6

    def test_cpu_count_otherwise(self, monkeypatch):
        monkeypatch.delattr(resources.os, "sched_getaffinity", raising=False)
        monkeypatch.setattr(resources.os, "cpu_count", lambda: 12)
        assert resources.available_cpus() == 12

    def test_one_when_nothing_is_known(self, monkeypatch):
        monkeypatch.delattr(resources.os, "sched_getaffinity", raising=False)
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "not a number")
        monkeypatch.setattr(resources.os, "cpu_count", lambda: None)
        assert resources.available_cpus() == 1

    def test_unreadable_affinity_falls_back(self, monkeypatch):
        def broken(pid):
            raise OSError("not permitted")

        monkeypatch.setattr(resources.os, "sched_getaffinity", broken, raising=False)
        monkeypatch.setattr(resources.os, "cpu_count", lambda: 5)
        assert resources.available_cpus() == 5


class TestMemoryLimit:
    def test_cgroup_v2_takes_the_smallest_limit_on_the_way_to_the_root(self, tmp_path):
        proc, sys_root = _cgroup_v2(
            tmp_path,
            "/slurm/job_1/step_0",
            {"/slurm": str(64 * GIB), "/slurm/job_1": str(8 * GIB), "/slurm/job_1/step_0": "max"},
        )
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=sys_root) == 8 * GIB

    def test_cgroup_v2_all_max_falls_through_to_none(self, tmp_path):
        proc, sys_root = _cgroup_v2(tmp_path, "/user.slice", {"/user.slice": "max"})
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=sys_root) is None

    def test_cgroup_v2_root_group(self, tmp_path):
        proc, sys_root = _cgroup_v2(tmp_path, "/", {"/": str(4 * GIB)})
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=sys_root) == 4 * GIB

    def test_cgroup_v1(self, tmp_path):
        proc = tmp_path / "proc"
        sys_root = tmp_path / "sys"
        _write(proc / "self" / "cgroup", "4:memory:/job\n")
        _write(sys_root / "fs" / "cgroup" / "memory" / "memory.limit_in_bytes", f"{2 * GIB}\n")
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=sys_root) == 2 * GIB

    def test_cgroup_v1_unlimited_is_ignored(self, tmp_path):
        sys_root = tmp_path / "sys"
        _write(sys_root / "fs" / "cgroup" / "memory" / "memory.limit_in_bytes", "9223372036854771712\n")
        assert resources.memory_limit_bytes(proc_root=tmp_path / "proc", sys_root=sys_root) is None

    def test_cgroup_v1_slurm_job_group_below_an_unlimited_root(self, tmp_path):
        # A SLURM job on a cgroup v1 host without a cgroup namespace: the root says unlimited,
        # the job's own group holds the limit.
        proc = tmp_path / "proc"
        memory = tmp_path / "sys" / "fs" / "cgroup" / "memory"
        _write(proc / "self" / "cgroup", "11:cpuset:/slurm/uid_1/job_7\n4:memory:/slurm/uid_1/job_7/step_0\n")
        _write(memory / "memory.limit_in_bytes", "9223372036854771712\n")
        _write(memory / "slurm" / "uid_1" / "job_7" / "memory.limit_in_bytes", f"{16 * GIB}\n")
        _write(memory / "slurm" / "uid_1" / "job_7" / "step_0" / "memory.limit_in_bytes", "9223372036854771712\n")
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=tmp_path / "sys") == 16 * GIB

    def test_cgroup_v1_combined_controller_entry(self, tmp_path):
        proc = tmp_path / "proc"
        memory = tmp_path / "sys" / "fs" / "cgroup" / "memory"
        _write(proc / "self" / "cgroup", "0::/ignored\n6:cpuacct,memory:/job\n")
        _write(memory / "job" / "memory.limit_in_bytes", f"{3 * GIB}\n")
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=tmp_path / "sys") == 3 * GIB

    def test_slurm_mem_per_cpu_times_available_cpus(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SLURM_MEM_PER_CPU", "4000")
        monkeypatch.setattr(resources, "available_cpus", lambda: 4)
        assert resources.memory_limit_bytes(proc_root=tmp_path / "p", sys_root=tmp_path / "s") == 16000 * 1024**2

    def test_slurm_mem_per_node_wins_over_per_cpu(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "1000")
        monkeypatch.setenv("SLURM_MEM_PER_CPU", "4000")
        assert resources.memory_limit_bytes(proc_root=tmp_path / "p", sys_root=tmp_path / "s") == 1000 * 1024**2

    def test_slurm_fallback_in_megabytes(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "16000")
        assert resources.memory_limit_bytes(proc_root=tmp_path / "p", sys_root=tmp_path / "s") == 16000 * 1024**2

    def test_cgroup_wins_over_slurm(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "16000")
        proc, sys_root = _cgroup_v2(tmp_path, "/job", {"/job": str(GIB)})
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=sys_root) == GIB

    def test_none_without_cgroups_or_slurm(self, tmp_path):
        # macOS, or a Linux host with no limit: nothing to read.
        assert resources.memory_limit_bytes(proc_root=tmp_path / "p", sys_root=tmp_path / "s") is None

    def test_garbage_is_ignored(self, tmp_path, monkeypatch):
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "lots")
        proc, sys_root = _cgroup_v2(tmp_path, "/job", {"/job": "unlimited"})
        assert resources.memory_limit_bytes(proc_root=proc, sys_root=sys_root) is None


class TestParseMemory:
    @pytest.mark.parametrize(
        "value, limit, expected",
        [
            ("auto", 10 * GIB, int(0.9 * 10 * GIB)),
            ("AUTO", 1000, 900),
            ("auto", None, None),
            ("0.5", None, 0.5),
            ("0.5", 10 * GIB, 0.5),
            (".25", None, 0.25),
            ("1.0", None, 1.0),
            ("1", None, 1.0),
            ("32G", None, 32 * GIB),
            ("32g", None, 32 * GIB),
            ("32000M", None, 32000 * 1024**2),
            ("512K", None, 512 * 1024),
            ("1T", None, 1024**4),
            ("1.5G", None, int(1.5 * GIB)),
            ("123456789", None, 123456789),
        ],
    )
    def test_table(self, value, limit, expected):
        result = resources.parse_memory(value, limit)
        assert result == expected
        assert type(result) is type(expected)

    @pytest.mark.parametrize("value", ["", "lots", "-1G", "2.5", "1.5x", "2", "64000"])
    def test_rejects(self, value):
        with pytest.raises(ValueError):
            resources.parse_memory(value, None)


class TestMegahitMemory:
    def _run(self, tmp_path, memory):
        reads = [tmp_path / "r_1.fq", tmp_path / "r_2.fq"]
        with patch.object(assembly_mod.SecureSubprocess, "run_secure", side_effect=_fake_tools({})) as run:
            assembly_mod.assemble_extracted_reads(reads, tmp_path / "asm", preset=None, memory=memory)
        return run.call_args.args[1]

    def test_memory_flag_when_resolved(self, tmp_path):
        argv = self._run(tmp_path, 8 * GIB)
        assert argv[argv.index("--memory") + 1] == str(8 * GIB)

    def test_fraction_passes_through(self, tmp_path):
        argv = self._run(tmp_path, 0.5)
        assert argv[argv.index("--memory") + 1] == "0.5"

    def test_no_memory_flag_when_unresolved(self, tmp_path):
        assert "--memory" not in self._run(tmp_path, None)

    def test_assembly_memory_resolves_auto_from_the_detected_limit(self, monkeypatch):
        monkeypatch.setattr(resources, "memory_limit_bytes", lambda: 10 * GIB)
        assert assembly_mod.assembly_memory("auto") == int(0.9 * 10 * GIB)

    def test_assembly_memory_auto_without_a_limit_is_none(self, monkeypatch):
        monkeypatch.setattr(resources, "memory_limit_bytes", lambda: None)
        assert assembly_mod.assembly_memory("auto") is None


class TestAssemblyMemoryFlag:
    @staticmethod
    def _memory(args):
        return assembly_cmd.resolve_assembly_memory(args, ExtractTargetReadsCommand().logger)

    def _parse(self, *argv):
        parser = argparse.ArgumentParser()
        ExtractTargetReadsCommand().configure_parser(parser)
        return parser.parse_args(["--parsed-containment", "p.tsv", "--genome-id", "g", "--genome-fasta", "g.fa", *argv])

    def test_flag_value_in_bytes(self):
        assert self._memory(self._parse("--assembly-memory", "16G")) == 16 * GIB

    def test_setting_decides_without_the_flag(self, monkeypatch):
        monkeypatch.setenv("PREPARE_ASSEMBLY_MEMORY", "0.5")
        assert self._memory(self._parse()) == 0.5

    def test_auto_without_a_limit_leaves_the_flag_out(self, monkeypatch):
        monkeypatch.setattr(resources, "memory_limit_bytes", lambda: None)
        assert self._memory(self._parse()) is None

    def test_bad_value_is_a_configuration_error(self):
        with pytest.raises(ConfigurationError, match="--assembly-memory"):
            self._memory(self._parse("--assembly-memory", "2.5"))

    def test_resolved_memory_reaches_megahit(self, tmp_path, monkeypatch):
        """The command passes the resolved value to assemble_extracted_reads."""
        monkeypatch.setattr(resources, "memory_limit_bytes", lambda: 10 * GIB)
        args = self._parse("--output-folder", str(tmp_path / "out"), "--registry", str(tmp_path / "reg.json"))
        genome = tmp_path / "g.fa"
        genome.write_text(">g\nACGT\n")
        args.genome_fasta = str(genome)
        with patch.object(assembly_cmd, "assemble_extracted_reads", side_effect=RuntimeError("stop here")) as asm:
            with patch.object(assembly_cmd, "megahit_version", return_value="1.2.9"):
                with pytest.raises(RuntimeError):
                    assembly_cmd.assemble_samples(ExtractTargetReadsCommand(), args, {"SRR1": [tmp_path / "r.fq"]}, {})
        assert asm.call_args.kwargs["memory"] == int(0.9 * 10 * GIB)
