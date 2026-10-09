"""Assembly identity marker and staged megahit output (``prepare.data.assembly_identity``).

megahit is replaced by the fake from ``tests/helpers_extraction.py``; no real tool runs.
"""

import json
import logging
import os
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.core.exceptions import ProcessingError
from prepare.data import assembly_identity as ident
from prepare.data import registry as reg
from prepare.data.assembly import assemble_extracted_reads
from prepare.data.assembly_identity import (
    MARKER_NAME,
    AssemblyInputs,
    assembly_inputs,
    assembly_state,
    publish_assembly,
    read_marker,
    sweep_staging,
    write_marker,
)
from helpers_extraction import _fake_tools

RUN = "prepare.data.assembly.SecureSubprocess.run_secure"


def _reads(tmp_path, sizes=(10, 12)):
    paths = []
    for index, size in enumerate(sizes, start=1):
        path = tmp_path / f"GCF_1_{index}.fastq.gz"
        path.write_bytes(b"x" * size)
        paths.append(path)
    return paths


def _inputs(tmp_path, **changes):
    values = dict(extraction_date="2026-10-01T10:00:00", preset="meta-sensitive", min_contig_len=None)
    values.update(changes)
    return assembly_inputs(_reads(tmp_path), **values)


def _old_assembly(out_dir, text=">old len=4\nACGT\n"):
    out_dir.mkdir(parents=True)
    (out_dir / "final.contigs.fa").write_text(text)


def _staging_left(out_dir):
    return sorted(p.name for p in out_dir.parent.iterdir() if p.name.startswith(f".{out_dir.name}."))


def _megahit_calls(state):
    return [args for exe, args in state.get("calls", []) if exe == "megahit" and "-o" in args]


class TestAssemblyInputs:
    def test_reads_are_names_and_sizes(self, tmp_path):
        inputs = _inputs(tmp_path)
        assert inputs.reads == (("GCF_1_1.fastq.gz", 10), ("GCF_1_2.fastq.gz", 12))

    def test_mtime_is_not_part_of_the_identity(self, tmp_path):
        first = _inputs(tmp_path)
        for path in _reads(tmp_path):
            os.utime(path, (1, 1))
        assert _inputs(tmp_path) == first

    def test_default_preset_is_normalised_to_none(self, tmp_path):
        assert _inputs(tmp_path, preset="default") == _inputs(tmp_path, preset=None)

    def test_k_flags_are_sorted_pairs_without_dashes(self, tmp_path):
        inputs = _inputs(tmp_path, preset=None, k_flags={"k-max": 99, "--k-min": 21})
        assert inputs.k_flags == (("k-max", 99), ("k-min", 21))

    def test_round_trip(self, tmp_path):
        inputs = _inputs(tmp_path, min_contig_len=500, preset=None, k_flags={"k-min": 21})
        data = json.loads(json.dumps(inputs.to_dict()))
        assert AssemblyInputs.from_dict(data) == inputs

    @pytest.mark.parametrize("data", [None, [], {"reads": "x"}, {"reads": [["a"]]}, {"reads": [], "k_flags": [[1]]}])
    def test_from_dict_rejects_malformed(self, data):
        assert AssemblyInputs.from_dict(data) is None

    def test_differences_name_each_changed_field(self, tmp_path):
        base = _inputs(tmp_path)
        other = AssemblyInputs(
            reads=(("GCF_1_1.fastq.gz", 99),),
            extraction_date="2026-10-02",
            preset=None,
            min_contig_len=1000,
            k_flags=(("k-min", 21),),
        )
        changed = base.differences(other)
        assert [line.split(":")[0] for line in changed] == [
            "reads",
            "extraction_date",
            "preset",
            "min_contig_len",
            "k_flags",
        ]
        assert base.differences(base) == []


class TestMarker:
    def test_write_and_read(self, tmp_path):
        inputs = _inputs(tmp_path)
        out = tmp_path / "asm"
        out.mkdir()
        path = write_marker(out, inputs, "MEGAHIT v1.2.9", {"threads": 4})
        assert path == out / MARKER_NAME
        marker = read_marker(out)
        assert AssemblyInputs.from_dict(marker["inputs"]) == inputs
        assert marker["megahit_version"] == "MEGAHIT v1.2.9"
        assert marker["params"] == {"threads": 4}

    def test_missing_or_unreadable_marker_reads_as_none(self, tmp_path):
        assert read_marker(tmp_path) is None
        (tmp_path / MARKER_NAME).write_text("{not json")
        assert read_marker(tmp_path) is None


class TestAssemblyState:
    def test_absent(self, tmp_path):
        assert assembly_state(tmp_path / "asm", _inputs(tmp_path), True)[0] == "absent"

    def test_incomplete_without_contigs(self, tmp_path):
        (tmp_path / "asm").mkdir()
        assert assembly_state(tmp_path / "asm", _inputs(tmp_path), True)[0] == "incomplete"

    def test_marker_matches(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        write_marker(out, _inputs(tmp_path), "v1", {})
        assert assembly_state(out, _inputs(tmp_path), False)[0] == "current"

    def test_megahit_version_is_not_compared(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        write_marker(out, _inputs(tmp_path), "MEGAHIT v1.0", {})
        assert assembly_state(out, _inputs(tmp_path), False)[0] == "current"

    def test_marker_differs(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        write_marker(out, _inputs(tmp_path), "v1", {})
        state, reason = assembly_state(out, _inputs(tmp_path, min_contig_len=1000), True)
        assert state == "stale" and "min_contig_len" in reason

    def test_unmarked_folder_per_accept_unmarked(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        assert assembly_state(out, _inputs(tmp_path), True)[0] == "current"
        assert assembly_state(out, _inputs(tmp_path), False) == ("stale", "no marker")

    def test_unreadable_marker_is_stale(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        (out / MARKER_NAME).write_text("{not json")
        assert assembly_state(out, _inputs(tmp_path), True)[0] == "stale"

    def test_without_expected_contigs_are_current(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        assert assembly_state(out, None, False)[0] == "current"


class TestPublishAndSweep:
    def test_publish_replaces_existing_folder(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        staging = tmp_path / ".asm.new.tmp"
        _old_assembly(staging, ">new len=4\nACGT\n")
        publish_assembly(staging, out)
        assert ">new" in (out / "final.contigs.fa").read_text()
        assert not staging.exists()
        assert _staging_left(out) == []

    def test_publish_restores_old_folder_when_the_rename_fails(self, tmp_path):
        out = tmp_path / "asm"
        _old_assembly(out)
        staging = tmp_path / ".asm.new.tmp"
        _old_assembly(staging, ">new len=4\nACGT\n")
        real_rename = os.rename

        def failing(src, dst):
            if Path(src) == staging:
                raise OSError("rename failed")
            return real_rename(src, dst)

        with patch.object(ident.os, "rename", side_effect=failing):
            with pytest.raises(OSError):
                publish_assembly(staging, out)
        assert ">old" in (out / "final.contigs.fa").read_text()
        assert _staging_left(out) == [staging.name]

    def test_sweep_removes_staging_siblings_only(self, tmp_path):
        out = tmp_path / "GCF_1_assembly"
        _old_assembly(out)
        (tmp_path / ".GCF_1_assembly.host.1.ab.tmp").mkdir()
        (tmp_path / ".GCF_1.1_assembly.host.1.ab.tmp").mkdir()
        removed = sweep_staging(out)
        assert [p.name for p in removed] == [".GCF_1_assembly.host.1.ab.tmp"]
        assert (tmp_path / ".GCF_1.1_assembly.host.1.ab.tmp").exists()
        assert out.exists()

    def test_sweep_restores_an_aside_copy_when_the_folder_is_missing(self, tmp_path):
        out = tmp_path / "GCF_1_assembly"
        aside = ident._aside_path(out)
        _old_assembly(aside)
        (tmp_path / ".GCF_1_assembly.host.1.ab.tmp").mkdir()
        sweep_staging(out)
        assert (out / "final.contigs.fa").exists()
        assert _staging_left(out) == []

    def test_sweep_without_parent_folder(self, tmp_path):
        assert sweep_staging(tmp_path / "missing" / "asm") == []

    def test_sweep_does_not_bring_back_a_folder_removed_after_a_finished_publish(self, tmp_path):
        """An old copy with no staging folder beside it is a leftover of a finished publish."""
        out = tmp_path / "GCF_1_assembly"
        aside = ident._aside_path(out)
        _old_assembly(aside)

        removed = sweep_staging(out)

        assert not out.exists()
        assert removed == [aside]
        assert _staging_left(out) == []

    def test_sweep_publishes_a_finished_staging_folder_over_the_old_copy(self, tmp_path):
        """A kill between the two renames, after the new assembly and its marker were written."""
        out = tmp_path / "GCF_1_assembly"
        aside = ident._aside_path(out)
        _old_assembly(aside)
        staging = tmp_path / ".GCF_1_assembly.host.1.ab.tmp"
        _old_assembly(staging, ">new len=4\nACGT\n")
        write_marker(staging, _inputs(tmp_path), "v1", {})

        sweep_staging(out)

        assert ">new" in (out / "final.contigs.fa").read_text()
        assert (out / MARKER_NAME).is_file()
        assert _staging_left(out) == []

    def test_publish_keeps_the_first_error_when_the_restore_fails_too(self, tmp_path, caplog):
        out = tmp_path / "asm"
        _old_assembly(out)
        staging = tmp_path / ".asm.new.tmp"
        _old_assembly(staging, ">new len=4\nACGT\n")
        real_rename = os.rename

        def failing(src, dst):
            if Path(src) == staging:
                raise OSError("first rename failed")
            if Path(dst) == out:
                raise OSError("restore failed")
            return real_rename(src, dst)

        with patch.object(ident.os, "rename", side_effect=failing), caplog.at_level(logging.ERROR):
            with pytest.raises(OSError, match="first rename failed"):
                publish_assembly(staging, out)
        assert "restore failed" in caplog.text
        # The old copy and the staging folder are both left beside out; the next sweep puts the old copy back.
        assert not out.exists()
        sweep_staging(out)
        assert ">old" in (out / "final.contigs.fa").read_text()


class TestAssembleWithIdentity:
    def test_matching_marker_is_skipped(self, tmp_path):
        state = {}
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        write_marker(out, _inputs(tmp_path), "v1", {})
        with patch(RUN, side_effect=_fake_tools(state)):
            _, ran = assemble_extracted_reads(_reads(tmp_path), out, expected=_inputs(tmp_path))
        assert ran is False and _megahit_calls(state) == []

    def test_stale_marker_is_redone(self, tmp_path, caplog):
        state = {}
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        write_marker(out, _inputs(tmp_path, min_contig_len=1000), "v1", {})
        expected = _inputs(tmp_path)
        with patch(RUN, side_effect=_fake_tools(state)), caplog.at_level(logging.INFO):
            _, ran = assemble_extracted_reads(
                _reads(tmp_path), out, expected=expected, version="MEGAHIT v1.2.9", params={"threads": 4}
            )
        assert ran is True and len(_megahit_calls(state)) == 1
        assert ">c1" in (out / "final.contigs.fa").read_text()
        marker = read_marker(out)
        assert AssemblyInputs.from_dict(marker["inputs"]) == expected
        assert marker["megahit_version"] == "MEGAHIT v1.2.9" and marker["params"] == {"threads": 4}
        assert "min_contig_len" in caplog.text
        assert _staging_left(out) == []

    def test_megahit_output_folder_does_not_exist_when_it_starts(self, tmp_path):
        """Real megahit refuses an existing -o folder; the staging folder is only a name until it runs."""
        state = {}
        fake = _fake_tools(state)
        existed = []

        def _checking(executable, args, **kwargs):
            command = [str(part) for part in args]
            if executable == "megahit" and "-o" in command:
                existed.append(Path(command[command.index("-o") + 1]).exists())
            return fake(executable, args, **kwargs)

        out = tmp_path / "SRR1" / "GCF_1_assembly"
        with patch(RUN, side_effect=_checking):
            _, ran = assemble_extracted_reads(_reads(tmp_path), out, expected=_inputs(tmp_path))
        assert ran is True
        assert existed == [False]

    def test_megahit_writes_into_hidden_staging(self, tmp_path):
        state = {}
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        with patch(RUN, side_effect=_fake_tools(state)):
            assemble_extracted_reads(_reads(tmp_path), out, expected=_inputs(tmp_path))
        target = Path(_megahit_calls(state)[0][_megahit_calls(state)[0].index("-o") + 1])
        assert target.parent == out.parent and target.name.startswith(".GCF_1_assembly.")
        assert target.name.endswith(".tmp")
        assert not (out / "intermediate_contigs").exists()

    def test_unmarked_folder_accepted_gets_a_marker(self, tmp_path):
        state = {}
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        expected = _inputs(tmp_path)
        with patch(RUN, side_effect=_fake_tools(state)):
            _, ran = assemble_extracted_reads(_reads(tmp_path), out, expected=expected, version="v9")
        assert ran is False and _megahit_calls(state) == []
        assert ">old" in (out / "final.contigs.fa").read_text()
        assert AssemblyInputs.from_dict(read_marker(out)["inputs"]) == expected

    def test_unmarked_folder_redone_when_not_accepted(self, tmp_path):
        state = {}
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        with patch(RUN, side_effect=_fake_tools(state)):
            _, ran = assemble_extracted_reads(_reads(tmp_path), out, expected=_inputs(tmp_path), accept_unmarked=False)
        assert ran is True and ">c1" in (out / "final.contigs.fa").read_text()
        assert read_marker(out) is not None

    def test_force_redoes_a_current_assembly(self, tmp_path):
        state = {}
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        write_marker(out, _inputs(tmp_path), "v1", {})
        with patch(RUN, side_effect=_fake_tools(state)):
            _, ran = assemble_extracted_reads(_reads(tmp_path), out, expected=_inputs(tmp_path), force=True)
        assert ran is True and ">c1" in (out / "final.contigs.fa").read_text()

    def test_megahit_failure_leaves_old_folder_and_no_staging(self, tmp_path):
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        fake = _fake_tools({})

        def failing(executable, args, **kwargs):
            if executable == "megahit" and "-o" in args:
                fake(executable, args, **kwargs)  # leaves a half-written staging folder behind
                raise subprocess.CalledProcessError(1, ["megahit"], stderr="boom")
            return fake(executable, args, **kwargs)

        with patch(RUN, side_effect=failing):
            with pytest.raises(ProcessingError, match="megahit failed"):
                assemble_extracted_reads(_reads(tmp_path), out, expected=_inputs(tmp_path), force=True)
        assert ">old" in (out / "final.contigs.fa").read_text()
        assert _staging_left(out) == []

    def test_megahit_without_contigs_is_an_error_and_keeps_the_old_folder(self, tmp_path):
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        with patch(RUN):
            with pytest.raises(ProcessingError, match="final.contigs.fa"):
                assemble_extracted_reads(_reads(tmp_path), out, force=True)
        assert ">old" in (out / "final.contigs.fa").read_text()
        assert _staging_left(out) == []

    def test_tmp_dir_inside_out_dir_is_refused(self, tmp_path):
        out = tmp_path / "SRR1" / "GCF_1_assembly"
        with patch(RUN, side_effect=_fake_tools({})) as run:
            with pytest.raises(ProcessingError, match="must not be output_dir"):
                assemble_extracted_reads(_reads(tmp_path), out, tmp_dir=out / "tmp")
        run.assert_not_called()
        assert not out.exists()

    def test_scan_assemblies_ignores_staging(self, tmp_path):
        targeted = tmp_path / "targeted"
        out = targeted / "SRR1" / "GCF_1_assembly"
        _old_assembly(out)
        staging = ident.unique_temp_path(out)
        staging.mkdir()
        aside = ident._aside_path(out)
        aside.mkdir()
        assert reg.scan_assemblies(targeted, []) == {"SRR1": {"GCF_1": out}}
