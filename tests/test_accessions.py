import pytest
from prepare.core.exceptions import ValidationError
from prepare.data.accessions import read_accessions_file, require_accession


def test_reader_skips_comments_blanks_and_duplicates(tmp_path):
    f = tmp_path / "acc.txt"
    f.write_text("# header\nSRR1\n\n  SRR2  \nSRR1\n")
    assert read_accessions_file(f) == ["SRR1", "SRR2"]


@pytest.mark.parametrize("bad", ["../../victim", "..", "/tmp/x", "SRR1/..", "SRR 1", "srr1"])
def test_reader_refuses_unsafe_lines(tmp_path, bad):
    f = tmp_path / "acc.txt"
    f.write_text(f"SRR1\n{bad}\n")
    with pytest.raises(ValidationError) as info:
        read_accessions_file(f)
    assert "line 2" in str(info.value) and bad in str(info.value)


def test_require_accession():
    assert require_accession("ERR123") == "ERR123"
    with pytest.raises(ValidationError):
        require_accession("../x")


def test_reader_accepts_a_pattern_and_can_skip_validation(tmp_path):
    f = tmp_path / "genomes.txt"
    f.write_text("# genomes\nGCF_000006945.2\nGCA_000001405.29\nGCF_000006945.2\n")
    assert read_accessions_file(f, validate=False) == ["GCF_000006945.2", "GCA_000001405.29"]
    assert read_accessions_file(f, pattern=r"GC[FA]_[0-9]+\.[0-9]+") == ["GCF_000006945.2", "GCA_000001405.29"]
    with pytest.raises(ValidationError):
        read_accessions_file(f)
    with pytest.raises(ValidationError, match="line 3"):
        read_accessions_file(f, pattern=r"GCF_[0-9]+\.[0-9]+")
