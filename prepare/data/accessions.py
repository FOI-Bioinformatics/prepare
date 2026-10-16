"""The accessions file: one SRA run accession per line."""

import re
from pathlib import Path
from typing import List, Union

from prepare.core.constants import SRA_ACCESSION_PATTERN
from prepare.core.exceptions import DataAccessError, ValidationError

_ACCESSION_RE = re.compile(SRA_ACCESSION_PATTERN)


def require_accession(accession: str) -> str:
    """``accession`` when it is a well-formed SRA run accession; ``ValidationError`` otherwise.

    The pattern admits only upper-case letters and digits, so a value that passes is safe
    to use as a folder name.
    """
    if not _ACCESSION_RE.fullmatch(accession or ""):
        raise ValidationError(f"Invalid SRA accession format: {accession!r}")
    return accession


def read_accessions_file(
    path: Union[str, Path], validate: bool = True, pattern: str = SRA_ACCESSION_PATTERN
) -> List[str]:
    """Accessions listed in ``path`` in file order, without blanks, comments or repeats.

    With ``validate`` each line must match ``pattern`` (by default a well-formed SRA run
    accession); the error names the file, the line number and the offending text. Files of other
    kinds, such as genome assembly accessions, are read with ``validate=False`` or with their own
    ``pattern``. An unreadable file raises ``DataAccessError``.
    """
    try:
        with open(path, "r") as handle:
            lines = handle.read().splitlines()
    except OSError as e:
        raise DataAccessError(f"Cannot read accessions file {path}: {e}") from e
    accession_re = _ACCESSION_RE if pattern == SRA_ACCESSION_PATTERN else re.compile(pattern)
    seen: dict = {}
    for number, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if validate and not accession_re.fullmatch(line):
            what = "an SRA accession" if accession_re is _ACCESSION_RE else "an accepted accession"
            raise ValidationError(f"{path}, line {number}: {line!r} is not {what}")
        seen.setdefault(line, None)
    return list(seen)
