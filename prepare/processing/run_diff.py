"""
Comparisons between recorded runs of the project's run log (``prepare.data.run_log``).

The ``runs`` command uses these to diff two runs and to trace one accession through the log.
A run's detail is whatever its command passed to ``note_run``; the rows compared here are the
mappings of plain values held in its sections, at any depth, keyed by their own key (an accession, or
``accession/genome``). For example, ``sra_profile`` notes ``{"analyses": {"profile": {acc: {...}}}}``
and each accession's mapping is one row. When one row key appears in more than one section, each of
its field names is prefixed with the full path of its section (``analyses.profile.gc_percent``) so the
sections stay apart. A flat mapping beside the sections (settings, counts) is not a row.
"""

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from prepare.data import run_log
from prepare.data.run_log import RunRecord

PathLike = Union[str, Path]
DELTA_DECIMALS = 6


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _flatten(values: Mapping[str, Any], prefix: str = "") -> Dict[str, Any]:
    """``values`` with nested mappings flattened into dotted keys."""
    flat: Dict[str, Any] = {}
    for key, value in values.items():
        name = f"{prefix}{key}"
        if isinstance(value, Mapping) and value:
            flat.update(_flatten(value, f"{name}."))
        else:
            flat[name] = value
    return flat


def _delta(before: Any, after: Any) -> Optional[float]:
    if _is_number(before) and _is_number(after):
        delta = after - before
        return delta if isinstance(delta, int) else round(delta, DELTA_DECIMALS)
    return None


def diff_summaries(a: Mapping[str, Any], b: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Every summary key of two runs, sorted: ``key``, ``before`` (a), ``after`` (b), ``delta``, ``changed``.

    Nested mappings are flattened into dotted keys; a key missing on one side has None there.
    ``delta`` is ``after - before`` when both are numbers (not booleans), rounded to
    ``DELTA_DECIMALS`` places for floats, and None otherwise.
    """
    flat_a, flat_b = _flatten(a or {}), _flatten(b or {})
    rows = []
    for key in sorted(set(flat_a) | set(flat_b)):
        before, after = flat_a.get(key), flat_b.get(key)
        rows.append(
            {"key": key, "before": before, "after": after, "delta": _delta(before, after), "changed": before != after}
        )
    return rows


def _is_row(value: Any) -> bool:
    """A non-empty mapping of plain values (no nested mapping)."""
    return isinstance(value, Mapping) and bool(value) and not any(isinstance(item, Mapping) for item in value.values())


def _row_sections(detail: Mapping[str, Any], path: Tuple[str, ...] = ()) -> List[Tuple[Tuple[str, ...], str, Dict]]:
    """Every row in ``detail`` as (section path, row key, row).

    A section is a mapping whose mapping values are all rows (plain values beside them are
    ignored); its values are the rows. A mapping that also holds a non-row mapping is not a
    section: its non-row mappings are searched in turn and its flat mappings (settings, counts)
    are not rows and are ignored.
    """
    nested = [value for value in detail.values() if isinstance(value, Mapping)]
    if nested and all(_is_row(value) for value in nested):
        return [(path, str(key), dict(value)) for key, value in detail.items() if isinstance(value, Mapping)]
    found: List[Tuple[Tuple[str, ...], str, Dict]] = []
    for key, value in detail.items():
        if isinstance(value, Mapping) and not _is_row(value):
            found.extend(_row_sections(value, path + (str(key),)))
    return found


def detail_rows(detail: Optional[Mapping[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """The rows of one run's detail, keyed by row key (an accession, or ``accession/genome``).

    A row is a mapping of plain values inside a section (see ``_row_sections``). When the same
    row key appears in more than one section, each of its fields is prefixed with the full path of
    the section holding it (``analyses.profile.gc_percent``); a key found in one section only
    keeps its field names, so they do not change with the other sections a detail holds.
    Empty for None or a detail without rows.
    """
    if not detail:
        return {}
    found = _row_sections(detail)
    sections_of: Dict[str, set] = {}
    for path, key, _ in found:
        sections_of.setdefault(key, set()).add(path)
    rows: Dict[str, Dict[str, Any]] = {}
    for path, key, row in found:
        if len(sections_of[key]) > 1 and path:
            prefix = ".".join(path)
            row = {f"{prefix}.{field}": value for field, value in row.items()}
        rows.setdefault(key, {}).update(row)
    return rows


def diff_details(a: Optional[Mapping[str, Any]], b: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    """Compare the rows of two runs' details (see ``detail_rows``).

    Returns ``added`` (row keys only in b), ``removed`` (only in a), ``changed``
    (``{key: {field: [value in a, value in b]}}``, a field missing on one side as None) and
    ``unchanged`` (the number of rows identical in both).
    """
    rows_a, rows_b = detail_rows(a), detail_rows(b)
    changed: Dict[str, Dict[str, List[Any]]] = {}
    unchanged = 0
    for key in sorted(set(rows_a) & set(rows_b)):
        row_a, row_b = rows_a[key], rows_b[key]
        fields = {
            field: [row_a.get(field), row_b.get(field)]
            for field in sorted(set(row_a) | set(row_b))
            if row_a.get(field) != row_b.get(field)
        }
        if fields:
            changed[key] = fields
        else:
            unchanged += 1
    return {
        "added": sorted(set(rows_b) - set(rows_a)),
        "removed": sorted(set(rows_a) - set(rows_b)),
        "changed": changed,
        "unchanged": unchanged,
    }


def detail_kept(project: PathLike, record: RunRecord) -> bool:
    """Whether ``record``'s detail file is still on disk (False when it had none or it was pruned)."""
    if not record.detail:
        return False
    return (run_log.runs_dir(project) / Path(record.detail).name).is_file()


def _matches(key: str, accession: str) -> bool:
    """Whether row ``key`` belongs to ``accession``: the key itself or ``accession/<genome>``."""
    return key == accession or key.startswith(f"{accession}/")


def _names_accession(argv: Sequence[str], accession: str) -> bool:
    """Whether a command line names ``accession`` as an argument of its own or as ``--option=ACC``.

    The contents of an accessions file named on the command line are not read: a pruned run that
    took its accessions from ``--accessions-file`` cannot be matched.
    """
    return any(item == accession or item.partition("=")[2] == accession for item in argv)


def accession_history(project: PathLike, records: Sequence[RunRecord], accession: str) -> List[Dict[str, Any]]:
    """The runs, oldest first, that recorded values for ``accession``.

    Each entry holds ``run_id``, ``command``, ``started``, ``exit_code``, ``detail_kept`` and
    ``values`` (``{row key: row}`` for the rows of that accession). A run whose detail is no
    longer kept is included, with ``values`` None, only when its command line names the
    accession (``_names_accession``); other runs without a detail cannot be searched and are left out.
    """
    history = []
    for record in records:
        kept = detail_kept(project, record)
        values: Optional[Dict[str, Dict[str, Any]]] = None
        if kept:
            rows = detail_rows(run_log.read_detail(project, record))
            values = {key: row for key, row in rows.items() if _matches(key, accession)} or None
            if values is None:
                continue
        elif not _names_accession(record.argv, accession):
            continue
        history.append(
            {
                "run_id": record.run_id,
                "command": record.command,
                "started": record.started,
                "exit_code": record.exit_code,
                "detail_kept": kept,
                "values": values,
            }
        )
    return history
