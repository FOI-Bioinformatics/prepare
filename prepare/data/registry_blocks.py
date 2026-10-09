"""Typed forms of the blocks a project registry records, per dataset and per project.

The registry file stays plain JSON (layout and ``SCHEMA_VERSION`` unchanged); these dataclasses
are the in-memory view writers build and readers consume instead of chains of ``dict.get``.

Every block converts with ``from_dict`` (tolerant of missing keys) and ``to_dict``, and the
conversion is faithful in both directions:

- A key the block does not declare is kept in ``extra`` and written back unchanged, so a key
  added by a newer version, or by hand, is never lost by a rewrite.
- A block loaded with ``from_dict`` writes back exactly the keys it was loaded with, plus any
  field assigned or changed in place since (``discard`` removes one again). A block that some
  writer only ever wrote in part, e.g. a download record holding nothing but a verdict, keeps
  that shape.
- A block constructed directly writes every field, except the optional ones (``ranked``,
  ``inferred``, a download's ``source``, ...) while they hold their default, which matches
  what the writers wrote before these classes existed.

Values are passed through without coercion, so an integer stays an integer and a JSON
``null`` stays ``None``.
"""

from __future__ import annotations

from dataclasses import MISSING, Field, dataclass, field, fields
from datetime import datetime
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Self,
    Set,
    Tuple,
    Type,
    TypeVar,
)

if TYPE_CHECKING:
    # Only for annotations: prepare.data.registry imports this module.
    from prepare.data.registry import Registry

# Field metadata keys: "omit" marks a field left out of a fresh block's dict while it holds its
# default; "load"/"dump" convert a nested value (skipped for None).
_OMIT = "omit"
_LOAD = "load"
_DUMP = "dump"
# Attributes of the base class itself, never written as keys.
_INTERNAL = ("extra", "_keys")


def _optional(default: Any = None) -> Any:
    """A field whose key a fresh block leaves out while it holds ``default``."""
    return field(default=default, metadata={_OMIT: True})


# Metadata of the optional timing fields (``started``, ``seconds``), left out while None. They are
# declared with ``init=False`` (inline, so mypy sees it): timing is set only by
# ``prepare.data.registry_timing`` after a step is recorded, never when a block is built.
_TIMING = {_OMIT: True}


def _copy(value: Any) -> Any:
    """A shallow copy of a list or dict, so a block never aliases the registry's own containers."""
    if isinstance(value, dict):
        return dict(value)
    if isinstance(value, list):
        return list(value)
    return value


class _FieldSpec(NamedTuple):
    """What the conversion needs to know about one declared field, worked out once per class."""

    name: str
    default: Any
    factory: Optional[Callable[[], Any]]
    load: Optional[Callable[[Any], Any]]
    dump: Optional[Callable[[Any], Any]]
    omit: bool

    def default_value(self) -> Any:
        """The value the field takes when its key is absent."""
        return self.factory() if self.factory is not None else self.default


_SPECS: Dict[type, Tuple[_FieldSpec, ...]] = {}


def _field_specs(cls: type) -> Tuple[_FieldSpec, ...]:
    """The declared (key-carrying) fields of a block class, cached per class."""
    found = _SPECS.get(cls)
    if found is None:
        found = _SPECS[cls] = tuple(_spec_of(f) for f in fields(cls) if f.name not in _INTERNAL)
    return found


def _spec_of(f: Field[Any]) -> _FieldSpec:
    factory = f.default_factory if f.default_factory is not MISSING else None
    default = f.default if f.default is not MISSING else None
    meta = f.metadata
    return _FieldSpec(f.name, default, factory, meta.get(_LOAD), meta.get(_DUMP), bool(meta.get(_OMIT)))


@dataclass
class RegistryBlock:
    """Shared conversion logic; see the module docstring for the round-trip rules."""

    extra: Dict[str, Any] = field(default_factory=dict, kw_only=True)
    # Declared keys the source dict held, plus fields assigned since; None for a fresh block.
    _keys: Optional[Set[str]] = field(default=None, init=False, repr=False, compare=False)

    def __setattr__(self, name: str, value: Any) -> None:
        keys = self.__dict__.get("_keys")
        if keys is not None and name not in _INTERNAL:
            keys.add(name)
        object.__setattr__(self, name, value)

    @classmethod
    def from_dict(cls, data: Optional[Mapping[str, Any]]) -> Self:
        """Build the block from its JSON form; missing keys take the field defaults."""
        data = data or {}
        state: Dict[str, Any] = {}
        present: Set[str] = set()
        for spec in _field_specs(cls):
            if spec.name in data:
                value = data[spec.name]
                state[spec.name] = spec.load(value) if spec.load is not None and value is not None else _copy(value)
                present.add(spec.name)
            else:
                state[spec.name] = spec.default_value()
        # Filled in directly rather than through __init__: a registry can hold tens of thousands
        # of datasets, and status converts each of their blocks several times over.
        block = cls.__new__(cls)
        block.__dict__.update(state)
        block.__dict__["extra"] = {k: _copy(v) for k, v in data.items() if k not in state}
        block.__dict__["_keys"] = present
        return block

    def to_dict(self) -> Dict[str, Any]:
        """The JSON form: the keys this block holds (see the module docstring), then ``extra``."""
        out: Dict[str, Any] = {}
        keys = self._keys
        for spec in _field_specs(type(self)):
            value = getattr(self, spec.name)
            if (keys is None and not spec.omit) or (keys is not None and spec.name in keys):
                pass
            elif value == spec.default_value():
                continue
            out[spec.name] = spec.dump(value) if spec.dump is not None and value is not None else _copy(value)
        out.update(self.extra)
        return out

    def discard(self, name: str) -> None:
        """Reset field ``name`` to its default and leave its key out of ``to_dict``."""
        spec = next(spec for spec in _field_specs(type(self)) if spec.name == name)
        object.__setattr__(self, name, spec.default_value())
        if self._keys is not None:
            self._keys.discard(name)


def _load_or_none(cls: type[RegistryBlock], value: Any) -> Any:
    """One nested block from its dict; a JSON null stays None so it is written back as null."""
    return cls.from_dict(value) if value is not None else None


def _dump_or_none(block: Optional[RegistryBlock]) -> Any:
    return block.to_dict() if block is not None else None


def _nested(cls: type[RegistryBlock], default: Any = None, omit: bool = False) -> Any:
    """A field holding one nested block (or None)."""
    metadata = {_LOAD: cls.from_dict, _DUMP: lambda block: block.to_dict(), _OMIT: omit}
    return field(default=default, metadata=metadata)


def _nested_list(cls: type[RegistryBlock]) -> Any:
    """A field holding a list of nested blocks."""
    metadata = {
        _LOAD: lambda items: [_load_or_none(cls, item) for item in items],
        _DUMP: lambda blocks: [_dump_or_none(block) for block in blocks],
    }
    return field(default_factory=list, metadata=metadata)


def _nested_map(cls: type[RegistryBlock], omit: bool = False) -> Any:
    """A field holding a mapping of name to nested block."""
    metadata = {
        _LOAD: lambda items: {k: _load_or_none(cls, v) for k, v in items.items()},
        _DUMP: lambda blocks: {k: _dump_or_none(block) for k, block in blocks.items()},
        _OMIT: omit,
    }
    return field(default_factory=dict, metadata=metadata)


# ------------------------------------------------------------------ screening


@dataclass
class ScreeningEntry(RegistryBlock):
    """One genome's screening result for one accession."""

    containment: Optional[float] = None
    cani: Optional[float] = None
    csv: Optional[str] = None
    source: str = ""
    query_threshold: float = 0.0


@dataclass
class ScreeningBlock(RegistryBlock):
    """An accession's screening results, keyed by genome id."""

    date: str = ""
    genomes: Dict[str, ScreeningEntry] = _nested_map(ScreeningEntry)
    inferred: bool = _optional(False)


# ---------------------------------------------------------- selection, exclusion


@dataclass
class SelectionBlock(RegistryBlock):
    """Whether an accession is selected, and by which criteria and output file."""

    selected: bool = False
    date: str = ""
    criteria: Dict[str, Any] = field(default_factory=dict)
    output: str = ""
    ranked: Optional[List[Dict[str, Any]]] = _optional()
    inferred: bool = _optional(False)


@dataclass
class ExclusionBlock(RegistryBlock):
    """An accession's exclusion (or a reversed one, with ``excluded`` False)."""

    excluded: bool = False
    reason: str = ""
    source: str = ""
    date: str = ""


# -------------------------------------------------------------------- download


@dataclass
class FileEntry(RegistryBlock):
    """One downloaded read file: project-relative path, size and modification time."""

    path: str = ""
    bytes: int = 0
    mtime: str = ""


@dataclass
class Verdict(RegistryBlock):
    """A download's completeness verdict; producers differ in which keys they fill in."""

    method: Optional[str] = _optional()
    ratio: Optional[float] = _optional()
    verdict: Optional[str] = _optional()
    expected_spots: Optional[int] = _optional()
    reads_r1: Optional[int] = _optional()

    def carrying_counts_from(self, previous: Optional[Verdict]) -> Verdict:
        """This verdict, about to replace ``previous``, with the read counts it may inherit filled in.

        A count this verdict leaves as ``None`` inherits the previous count only when both carry
        the same verdict (a complete-to-complete relink describes the same reads). Across a
        verdict change the previous count describes other files, e.g. a truncated project copy
        relinked to a complete store copy, so it stays ``None``; ``store_verify --rescan`` can
        fill in a real count. Only ``reads_r1`` is carried: ``ratio`` and ``expected_spots``
        always come from the new verdict, so the result never mixes two downloads' verdicts.
        """
        if previous is not None and previous.verdict == self.verdict:
            if self.reads_r1 is None and previous.reads_r1 is not None:
                self.reads_r1 = previous.reads_r1
        return self


@dataclass
class DownloadBlock(RegistryBlock):
    """An accession's download outcome, completeness verdict and cached mate read counts."""

    attempts: int = 0
    state: str = ""
    date: str = ""
    files: List[FileEntry] = _nested_list(FileEntry)
    bytes_total: int = 0
    message: str = ""
    complete: Optional[Verdict] = _nested(Verdict, omit=True)
    source: Optional[str] = _optional()
    store_name: Optional[str] = _optional()
    mate_reads: Optional[List[int]] = _optional()
    # ``[[file name, size bytes, mtime], ...]`` of the mate files ``mate_reads`` was counted from.
    mate_reads_signature: Optional[List[List[Any]]] = _optional()
    inferred: bool = _optional(False)
    # When the last download attempt started (ISO 8601 UTC) and how many seconds it took.
    started: Optional[str] = field(default=None, init=False, metadata=_TIMING)
    seconds: Optional[float] = field(default=None, init=False, metadata=_TIMING)


# -------------------------------------------------------------------- metadata


@dataclass
class MetadataBlock(RegistryBlock):
    """The NCBI run metadata recorded for an accession."""

    xml: str = ""
    date: str = ""
    run_size: Any = None
    run_md5: Any = None
    assay_type: Any = None
    organism: Any = None
    collection_date: Any = None
    library_layout: Any = None
    platform: Any = None
    library_strategy: Any = None
    run_total_spots: Optional[int] = None
    run_total_bases: Optional[int] = None
    inferred: bool = _optional(False)


# ---------------------------------------------------------- extraction, assembly


@dataclass
class AssemblyBlock(RegistryBlock):
    """One assembly's provenance and stats; stats beyond the four required ones live in ``extra``."""

    date: str = ""
    dir: str = ""
    tool: str = ""
    version: str = ""
    params: Dict[str, Any] = field(default_factory=dict)
    contigs: int = 0
    total_bp: int = 0
    n50: int = 0
    largest: int = 0
    # What the assembly was built from, left out while None. Like the timing fields it is set only
    # after the assembly is recorded (``prepare.data.registry_assembly``), never when it is built.
    inputs: Optional[Dict[str, Any]] = field(default=None, init=False, metadata={_OMIT: True})
    # When megahit started (ISO 8601 UTC) and how many seconds it ran.
    started: Optional[str] = field(default=None, init=False, metadata=_TIMING)
    seconds: Optional[float] = field(default=None, init=False, metadata=_TIMING)


@dataclass
class ExtractionBlock(RegistryBlock):
    """One accession's targeted read extraction against one genome, with its assembly if any."""

    date: str = ""
    genome_fasta: Optional[str] = None
    preset: Any = None
    threshold: Any = None
    filter_flags: Any = None
    min_mapq: Any = None
    index: Any = None
    mapped_reads: Optional[int] = None
    mapped_total: Optional[int] = None
    unequal_mates: bool = False
    files: List[str] = field(default_factory=list)
    breadth: Optional[float] = None
    mean_depth: Optional[float] = None
    coverage_tsv: Optional[str] = None
    assembly: Optional[AssemblyBlock] = _nested(AssemblyBlock)
    inferred: bool = _optional(False)
    # When the extraction of this sample started (ISO 8601 UTC) and how many seconds it took.
    started: Optional[str] = field(default=None, init=False, metadata=_TIMING)
    seconds: Optional[float] = field(default=None, init=False, metadata=_TIMING)


# ------------------------------------------------------------ analyses, exports


@dataclass
class AnalysisEntry(RegistryBlock):
    """The latest run of one named analysis for an accession."""

    date: str = ""
    output: str = ""
    summary: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ExportEntry(RegistryBlock):
    """The latest run of one project-level export, e.g. the results table."""

    date: str = ""
    output: str = ""
    summary: Dict[str, Any] = field(default_factory=dict)


# -------------------------------------------------------------- project, store


@dataclass
class ProjectBlock(RegistryBlock):
    """This project's identity (bound by ``store_init`` or minted on first store use) and its exports."""

    id: Optional[str] = _optional()
    name: Optional[str] = _optional()
    path: Optional[str] = _optional()
    created: Optional[str] = _optional()
    exports: Dict[str, ExportEntry] = _nested_map(ExportEntry, omit=True)


@dataclass
class StoreBlock(RegistryBlock):
    """The shared data store this project is bound to, and the accessions it links from it."""

    root: Optional[str] = None
    mode: str = "symlink"
    linked: List[str] = field(default_factory=list)


# ------------------------------------------------------------ typed block access
#
# Readers take a typed copy of a block with ``<name>_block``; writers change the copy and store
# it with ``set_<name>_block``. The copy never aliases the registry's dicts, so a block is only
# written back when a setter is called. The module docstring says how unknown
# keys and partly written blocks survive the round trip.

_B = TypeVar("_B", bound=RegistryBlock)


def _dataset_block(registry: Registry, accession: str, key: str, cls: Type[_B]) -> Optional[_B]:
    raw = registry.datasets.get(accession, {}).get(key)
    return cls.from_dict(raw) if isinstance(raw, dict) else None


def raw(registry: Registry, accession: str, block: str, key: str, default: Any = None) -> Any:
    """One value from ``accession``'s ``block`` read straight from the registry dict, or ``default``.

    For stage checks that look at a single field of every dataset (``query``, ``stage_counts``):
    building a full typed block per dataset just to read one value costs about ten times more
    on a large registry. The value is the registry's own, not a copy; do not change it. For
    ``block="extractions"`` the ``key`` is a genome id and the value is that extraction's dict.
    """
    found = registry.datasets.get(accession, {}).get(block)
    if not isinstance(found, dict):
        return default
    return found.get(key, default)


def screening_block(registry: Registry, accession: str) -> Optional[ScreeningBlock]:
    """``accession``'s screening block, or None when it was never screened."""
    return _dataset_block(registry, accession, "screening", ScreeningBlock)


def selection_block(registry: Registry, accession: str) -> Optional[SelectionBlock]:
    """``accession``'s selection block, or None when no selection ever named it."""
    return _dataset_block(registry, accession, "selection", SelectionBlock)


def exclusion_block(registry: Registry, accession: str) -> Optional[ExclusionBlock]:
    """``accession``'s exclusion block, or None when it was never excluded."""
    return _dataset_block(registry, accession, "exclusion", ExclusionBlock)


def download_block(registry: Registry, accession: str) -> Optional[DownloadBlock]:
    """``accession``'s download block, or None when nothing about its download was recorded."""
    return _dataset_block(registry, accession, "download", DownloadBlock)


def download_verdict(registry: Registry, accession: str) -> Optional[Verdict]:
    """``accession``'s recorded completeness verdict, or None when it has none."""
    download = download_block(registry, accession)
    return download.complete if download is not None else None


def download_block_for_write(registry: Registry, accession: str) -> DownloadBlock:
    """``accession``'s download block, or the ``{"attempts": 0}`` block a first write starts from."""
    return download_block(registry, accession) or DownloadBlock.from_dict({"attempts": 0})


def set_mate_reads(registry: Registry, accession: str, mate_reads: List[int], signature: List[List[Any]]) -> None:
    """Cache ``accession``'s mate read counts, with the ``[name, size, mtime]`` signature of the files counted."""
    download = download_block_for_write(registry, accession)
    download.mate_reads = mate_reads
    download.mate_reads_signature = signature
    set_download_block(registry, accession, download)


def metadata_block(registry: Registry, accession: str) -> Optional[MetadataBlock]:
    """``accession``'s NCBI metadata block, or None when none was recorded."""
    return _dataset_block(registry, accession, "metadata", MetadataBlock)


def extraction_block(registry: Registry, accession: str, genome_id: str) -> Optional[ExtractionBlock]:
    """``accession``'s extraction against ``genome_id``, or None when there is none."""
    raw = registry.datasets.get(accession, {}).get("extractions", {}).get(genome_id)
    return ExtractionBlock.from_dict(raw) if isinstance(raw, dict) else None


def extraction_blocks(registry: Registry, accession: str) -> Dict[str, ExtractionBlock]:
    """Every recorded extraction of ``accession``, keyed by genome id (entries set to null are skipped)."""
    extractions = registry.datasets.get(accession, {}).get("extractions") or {}
    return {g: ExtractionBlock.from_dict(raw) for g, raw in extractions.items() if isinstance(raw, dict)}


def _analysis_summary(registry: Registry, accession: str, name: str) -> Dict[str, Any]:
    """The summary of ``accession``'s analysis ``name``, or ``{}`` when it was never recorded."""
    raw = (registry.datasets.get(accession, {}).get("analyses") or {}).get(name)
    return AnalysisEntry.from_dict(raw).summary if isinstance(raw, dict) else {}


def _analysis_entry(registry: Registry, accession: str, name: str) -> Optional[AnalysisEntry]:
    """``accession``'s analysis ``name`` as a typed entry (date and summary), or None if unrecorded."""
    raw = (registry.datasets.get(accession, {}).get("analyses") or {}).get(name)
    return AnalysisEntry.from_dict(raw) if isinstance(raw, dict) else None


def _recency(date: str) -> float:
    """A sortable recency score for ``date`` (``datetime.fromisoformat``); unparsable text is oldest."""
    try:
        return datetime.fromisoformat(date).timestamp()
    except (TypeError, ValueError):
        return float("-inf")


def _quality_fields(entry: AnalysisEntry) -> Dict[str, Any]:
    """``total_reads``, ``gc_percent`` and ``quality_grade`` out of one analysis entry's summary."""
    return {key: entry.summary.get(key) for key in ("total_reads", "gc_percent", "quality_grade")}


def _legacy_quality_summary(registry: Registry, accession: str) -> Tuple[Dict[str, Any], bool]:
    """``total_reads``, ``gc_percent`` and ``quality_grade`` from a pre-0.5.0 registry, and whether either existed.

    A registry written before 0.5.0 holds the ``"sra_stats"`` analysis (GC already in percent) and
    the ``"quality"`` one of ``sra_profile_quality`` (GC as a 0-1 fraction, the grade under
    ``"grade"``); ``sra_stats`` is read first for totals and GC.
    """
    stats = _analysis_summary(registry, accession, "sra_stats")
    quality = _analysis_summary(registry, accession, "quality")
    quality_gc = quality.get("gc_content")
    gc_percent = stats.get("gc_content")
    if gc_percent is None and quality_gc is not None:
        gc_percent = quality_gc * 100
    total_reads = stats.get("total_reads")
    summary = {
        "total_reads": total_reads if total_reads is not None else quality.get("total_reads"),
        "gc_percent": gc_percent,
        "quality_grade": quality.get("grade"),
    }
    return summary, bool(stats) or bool(quality)


def quality_summary(registry: Registry, accession: str) -> Tuple[Dict[str, Any], Optional[str]]:
    """``total_reads``, ``gc_percent`` and ``quality_grade`` for ``accession``, and where they came from.

    ``sra_profile`` records the ``"profile"`` analysis; ``sra_report`` records ``"report"``, which
    also carries these three fields. Whichever is newer by ``datetime.fromisoformat`` (an
    unparsable date counts as the oldest; equal dates favour ``"profile"``) supplies the result,
    with any field it leaves ``None`` filled in from the other one when both exist, and a field
    still ``None`` after that filled in from the pre-0.5.0 analyses when the registry holds them
    (the source stays ``"profile"`` or ``"report"``). When neither exists, the pre-0.5.0
    ``"sra_stats"``/``"quality"`` analyses are read instead (see ``_legacy_quality_summary``).
    The second element of the pair names the source used: ``"profile"``, ``"report"``,
    ``"legacy"``, or None when nothing was ever recorded (every field is then None too).
    """
    profile_entry = _analysis_entry(registry, accession, "profile")
    report_entry = _analysis_entry(registry, accession, "report")
    if profile_entry is not None or report_entry is not None:
        primary: AnalysisEntry
        secondary: Optional[AnalysisEntry]
        source: str
        if profile_entry is not None and report_entry is not None:
            if _recency(profile_entry.date) >= _recency(report_entry.date):
                primary, secondary, source = profile_entry, report_entry, "profile"
            else:
                primary, secondary, source = report_entry, profile_entry, "report"
        elif profile_entry is not None:
            primary, secondary, source = profile_entry, None, "profile"
        else:
            assert report_entry is not None  # the outer `or` guarantees this
            primary, secondary, source = report_entry, None, "report"
        primary_fields = _quality_fields(primary)
        secondary_fields = _quality_fields(secondary) if secondary is not None else {}
        summary = {
            key: value if value is not None else secondary_fields.get(key) for key, value in primary_fields.items()
        }
        if any(value is None for value in summary.values()):
            legacy, _ = _legacy_quality_summary(registry, accession)
            summary = {key: value if value is not None else legacy.get(key) for key, value in summary.items()}
        return summary, source
    legacy, found = _legacy_quality_summary(registry, accession)
    return legacy, ("legacy" if found else None)


def profile_summary(registry: Registry, accession: str) -> Dict[str, Any]:
    """``total_reads``, ``gc_percent`` and ``quality_grade`` for ``accession``; see ``quality_summary``."""
    return quality_summary(registry, accession)[0]


def project_block(registry: Registry) -> ProjectBlock:
    """This project's identity and exports (empty fields when ``store_init`` never ran)."""
    return ProjectBlock.from_dict(registry.project)


def store_block(registry: Registry) -> StoreBlock:
    """The shared store this project is bound to (``root`` None when it is bound to none)."""
    return StoreBlock.from_dict(registry.store)


def set_screening_block(registry: Registry, accession: str, block: ScreeningBlock) -> None:
    """Store ``block`` as ``accession``'s screening block."""
    registry.datasets.setdefault(accession, {})["screening"] = block.to_dict()


def set_selection_block(registry: Registry, accession: str, block: SelectionBlock) -> None:
    """Store ``block`` as ``accession``'s selection block."""
    registry.datasets.setdefault(accession, {})["selection"] = block.to_dict()


def set_exclusion_block(registry: Registry, accession: str, block: ExclusionBlock) -> None:
    """Store ``block`` as ``accession``'s exclusion block."""
    registry.datasets.setdefault(accession, {})["exclusion"] = block.to_dict()


def set_download_block(registry: Registry, accession: str, block: DownloadBlock) -> None:
    """Store ``block`` as ``accession``'s download block."""
    registry.datasets.setdefault(accession, {})["download"] = block.to_dict()


def set_metadata_block(registry: Registry, accession: str, block: MetadataBlock) -> None:
    """Store ``block`` as ``accession``'s metadata block."""
    registry.datasets.setdefault(accession, {})["metadata"] = block.to_dict()


def set_extraction_block(registry: Registry, accession: str, genome_id: str, block: ExtractionBlock) -> None:
    """Store ``block`` as ``accession``'s extraction against ``genome_id``."""
    registry.datasets.setdefault(accession, {}).setdefault("extractions", {})[genome_id] = block.to_dict()


def set_project_block(registry: Registry, block: ProjectBlock) -> None:
    """Store ``block`` as this project's identity and exports."""
    registry.project = block.to_dict()


def set_store_block(registry: Registry, block: StoreBlock) -> None:
    """Store ``block`` as this project's store binding."""
    registry.store = block.to_dict()


# The dataset blocks bootstrap and reconcile can mark as reconstructed from disk.
_INFERABLE: Dict[str, Type[RegistryBlock]] = {
    "screening": ScreeningBlock,
    "selection": SelectionBlock,
    "download": DownloadBlock,
    "metadata": MetadataBlock,
    "extractions": ExtractionBlock,
}


def mark_inferred(registry: Registry, accession: str, key: str, genome_id: Optional[str] = None, **values: Any) -> None:
    """Flag a block just recorded for ``accession`` as reconstructed from disk, and set ``values`` on it.

    ``key`` names the dataset block (``"screening"``, ``"selection"``, ``"download"``,
    ``"metadata"``), or ``"extractions"`` together with ``genome_id`` for one extraction. A block
    that was never recorded is left alone.
    """
    container = registry.datasets.get(accession, {})
    name = key
    if genome_id is not None:
        container, name = container.get(key) or {}, genome_id
    raw = container.get(name)
    if not isinstance(raw, dict):
        return
    block = _INFERABLE[key].from_dict(raw)
    for field_name, value in {"inferred": True, **values}.items():
        setattr(block, field_name, value)
    container[name] = block.to_dict()
