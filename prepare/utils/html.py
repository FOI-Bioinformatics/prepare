"""Shared HTML assets for PREPARE report and dashboard generators.

Centralizes the inline Plotly ``<script>`` tag, the shared "assay readout"
stylesheet and Plotly styling, and the sortable/filterable table JavaScript so
the several HTML generators do not each carry their own copy.
"""

from typing import Any, Dict, List

from prepare.core.optional import require

# --- Design system palette ----------------------------------------------------
# Categorical hues (validated as a set, worst adjacent CVD dE 24.2), assigned in
# fixed order; used for taxonomy families, dataset groups, and quality grades.
CATEGORICAL_COLORS: List[str] = [
    "#2a78d6",  # blue
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#008300",  # green
    "#4a3aa7",  # violet
    "#e34948",  # red
    "#e87ba4",  # magenta
    "#eb6834",  # orange
]

# Single-hue sequential ramp (teal-green, light -> dark) for containment, the one
# magnitude the reports encode. Used identically in the sunburst, the heatmap,
# and the table heat-bars.
CONTAINMENT_SCALE: List[str] = [
    "#e8f5ef",
    "#c2e6d7",
    "#93d3ba",
    "#5ebd9a",
    "#2ea67d",
    "#1b8663",
    "#12664a",
    "#0b4735",
]

_SANS = 'system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'
_MONO = 'ui-monospace, "SF Mono", "JetBrains Mono", Menlo, Consolas, monospace'

# Ink used on the (always light) chart panels; see REPORT_CSS for the token set.
_PANEL_INK = "#172420"
_PANEL_MUTED = "#5f6b66"
_PANEL_LINE = "#e7ebe8"


def plotly_layout(**overrides: Any) -> Dict[str, Any]:
    """Return shared Plotly ``update_layout`` kwargs for the report figures.

    Figures are transparent so they sit on the light chart panel, use the
    report's sans/mono fonts, recede their axes and gridlines, and share the
    validated categorical colourway. Pass overrides (e.g. ``height=500``) to
    extend or replace individual keys.
    """
    layout: Dict[str, Any] = dict(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        colorway=CATEGORICAL_COLORS,
        font=dict(family=_SANS, size=13, color=_PANEL_INK),
        title=dict(font=dict(family=_SANS, size=15, color=_PANEL_INK), x=0.01, xanchor="left"),
        margin=dict(t=46, l=10, r=10, b=10),
        legend=dict(font=dict(family=_SANS, size=11, color=_PANEL_INK)),
        xaxis=dict(
            gridcolor=_PANEL_LINE,
            zeroline=False,
            linecolor="#dfe4e0",
            tickfont=dict(family=_MONO, size=11, color=_PANEL_MUTED),
            title=dict(font=dict(family=_SANS, size=12, color=_PANEL_MUTED)),
        ),
        yaxis=dict(
            gridcolor=_PANEL_LINE,
            zeroline=False,
            linecolor="#dfe4e0",
            tickfont=dict(family=_MONO, size=11, color=_PANEL_MUTED),
            title=dict(font=dict(family=_SANS, size=12, color=_PANEL_MUTED)),
        ),
        coloraxis=dict(
            colorbar=dict(
                tickfont=dict(family=_MONO, size=10, color=_PANEL_MUTED),
                outlinewidth=0,
                thickness=12,
            )
        ),
    )
    layout.update(overrides)
    return layout


def plotly_js_script() -> str:
    """Return a ``<script>`` tag with plotly.js embedded inline.

    Uses the plotly.js bundled with the installed plotly.py, so the runtime
    always matches the figures it renders and the report renders with no
    network access (the individual figures are emitted with
    ``include_plotlyjs=False`` and share this single embedded copy).

    Raises:
        ConfigurationError: If plotly is not installed; the callers check for it
            before drawing any figure, so this is reached only by a direct call.
    """
    offline = require("plotly.offline", "interactive", "Embedding plotly.js in a report")
    return f'<script type="text/javascript">{offline.get_plotlyjs()}</script>'


# The shared "assay readout" stylesheet. Data (accessions, containment values,
# counts) is set in a monospace face; containment gets one teal ramp everywhere;
# in dark mode the page darkens but the chart/table panels stay light (a
# "lightbox"), which also lets the static Plotly figures read in both themes.
# Literal CSS braces are fine: this is a plain string dropped in via a Jinja
# ``|safe`` variable or substituted into an f-string, never re-parsed.
REPORT_CSS = """
:root {
  --paper-bg: #f4f6f4;
  --paper-ink: #172420;
  --paper-muted: #5f6b66;
  --panel-bg: #ffffff;
  --panel-ink: #172420;
  --panel-muted: #5f6b66;
  --line: #dfe4e0;
  --line-soft: #e9ede9;
  --accent: #b26a2b;
  --seq-track: #e8f0ec;
  --seq-fill: #1b8663;
  --seq-strong: #12664a;
  --radius: 10px;
  --sans: system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  --mono: ui-monospace, "SF Mono", "JetBrains Mono", Menlo, Consolas, monospace;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--paper-bg);
  color: var(--paper-ink);
  font-family: var(--sans);
  font-size: 15px;
  line-height: 1.5;
  -webkit-font-smoothing: antialiased;
}
.mq-skip {
  position: absolute; left: -999px; top: 0;
  background: var(--panel-bg); color: var(--panel-ink);
  padding: 8px 14px; border-radius: 6px; z-index: 20;
}
.mq-skip:focus { left: 12px; top: 12px; }
.mq-wrap { max-width: 1240px; margin: 0 auto; padding: 0 24px 64px; }

.mq-header { padding: 30px 0 20px; border-bottom: 1px solid var(--line); margin-bottom: 28px; }
.mq-header .mq-wrap { padding-bottom: 0; padding-top: 0; }
.mq-eyebrow {
  font-family: var(--mono); font-size: 11px; letter-spacing: 0.14em;
  text-transform: uppercase; color: var(--paper-muted); margin: 0 0 6px;
}
.mq-title {
  font-size: 1.9rem; font-weight: 650; letter-spacing: -0.015em; margin: 0;
  line-height: 1.1; display: flex; align-items: baseline; gap: 12px;
}
.mq-title::before {
  content: ""; width: 10px; height: 26px; border-radius: 2px;
  background: var(--accent); transform: translateY(3px); flex: none;
}
.mq-readout {
  font-family: var(--mono); font-size: 12.5px; color: var(--paper-muted);
  margin: 12px 0 0; display: flex; flex-wrap: wrap; gap: 6px 14px;
}
.mq-readout b { color: var(--paper-ink); font-weight: 600; }

.mq-stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 14px; margin: 0 0 34px; }
.mq-stat {
  background: var(--panel-bg); border: 1px solid var(--line);
  border-radius: var(--radius); padding: 15px 17px;
}
.mq-stat .k {
  font-family: var(--mono); font-size: 10.5px; letter-spacing: 0.08em;
  text-transform: uppercase; color: var(--panel-muted); margin: 0 0 8px;
}
.mq-stat .v {
  font-family: var(--mono); font-size: 1.7rem; font-weight: 600;
  color: var(--panel-ink); line-height: 1; font-variant-numeric: tabular-nums;
}
.mq-stat.wide { grid-column: span 2; }
.mq-stat .sub { margin-top: 8px; font-size: 12.5px; color: var(--panel-muted); line-height: 1.45; }

.mq-section { margin: 0 0 34px; }
.mq-section > h2 {
  font-size: 0.82rem; font-weight: 600; letter-spacing: 0.12em; text-transform: uppercase;
  color: var(--paper-muted); margin: 0 0 14px; padding-bottom: 8px;
  border-bottom: 1px solid var(--line);
}
.mq-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 20px; }
.mq-panel {
  background: var(--panel-bg); border: 1px solid var(--line);
  border-radius: var(--radius); padding: 10px 12px; overflow: hidden;
}

.mq-toolbar { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-bottom: 12px; }
.mq-toolbar input[type="text"] {
  flex: 1 1 260px; min-width: 200px; padding: 9px 13px;
  border: 1px solid var(--line); border-radius: 8px; background: var(--panel-bg);
  color: var(--panel-ink); font-family: var(--mono); font-size: 13px;
}
.mq-toolbar input::placeholder { color: var(--panel-muted); }
.mq-btn {
  padding: 9px 15px; border: 1px solid var(--line); border-radius: 8px;
  background: var(--panel-bg); color: var(--panel-ink); font-family: var(--sans);
  font-size: 13px; font-weight: 550; cursor: pointer;
}
.mq-btn:hover { border-color: var(--seq-fill); color: var(--seq-strong); }

.mq-table-wrap {
  background: var(--panel-bg); border: 1px solid var(--line); border-radius: var(--radius);
  max-height: 620px; overflow: auto;
}
table.mq-table { width: 100%; border-collapse: collapse; font-size: 13px; }
table.mq-table th {
  position: sticky; top: 0; z-index: 1; text-align: left; cursor: pointer;
  background: var(--panel-bg); color: var(--panel-muted);
  font-family: var(--mono); font-size: 10.5px; letter-spacing: 0.06em; text-transform: uppercase;
  font-weight: 600; padding: 11px 14px; border-bottom: 1px solid var(--line);
  white-space: nowrap; user-select: none;
}
table.mq-table th::after { content: "\\2195"; opacity: 0.25; margin-left: 6px; }
table.mq-table th.is-sorted { color: var(--seq-strong); }
table.mq-table th.is-sorted::after { content: "\\2191"; opacity: 0.9; }
table.mq-table th.is-sorted.is-desc::after { content: "\\2193"; }
table.mq-table td {
  padding: 9px 14px; border-bottom: 1px solid var(--line-soft);
  font-family: var(--mono); color: var(--panel-ink); white-space: nowrap;
  font-variant-numeric: tabular-nums;
}
table.mq-table tbody tr:nth-child(even) { background: rgba(27,134,99,0.035); }
table.mq-table tbody tr:hover { background: rgba(27,134,99,0.09); }
td.mq-num { text-align: right; }
.mq-cell { display: inline-flex; align-items: center; gap: 9px; }
.mq-heatbar {
  --v: 0; position: relative; display: inline-block; width: 54px; height: 9px;
  border-radius: 2px; background: var(--seq-track); flex: none;
}
.mq-heatbar::before {
  content: ""; position: absolute; inset: 0 auto 0 0; min-width: 2px;
  width: calc(var(--v) * 100%); background: var(--seq-fill); border-radius: 2px;
}

.mq-ok { color: #0b6b45; font-weight: 600; }
.mq-bad { color: #c0322f; font-weight: 600; }
.mq-note {
  background: var(--panel-bg); border: 1px solid var(--line);
  border-left: 3px solid var(--seq-fill); border-radius: var(--radius);
  padding: 14px 18px; margin: 0 0 18px; color: var(--panel-ink);
}
.mq-note.warn { border-left-color: var(--accent); }
.mq-note h3 { margin: 0 0 8px; font-size: 0.95rem; }
.mq-note ul { margin: 6px 0 0; padding-left: 20px; }
.mq-note li { margin: 3px 0; }

.mq-footer {
  margin-top: 40px; padding-top: 18px; border-top: 1px solid var(--line);
  font-family: var(--mono); font-size: 11.5px; color: var(--paper-muted);
}

:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; border-radius: 3px; }

@media (max-width: 720px) {
  .mq-wrap { padding: 0 16px 48px; }
  .mq-title { font-size: 1.5rem; }
  .mq-grid { grid-template-columns: 1fr; }
  .mq-stat.wide { grid-column: span 1; }
}

@media (prefers-reduced-motion: reduce) {
  * { transition: none !important; animation: none !important; }
}

@media (prefers-color-scheme: dark) {
  :root {
    --paper-bg: #0e1512;
    --paper-ink: #e8ede9;
    --paper-muted: #93a49c;
    --line: #26312c;
  }
  body { -webkit-font-smoothing: auto; }
  /* Panels (stats, charts, table) intentionally stay light: a lightbox for the
     data, and it keeps the static Plotly figures legible in dark mode. */
}
"""


# Client-side helpers for the interactive results table: full-text filtering,
# column sorting (numeric-aware, keyboard-operable, with aria-sort), and CSV
# export of the visible rows. Uses single braces so it can be dropped into
# either a Jinja template (via a ``|safe`` variable) or an f-string.
TABLE_SCRIPT = """<script>
function filterTable() {
    var input = document.getElementById("searchInput").value.toUpperCase();
    var table = document.getElementById("resultsTable");
    var tr = table.getElementsByTagName("tr");
    for (var i = 1; i < tr.length; i++) {
        var cells = tr[i].getElementsByTagName("td");
        var match = false;
        for (var j = 0; j < cells.length; j++) {
            if (cells[j].textContent.toUpperCase().indexOf(input) > -1) {
                match = true; break;
            }
        }
        tr[i].style.display = match ? "" : "none";
    }
}

function sortTable(col) {
    var table = document.getElementById("resultsTable");
    var rows = Array.from(table.tBodies[0].rows);
    var asc = table.getAttribute("data-sort-col") == col
              && table.getAttribute("data-sort-dir") != "asc";
    rows.sort(function(a, b) {
        var va = a.cells[col].textContent.trim();
        var vb = b.cells[col].textContent.trim();
        var na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
        return asc ? va.localeCompare(vb) : vb.localeCompare(va);
    });
    rows.forEach(function(r) { table.tBodies[0].appendChild(r); });
    table.setAttribute("data-sort-col", col);
    table.setAttribute("data-sort-dir", asc ? "asc" : "desc");
    var ths = table.tHead ? table.tHead.rows[0].cells : [];
    for (var k = 0; k < ths.length; k++) {
        if (k === col) {
            ths[k].setAttribute("aria-sort", asc ? "ascending" : "descending");
            ths[k].classList.add("is-sorted");
            ths[k].classList.toggle("is-desc", !asc);
        } else {
            ths[k].removeAttribute("aria-sort");
            ths[k].classList.remove("is-sorted", "is-desc");
        }
    }
}

function sortKey(e, col) {
    if (e.key === "Enter" || e.key === " ") { e.preventDefault(); sortTable(col); }
}

function exportCSV() {
    var table = document.getElementById("resultsTable");
    var rows = table.querySelectorAll("tr");
    var csv = [];
    rows.forEach(function(row) {
        if (row.style.display === "none") return;
        var cols = row.querySelectorAll("td, th");
        var line = [];
        cols.forEach(function(col) {
            var val = col.textContent.replace(/"/g, '""');
            line.push('"' + val + '"');
        });
        csv.push(line.join(","));
    });
    var blob = new Blob([csv.join("\\n")], {type: "text/csv"});
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "containment_results.csv";
    a.click();
}
</script>"""
