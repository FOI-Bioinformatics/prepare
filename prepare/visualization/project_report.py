"""The project report as one self-contained HTML page.

``render_html`` draws the same section blocks as the Markdown report
(``prepare.processing.project_report_markdown.section_blocks``), each section in a
``<section id="<name>">``, plus two Plotly figures: the cross-stage funnel and the extracted and
assembled counts per genome. plotly and jinja2 come from the ``interactive`` extra and are
imported through ``require``; ``require_html`` checks for both before any file is written. The
template autoescapes every value, plotly.js is embedded inline (``plotly_js_script``), and the
page uses the shared ``REPORT_CSS``.
"""

from typing import Any, Dict, List

from prepare.core.optional import require
from prepare.processing.project_report import SECTION_TITLES
from prepare.processing.project_report_markdown import format_cell, section_blocks, subtitle
from prepare.utils.html import REPORT_CSS, plotly_js_script, plotly_layout

PURPOSE = "The project report's HTML page"
FUNNEL_STAGES = ("screened", "selected", "downloaded", "analysed", "extracted", "assembled")


def require_html() -> None:
    """Raise a ``ConfigurationError`` naming the ``interactive`` extra unless plotly and jinja2 import."""
    require("plotly.graph_objects", "interactive", PURPOSE)
    require("plotly.offline", "interactive", PURPOSE)
    require("jinja2", "interactive", PURPOSE)


def _figure_div(figure: Any) -> str:
    offline = require("plotly.offline", "interactive", PURPOSE)
    return str(offline.plot(figure, output_type="div", include_plotlyjs=False))


def _funnel_figure(funnel: Dict[str, Any]) -> str:
    go = require("plotly.graph_objects", "interactive", PURPOSE)
    stages = [stage for stage in FUNNEL_STAGES if stage in funnel]
    figure = go.Figure(go.Funnel(y=stages, x=[funnel[stage].get("accessions", 0) for stage in stages]))
    figure.update_layout(**plotly_layout(title="Datasets reaching each stage", height=360))
    return _figure_div(figure)


def _genome_figure(genomes: Dict[str, Any]) -> str:
    go = require("plotly.graph_objects", "interactive", PURPOSE)
    names = list(genomes)
    figure = go.Figure(
        [
            go.Bar(name="extracted", x=names, y=[genomes[g]["extracted"] for g in names]),
            go.Bar(name="assembled", x=names, y=[genomes[g]["assembled"] for g in names]),
        ]
    )
    figure.update_layout(**plotly_layout(title="Datasets per genome", barmode="group", height=360))
    return _figure_div(figure)


def _sections(report: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Template data per section: id, title, blocks with formatted cells and an optional figure."""
    figures = {"funnel": _funnel_figure(report["funnel"])}
    if report["genomes"]:
        figures["genomes"] = _genome_figure(report["genomes"])
    sections = []
    for name, blocks in section_blocks(report).items():
        sections.append(
            {
                "id": name,
                "title": SECTION_TITLES[name],
                "figure": figures.get(name),
                "blocks": [
                    {
                        "text": block.text,
                        "header": list(block.header),
                        "rows": [[format_cell(cell) for cell in row] for row in block.rows],
                    }
                    for block in blocks
                ],
            }
        )
    return sections


def render_html(report: Dict[str, Any]) -> str:
    """The report as one HTML page with every section of ``SECTION_TITLES``; needs plotly and jinja2."""
    require_html()
    jinja2 = require("jinja2", "interactive", PURPOSE)
    template = jinja2.Environment(loader=jinja2.BaseLoader(), autoescape=True).from_string(_PAGE_TEMPLATE)
    project = report["project"]
    return str(
        template.render(
            title="Project report",
            subtitle=subtitle(report),
            project=project,
            funnel=report["funnel"],
            sections=_sections(report),
            plotly_js=plotly_js_script(),
            report_css=REPORT_CSS,
        )
    )


_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{{ title }}</title>
    {{ plotly_js|safe }}
    <style>{{ report_css|safe }}</style>
    <style>
        /* Static tables: no sort control, and space between the blocks of a section. */
        table.pre-table th { cursor: default; }
        table.pre-table th::after { content: none; }
        .pre-section .pre-table-wrap, .pre-section .pre-panel { margin: 0 0 14px; }
    </style>
</head>
<body>
    <a class="pre-skip" href="#project">Skip to the report</a>
    <header class="pre-header"><div class="pre-wrap">
        <p class="pre-eyebrow">PREPARE &middot; project report</p>
        <h1 class="pre-title">{{ project.name or title }}</h1>
        <p class="pre-readout"><span>{{ subtitle }}</span></p>
    </div></header>
    <main class="pre-wrap">
        <section class="pre-stats" aria-label="Funnel summary">
        {% for stage in ["screened", "selected", "downloaded", "analysed", "extracted", "assembled"] %}
            {% if stage in funnel %}
            <div class="pre-stat"><p class="k">{{ stage }}</p><div class="v">{{ funnel[stage].accessions }}</div></div>
            {% endif %}
        {% endfor %}
        </section>
        {% for section in sections %}
        <section class="pre-section" id="{{ section.id }}">
            <h2>{{ section.title }}</h2>
            {% if section.figure %}<div class="pre-panel">{{ section.figure|safe }}</div>{% endif %}
            {% for block in section.blocks %}
                {% if block.header %}
                <div class="pre-table-wrap"><table class="pre-table">
                    <thead><tr>{% for cell in block.header %}<th scope="col">{{ cell }}</th>{% endfor %}</tr></thead>
                    <tbody>
                    {% for row in block.rows %}
                        <tr>{% for cell in row %}<td>{{ cell }}</td>{% endfor %}</tr>
                    {% endfor %}
                    </tbody>
                </table></div>
                {% else %}
                <p>{{ block.text }}</p>
                {% endif %}
            {% endfor %}
        </section>
        {% endfor %}
        <footer class="pre-footer">Written by prepare project_report; the same content is in
            project_report.md and project_report.json.</footer>
    </main>
</body>
</html>
"""
