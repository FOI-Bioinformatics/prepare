"""
Advanced SRA Reporting and Dashboard Generation.

This module provides comprehensive reporting capabilities including:
- Interactive HTML dashboards with Plotly visualizations
- Quality control reports with recommendations
- Comparative analysis reports across datasets

Plotly and jinja2 come from the ``interactive`` extra. Each report checks for both before it
profiles anything or writes a file, so a missing package is an error naming the extra rather
than a report without charts.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Union, Any

import pandas as pd
import numpy as np

from prepare.sra.analytics import (
    HIGH_ADAPTER_FRACTION,
    NUMERIC_COLUMNS,
    AnomalyReport,
    QualityProfile,
    ComparativeAnalysis,
    SRADatasetAnalyzer,
)
from prepare.core.optional import require
from prepare.data.file_io import write_text_atomic
from prepare.utils.html import (
    CATEGORICAL_COLORS,
    chart_panel,
    figure_title,
    page_footer,
    page_head,
    plotly_layout,
    utc_stamp,
)

# Quality-grade colours drawn from the validated categorical palette (ordinal:
# excellent -> poor).
_GRADE_COLORS = {"excellent": "#008300", "good": "#2a78d6", "fair": "#eda100", "poor": "#e34948"}

logger = logging.getLogger(__name__)

# How each numeric column is named, with its unit, in figures and test headings.
METRIC_LABELS = {
    "avg_read_length": "Mean read length (bp)",
    "gc_percent": "GC content (%)",
    "total_reads": "Reads",
    "complexity_score": "Complexity score",
}

_PURPOSE = "An SRA HTML report"


def _require_report_packages() -> None:
    """Raise a ConfigurationError naming the interactive extra unless plotly and jinja2 import."""
    require("plotly.graph_objects", "interactive", _PURPOSE)
    require("jinja2", "interactive", _PURPOSE)


def _plotly() -> tuple:
    """Return (plotly.graph_objects, plotly.offline)."""
    return (
        require("plotly.graph_objects", "interactive", _PURPOSE),
        require("plotly.offline", "interactive", _PURPOSE),
    )


def _jinja_template(template_str: str) -> Any:
    """Compile an autoescaping jinja2 template."""
    jinja2 = require("jinja2", "interactive", _PURPOSE)
    return jinja2.Environment(loader=jinja2.BaseLoader(), autoescape=True).from_string(template_str)


class SRAReportGenerator:
    """Generate comprehensive reports for SRA download and analysis sessions."""

    def __init__(self, output_dir: Union[str, Path], fastq_dir: Optional[Union[str, Path]] = None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.analyzer = SRADatasetAnalyzer(fastq_dir=fastq_dir)

    def quality_section(
        self, profiles: Dict[str, QualityProfile], anomalies: Optional[AnomalyReport] = None
    ) -> Dict[str, Any]:
        """The data of a report's quality section for already computed ``profiles``.

        Nothing is profiled here. ``anomalies`` is computed from ``profiles`` when not given.
        """
        return {
            "total_datasets": len(profiles),
            "summary_stats": self._calculate_quality_summary(profiles),
            "anomaly_report": anomalies or self.analyzer.detect_dataset_anomalies(list(profiles), profiles=profiles),
            "plots": self._create_quality_plots(profiles),
        }

    def generate_report(
        self,
        profiles: Dict[str, QualityProfile],
        title: str = "SRA Report",
        comparison: Optional[ComparativeAnalysis] = None,
        anomalies: Optional[AnomalyReport] = None,
        filename: str = "sra_report.html",
    ) -> Path:
        """Write one HTML report: a quality section and, with ``comparison``, a comparative one.

        Both sections are built from the ``profiles`` (and the ``comparison`` computed from
        them) that the caller supplies, so no dataset is profiled here.

        Raises:
            ConfigurationError: If plotly or jinja2 is not installed
            ValueError: If ``profiles`` is empty
        """
        _require_report_packages()
        if not profiles:
            raise ValueError("No datasets could be profiled")
        comparative = None
        if comparison is not None:
            comparative = {
                "comparison": comparison,
                "group_counts": {name: len(accessions) for name, accessions in comparison.dataset_groups.items()},
                "plots": self._create_comparative_plots(comparison),
            }
        timestamp = utc_stamp()
        html = _render_page(
            "SRA report", title, timestamp, quality=self.quality_section(profiles, anomalies), comparative=comparative
        )
        report_path = self.output_dir / filename
        write_text_atomic(report_path, html)
        logger.info(f"SRA report saved to {report_path}")
        return report_path

    def _create_quality_plots(self, profiles: Dict[str, QualityProfile]) -> Dict[str, str]:
        """Create interactive plots for quality dashboard."""
        plots: dict = {}
        go, pyo = _plotly()

        # Quality grade distribution
        grades = [p.quality_grade for p in profiles.values()]
        grade_counts = pd.Series(grades).value_counts()

        fig_grades = go.Figure(
            data=[
                go.Bar(
                    x=grade_counts.index,
                    y=grade_counts.values,
                    marker_color=[_GRADE_COLORS.get(str(g), CATEGORICAL_COLORS[0]) for g in grade_counts.index],
                )
            ]
        )
        fig_grades.update_layout(**plotly_layout())
        fig_grades.update_layout(title_text="Quality grade distribution")
        plots["quality_grades"] = chart_panel(
            pyo.plot(fig_grades, output_type="div", include_plotlyjs=False), figure_title(fig_grades)
        )

        # GC content distribution
        gc_percents = [p.gc_percent for p in profiles.values()]
        fig_gc = go.Figure(data=[go.Histogram(x=gc_percents, nbinsx=25)])
        fig_gc.update_layout(**plotly_layout())
        fig_gc.update_layout(title_text="GC content distribution", xaxis_title="GC content (%)", yaxis_title="Count")
        plots["gc_distribution"] = chart_panel(
            pyo.plot(fig_gc, output_type="div", include_plotlyjs=False), figure_title(fig_gc)
        )

        # Read length vs complexity scatter
        read_lengths = [p.avg_read_length for p in profiles.values()]
        complexities = [p.complexity_score for p in profiles.values()]
        accessions = list(profiles.keys())

        fig_complexity = go.Figure(
            data=[
                go.Scatter(
                    x=read_lengths,
                    y=complexities,
                    mode="markers",
                    text=accessions,
                    hovertemplate="<b>%{text}</b><br>Length: %{x:.0f} bp<br>Complexity: %{y:.2f}",
                )
            ]
        )
        fig_complexity.update_layout(**plotly_layout())
        fig_complexity.update_layout(
            title_text="Read length vs sequence complexity",
            xaxis_title="Average read length (bp)",
            yaxis_title="Complexity score",
        )
        plots["complexity_scatter"] = chart_panel(
            pyo.plot(fig_complexity, output_type="div", include_plotlyjs=False), figure_title(fig_complexity)
        )

        return plots

    def _create_comparative_plots(self, comparison: ComparativeAnalysis) -> Dict[str, str]:
        """Create interactive plots for comparative analysis."""
        plots: dict = {}
        if not comparison.visualization_data:
            return plots
        go, pyo = _plotly()

        # Box plots for numeric variables
        boxplot_data = comparison.visualization_data.get("boxplot_data", [])
        if boxplot_data:
            df = pd.DataFrame(boxplot_data)

            for col in NUMERIC_COLUMNS:
                if col in df.columns:
                    fig_box = go.Figure()

                    for group in df["group"].unique():
                        group_data = df[df["group"] == group][col]
                        fig_box.add_trace(go.Box(y=group_data, name=group, boxpoints="outliers"))

                    fig_box.update_layout(**plotly_layout())
                    label = METRIC_LABELS.get(col, col)
                    fig_box.update_layout(title_text=f"{label} by group", yaxis_title=label)
                    plots[f"{col}_boxplot"] = chart_panel(
                        pyo.plot(fig_box, output_type="div", include_plotlyjs=False), figure_title(fig_box)
                    )

        return plots

    def _calculate_quality_summary(self, profiles: Dict[str, QualityProfile]) -> Dict[str, Any]:
        """Calculate summary statistics for quality profiles."""
        if not profiles:
            return {}

        # Totals come from the datasets with an exact count (both mates counted, as the
        # statistics record counts them); a profile that only has its sample's figures is
        # left out and reported as sampled, so a sample size is never added to a total.
        exact = [p for p in profiles.values() if not getattr(p, "sampled", False)]
        total_reads = sum(p.total_reads for p in exact)
        total_bases = sum(p.total_bases for p in exact)
        avg_gc = np.mean([p.gc_percent for p in profiles.values()])
        avg_complexity = np.mean([p.complexity_score for p in profiles.values()])

        # Quality grade distribution
        grades = [p.quality_grade for p in profiles.values()]
        grade_dist = pd.Series(grades).value_counts().to_dict()

        # Contamination statistics
        adapter_contamination = [p.contamination_indicators.get("adapter_contamination", 0) for p in profiles.values()]
        avg_contamination = np.mean(adapter_contamination)
        high_contamination = sum(1 for c in adapter_contamination if c > HIGH_ADAPTER_FRACTION)

        return {
            "total_datasets": len(profiles),
            "total_reads": total_reads,
            "total_bases": total_bases,
            "sampled_datasets": len(profiles) - len(exact),
            "average_gc_percent": avg_gc,
            "average_complexity": avg_complexity,
            "quality_grade_distribution": grade_dist,
            "average_contamination": avg_contamination,
            "high_contamination_count": high_contamination,
        }


_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    {{ head_html|safe }}
</head>
<body>
    <a class="pre-skip" href="#report">Skip to the report</a>
    <header class="pre-header"><div class="pre-wrap">
        <p class="pre-eyebrow">PREPARE &middot; {{ eyebrow }}</p>
        <h1 class="pre-title">{{ title }}</h1>
        <p class="pre-readout"><span>generated <b>{{ timestamp }}</b></span>
        {% if quality %}<span><b>{{ quality.total_datasets }}</b> datasets</span>{% endif %}</p>
    </div></header>
    <main class="pre-wrap" id="report">
    {% if quality %}{% with summary_stats=quality.summary_stats, anomaly_report=quality.anomaly_report,
        plots=quality.plots %}
        <section class="pre-stats" aria-label="Quality summary">
            <div class="pre-stat"><p class="k">Total reads (mates counted)</p>
                <div class="v">{{ "{:,.0f}".format(summary_stats.total_reads) }}</div>
                {% if summary_stats.sampled_datasets %}<p class="k">{{ summary_stats.sampled_datasets }}
                    dataset(s) sampled only, not counted</p>{% endif %}</div>
            <div class="pre-stat"><p class="k">Average GC content</p>
                <div class="v">{{ "%.1f"|format(summary_stats.average_gc_percent) }}%</div></div>
            <div class="pre-stat"><p class="k">High quality</p>
                <div class="v">{{ summary_stats.quality_grade_distribution.get('excellent', 0)
                    + summary_stats.quality_grade_distribution.get('good', 0) }}</div></div>
            <div class="pre-stat"><p class="k">Contamination issues</p>
                <div class="v">{{ summary_stats.high_contamination_count }}</div></div>
        </section>

        {% if anomaly_report.anomalous_datasets %}
        <div class="pre-note warn">
            <h3>Anomalies detected</h3>
            <p><b>{{ anomaly_report.anomalous_datasets|length }}</b> datasets flagged for review:</p>
            <ul>
                {% for dataset in anomaly_report.anomalous_datasets[:10] %}
                <li><b>{{ dataset }}</b>: {{ anomaly_report.explanations.get(dataset, 'Multiple issues') }}</li>
                {% endfor %}
                {% if anomaly_report.anomalous_datasets|length > 10 %}
                <li>&hellip; and {{ anomaly_report.anomalous_datasets|length - 10 }} more</li>
                {% endif %}
            </ul>
        </div>
        {% endif %}

        {% if plots %}
        <section class="pre-section">
            <h2>Quality metrics</h2>
            <div class="pre-grid">
            {% for plot_name, plot_html in plots.items() %}
                {{ plot_html|safe }}
            {% endfor %}
            </div>
        </section>
        {% endif %}
    {% endwith %}{% endif %}

    {% if comparative %}{% with comparison=comparative.comparison, group_counts=comparative.group_counts,
        plots=comparative.plots %}
        <section class="pre-stats" aria-label="Group summary">
            {% for group_name, count in group_counts.items() %}
            <div class="pre-stat"><p class="k">{{ group_name }}</p>
                <div class="v">{{ count }}<span style="font-size:0.9rem"> datasets</span></div></div>
            {% endfor %}
        </section>

        {% if comparison.statistical_tests %}
        <section class="pre-section">
            <h2>Statistical tests</h2>
            {% for test_name, test_result in comparison.statistical_tests.items() %}
            <div class="pre-note{% if test_result.significant %} warn{% endif %}">
                <h3>{{ metric_labels.get(test_name, test_name) }}</h3>
                {# A skipped test carries no p-value (None) or NaN, which is not equal to itself. #}
                {% if test_result.p_value is none or test_result.p_value != test_result.p_value %}
                <p>{{ test_result.test }}{% if test_result.note %} &middot; {{ test_result.note }}{% endif %}</p>
                {% else %}
                <p>{{ test_result.test }} &middot; p {{ "%.4f"|format(test_result.p_value) }}
                    {%- if test_result.p_adjusted is defined -%}
                    , adjusted {{ "%.4f"|format(test_result.p_adjusted) }} (Holm)
                    {%- endif %} &middot;
                    <b>{{ "significant" if test_result.significant else "not significant" }}</b></p>
                {% endif %}
            </div>
            {% endfor %}
        </section>
        {% endif %}

        {% if plots %}
        <section class="pre-section">
            <h2>Comparative visualizations</h2>
            <div class="pre-grid">
            {% for plot_name, plot_html in plots.items() %}
                {{ plot_html|safe }}
            {% endfor %}
            </div>
        </section>
        {% endif %}

        {% if comparison.recommendations %}
        <div class="pre-note">
            <h3>Recommendations</h3>
            <ul>
                {% for rec in comparison.recommendations %}<li>{{ rec }}</li>{% endfor %}
            </ul>
        </div>
        {% endif %}
    {% endwith %}{% endif %}
        {{ footer_html|safe }}
    </main>
</body>
</html>
"""


def _render_page(
    eyebrow: str,
    title: str,
    timestamp: str,
    quality: Optional[Dict[str, Any]],
    comparative: Optional[Dict[str, Any]],
) -> str:
    """Render one report page with a quality section, a comparative section, or both."""
    template = _jinja_template(_PAGE_TEMPLATE)
    return template.render(
        head_html=page_head(title),
        footer_html=page_footer(),
        metric_labels=METRIC_LABELS,
        eyebrow=eyebrow,
        title=title,
        timestamp=timestamp,
        quality=quality,
        comparative=comparative,
    )
