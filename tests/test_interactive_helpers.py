"""Tests for the shared figure writer and matrix preparation in the interactive plots."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("plotly")

from prepare.core.exceptions import VisualizationError  # noqa: E402
from prepare.visualization.interactive import _prepare_matrix, _write_figure  # noqa: E402


class _FailingFigure:
    def write_html(self, path):
        with open(path, "w") as handle:
            handle.write("partial")
        raise RuntimeError("disk full")


class _GoodFigure:
    def write_html(self, path):
        with open(path, "w") as handle:
            handle.write("<html></html>")


def test_write_figure_leaves_no_temp_file_on_failure(tmp_path):
    target = tmp_path / "sub" / "plot.html"
    with pytest.raises(RuntimeError):
        _write_figure(_FailingFigure(), target)
    assert not target.exists()
    assert list(target.parent.iterdir()) == []


def test_write_figure_replaces_target_and_creates_folder(tmp_path):
    target = tmp_path / "new" / "plot.html"
    assert _write_figure(_GoodFigure(), target) == target
    assert target.read_text() == "<html></html>"
    assert [p.name for p in target.parent.iterdir()] == ["plot.html"]


def test_prepare_matrix_names_and_shape_checks():
    frame = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], index=["a", "b"])
    X, names = _prepare_matrix(frame)
    assert names == ["a", "b"] and X.shape == (2, 2)
    _, generated = _prepare_matrix(np.ones((3, 2)))
    assert generated == ["Sample_0", "Sample_1", "Sample_2"]
    with pytest.raises(VisualizationError, match="t-SNE requires at least 2 samples"):
        _prepare_matrix(np.ones((1, 2)), method="t-SNE")
    with pytest.raises(VisualizationError, match="PCA requires at least 1 feature"):
        _prepare_matrix(np.ones((3, 0)))
