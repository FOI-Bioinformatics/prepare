"""Tests for the plot_map command."""

import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd

from prepare.cli.commands.metadata import PlotMapCommand
from prepare.core.exceptions import ExitCode
from prepare.core.optional import install_hint
from prepare.plugins.base import visualizer_registry


def _args(tmp_path, **overrides):
    table = tmp_path / "metadata.tsv"
    pd.DataFrame({"acc": ["SRR1", "SRR2"], "lat_lon": ["10 N 20 E", "5 S 3 W"]}).to_csv(table, sep="\t", index=False)
    values = {"metadata_file": str(table), "output": str(tmp_path / "map.png"), "lat_lon_column": "lat_lon"}
    values.update(overrides)
    return argparse.Namespace(**values)


def test_plot_map_reports_missing_maps_extra(tmp_path, monkeypatch, caplog):
    monkeypatch.setitem(sys.modules, "cartopy", None)
    assert PlotMapCommand().execute(_args(tmp_path)) == ExitCode.CONFIGURATION
    assert install_hint("maps") in caplog.text
    assert not (tmp_path / "map.png").exists()


def test_plot_map_writes_the_figure_and_closes_it(tmp_path, monkeypatch):
    seen = {}

    def create_plot(data, **kwargs):
        seen["kwargs"] = kwargs
        seen["rows"] = len(data)
        seen["fig"] = plt.figure()
        seen["fig"].savefig(kwargs["output_file"], format=kwargs["output_format"])
        return seen["fig"]

    plugin = visualizer_registry.get("map")
    monkeypatch.setattr(plugin, "create_plot", staticmethod(create_plot))
    assert PlotMapCommand().execute(_args(tmp_path, lat_lon_column="coords")) == 0
    assert (tmp_path / "map.png").stat().st_size > 0
    assert seen["rows"] == 2 and seen["kwargs"]["lat_lon_column"] == "coords"
    assert not plt.fignum_exists(seen["fig"].number)


def test_plot_map_missing_input_fails_cleanly(tmp_path):
    assert PlotMapCommand().execute(_args(tmp_path, metadata_file=str(tmp_path / "none.tsv"))) == ExitCode.FAILURE
