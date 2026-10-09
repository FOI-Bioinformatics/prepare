"""Tests for the browser-opening helper in prepare.utils.browser."""

from pathlib import Path
from unittest.mock import patch

from prepare.utils import browser


class TestOpenInBrowser:
    """Tests for open_in_browser across platforms."""

    def test_macos_uses_open_command(self, tmp_path):
        """On macOS the `open` command is invoked with the resolved path."""
        html = tmp_path / "dashboard.html"
        html.write_text("<html></html>")

        with patch.object(browser.platform, "system", return_value="Darwin"):
            with patch.object(browser.subprocess, "run") as mock_run:
                result = browser.open_in_browser(html)

        assert result is True
        args, _ = mock_run.call_args
        assert args[0] == ["open", str(html.resolve())]

    def test_macos_returns_false_when_open_fails(self, tmp_path):
        """A failing `open` command yields False rather than raising."""
        html = tmp_path / "dashboard.html"
        html.write_text("<html></html>")

        with patch.object(browser.platform, "system", return_value="Darwin"):
            with patch.object(browser.subprocess, "run", side_effect=OSError("boom")):
                assert browser.open_in_browser(html) is False

    def test_non_macos_uses_webbrowser_with_encoded_uri(self, tmp_path):
        """Other platforms hand webbrowser a percent-encoded file URI."""
        # A space in the path must be encoded so the browser can open it.
        html = tmp_path / "my report.html"
        html.write_text("<html></html>")

        with patch.object(browser.platform, "system", return_value="Linux"):
            with patch.object(browser.webbrowser, "open", return_value=True) as mock_open:
                result = browser.open_in_browser(html)

        assert result is True
        (uri,) = mock_open.call_args[0]
        assert uri.startswith("file://")
        assert "my%20report.html" in uri
        assert uri == Path(html).resolve().as_uri()

    def test_non_macos_returns_false_on_webbrowser_error(self, tmp_path):
        """A webbrowser failure yields False rather than raising."""
        html = tmp_path / "dashboard.html"
        html.write_text("<html></html>")

        with patch.object(browser.platform, "system", return_value="Linux"):
            with patch.object(browser.webbrowser, "open", side_effect=RuntimeError("no display")):
                assert browser.open_in_browser(html) is False
