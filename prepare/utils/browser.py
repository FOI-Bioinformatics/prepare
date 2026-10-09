"""Helpers for opening generated HTML output in a web browser.

Used by the commands that produce interactive HTML (dashboards, the
containment explorer) so a developer can view the result immediately.
"""

import logging
import platform
import subprocess
import webbrowser
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


def open_in_browser(path: Union[str, Path]) -> bool:
    """Open a local HTML file in the default browser.

    The path is turned into a percent-encoded ``file://`` URI so locations
    containing spaces work. On macOS the ``open`` command is used directly:
    Python's ``webbrowser`` backend there (``MacOSXOSAScript``) can report
    success without actually opening anything. Other platforms use the
    ``webbrowser`` module.

    Args:
        path: Path to the HTML file to open.

    Returns:
        True if a browser was launched, False otherwise.
    """
    resolved = Path(path).resolve()
    uri = resolved.as_uri()

    if platform.system() == "Darwin":
        try:
            subprocess.run(["open", str(resolved)], check=True)
            return True
        except Exception as e:
            logger.warning("Could not open browser with `open`: %s", e)
            return False

    try:
        return webbrowser.open(uri)
    except Exception as e:
        logger.warning("Could not open browser: %s", e)
        return False
