"""
Optional dependencies, imported only when a command needs them.

The core install (pandas, numpy, matplotlib, biopython, lxml, requests) is enough to load the
CLI and run the download and containment steps. Packages behind an extra are imported at the
point of use through ``require``, which turns a missing package into a ``ConfigurationError``
that names the extra to install and the interpreter to install it into.
"""

import importlib
import sys
from types import ModuleType

from prepare.core.exceptions import ConfigurationError

# Import names whose PyPI distribution is named differently, for the error message.
_DISTRIBUTION_NAMES = {"sklearn": "scikit-learn"}


def install_hint(extra: str) -> str:
    """Install command for one extra, naming the interpreter this process runs under.

    A user with several Python environments then installs into the one PREPARE actually uses.
    """
    return f"{sys.executable} -m pip install 'prepare[{extra}]'"


def require(module: str, extra: str, purpose: str) -> ModuleType:
    """Import ``module`` or raise a ``ConfigurationError`` naming the extra that provides it.

    Args:
        module: Module to import, for example ``"sklearn"`` or ``"plotly.graph_objects"``.
        extra: The PREPARE extra that installs it, for example ``"analysis"``.
        purpose: What needs it, used as the start of the message ("Diversity analysis").

    Returns:
        The imported module.

    Raises:
        ConfigurationError: If the module or its top-level package cannot be imported.
    """
    package = module.partition(".")[0]
    try:
        # Import the top-level package first: a submodule already in sys.modules would
        # otherwise be returned even when the package itself is blocked or absent.
        importlib.import_module(package)
        return importlib.import_module(module)
    except ImportError as e:
        name = _DISTRIBUTION_NAMES.get(package, package)
        raise ConfigurationError(
            f"{purpose} needs the '{name}' package. Install it into this interpreter with: {install_hint(extra)}"
        ) from e
