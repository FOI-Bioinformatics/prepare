"""
Basic unit tests for PREPARE.
"""

import prepare


def test_version():
    """Test that version is defined."""
    assert prepare.__version__ is not None
