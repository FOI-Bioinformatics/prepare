"""Tests for prepare.core.validation.validate_folder."""

import pytest

from prepare.core.exceptions import DataAccessError, ValidationError, exit_code_for
from prepare.core.validation import validate_folder


def test_a_folder_that_cannot_be_created_is_a_data_access_error(tmp_path):
    blocker = tmp_path / "plain_file"
    blocker.write_text("x")

    with pytest.raises(DataAccessError, match="Failed to create directory") as caught:
        validate_folder(blocker / "sub", create=True)

    assert not isinstance(caught.value, ValidationError)
    assert exit_code_for(caught.value) == 1


def test_a_missing_folder_without_create_is_invalid_input(tmp_path):
    with pytest.raises(ValidationError, match="Directory does not exist"):
        validate_folder(tmp_path / "absent")
