"""``prepare.utils.http.retrying_session``: the retry policy shared by every HTTP client."""

import pytest

from prepare.utils.http import RETRY_STATUS_FORCELIST, RETRY_TOTAL, retrying_session
from tests.fake_http import serve


def test_session_is_mounted_on_both_schemes():
    session = retrying_session()
    assert "https://" in session.adapters
    assert "http://" in session.adapters


def test_default_retry_does_not_raise_on_exhausted_status():
    """``raise_on_status=False``: a session built with defaults returns the last response, never raises."""
    session = retrying_session()
    retry = session.adapters["https://"].max_retries
    assert retry.raise_on_status is False
    assert retry.total == RETRY_TOTAL
    assert set(retry.status_forcelist) == set(RETRY_STATUS_FORCELIST)


def test_retries_a_persistent_503_then_gives_the_last_response(monkeypatch):
    requested = serve(monkeypatch, lambda path: (503, b"unavailable"))
    session = retrying_session()

    response = session.get("https://example.invalid/thing")

    assert response.status_code == 503
    assert len(requested) == RETRY_TOTAL + 1  # one first attempt plus RETRY_TOTAL retries


def test_a_transient_failure_then_success_is_retried_through(monkeypatch):
    attempts = {"n": 0}

    def handler(path):
        attempts["n"] += 1
        if attempts["n"] < 2:
            return (503, b"unavailable")
        return (200, b"ok")

    serve(monkeypatch, handler)
    session = retrying_session()

    response = session.get("https://example.invalid/thing")

    assert response.status_code == 200
    assert response.content == b"ok"


def test_custom_policy_overrides_the_module_defaults(monkeypatch):
    requested = serve(monkeypatch, lambda path: (500, b"err"))
    session = retrying_session(total=1, status_forcelist=(500,))

    response = session.get("https://example.invalid/thing")

    assert response.status_code == 500
    assert len(requested) == 2  # one first attempt plus 1 retry


if __name__ == "__main__":
    pytest.main([__file__])
