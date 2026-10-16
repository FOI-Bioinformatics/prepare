"""Tests for prepare.data.ncbi_pacing.RequestPacer, the shared NCBI request pacer."""

import threading
from unittest.mock import patch

import pytest

from prepare.data.ncbi_pacing import RequestPacer


class _FakeClock:
    """A monotonic clock that advances only when something sleeps on it."""

    def __init__(self, start: float = 100.0):
        self.now = start
        self.sleeps = []

    def monotonic(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.sleeps.append(seconds)
        self.now += seconds


@pytest.fixture
def clock():
    fake = _FakeClock()
    with patch("prepare.data.ncbi_pacing.time.monotonic", side_effect=fake.monotonic):
        with patch("prepare.data.ncbi_pacing.time.sleep", side_effect=fake.sleep):
            yield fake


def test_delay_without_key_is_one_third_second():
    assert RequestPacer(None).delay == pytest.approx(1 / 3)


def test_delay_with_key_is_one_tenth_second():
    assert RequestPacer("my-key").delay == pytest.approx(1 / 10)


def test_custom_rates():
    assert RequestPacer(None, per_second_without_key=2).delay == pytest.approx(0.5)
    assert RequestPacer("k", per_second_with_key=5).delay == pytest.approx(0.2)


def test_first_wait_does_not_sleep(clock):
    RequestPacer(None).wait()
    assert clock.sleeps == []


def test_two_waits_are_at_least_delay_apart(clock):
    pacer = RequestPacer(None)
    pacer.wait()
    first = clock.now
    pacer.wait()
    assert clock.now - first >= pacer.delay - 1e-9
    assert clock.sleeps == [pytest.approx(1 / 3)]


def test_sleeps_only_the_remaining_gap(clock):
    pacer = RequestPacer("my-key")
    pacer.wait()
    clock.now += 0.04
    pacer.wait()
    assert clock.sleeps == [pytest.approx(0.06)]


def test_no_sleep_when_delay_already_elapsed(clock):
    pacer = RequestPacer(None)
    pacer.wait()
    clock.now += 1.0
    pacer.wait()
    assert clock.sleeps == []


def test_concurrent_waits_are_each_spaced_by_delay(clock):
    pacer = RequestPacer("my-key")
    threads = [threading.Thread(target=pacer.wait) for _ in range(5)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    # The fake clock moves only on sleep, so every wait after the first sleeps one full gap.
    assert clock.sleeps == [pytest.approx(pacer.delay)] * 4
