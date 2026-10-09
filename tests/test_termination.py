"""Tests for prepare.utils.termination and the BaseCommand.run wrapper."""

import argparse
import logging
import os
import signal
import threading
from unittest.mock import patch

import pytest

import prepare.utils.termination as term_mod
from prepare.cli.base import BaseCommand, CommandRegistry
from prepare.data import registry_blocks as rb
from prepare.data.registry import load_registry, record_download
from prepare.data.registry_batch import RegistryBatch, registry_batch
from prepare.utils.termination import ABANDON_AFTER_SIGNALS, Termination, graceful_termination

HANDLED = [s for s in (signal.SIGINT, signal.SIGTERM, getattr(signal, "SIGHUP", None)) if s is not None]


def _current_handlers():
    return {s: signal.getsignal(s) for s in HANDLED}


class TestGracefulTermination:
    """The signal state machine: first raises, later ones are counted, the last abandons."""

    def test_first_signal_raises_second_is_logged_third_abandons(self, caplog):
        before = _current_handlers()
        stop = threading.Event()
        with (
            patch.object(term_mod.os, "_exit") as fake_exit,
            patch.object(term_mod.SecureSubprocess, "terminate_children", return_value=0) as terminate,
        ):
            with caplog.at_level(logging.WARNING, logger=term_mod.__name__):
                with graceful_termination(stop) as term:
                    with pytest.raises(KeyboardInterrupt):
                        os.kill(os.getpid(), signal.SIGTERM)
                    assert term.signum == signal.SIGTERM and term.requested and stop.is_set()

                    os.kill(os.getpid(), signal.SIGTERM)
                    assert term.repeats == 1
                    fake_exit.assert_not_called()
                    assert any("again; 1 more will abandon the write" in r.getMessage() for r in caplog.records)

                    os.kill(os.getpid(), signal.SIGTERM)
                    assert term.repeats == 2
        terminate.assert_called_once_with(grace=0)
        fake_exit.assert_called_once_with(130)
        assert _current_handlers() == before

    def test_a_signal_as_the_first_handler_goes_in_still_restores_it(self):
        """The handler can run as ``signal.signal`` returns; the old handler must already be recorded."""
        before = _current_handlers()
        real_signal = signal.signal
        calls = []

        def install_then_interrupt(signum, handler):
            old = real_signal(signum, handler)
            calls.append(signum)
            if len(calls) == 1:
                raise KeyboardInterrupt  # what the new handler does when a signal is pending
            return old

        with patch.object(term_mod.signal, "signal", side_effect=install_then_interrupt):
            with pytest.raises(KeyboardInterrupt):
                with graceful_termination():
                    pass  # pragma: no cover - the install is interrupted
        assert _current_handlers() == before

    def test_abandon_after_one_repeat(self):
        with (
            patch.object(term_mod.os, "_exit") as fake_exit,
            patch.object(term_mod.SecureSubprocess, "terminate_children", return_value=0),
        ):
            with graceful_termination(abandon_after=2):
                with pytest.raises(KeyboardInterrupt):
                    os.kill(os.getpid(), signal.SIGINT)
                os.kill(os.getpid(), signal.SIGINT)
        fake_exit.assert_called_once_with(130)

    def test_default_abandon_count(self):
        assert ABANDON_AFTER_SIGNALS == 3

    def test_handlers_restored_when_the_block_raises(self):
        before = _current_handlers()
        with pytest.raises(ValueError):
            with graceful_termination():
                assert all(signal.getsignal(s) is not before[s] for s in HANDLED)
                raise ValueError("boom")
        assert _current_handlers() == before

    def test_no_signal_means_not_requested(self):
        with graceful_termination() as term:
            pass
        assert not term.requested and term.signum is None and term.repeats == 0

    def test_an_ignored_signal_stays_ignored(self):
        if not hasattr(signal, "SIGHUP"):
            pytest.skip("no SIGHUP on this platform")
        previous = signal.signal(signal.SIGHUP, signal.SIG_IGN)
        try:
            with graceful_termination():
                # A run started under nohup keeps ignoring hangups.
                assert signal.getsignal(signal.SIGHUP) is signal.SIG_IGN
            assert signal.getsignal(signal.SIGHUP) is signal.SIG_IGN
        finally:
            signal.signal(signal.SIGHUP, previous)

    def test_nothing_installed_off_the_main_thread(self):
        seen = {}

        def work():
            try:
                with patch.object(term_mod.signal, "signal") as spy:
                    own_stop = threading.Event()
                    with graceful_termination(own_stop) as term:
                        seen["term"] = term
                    seen["calls"] = spy.call_count
                    seen["stop"] = term.stop is own_stop
                with graceful_termination() as fresh:
                    seen["fresh"] = isinstance(fresh.stop, threading.Event) and not fresh.requested
            except BaseException as e:  # noqa: B902 - reported to the main thread below
                seen["error"] = e

        worker = threading.Thread(target=work)
        worker.start()
        worker.join()
        assert "error" not in seen
        assert seen["calls"] == 0 and seen["stop"] and seen["fresh"]
        assert isinstance(seen["term"], Termination)

    def test_a_nested_context_shares_the_outer_one(self):
        outer_stop, inner_stop = threading.Event(), threading.Event()
        with graceful_termination(outer_stop) as outer:
            installed = _current_handlers()
            with patch.object(term_mod.signal, "signal") as spy:
                with graceful_termination(inner_stop) as inner:
                    assert inner is outer
                    with pytest.raises(KeyboardInterrupt):
                        os.kill(os.getpid(), signal.SIGTERM)
            assert spy.call_count == 0
            assert outer_stop.is_set() and inner_stop.is_set()
            assert _current_handlers() == installed

    def test_second_sigint_during_the_final_flush_lets_it_complete(self, tmp_path, caplog):
        """A repeated Ctrl-C while a batch writes on the way out is logged; the outcome still reaches the file."""
        registry_file = tmp_path / "prepare_registry.json"
        real_apply_all = RegistryBatch._apply_all
        sent = []

        def apply_all_with_signal(self, registry, pending):
            sent.append(True)
            os.kill(os.getpid(), signal.SIGINT)
            real_apply_all(self, registry, pending)

        with caplog.at_level(logging.WARNING, logger=term_mod.__name__):
            with patch.object(RegistryBatch, "_apply_all", apply_all_with_signal):
                with pytest.raises(KeyboardInterrupt):
                    with graceful_termination() as term, registry_batch(registry_file, flush_every=None) as batch:
                        batch.apply(lambda r: record_download(r, "SRR1", "failed", tmp_path / "fastq", "t"), "SRR1")
                        os.kill(os.getpid(), signal.SIGINT)
        assert sent == [True]
        assert term.repeats == 1
        assert rb.download_block(load_registry(registry_file), "SRR1").state == "failed"
        assert any("SIGINT again" in r.getMessage() for r in caplog.records)


class _FakeCommand(BaseCommand):
    """A command whose execute does whatever the test sets on it."""

    def __init__(self, action, graceful=True):
        super().__init__()
        self._action = action
        self.graceful_shutdown = graceful

    @property
    def name(self) -> str:
        return "fake"

    @property
    def help(self) -> str:
        return "A fake command"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument("--value", type=int, default=0)

    def execute(self, args: argparse.Namespace) -> int:
        return self._action(args)


class TestBaseCommandRun:
    """BaseCommand.run wraps execute in graceful_termination."""

    def test_return_value_passes_through_and_termination_is_on_args(self):
        seen = {}

        def action(args):
            seen["term"] = args._termination
            return 7

        args = argparse.Namespace()
        assert _FakeCommand(action).run(args) == 7
        assert isinstance(seen["term"], Termination)

    def test_keyboard_interrupt_returns_130_logs_and_stops_children(self, caplog):
        def action(args):
            raise KeyboardInterrupt

        args = argparse.Namespace()
        with patch("prepare.cli.base.SecureSubprocess.terminate_children", return_value=0) as terminate:
            with caplog.at_level(logging.ERROR):
                assert _FakeCommand(action).run(args) == 130
        terminate.assert_called_once_with(stop=args._termination.stop)
        assert any(r.getMessage() == "Interrupted (Ctrl-C)" for r in caplog.records)

    def test_a_signal_names_itself_in_the_log(self, caplog):
        def action(args):
            os.kill(os.getpid(), signal.SIGTERM)
            return 0

        with patch("prepare.cli.base.SecureSubprocess.terminate_children", return_value=0):
            with caplog.at_level(logging.ERROR):
                assert _FakeCommand(action).run(argparse.Namespace()) == 130
        assert any(r.getMessage() == "Interrupted (SIGTERM)" for r in caplog.records)

    def test_an_interrupt_while_handlers_are_installed_returns_130_without_a_traceback(self, caplog):
        from contextlib import contextmanager

        @contextmanager
        def interrupted_install():
            raise KeyboardInterrupt
            yield  # pragma: no cover - never reached

        with patch("prepare.cli.base.graceful_termination", interrupted_install):
            with caplog.at_level(logging.ERROR):
                assert _FakeCommand(lambda args: 0).run(argparse.Namespace()) == 130
        assert any(r.getMessage() == "Interrupted" for r in caplog.records)

    def test_graceful_shutdown_false_installs_nothing(self):
        with patch.object(term_mod.signal, "signal") as spy:
            assert _FakeCommand(lambda args: 3, graceful=False).run(argparse.Namespace()) == 3
        assert spy.call_count == 0

    def test_graceful_shutdown_defaults_to_true(self):
        assert _FakeCommand(lambda args: 0).graceful_shutdown is True
        assert BaseCommand.graceful_shutdown is True

    def test_parser_binds_run(self):
        registry = CommandRegistry()
        command = _FakeCommand(lambda args: 5)
        registry.register(command)
        parser = argparse.ArgumentParser()
        registry.setup_parsers(parser)
        args = parser.parse_args(["fake", "--value", "2"])
        assert args.func == command.run
        assert args.func(args) == 5


def test_download_sra_alias_delegates_to_graceful_termination():
    import prepare.cli.commands.sra as sra_mod

    with sra_mod._termination_raises_interrupt() as term:
        assert isinstance(term, Termination)
        with pytest.raises(KeyboardInterrupt):
            os.kill(os.getpid(), signal.SIGTERM)


def test_terminate_children_does_not_deadlock_when_its_own_thread_holds_the_table_lock():
    """The third-signal handler runs on the main thread, possibly inside run_secure's locked section."""
    import time
    from unittest.mock import MagicMock

    from prepare.utils.security import SecureSubprocess

    child = MagicMock()
    child.poll.return_value = None
    token = threading.Event()
    SecureSubprocess._children[child] = token
    try:
        with SecureSubprocess._children_lock:
            started = time.monotonic()
            stopped = SecureSubprocess.terminate_children(grace=0.1, stop=token)
            elapsed = time.monotonic() - started
    finally:
        SecureSubprocess._children.pop(child, None)
    assert stopped == 1
    child.terminate.assert_called_once()
    assert token.is_set()
    assert elapsed < 2.0
