"""
Request pacing for the NCBI E-utilities clients in PREPARE.

NCBI allows roughly 3 requests per second without an API key and 10 per second with one.
The taxonomy client, the SRA metadata client and the Entrez metadata download each hold a
RequestPacer and call ``wait()`` before every request, so all three follow the same rule.
"""

import threading
import time
from typing import Optional


class RequestPacer:
    """Keep successive requests at least ``delay`` seconds apart.

    The pacer is thread-safe: the lock is held across the sleep, so two threads racing on a
    request each wait out a full gap instead of both sleeping the same remainder. It uses the
    monotonic clock, so a wall-clock adjustment cannot shorten or lengthen a gap.
    """

    def __init__(self, api_key: Optional[str], per_second_without_key: int = 3, per_second_with_key: int = 10):
        """
        Initialize the pacer.

        Args:
            api_key: NCBI API key, or None; only its presence matters, for the request rate
            per_second_without_key: Requests per second allowed without a key
            per_second_with_key: Requests per second allowed with a key
        """
        self._delay = 1.0 / (per_second_with_key if api_key else per_second_without_key)
        self._last_request = float("-inf")
        self._lock = threading.Lock()

    @property
    def delay(self) -> float:
        """Minimum number of seconds between two successive requests."""
        return self._delay

    def wait(self) -> None:
        """Sleep only as long as needed since the previous request, then record this one."""
        with self._lock:
            elapsed = time.monotonic() - self._last_request
            if elapsed < self._delay:
                time.sleep(self._delay - elapsed)
            self._last_request = time.monotonic()
