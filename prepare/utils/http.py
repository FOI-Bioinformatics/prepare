"""One retrying HTTP session builder, shared by the NCBI taxonomy, GTDB and SRA metadata clients.

(The Branchwater client and the Entrez calls of ``data/metadata.py`` handle their own retries.)
NCBI and GTDB are flaky in the same way: an occasional connection failure, a 429
(rate limited) or a 5xx while the service recovers. ``retrying_session`` mounts one urllib3
``Retry`` policy on both ``http://`` and ``https://`` so a caller's ``session.get(...)`` retries
those failures with backoff automatically, without each client re-implementing the same few lines
of ``requests.adapters``/``urllib3.util.retry`` wiring. ``raise_on_status=False`` means a session
built this way never raises urllib3's own ``MaxRetryError`` when every retry is exhausted; the
last response is returned instead, so the caller's own ``response.raise_for_status()`` reports the
plain HTTP error it would have reported without retries at all.
"""

from typing import Sequence

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Defaults used by every current caller (NCBI taxonomy, GTDB); a caller with different needs
# passes its own values rather than changing these.
RETRY_TOTAL = 3
RETRY_BACKOFF_FACTOR = 0.5
RETRY_STATUS_FORCELIST = (429, 500, 502, 503, 504)
RETRY_ALLOWED_METHODS = ("GET",)


def retrying_session(
    total: int = RETRY_TOTAL,
    backoff_factor: float = RETRY_BACKOFF_FACTOR,
    status_forcelist: Sequence[int] = RETRY_STATUS_FORCELIST,
    allowed_methods: Sequence[str] = RETRY_ALLOWED_METHODS,
) -> requests.Session:
    """Build a ``requests.Session`` that retries a connection failure or a 429/5xx response.

    Every retry sleeps ``backoff_factor * (2 ** (attempt - 1))`` seconds before the next
    attempt, up to ``total`` retries. When every attempt fails, the *last* response is returned
    to the caller rather than raised as a ``urllib3.exceptions.MaxRetryError``, so
    ``response.raise_for_status()`` reports the plain HTTP error.
    """
    session = requests.Session()
    retry = Retry(
        total=total,
        backoff_factor=backoff_factor,
        status_forcelist=list(status_forcelist),
        allowed_methods=list(allowed_methods),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session
