"""
Shared data store CLI commands.

`store_init` creates (or reuses) a shared data store folder and records it, along with this
project's identity, in the project registry. `store_status` and `store_reindex` operate
against whichever store root resolves for the current project (an explicit `--data-root`, the
`PREPARE_DATA` environment variable, the registry's recorded `store.root`, or the user's
default config), via `prepare.store.resolve.resolve_store_root`.

Every command here exists to operate on the store, so an unreachable one is an error with
exit 1, not something to work around; the analysis and reporting commands degrade instead
(see `prepare.store.resolve.resolve_optional_store`). A project that links from the store
without ever running `store_init` has its identity minted on the spot
(`prepare.store.usage.ensure_project_identity`), so `store_gc` can always see who uses what.

Each command lives in its own module; helpers used by more than one live in `_shared`.
"""

from prepare.cli.commands.store.adopt import StoreAdoptCommand
from prepare.cli.commands.store.gc import StoreGcCommand
from prepare.cli.commands.store.init import StoreInitCommand
from prepare.cli.commands.store.link import StoreLinkCommand, StoreUnlinkCommand
from prepare.cli.commands.store.reindex import StoreReindexCommand
from prepare.cli.commands.store.status import StoreStatusCommand
from prepare.cli.commands.store.usage import StoreUsageCommand
from prepare.cli.commands.store.verify import StoreVerifyCommand

__all__ = [
    "StoreAdoptCommand",
    "StoreGcCommand",
    "StoreInitCommand",
    "StoreLinkCommand",
    "StoreReindexCommand",
    "StoreStatusCommand",
    "StoreUnlinkCommand",
    "StoreUsageCommand",
    "StoreVerifyCommand",
]
