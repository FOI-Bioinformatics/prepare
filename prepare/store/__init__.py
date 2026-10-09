"""
Shared data store package.

The store is one folder holding a single copy of every downloaded
metagenome, used by several organism projects. This package resolves where
that folder lives and defines its on-disk layout. Nothing outside this
package uses it yet.
"""

from prepare.store.layout import (
    StorePaths,
    init_store,
    lock_path,
    read_marker,
    sidecar_path,
    sra_dir,
    store_paths,
)
from prepare.store.resolve import (
    config_path,
    read_config,
    resolve_store_root,
    write_config_data_root,
)
from prepare.store.sidecar import (
    SIDECAR_SCHEMA,
    Sidecar,
    build_sidecar,
    ncbi_from_metadata_xml,
    read_sidecar,
    write_sidecar,
)
from prepare.store.catalog import Catalog, catalog_write
from prepare.store.link import (
    LINK_MODES,
    dangling_links,
    is_store_link,
    link_dataset,
    unlink_dataset,
)
from prepare.store.adopt import AdoptReport, adopt

__all__ = [
    "StorePaths",
    "init_store",
    "lock_path",
    "read_marker",
    "sidecar_path",
    "sra_dir",
    "store_paths",
    "config_path",
    "read_config",
    "resolve_store_root",
    "write_config_data_root",
    "SIDECAR_SCHEMA",
    "Sidecar",
    "build_sidecar",
    "ncbi_from_metadata_xml",
    "read_sidecar",
    "write_sidecar",
    "Catalog",
    "catalog_write",
    "LINK_MODES",
    "dangling_links",
    "is_store_link",
    "link_dataset",
    "unlink_dataset",
    "AdoptReport",
    "adopt",
]
