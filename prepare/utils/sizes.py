"""Disk-size helpers shared by the store, cleanup and link code."""

import os
import stat
from pathlib import Path
from typing import Union


def tree_bytes(path: Union[str, Path]) -> int:
    """Total bytes of the regular files at or under ``path``.

    A regular file counts its own size and a directory counts the sum of the files beneath it.
    Symbolic links are never followed, so a link that points outside the tree adds nothing and a
    link to a directory is not entered. An entry that cannot be examined, or that disappears
    during the walk, counts as 0, and a missing ``path`` returns 0.
    """
    try:
        info = os.lstat(path)
    except OSError:
        return 0
    if stat.S_ISREG(info.st_mode):
        return info.st_size
    if not stat.S_ISDIR(info.st_mode):
        return 0
    total = 0
    pending = [os.fspath(path)]
    while pending:
        current = pending.pop()
        try:
            with os.scandir(current) as entries:
                for entry in entries:
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            pending.append(entry.path)
                        elif entry.is_file(follow_symlinks=False):
                            total += entry.stat(follow_symlinks=False).st_size
                    except OSError:
                        continue
        except OSError:
            continue
    return total
