#!/usr/bin/env python
"""pathpatch.py -- the reversible lab-path substitution used for the shipped copies (see PATCHES.md).

The lab's files name lab folders in a few places (sys.path inserts, data file names, usage lines, recorded logs and
ledger rows). The shipped copies replace each such prefix with a placeholder from one fixed table
(tools/path_patch.json, applied in order; the original strings are stored there as hex so that the table
itself carries no path text). The substitution is exactly reversible: restore() maps every
placeholder back, and the restored bytes have the original sha256 (the ledger pin for sealed files). Nothing else in
any file is changed.

Library use: table(), patch_bytes(data), restore_bytes(data), patched_files(), restore_file(path).
Command line:  python tools/pathpatch.py restore DIR     restore every file under DIR that carries a placeholder
               (in place); used by tools/reproduce.py on its work folder.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load():
    return json.load(open(os.path.join(ROOT, 'tools', 'path_patch.json'), encoding='utf-8'))


def table():
    """[(original bytes, placeholder bytes)], in the order they are applied (restore runs them in reverse)."""
    return [(bytes.fromhex(e["original_hex"]), e["placeholder"].encode("utf-8")) for e in load()["table"]
            if not e.get("redacted")]  # redacted entries (privacy, 2026-11-17) are never restored


def patch_bytes(data):
    for o, p in table():
        data = data.replace(o, p)
    return data


def restore_bytes(data):
    for o, p in reversed(table()):
        data = data.replace(p, o)
    return data


def patched_files():
    return load()['files']


def restore_file(path):
    with open(path, 'rb') as f:
        data = f.read()
    with open(path, 'wb') as f:
        f.write(restore_bytes(data))


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] != 'restore':
        sys.exit(__doc__)
    n = 0
    for d, _, fs in os.walk(sys.argv[2]):
        for f in fs:
            p = os.path.join(d, f)
            with open(p, 'rb') as fh:
                data = fh.read()
            new = restore_bytes(data)
            if new != data:
                with open(p, 'wb') as fh:
                    fh.write(new)
                n += 1
    print('pathpatch: restored %d file(s) under %s' % (n, sys.argv[2]))


def redacted_ids():
    """Ids of table entries whose original was removed for privacy (2026-11-17)."""
    return {e["id"] for e in load()["table"] if e.get("redacted")}
