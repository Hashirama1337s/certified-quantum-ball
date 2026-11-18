#!/usr/bin/env python
"""make_manifest.py -- write MANIFEST.sha256 (sha256 of every file in the repository except the manifest itself,
version-control metadata and caches). Re-run after any file changes; tools/check_pins.py verifies it.
Run from anywhere:  python tools/make_manifest.py"""
import hashlib
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIP_DIRS = {'.git', '__pycache__', 'work'}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main():
    rows = []
    for d, dirs, files in os.walk(ROOT):
        dirs[:] = sorted(x for x in dirs if x not in SKIP_DIRS)
        for f in files:
            p = os.path.join(d, f)
            rel = os.path.relpath(p, ROOT).replace('\\', '/')
            if rel == 'MANIFEST.sha256':
                continue
            rows.append((rel, sha(p)))
    rows.sort()
    with open(os.path.join(ROOT, 'MANIFEST.sha256'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(''.join('%s  %s\n' % (h, rel) for rel, h in rows))
    print('MANIFEST.sha256: %d files' % len(rows))


if __name__ == '__main__':
    main()
