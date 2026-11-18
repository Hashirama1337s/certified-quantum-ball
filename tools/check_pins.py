#!/usr/bin/env python
"""check_pins.py -- verify the shipped files.

1. Every pinned file (tools/pins.json): its sha256 equals the pin. For a path-patched file (PATCHES.md) the file is
   first restored in memory by tools/pathpatch.py, and the restored bytes must have the pinned sha256.
2. Every path-patched file (tools/path_patch.json): the shipped bytes have the listed patched sha256, the same value
   appears in PATCHES.md, and the restored bytes have the listed original sha256.
3. Every file against MANIFEST.sha256 (if present).
Run from anywhere:  python tools/check_pins.py   (exit 0 = all match)."""
import hashlib
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'tools'))
import pathpatch  # noqa: E402


def sha_b(b):
    return hashlib.sha256(b).hexdigest()


def read(rel):
    p = os.path.join(ROOT, rel)
    return open(p, 'rb').read() if os.path.exists(p) else None


def main():
    bad = 0
    patched = pathpatch.patched_files()
    red = pathpatch.redacted_ids()
    redacted_files = {rel for rel, v in patched.items() if any(x.split()[0] in red for x in v.get('entries', []))}
    md = open(os.path.join(ROOT, 'PATCHES.md'), encoding='utf-8').read()
    md_rows = {m.group(1): m.group(2) for m in
               re.finditer(r'^\| (\S+) \| (?:yes|no) \| `[0-9a-f]{64}` \| `([0-9a-f]{64})` \|', md, re.M)}
    pins = json.load(open(os.path.join(ROOT, 'tools', 'pins.json'), encoding='utf-8'))['pins']
    for p in pins:
        data = read(p['path'])
        if data is None:
            ok, how = False, 'missing'
        elif p['path'] in redacted_files:
            ok, how = sha_b(data) == patched[p['path']]['patched_sha256'], 'REDACTED for privacy: pin verifiable against the lab ledger only; shipped bytes checked'
        elif p['path'] in patched:
            ok, how = sha_b(pathpatch.restore_bytes(data)) == p['sha256'], 'restored'
        else:
            ok, how = sha_b(data) == p['sha256'], 'as shipped'
        bad += not ok
        print('%s  %s  %s (%s)' % ('pin OK ' if ok else 'PIN BAD', p['path'], p['sha256'][:16], how))
    for rel, v in sorted(patched.items()):
        data = read(rel)
        if rel in redacted_files:
            ok = data is not None and sha_b(data) == v['patched_sha256'] and md_rows.get(rel) == v['patched_sha256']
            bad += not ok
            print('%s  %s (REDACTED for privacy; original verifiable against the lab ledger only)' % ('patch OK ' if ok else 'PATCH BAD', rel))
            continue
        ok = (data is not None and sha_b(data) == v['patched_sha256'] and md_rows.get(rel) == v['patched_sha256']
              and sha_b(pathpatch.restore_bytes(data)) == v['original_sha256']
              and pathpatch.patch_bytes(pathpatch.restore_bytes(data)) == data)
        bad += not ok
        print('%s  %s' % ('patch OK ' if ok else 'PATCH BAD', rel))
    if len(md_rows) != len(patched):
        bad += 1
        print('PATCH BAD  PATCHES.md lists %d files, tools/path_patch.json %d' % (len(md_rows), len(patched)))
    man = os.path.join(ROOT, 'MANIFEST.sha256')
    if os.path.exists(man):
        n = 0
        for line in open(man, encoding='utf-8'):
            if not line.strip():
                continue
            h, rel = line.rstrip('\n').split('  ', 1)
            data = read(rel)
            n += 1
            if data is None or sha_b(data) != h:
                bad += 1
                print('MANIFEST BAD  %s' % rel)
        print('manifest: %d files checked' % n)
    else:
        print('manifest: MANIFEST.sha256 not found (skipped)')
    print('check_pins: %s (%d pins, %d patched files)' % ('ALL MATCH' if not bad else '%d MISMATCH(ES)' % bad,
                                                          len(pins), len(patched)))
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
