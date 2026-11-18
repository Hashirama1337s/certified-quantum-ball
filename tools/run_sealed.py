#!/usr/bin/env python
"""run_sealed.py -- run a sealed instrument byte-for-byte unchanged, with the lab folder it names pointed at a work folder.

The instruments were written inside one lab folder and name it in string literals: sys.path inserts
(`sys.path.insert(0, '<LAB>/engine/cert')`), input and output files (`'<LAB>/engine/cert/lin_fw4_bg3.pkl'`), a
default cache location, a runpy target, and one recorded config file (b5/cfg_ec.json names its input by its lab
path). Editing those literals would change the files' sha256, which the lab ledger pinned before the certified runs
and which some instruments record or check themselves (count_b4v3_seal.py writes the sha256 of every engine file into
its ledger header; fw_tail2.py logs its own sha256; close_b5.py keys its cache on the sha256 of its producers). So the
files on disk are never edited. Instead this launcher, for every Python file loaded from the work folder (the main
script, every module it imports, and every runpy.run_path target), reads the source bytes, replaces ONLY the lab-root
string (entry R3 of tools/path_patch.json, stored there as hex) with the work folder, and compiles the result under the
file's own name. Everything that reads a file itself (sha256 of a module, a recorded engine hash) still reads the
untouched sealed bytes. Data files are not rewritten either: an open() of a path that starts with the lab root (any of
the root forms R1-R4, any case or slash direction) is redirected to the same relative path under the work folder.
The lab's tool folder (entry R6; the BG2 gate imports the shared count engine gravity_engine from it) is treated the
same way, pointed at WORK/_pocketlab.

Every substitution and redirection is printed on stderr, one line each (`[run_sealed] ...`).
Guard: an audit hook refuses (PermissionError) any file open that still resolves under the lab root after these
rewrites, so a missed literal fails loudly instead of reading a file outside the work folder. On a machine without
the lab folder this cannot happen anyway; on the lab's own machine it proves that the run used only the work folder.

Usage:  python tools/run_sealed.py --root WORK [--trace FILE] SCRIPT.py [script arguments ...]
        SCRIPT.py is a path inside WORK. --trace appends every file opened under WORK (read or write) to FILE.
"""
import builtins
import importlib.machinery
import io
import json
import os
import runpy
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _norm(p):
    return os.path.normcase(os.path.abspath(p)).replace('\\', '/')


def main(argv):
    trace = None
    if len(argv) >= 2 and argv[0] == '--root':
        work = os.path.abspath(argv[1]).replace('\\', '/')
        argv = argv[2:]
    else:
        sys.exit(__doc__)
    if len(argv) >= 2 and argv[0] == '--trace':
        trace = argv[1]
        argv = argv[2:]
    if not argv:
        sys.exit(__doc__)
    tab = json.load(open(os.path.join(HERE, 'path_patch.json'), encoding='utf-8'))
    ent = {e['id']: bytes.fromhex(e['original_hex']).decode('utf-8') for e in tab['table']}
    root = ent['R3']                                   # the literal as it appears in the sources
    root_b = root.encode('utf-8')
    work_b = work.encode('utf-8')
    forms = sorted({_norm(ent[i].replace('\\\\', '\\')) for i in tab['root_ids']} | {_norm(root)}, key=len, reverse=True)
    lab_n = _norm(root)
    tools = ent['R6']                                  # the lab's tool folder (gravity_engine)
    tools_b = tools.encode('utf-8')
    wtools = work + '/_pocketlab'
    tforms = sorted({_norm(ent[i].replace('\\\\', '\\')) for i in tab.get('tools_ids', [])} | {_norm(tools)},
                    key=len, reverse=True)
    work_n = _norm(work)
    if work_n == lab_n:
        sys.exit('run_sealed: the work folder must not be the lab folder itself')

    def say(msg):
        print('[run_sealed] ' + msg, file=sys.stderr, flush=True)

    def rel(p):
        n = _norm(p)
        return os.path.abspath(p).replace(os.sep, '/')[len(work_n) + 1:] if n.startswith(work_n + '/') else p

    def rewrite(data, path):
        if not _norm(path).startswith(work_n + '/'):
            return data
        n = data.count(root_b)
        if n:
            say('%s: %d lab-root literal(s) -> work folder (file on disk unchanged)' % (rel(path), n))
            data = data.replace(root_b, work_b)
        n = data.count(tools_b)
        if n:
            say('%s: %d lab-tool-folder literal(s) -> work folder/_pocketlab (file on disk unchanged)' % (rel(path), n))
            data = data.replace(tools_b, wtools.encode('utf-8'))
        return data

    # 1. import hook: every .py loaded from the work folder is compiled from rewritten bytes (no bytecode cache)
    class Loader(importlib.machinery.SourceFileLoader):
        def get_code(self, fullname):
            path = self.get_filename(fullname)
            with io.open_code(path) as f:
                data = f.read()
            return compile(rewrite(data, path), path, 'exec', dont_inherit=True)

    sys.dont_write_bytecode = True
    hook = importlib.machinery.FileFinder.path_hook(
        (importlib.machinery.ExtensionFileLoader, importlib.machinery.EXTENSION_SUFFIXES),
        (Loader, importlib.machinery.SOURCE_SUFFIXES),
        (importlib.machinery.SourcelessFileLoader, importlib.machinery.BYTECODE_SUFFIXES))
    sys.path_hooks.insert(0, hook)
    sys.path_importer_cache.clear()

    # 2. runpy.run_path targets (count_b4v3.py runs lin_b4_smoke3.py this way)
    _orig_get = runpy._get_code_from_file

    def _get_code_from_file(*a):
        fname = a[-1]
        with io.open_code(os.path.abspath(fname)) as f:
            data = f.read()
        code = compile(rewrite(data, fname), fname, 'exec')
        return (code, fname) if len(a) == 2 else code
    runpy._get_code_from_file = _get_code_from_file

    # 3. data paths that start with the lab root (e.g. a path recorded inside a config file)
    seen = set()

    def redirect(p):
        if isinstance(p, bytes) or not isinstance(p, (str, os.PathLike)):
            return p
        s = os.fspath(p)
        if not isinstance(s, str):
            return p
        n = _norm(s)
        if n == work_n or n.startswith(work_n + '/'):
            return p
        for fms, base in ((forms, work), (tforms, wtools)):
            for fm in fms:
                if n == fm or n.startswith(fm + '/'):
                    new = base + n[len(fm):]
                    if s not in seen:
                        seen.add(s)
                        say('open: lab path %s -> work folder %s' % (n[len(fm):] or '/', rel(new)))
                    return new
        return p
    _open = builtins.open

    def open_(file, *a, **k):
        return _open(redirect(file), *a, **k)
    builtins.open = open_
    io.open = open_

    # 4. guard + trace
    tf = _open(trace, 'a', encoding='utf-8') if trace else None

    def audit(event, args):
        if event == 'open' and args and isinstance(args[0], (str, bytes, os.PathLike)):
            try:
                s = os.fsdecode(args[0])
            except Exception:
                return
            n = _norm(s)
            for fm in forms + tforms:
                if n == fm or n.startswith(fm + '/'):
                    if not (n == work_n or n.startswith(work_n + '/')):
                        raise PermissionError('run_sealed guard: refused to open %s (outside the work folder)' % s)
            if tf is not None and n.startswith(work_n + '/'):
                mode = args[1] if len(args) > 1 else 'r'
                tf.write('%s\t%s\n' % (mode, rel(s)))
                tf.flush()
    sys.addaudithook(audit)

    path = os.path.abspath(argv[0])
    if not _norm(path).startswith(work_n + '/'):
        sys.exit('run_sealed: %s is not inside the work folder %s' % (path, work))
    with io.open_code(path) as f:
        src = f.read()
    sys.argv = [path] + argv[1:]
    sys.path.insert(0, os.path.dirname(path))
    code = compile(rewrite(src, path), path, 'exec', dont_inherit=True)
    g = {'__name__': '__main__', '__file__': path, '__builtins__': builtins, '__spec__': None}
    exec(code, g)


if __name__ == '__main__':
    main(sys.argv[1:])
