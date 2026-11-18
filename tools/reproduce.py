#!/usr/bin/env python
"""reproduce.py -- re-run the sealed instruments of this repository in a clean work folder and compare with the
shipped outputs.

    python tools/reproduce.py --tier quick     pins + the short checks (symbolic derivations, QB1 main run, QB1b
                                               tasks and one sign-map region, centre closed form, the QB2c replay of
                                               the matching / readout and the replication row)
    python tools/reproduce.py --tier full      everything, every data file regenerated (see tools/reproduce.md)
    python tools/reproduce.py --tier bg2       the BG2 background record (sonic box + outward leg + gate), regenerated
    python tools/reproduce.py --only NAME ...  named steps (the data files they need are prepared first)
Options: --work DIR (default ./work; must be absent or empty); --python "CMD" (default: this interpreter);
         --regen (regenerate data files instead of using the shipped, sha256-checked copies); --list.

What it does
  1. tools/check_pins.py (every pinned file; every path-patched copy; MANIFEST.sha256).
  2. WORK = a mirror of the lab tree: code/ is copied to WORK and every path-patched copy is restored to its original
     (sealed) bytes. The lab's tool folder (code/_pocketlab) sits at WORK/_pocketlab.
  3. Data files (pickles). Loading a pickle can execute code. The quick tier copies a shipped pickle into WORK only
     after its sha256 matches tools/layout.json (and the ledger pin where one exists); pickles that are not shipped
     are regenerated. The full tier regenerates every pickle and checks it byte for byte.
  4. Every step runs as `python tools/run_sealed.py --root WORK <instrument> <args>` from the folder the lab ran it
     in, with output files the step writes deleted first. Logs go to WORK/_out/<step>.log.
  5. Outputs are compared with the shipped copies by canonical content (tools/canon.py: elapsed times, peak memory,
     the project folder and run-time keys removed; every other character must agree).
Last line: [reproduce] ALL CHECKS PASS (or the failures). JSON report: WORK/_out/reproduce_report.json.
"""
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOLS = os.path.join(ROOT, 'tools')
sys.path.insert(0, TOOLS)
import pathpatch  # noqa: E402
import canon  # noqa: E402

Q = 'quantum/seats/'
B1, B2, B3, B4, B5 = (Q + 'build-%d' % i for i in range(1, 6))
P2 = 'sideways/build-5/phase2'

# data files: key -> (lab path, cwd, argv, lab minutes, needs)
GENS = {
    'eqs': ('engine/eqs_gen.py', '', ['engine/derive.py'], 1, []),
    'fields': ('engine/cert/fields.pkl', '', ['engine/cert/fields.py'], 1, ['eqs']),
    'lin_fields': ('engine/cert/lin_fields.pkl', '', ['engine/cert/lin_sym.py'], 2, ['eqs']),
    'lin_fields2': ('engine/cert/lin_fields2.pkl', '', ['engine/cert/lin_mc.py'], 1, ['lin_fields']),
    'tube13': ('engine/cert/lin_tube2_k1_3.pkl', '', ['engine/cert/lin_tube_v2.py', '1/3'], 9, ['fields']),
    'cacheL': ('b4speed/L13_nsub2.pkl', '', ['b4speed/cacheL.py', 'build', '1/3', '2', '1e-20', 'b4speed/L13_nsub2.pkl'],
               8, ['fields', 'lin_fields2', 'tube13']),
    'eqs_l': ('sideways/build-1/eqs_l.pkl', 'sideways/build-1', ['sideways/build-1/eqs_l.py'], 1.2, ['fields']),
    'outward_rec': ('sideways/build-1/outward_rec.pkl', 'sideways/build-1', ['sideways/build-1/outward_rec.py'], 6,
                    ['fields', 'eqs_l']),
    'outward_C': ('sideways/build-5/phase2/outward_C.pkl', P2, None, 64, ['fields', 'outward_rec']),  # tier bg2
    'bg2_tube': (B4 + '/lin_tube_bg2_k1_3.pkl', B4, [B4 + '/tube_bg2.py'], 8.6, ['fields', 'outward_C']),
    'Lbg2': (B4 + '/Lbg2_nsub2.pkl', B4, None, 6, ['lin_fields2', 'bg2_tube']),        # written by step qb2b_stages
    'qb2_exprs': (B3 + '/qb2_exprs.pkl', B3, [B3 + '/qb2sym2.py'], 1, []),
    'null_in': (B1 + '/null_in.pkl', B1, [B1 + '/null_bg.py', 'in'], 8.5, ['fields']),
    'null_out': (B1 + '/null_out.pkl', B1, [B1 + '/null_bg.py', 'out'], 13.3, ['fields']),
}

# steps: (name, tier, cwd, argv, env, compares [(produced lab path, kind)], lab minutes, data needed)
#   kind: 'json' / 'jsonl' / 'text' / 'log:<lab log>' (stdout compared with that shipped log)
S = []


def step(name, tier, cwd, argv, outs, minutes, needs, env=None):
    S.append((name, tier, cwd, argv, env or {}, outs, minutes, needs))


# QB1 (54caa4012c), build-1
step('qb1_controls', 'full', B1, [B1 + '/qb1.py', 'controls'], [], 0.3, ['fields', 'tube13', 'outward_C'])
step('qb1_main', 'quick', B1, [B1 + '/qb1.py', 'main'], [(B1 + '/result_qb1.json', 'json'), ('', 'log:' + B1 + '/main.log')],
     2.3, ['fields', 'tube13', 'outward_C'])
step('qb1_dnew', 'full', B1, [B1 + '/dnew.py'], [(B1 + '/dnew_qb1.json', 'json'), ('', 'log:' + B1 + '/dnew.log')], 0.5,
     ['fields', 'lin_fields', 'lin_fields2', 'tube13', 'cacheL'])
step('qb1_null_eval', 'full', B1, [B1 + '/null_eval.py'], [(B1 + '/null_qb1.json', 'json')], 0.1,
     ['fields', 'tube13', 'outward_C', 'null_in', 'null_out'])
step('qb1_make_controls', 'full', B1, [B1 + '/make_controls.py'], [(B1 + '/controls_qb1.json', 'json')], 0.1,
     ['fields', 'tube13', 'outward_C'])
# QB1b (4d3b12ef39), build-2
REG = [('w', '-12', '-11/2', '13', '5', 'sm_w.json', 1.9), ('centre', '-11/2', '-9/2', '16', '5', 'sm_centre.json', 10),
       ('tube', '-9/2', '-1/16', '284', '5', 'sm_tube.json', 40), ('tube', '-1/16', '-13/256', '3', '5', 'sm_tube2.json', 5),
       ('sonic', '-1/16', '13/256', '29', '5', 'sm_sonic.json', 10), ('outw', '13/256', '3', '151', '5', 'sm_out1.json', 20),
       ('outw', '3', '9', '384', '7', 'sm_out2.json', 20), ('outw', '9', '21/2', '384', '7', 'sm_out3.json', 10),
       ('outw', '21/2', '12', '768', '7', 'sm_out4.json', 20)]
for r in REG:
    tag = r[5][:-5]
    step('qb1b_' + tag, 'quick' if tag == 'sm_w' else 'full', B2, [B2 + '/signmap2.py'] + list(r[:6]),
         [(B2 + '/' + r[5], 'json'), ('', 'log:%s/%s.log' % (B2, tag))], r[6], ['fields', 'tube13', 'outward_C'])
for t, out, lg in (('seam', 'pq8.json', 't2_seam.log'), ('null', 'pq6b.json', 't2_null.log'),
                   ('controls', 'controls_tm.json', 't2_controls.log')):
    step('qb1b_' + t, 'quick', B2, [B2 + '/tasks2.py', t], [(B2 + '/' + out, 'json'), ('', 'log:%s/%s' % (B2, lg))], 0.2,
         ['fields', 'tube13', 'outward_C'])
step('qb1b_dnew2', 'full', B2, [B2 + '/dnew2.py'], [(B2 + '/pq7b.json', 'json')], 0.5,
     ['fields', 'lin_fields', 'lin_fields2', 'tube13', 'cacheL'])
step('qb1b_finalize', 'quick', B2, [B2 + '/finalize.py'], [(B2 + '/result_qb1b.json', 'json')], 1,
     ['fields', 'tube13', 'outward_C'])
# QB2 (f2bcfadb3e), build-3
step('qb2_bachsym', 'quick', B3, [B3 + '/bachsym.py'], [(B3 + '/bach_exprs.txt', 'text'), ('', 'log:' + B3 + '/bachsym.log')],
     1.5, [])
step('qb2_bach_u', 'quick', B3, [B3 + '/bach_u.py'], [(B3 + '/bach_u.txt', 'text'), ('', 'log:' + B3 + '/bach_u.log')], 0.5, [])
step('qb2_n1_rows', 'quick', B3, [B3 + '/n1_rows.py'], [(B3 + '/n1_rows.txt', 'text'), ('', 'log:' + B3 + '/n1_rows.log')],
     0.4, [])
step('qb2_n0check', 'quick', B3, [B3 + '/n0check.py'], [(B3 + '/n0check.json', 'json')], 1, ['fields', 'tube13', 'qb2_exprs'])
step('qb2_stages', 'full', B3, [B3 + '/qb2stages.py'], [(B3 + '/stages.json', 'json')], 1,
     ['fields', 'lin_fields2', 'tube13', 'cacheL', 'qb2_exprs'])
step('qb2_pass', 'full', B3, ['@chunks', B3 + '/qb2pass.py'], [(B3 + '/pass_ckpt.json', 'json')], 80,
     ['fields', 'lin_fields2', 'tube13', 'cacheL', 'qb2_exprs'])
step('qb2_c2check', 'full', B3, [B3 + '/c2check.py'], [(B3 + '/c2.json', 'json')], 2,
     ['fields', 'lin_fields2', 'tube13', 'cacheL', 'qb2_exprs'])
step('qb2_final', 'full', B3, [B3 + '/qb2final.py'], [(B3 + '/final.json', 'json')], 3,
     ['fields', 'lin_fields2', 'tube13', 'cacheL', 'qb2_exprs'])
step('qb2_outputs', 'full', B3, [B3 + '/make_outputs.py'], [(B3 + '/result_qb2.json', 'json')], 0.2,
     ['fields', 'lin_fields2', 'tube13', 'cacheL', 'qb2_exprs'])
# QB2b (373f306e9f), build-4
step('qb2b_tube', 'full', B4, [B4 + '/tube_compare.py'], [(B4 + '/tube_compare.json', 'json')], 1, ['tube13', 'bg2_tube'])
step('qb2b_centre', 'quick', B4, [B4 + '/centre_exact.py'], [(B4 + '/centre_exact.json', 'json')], 3.4, [])
step('qb2b_thth_sym', 'quick', B4, [B4 + '/thth_sym.py'], [(B4 + '/thth_exprs.txt', 'text')], 2, [])
step('qb2b_stages', 'full', B4, [B4 + '/stages_b.py'], [(B4 + '/stages.json', 'json'), (B4 + '/Lbg2_nsub2.pkl', 'bytes')],
     6, ['fields', 'lin_fields2', 'tube13', 'bg2_tube', 'qb2_exprs', 'outward_C'], env={'QB2B_RUN1': '1'})
# The QB2b record is run 1 (plain homogeneous stage, cache Lbg2_nsub2.pkl): QB2B_RUN1=1 reproduces that configuration
# (build-4 RUN.md). Without it qb2b.py imports qb2b_tight and builds the tight cache of the stopped run 2 (> 40 min);
# the shipped run_pass_b.sh sets QB2B_NOPATCH=1, which is the run-2 configuration, not the record.
step('qb2b_pass', 'full', B4, ['@chunks', B4 + '/pass_b.py'], [(B4 + '/pass_ckpt.json', 'json')], 85,
     ['fields', 'lin_fields2', 'tube13', 'bg2_tube', 'Lbg2', 'qb2_exprs', 'outward_C'], env={'QB2B_RUN1': '1'})
step('qb2b_final', 'full', B4, [B4 + '/final_b.py'], [(B4 + '/final.json', 'json')], 3,
     ['fields', 'lin_fields2', 'tube13', 'bg2_tube', 'Lbg2', 'qb2_exprs', 'outward_C'], env={'QB2B_RUN1': '1'})
step('qb2b_thth', 'full', B4, [B4 + '/thth_eval.py'], [(B4 + '/thth.json', 'json')], 2,
     ['fields', 'lin_fields2', 'tube13', 'bg2_tube', 'Lbg2', 'qb2_exprs', 'outward_C'], env={'QB2B_RUN1': '1'})
step('qb2b_outputs', 'full', B4, [B4 + '/make_outputs_b.py'], [(B4 + '/result_qb2b.json', 'json')], 0.2,
     ['fields', 'lin_fields2', 'tube13', 'bg2_tube', 'Lbg2', 'qb2_exprs', 'outward_C'])
# QB2c (7aab7efe5a), build-5 (reads the QB2b run-1 record, read only)
D5 = ['fields', 'lin_fields2', 'tube13', 'bg2_tube', 'Lbg2', 'qb2_exprs', 'outward_C']
step('qb2c_final', 'quick', B5, [B5 + '/final_c.py'], [(B5 + '/final_c.json', 'json'), (B5 + '/conv_c.json', 'json')], 2, D5)
step('qb2c_thth', 'quick', B5, [B5 + '/thth_c.py'], [(B5 + '/thth_c.json', 'json')], 4, D5)
step('qb2c_deep', 'full', B5, [B5 + '/deep_ratio.py', '480'], [(B5 + '/deep_ratio.json', 'json')], 20, D5)
step('qb2c_outputs', 'quick', B5, [B5 + '/make_outputs_c.py'], [(B5 + '/result_qb2c.json', 'json')], 0.2, D5)
# BG2 (ed61d09268): sonic box, outward record, gate (tier bg2; regenerates outward_C.pkl)
step('bg2_box', 'bg2', P2, [P2 + '/b1_box.py', 'C', '256', '50', '24', '--nonewton'], [(P2 + '/b1_C.json', 'json')], 49,
     ['fields'], env={'LB_NGRID': '24'})
step('bg2_outward', 'bg2', P2, [P2 + '/b2_outward.py', 'C', 'b1_C.json', '256', '50'],
     [(P2 + '/b2_C.json', 'json'), (P2 + '/outward_C.pkl', 'bytes')], 15, ['fields', 'outward_rec'])
step('bg2_gate', 'bg2', P2, [P2 + '/b3_gate.py', 'C', 'outward_C.pkl'], [(P2 + '/b3_C.json', 'json')], 4.3,
     ['fields', 'tube13', 'eqs_l', 'outward_rec'])

GEN_STEP = {'Lbg2': 'qb2b_stages', 'outward_C': 'bg2_outward'}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load_json(rel):
    return json.load(open(os.path.join(ROOT, rel), encoding='utf-8'))


def roots():
    t = {e['id']: bytes.fromhex(e['original_hex']).decode() for e in load_json('tools/path_patch.json')['table']
         if e['original_hex']}
    return [t['R3']]


def setup(work):
    os.makedirs(work, exist_ok=True)
    if os.listdir(work):
        sys.exit('[reproduce] work folder %s is not empty' % work)
    shutil.copytree(os.path.join(ROOT, 'code'), work, dirs_exist_ok=True, ignore=shutil.ignore_patterns('__pycache__'))
    n = 0
    for rel, v in pathpatch.patched_files().items():
        if rel.startswith('code/'):
            p = os.path.join(work, rel[5:])
            pathpatch.restore_file(p)
            if sha(p) != v['original_sha256']:
                sys.exit('[reproduce] restore failed for %s' % rel)
            n += 1
    os.makedirs(os.path.join(work, '_out'), exist_ok=True)
    return n


def expected_sha():
    """lab path of each data file -> expected sha256 (tools/layout.json: shipped and regenerated-only data files)"""
    return {g['lab']: g['sha256'] for g in load_json('tools/layout.json')['pickles']}


def run(work, py, name, cwd, argv, env_extra):
    env = {k: v for k, v in os.environ.items() if not (k.startswith('LB_') or k.startswith('FW') or k.startswith('QB2B'))}
    env.update(env_extra)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['PYTHONIOENCODING'] = 'utf-8'
    raw = os.path.join(work, '_out', name + '.raw')
    t0 = time.time()
    chunks = argv[0] == '@chunks'
    if chunks:          # the lab's run_pass*.sh: 240-substep chunks, exit 3 = resume, exit 0 = done
        argv = argv[1:]
    k, rc = 0, None
    with open(raw, 'wb') as fo:
        while True:
            a = argv + ([str(k), str(k + 240)] if chunks else [])
            cmd = py + ['-u', os.path.join(TOOLS, 'run_sealed.py'), '--root', work, os.path.join(work, a[0])] + a[1:]
            rc = subprocess.run(cmd, cwd=os.path.join(work, cwd), env=env, stdout=fo, stderr=subprocess.STDOUT).returncode
            if not chunks or rc != 3:
                break
            k += 240
    lines = open(raw, 'rb').read().decode('utf-8', errors='replace').splitlines(True)
    with open(os.path.join(work, '_out', name + '.log'), 'w', encoding='utf-8', newline='') as f:
        f.writelines(l for l in lines if not l.startswith('[run_sealed]'))
    with open(os.path.join(work, '_out', name + '.launcher.log'), 'w', encoding='utf-8', newline='') as f:
        f.writelines(l for l in lines if l.startswith('[run_sealed]'))
    os.remove(raw)
    return rc, (time.time() - t0) / 60


def compare(work, name, outs):
    res = []
    for prod, kind in outs:
        if kind.startswith('log:'):
            ship = 'code/' + kind[4:]
            got = open(os.path.join(work, '_out', name + '.log'), 'rb').read()
            kind = 'log'
        else:
            ship = 'code/' + prod
            pp = os.path.join(work, prod)
            if not os.path.exists(pp):
                res.append(dict(produced=prod, status='MISSING'))
                continue
            got = open(pp, 'rb').read()
        if kind == 'bytes':
            want = expected_sha().get(prod)
            st = 'IDENTICAL (bytes)' if hashlib.sha256(got).hexdigest() == want else 'DIFFERENT (bytes)'
            res.append(dict(produced=prod, status=st))
            continue
        exp = pathpatch.restore_bytes(open(os.path.join(ROOT, ship), 'rb').read())
        k2 = 'log' if kind == 'text' else kind
        st, hp, he, diffs = canon.compare(k2, got, exp, roots() + [work])
        res.append(dict(produced=prod or '(stdout)', shipped=ship, status=st, sha_produced=hp, sha_expected=he,
                        differences=diffs))
    return res


def main(argv):
    work, py, tier, only, regen, lst = os.path.join(os.getcwd(), 'work'), [sys.executable], 'quick', None, False, False
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--work':
            work = os.path.abspath(argv[i + 1]); i += 2
        elif a == '--python':
            py = shlex.split(argv[i + 1]); i += 2
        elif a == '--tier':
            tier = argv[i + 1]; i += 2
        elif a == '--regen':
            regen = True; i += 1
        elif a == '--list':
            lst = True; i += 1
        elif a == '--only':
            only = []; i += 1
            while i < len(argv) and not argv[i].startswith('--'):
                only.append(argv[i]); i += 1
        else:
            sys.exit('[reproduce] unknown argument %s\n%s' % (a, __doc__))
    if tier not in ('quick', 'full', 'bg2'):
        sys.exit('[reproduce] --tier is quick, full or bg2')
    if lst:
        for kk, g in GENS.items():
            print('gen_%-13s data   %5.1f min  %s' % (kk, g[3], ' '.join(g[2]) if g[2] else 'written by ' + GEN_STEP[kk]))
        for s in S:
            print('%-17s %-6s %5.1f min  %s' % (s[0], s[1], s[6], ' '.join(s[3])))
        return 0
    if only is not None:
        todo = [s for s in S if s[0] in only]
        bad = set(only) - {s[0] for s in S}
        if bad:
            sys.exit('[reproduce] unknown step(s) %s' % ', '.join(sorted(bad)))
    elif tier == 'quick':
        todo = [s for s in S if s[1] == 'quick']
    elif tier == 'full':
        todo = [s for s in S if s[1] in ('quick', 'full')]
    else:
        todo = [s for s in S if s[1] == 'bg2']
    regen = regen or (tier == 'full' and only is None)

    def log(s):
        print('[reproduce] ' + s, flush=True)

    T0 = time.time()
    if subprocess.run(py + [os.path.join(TOOLS, 'check_pins.py')]).returncode != 0:
        log('pins FAILED: a shipped file differs from its pinned sha256; stopping')
        return 1
    log('pins: every pinned file matches')
    canon.selftest()       # run-time keys ignored at any depth; real value changes caught (else AssertionError)
    n = setup(work)
    log('work folder: %s (%d path-patched copies restored to their original bytes)' % (work, n))
    exp = expected_sha()
    shipped = {g['lab'] for g in load_json('tools/layout.json')['pickles'] if g.get('shipped', True)}
    need = []

    def want(kk):
        for d in GENS[kk][4]:
            want(d)
        if kk not in need:
            need.append(kk)
    for s in todo:
        for d in s[7]:
            want(d)
    results = []
    bg2_steps = {'bg2_box', 'bg2_outward'}
    for kk in need:
        lab, cwd, gargv, mins, _ = GENS[kk]
        dst = os.path.join(work, lab)
        if lab.endswith('.py'):           # eqs_gen.py is shipped as code; regenerated only with --regen
            if not regen:
                continue
        if not regen and lab in shipped and kk != 'Lbg2':
            if sha(dst) != exp.get(lab):
                log('data %s: shipped copy sha256 differs from tools/layout.json; not loading it' % lab)
                results.append(dict(name='data_' + kk, status='DIFFERENT'))
                break
            log('data %s: shipped copy sha256 matches' % lab)
            continue
        if gargv is None:                 # produced by a step (Lbg2 by qb2b_stages, outward_C by the bg2 tier)
            if kk == 'outward_C':
                if sha(dst) != exp.get(lab):
                    results.append(dict(name='data_' + kk, status='DIFFERENT'))
                    break
                log('data %s: shipped copy sha256 matches (regenerated only by --tier bg2)' % lab)
                continue
            st = [s for s in S if s[0] == GEN_STEP[kk]][0]
            if st not in todo:
                todo.insert(0, st)
            continue
        if os.path.exists(dst):
            os.remove(dst)
        log('regenerating %s (%s; about %s min in the lab)' % (lab, ' '.join(gargv), mins))
        rc, m = run(work, py, 'gen_' + kk, cwd, gargv, {'LB_NGRID': '24'} if kk == 'outward_rec' else {})
        got = sha(dst) if os.path.exists(dst) else None
        status = 'IDENTICAL (bytes)' if rc == 0 and got == exp.get(lab) else ('RUN FAILED (exit %d)' % rc if rc else 'DIFFERENT')
        results.append(dict(name='gen_' + kk, status=status, produced=got, expected=exp.get(lab), minutes=m))
        log('gen_%s: %s in %.1f min' % (kk, status, m))
        if not status.startswith('IDENTICAL'):
            break
    # Lbg2 is written by qb2b_stages: put that step first when a later step needs it
    order = sorted(todo, key=lambda s: 0 if s[0] in ({'qb2b_stages'} | bg2_steps) else 1)
    if all(r['status'].startswith('IDENTICAL') for r in results):
        for name, _, cwd, argv, env, outs, mins, _ in order:
            for prod, kind in outs:
                if prod and not kind == 'bytes':
                    p = os.path.join(work, prod)
                    if os.path.exists(p):
                        os.remove(p)
            log('running %s: %s (about %s min in the lab)' % (name, ' '.join(argv), mins))
            rc, m = run(work, py, name, cwd, argv, env)
            r = dict(name=name, rc=rc, minutes=m, files=[] if rc else compare(work, name, outs))
            r['status'] = ('RUN FAILED (exit %d)' % rc if rc else
                           ('IDENTICAL' if all(f['status'].startswith('IDENTICAL') for f in r['files']) else 'DIFFERENT'))
            results.append(r)
            log('%s: %s in %.1f min' % (name, r['status'], m))
            for f in r['files']:
                if not f['status'].startswith('IDENTICAL'):
                    log('  %s: %s' % (f['produced'], f['status']))
                    for d in f.get('differences') or []:
                        log('    ' + d)
    with open(os.path.join(work, '_out', 'reproduce_report.json'), 'w', encoding='utf-8') as f:
        json.dump(dict(tier=tier, only=only, regen=regen, minutes=(time.time() - T0) / 60, results=results), f, indent=1)
    badn = [r['name'] for r in results if not r['status'].startswith('IDENTICAL')]
    missing = [s[0] for s in todo if s[0] not in {r['name'] for r in results}]
    log('total %.1f min' % ((time.time() - T0) / 60))
    log('ALL CHECKS PASS' if not badn and not missing else 'FAILED: ' + ', '.join(badn + ['%s (not run)' % x for x in missing]))
    return 0 if not badn and not missing else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
