#!/usr/bin/env python
"""fill_balls.py -- every ball printed in the release text comes from a stored string, verbatim or rounded outward.

A ball is printed as `[mid ± rad]` (also `mid ± rad`, `+/-`, `+-`). For every printed ball in the text files below this
script finds the stored ball it stands for and then
  * leaves it alone if it is the stored string verbatim (same mid and radius strings), or
  * rewrites it as an OUTWARD rounding of the stored ball: the printed midpoint digits are kept, and the radius is the
    stored radius plus the midpoint shift, rounded UP to two significant digits. It then asserts, in exact decimal
    arithmetic, that the printed ball contains the stored ball: |printed mid - stored mid| + stored rad <= printed rad.
Stored balls are read from the shipped files only: every JSON file under results/ and code/ (pairs ["mid", "rad"] and
strings "[mid +/- rad]"), every log and Markdown file there ("[mid +/- rad]", "[mid ± rad]"), plus a few DERIVED balls
computed here with python-flint (arb) from stored strings (DERIVED below). A printed ball with no stored source is an
error.

Usage (from anywhere):
    python tools/fill_balls.py --check     report every printed ball: VERBATIM, OUTWARD, or NOT OUTWARD / NO SOURCE;
                                           exit 0 only if every print is verbatim or contains its stored ball
    python tools/fill_balls.py --write     rewrite the non-verbatim prints as outward roundings, then --check
Library: fix_text(text) -> (new text, report rows), used by tools/make_scorecard.py for SCORECARD.md.
"""
import json
import os
import re
import sys
from decimal import Decimal, ROUND_CEILING, getcontext

getcontext().prec = 80
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEXT = ['paper.md', 'README.md', 'CLAIMS.md', 'SCORECARD.md', 'CITATION.cff']
NUM = r'[−-]?\d+(?:\.\d+)?(?:e[−+-]?\d+)?'
BALL = re.compile(r'(?<![\w.])(' + NUM + r')(\s*)(±|\+/-|\+-)(\s*)(\d+(?:\.\d+)?(?:e[−+-]?\d+)?)(?!\w|\.\d)')
ARB = re.compile(r'\[\s*(' + NUM + r')\s*(?:\+/-|±)\s*(\d+(?:\.\d+)?(?:e[−+-]?\d+)?)\s*\]')
# derived balls: (name, function of the stored-ball index) -> (mid string, rad string)
PSI2 = ('3.4616358037141656915', '3.02e-12')   # results/result_qb2b.json PQ2b-5.N8.Psi2 (checked below)


def D(s):
    return Decimal(s.replace('−', '-'))


def _scan_json(x, out, src):
    if isinstance(x, list):
        if len(x) == 2 and all(isinstance(v, str) for v in x):
            try:
                m, r = D(x[0]), D(x[1])
                if r >= 0:
                    out.append((x[0], x[1], src))
                    return
            except Exception:
                pass
        for v in x:
            _scan_json(v, out, src)
    elif isinstance(x, dict):
        for v in x.values():
            _scan_json(v, out, src)
    elif isinstance(x, str):
        for m in ARB.finditer(x):
            out.append((m.group(1), m.group(2), src))


def stored_index():
    out = []
    for base in ('results', 'code'):
        for d, dirs, fs in os.walk(os.path.join(ROOT, base)):
            dirs[:] = [x for x in dirs if x != '__pycache__']
            for f in fs:
                p = os.path.join(d, f)
                rel = os.path.relpath(p, ROOT).replace('\\', '/')
                if f.endswith(('.json', '.jsonl')):
                    try:
                        txt = open(p, encoding='utf-8').read()
                        docs = [json.loads(txt)] if f.endswith('.json') else [json.loads(l) for l in txt.splitlines() if l.strip()]
                    except Exception:
                        continue
                    for doc in docs:
                        _scan_json(doc, out, rel)
                elif f.endswith(('.log', '.md', '.txt', '.out')):
                    for m in ARB.finditer(open(p, encoding='utf-8', errors='replace').read()):
                        out.append((m.group(1), m.group(2), rel))
    out += derived(out)
    return out


def derived(idx):
    """Balls computed from the stored Psi2 string with arb; the radius is enlarged and rounded up."""
    if not any(m == PSI2[0] and r == PSI2[1] and s.endswith('result_qb2b.json') for m, r, s in idx):
        sys.exit('fill_balls: the stored Psi2 ball was not found in results/result_qb2b.json')
    from flint import arb, ctx
    ctx.prec = 200
    p = arb(PSI2[0], PSI2[1])
    k = 1 / (30 * arb.pi())
    rows = []
    for name, v in (('K*Psi2 at N_c = 12 from result_qb2b.json PQ2b-5.N8.Psi2', 12 * p * k),
                    ('Psi2/(30 pi) from result_qb2b.json PQ2b-5.N8.Psi2', p * k)):
        mid = v.mid().str(30, radius=False)
        rad = fmt_rad(D(repr(float(v.rad()))) * Decimal('1.001'))
        assert arb(mid, rad).contains(v)
        rows.append((mid, rad, 'DERIVED: ' + name))
    return rows


def ceil2(x):
    if x == 0:
        return Decimal(0)
    e = x.adjusted() - 1
    q = Decimal(1).scaleb(e)
    return (x / q).to_integral_value(rounding=ROUND_CEILING) * q


def fmt_rad(x):
    s = '%.2g' % float(x)
    if D(s) < x:                      # float formatting must never round a radius down
        s = '%.2g' % float(ceil2(x * Decimal('1.000001')))
    s = s.replace('e-0', 'e-').replace('e+0', 'e').replace('e+', 'e')
    assert D(s) >= x
    return s


def ulp(s):
    s = s.replace('−', '-')
    mant, _, ex = s.partition('e')
    dec = len(mant.split('.')[1]) if '.' in mant else 0
    return Decimal(1).scaleb(-dec + (int(ex) if ex else 0))


def find(idx, pm, pr):
    pmD, prD = D(pm), D(pr)
    for m, r, s in idx:
        if m.replace('−', '-') == pm.replace('−', '-') and r == pr:
            return 'VERBATIM', (m, r, s)
    u = ulp(pm)
    best = None
    for m, r, s in idx:
        try:
            mD, rD = D(m), D(r)
        except Exception:
            continue
        if abs(mD - pmD) > u * Decimal('1.000001'):
            continue
        if rD == 0 and prD != 0:
            continue
        if prD and not (prD / 20 <= rD <= prD * 20):
            continue
        key = (abs(rD - prD) / (prD or 1), abs(mD - pmD))   # closest radius first, then closest midpoint
        if best is None or key < best[0]:
            best = (key, (m, r, s))
    if best is None:
        return 'NO SOURCE', None
    m, r, s = best[1]
    if abs(D(m) - pmD) + D(r) <= prD:
        return 'OUTWARD', best[1]
    return 'NOT OUTWARD', best[1]


def fix_text(text, idx=None, write=True):
    idx = idx if idx is not None else stored_index()
    rows = []

    def sub(mo):
        pm, sp1, pmsign, sp2, pr = mo.groups()
        st, src = find(idx, pm, pr)
        if st == 'NOT OUTWARD' and write:
            m, r, s = src
            new_r = fmt_rad(abs(D(m) - D(pm)) + D(r))
            assert abs(D(m) - D(pm)) + D(r) <= D(new_r), (pm, new_r, m, r)
            rows.append(('REWRITTEN OUTWARD', pm, pr, new_r, s))
            return pm + sp1 + pmsign + sp2 + new_r
        rows.append((st, pm, pr, src[1] if src else '', src[2] if src else ''))
        return mo.group(0)
    return BALL.sub(sub, text), rows


def main(argv):
    sys.stdout.reconfigure(encoding="utf-8")
    write = '--write' in argv
    idx = stored_index()
    bad = 0
    total = 0
    for f in TEXT:
        p = os.path.join(ROOT, f)
        t = open(p, encoding='utf-8').read()
        if write:
            t2, rows = fix_text(t, idx, write=True)
            if t2 != t:
                open(p, 'w', encoding='utf-8', newline='\n').write(t2)
            for r in rows:
                if r[0] == 'REWRITTEN OUTWARD':
                    print('%s: %s ± %s -> ± %s (stored %s)' % (f, r[1], r[2], r[3], r[4]))
            t = t2
        _, rows = fix_text(t, idx, write=False)
        for st, pm, pr, sr, src in rows:
            total += 1
            if st not in ('VERBATIM', 'OUTWARD'):
                bad += 1
                print('%-12s %s: %s ± %s  (closest stored rad %s, %s)' % (st, f, pm, pr, sr, src))
    print('fill_balls: %d printed balls checked in %s; %d not verbatim-or-outward' % (total, ', '.join(TEXT), bad))
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
