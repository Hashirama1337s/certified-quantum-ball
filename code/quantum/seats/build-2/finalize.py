"""build-2 QB1b: merge the sign-map parts, check exact coverage of [-12, 12], certify the zero (derivative TM on the
unresolved hull, bisection to <= 1e-8), and write result_qb1b.json + controls_qb1b.json.
usage (memcap 8): py -3.11 finalize.py"""
import sys, os, json, glob, hashlib, time
from fractions import Fraction as F
from flint import arb, ctx
ctx.cap = 64
import tm
from tm import TM
import tmcell, bg
HERE = os.path.dirname(os.path.abspath(__file__))


def A(fr):
    return arb(fr.numerator)/arb(fr.denominator)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def main():
    T0 = time.time()
    old = ctx.prec
    try:
        ctx.prec = 256
        parts = sorted(glob.glob(os.path.join(HERE, 'sm_*.json')))
        cells = []; unres = []; meta = []
        for p in parts:
            d = json.load(open(p))
            meta.append(dict(file=os.path.basename(p), region=d['region'], lo=d['lo'], hi=d['hi'], P=d['P'], evals=d['evals'],
                             cells=len(d['cells']), unresolved=len(d['unresolved']), time_s=d['time_s'], wpiece=d.get('wpiece')))
            for c in d['cells']:
                cells.append((F(c['a']), F(c['b']), c['sign'], d['region'], c['enc']))
            for u in d['unresolved']:
                unres.append((F(u['a']), F(u['b']), d['region'], u.get('err')))
        for a, b, *_ in cells + unres:
            assert a.denominator & (a.denominator - 1) == 0 and b.denominator & (b.denominator - 1) == 0
        # zero hull: unresolved cells away from the region seams
        z_unres = [u for u in unres if u[2] == 'tube' and F(-1) < u[0] < F(-1, 2)]
        def covered_by_cert(u):
            return any(c[0] <= u[0] and u[1] <= c[1] for c in cells)
        dropped = [u for u in unres if u not in z_unres and covered_by_cert(u)]
        other_unres = [u for u in unres if u not in z_unres and not covered_by_cert(u)]
        res = {}
        Bg = bg.Background()
        tm.P = 5
        zero = {}
        if z_unres:
            za = min(u[0] for u in z_unres); zb = max(u[1] for u in z_unres)
            zero['unresolved_hull'] = [str(za), str(zb)]
            # derivative on the hull, bisected until each piece is certified positive
            work = [(za, zb)]; dcells = []; ok = True
            while work:
                a, b = work.pop(0)
                try:
                    o = tmcell.cell(Bg, 'tube', A(a), A(b), deg=5, want_dx=True)
                    dv = tm.rng(o['db_nn_dx'])
                except ArithmeticError:
                    dv = None
                if dv is not None and dv > 0:
                    dcells.append(dict(a=str(a), b=str(b), db_nn_dx_lower=dv.lower().str(8, radius=False)))
                elif b - a > F(1, 10**9):
                    m = (a + b)/2; work[0:0] = [(a, m), (m, b)]
                else:
                    ok = False; break
            zero['derivative_cells'] = dcells
            zero['derivative_positive_on_hull'] = ok
            zero['db_nn_dx_lower_bound_on_hull'] = min((arb(c['db_nn_dx_lower']) for c in dcells), key=lambda v: float(v.mid())).str(8, radius=False) if dcells else None
            # neighbours
            left = [c for c in cells if c[1] == za]; right = [c for c in cells if c[0] == zb]
            zero['left_neighbour_sign'] = left[0][2] if left else None
            zero['right_neighbour_sign'] = right[0][2] if right else None
            # bisection (point TMs with H = 0, same s-chart path)
            def bval(x):
                o = tmcell.cell(Bg, 'tube', A(x), A(x)); return tm.rng(o['b_nn'])
            lo, hi = za, zb
            slo = -1
            it = 0
            while hi - lo > F(1, 2*10**8) and it < 100:
                m = (lo + hi)/2
                v = bval(m)
                s = 1 if v > 0 else (-1 if v < 0 else 0)
                if s == 0:
                    # midpoint within the point-evaluation resolution of the zero: bracket it by m -+ 2^-29
                    zero['bisection_midpoint_uncertified'] = str(m)
                    dlt = F(1, 2**29)
                    va, vb = bval(m - dlt), bval(m + dlt)
                    if va < 0 and vb > 0:
                        lo, hi = m - dlt, m + dlt
                    else:
                        zero['bisection_stopped_at'] = str(m)
                    break
                if s == slo: lo = m
                else: hi = m
                it += 1
            vlo, vhi = bval(lo), bval(hi)
            zero['zero_enclosure'] = [str(lo), str(hi)]
            zero['zero_enclosure_decimal'] = [A(lo).str(20, radius=False), A(hi).str(20, radius=False)]
            zero['zero_enclosure_width'] = float(hi - lo)
            zero['b_nn_at_enclosure_ends'] = [vlo.str(5), vhi.str(5)]
            zero['ends_signs_certified'] = bool(vlo < 0 and vhi > 0)
        # coverage
        cov = sorted([(a, b) for a, b, *_ in cells] + [(a, b) for a, b, *_ in z_unres])
        cur = F(-12); gaps = []
        for a, b in cov:
            if a > cur:
                gaps.append([str(cur), str(a)])
            cur = max(cur, b)
        covered = (not gaps) and cur >= 12 and cov[0][0] <= -12
        negs = [c for c in cells if c[2] < 0]; poss = [c for c in cells if c[2] > 0]
        unique = bool(z_unres) and zero.get('derivative_positive_on_hull') and zero.get('left_neighbour_sign') == -1 \
            and zero.get('right_neighbour_sign') == 1 and not other_unres and covered
        res['PQ2b'] = dict(
            method='Taylor models in eps = x - x_m (degree P, remainder), exact dyadic cells, regions w/centre/tube/sonic/outw; see RUN.md',
            parts=meta, n_certified_cells=len(cells), n_unresolved=len(unres), n_unresolved_inside_other_certified_cells=len(dropped), unresolved_other_than_zero=[[str(u[0]), str(u[1]), u[2], u[3]] for u in other_unres],
            coverage_minus12_to_12='complete, no gap' if covered else 'GAPS', gaps=gaps,
            negative_cells=len(negs), positive_cells=len(poss),
            max_negative_x=str(max(c[1] for c in negs)) if negs else None, min_positive_x=str(min(c[0] for c in poss)) if poss else None,
            min_cell_width=float(min(b - a for a, b, *_ in cells)), max_cell_width=float(max(b - a for a, b, *_ in cells)),
            zero=zero, unique_zero_on_minus12_to_12=unique,
            cells=[[str(a), str(b), s, r] for a, b, s, r, e in sorted(cells)])
        for k, fn in (('PQ6b', 'pq6b.json'), ('PQ7b', 'pq7b.json'), ('PQ8', 'pq8.json'), ('controls', 'controls_tm.json')):
            res[k] = json.load(open(os.path.join(HERE, fn)))
        if 'rows' in res['PQ6b']:
            pass
        res['time_s_finalize'] = round(time.time() - T0, 1)
        json.dump(res, open(os.path.join(HERE, 'result_qb1b.json'), 'w'), indent=1)
        files = ['tm.py', 'tmbg.py', 'tmcell.py', 'schart.py', 'signmap2.py', 'tasks2.py', 'dnew2.py', 'finalize.py']
        b1 = ['ajet.py', 'aser.py', 'charts.py', 'bg.py', 'ctl.py', 'qb1.py']
        ctl = dict(instrument_sha256={f: sha(os.path.join(HERE, f)) for f in files},
                   build1_imported_sha256={f: sha(os.path.join('<LAB>/quantum/seats/build-1', f)) for f in b1},
                   inputs_sha256={p: sha(p) for p in ['<LAB>/engine/cert/lin_tube2_k1_3.pkl',
                                                      '<LAB>/sideways/build-5/phase2/outward_C.pkl',
                                                      '<LAB>/engine/cert/fields.pkl', '<LAB>/engine/cert/taylor.py']},
                   parts_sha256={os.path.basename(p): sha(p) for p in parts},
                   result_qb1b_sha256=sha(os.path.join(HERE, 'result_qb1b.json')),
                   controls=res['controls'])
        json.dump(ctl, open(os.path.join(HERE, 'controls_qb1b.json'), 'w'), indent=1)
        print(json.dumps({k: v for k, v in res['PQ2b'].items() if k not in ('cells', 'parts')}, indent=1)[:3000])
    finally:
        ctx.prec = old


if __name__ == '__main__':
    main()
