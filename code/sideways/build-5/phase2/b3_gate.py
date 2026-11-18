"""b3_gate.py - build-5 phase 2, B3: the rigorous build-4 PL2b-3 16-point precision gate with a given outward record.

Everything as build-4 (l2b_ext.E5PathB + AxialB(lohner=True, match_k=k), l2b_cert.gate logic l.133-154, PREC 200, SU 5,
the pinned tube engine/cert/lin_tube2_k1_3.pkl); the ONLY change is the outward record file, set through l2_lib's module
attribute OUTW before the path is built (l2_lib.py:54, read at E5Path.__init__ l.366).  No edits, no shims, rigorous.
PROBE ONLY: no count.  Also records, per vertex, the relative radius of the light-cone Frobenius vector Y_lc (point value).

usage (memcap 8): py -3.11 b3_gate.py <tag> <outward pkl> [k,k,...]   -> b3_<tag>.json
"""
import os, sys, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
B4 = '<LAB>/sideways/build-4'
sys.path.insert(0, B4); sys.path.insert(0, '<TOOLS>')
import l2_cert as LC
import l2_lib as LL
import l2b_ext as X
from flint import arb, acb, ctx
ctx.cap = 64
assert ctx.cap >= 64


def relf(v):
    a = abs(v)
    if not a.mid() > 0:
        return float('inf')
    r = a.rad()/a.mid()
    return float(r.upper()) if r.is_finite() else float('inf')


def main():
    tag, pk = sys.argv[1], os.path.abspath(sys.argv[2])
    ks = [int(s) for s in sys.argv[3].split(',')] if len(sys.argv) > 3 else [21, 14, 7, 0]
    T0 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    old = ctx.prec
    try:
        ctx.prec = LC.PREC
        LL.OUTW = pk
        P2 = X.E5PathB(LL.Sys('axial', su=LC.SU), log=log)
        out = dict(tag=tag, outward=pk, outward_sha256=hashlib.sha256(open(pk, 'rb').read()).hexdigest(), prec=LC.PREC,
                   path_info=P2.info, gate_tol='1e-2', tried=[])
        _, V = LC.polygon_plan(complex(-1.8, 6.2), 0.5, 0.12)
        ylc = []
        for z in V:
            cx = LL.Ctx(acb(z.real, z.imag), arb(0), 0)
            Y, _ = LL.regsing_tm(cx, P2.light, None, P2.light['t'], S_set=[1])
            ylc.append(max(relf(cx.evalb(y, acb(0))) for y in Y))
        out['Ylc_rel_per_vertex'] = ['%.2e' % v for v in ylc]
        for k in ks:
            ax = X.AxialB(P2, lohner=True, match_k=k)
            rels = []
            for z in V:
                tm, cx = ax.D_tm(acb(z.real, z.imag), arb(0), 0)
                rels.append(relf(cx.evalb(tm, acb(0))))
            d = P2.near_pieces[k].d.str(5) if k < len(P2.near_pieces) else 'D_jF'
            out['tried'].append(dict(match_k=k, d_match=d, max_rel=max(rels), rel_radii=['%.2e' % x for x in rels],
                                     passed=bool(max(rels) <= 1e-2)))
            log('k=%d (d %s): max rel radius %.3e -> %s' % (k, d, max(rels), 'PASS' if max(rels) <= 1e-2 else 'fail'))
            json.dump(out, open(os.path.join(HERE, 'b3_%s.json' % tag), 'w'), indent=1)
        out['time_s'] = round(time.time() - T0, 1)
        json.dump(out, open(os.path.join(HERE, 'b3_%s.json' % tag), 'w'), indent=1)
        return 0
    finally:
        ctx.prec = old


if __name__ == '__main__':
    sys.exit(main())
