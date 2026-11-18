"""p_dover.py - control: point values D(kappa) at the 16 gate vertices with the OLD and the REBUILT outward record must
overlap (same function, two rigorous enclosures).  Same build-4 instrument (E5PathB + AxialB lohner, match_k given)."""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '<LAB>/sideways/build-4'); sys.path.insert(0, '<TOOLS>')
import l2_cert as LC
import l2_lib as LL
import l2b_ext as X
from flint import arb, acb, ctx
ctx.cap = 64
T0 = time.time()
k = int(sys.argv[1]) if len(sys.argv) > 1 else 0
old = ctx.prec
try:
    ctx.prec = LC.PREC
    _, V = LC.polygon_plan(complex(-1.8, 6.2), 0.5, 0.12)
    vals = {}
    for tag, pk in (('old', '<LAB>/sideways/build-1/outward_rec.pkl'), ('new', os.path.join(HERE, 'outward_C.pkl'))):
        LL.OUTW = pk
        P2 = X.E5PathB(LL.Sys('axial', su=LC.SU), log=lambda *a: None)
        ax = X.AxialB(P2, lohner=True, match_k=k)
        vals[tag] = []
        for z in V:
            tm, cx = ax.D_tm(acb(z.real, z.imag), arb(0), 0)
            vals[tag].append(cx.evalb(tm, acb(0)))
        print('%s built and evaluated (%.0fs)' % (tag, time.time() - T0), flush=True)
    ov = [bool(a.overlaps(b)) for a, b in zip(vals['old'], vals['new'])]
    out = dict(match_k=k, overlaps=ov, all_overlap=all(ov), D_new=[v.str(12, radius=True) for v in vals['new']],
               D_old=[v.str(6, radius=True) for v in vals['old']])
    json.dump(out, open(os.path.join(HERE, 'p_dover_k%d.json' % k), 'w'), indent=1)
    print('all overlap:', all(ov)); print('D_new v0, v8:', out['D_new'][0], out['D_new'][8])
finally:
    ctx.prec = old
