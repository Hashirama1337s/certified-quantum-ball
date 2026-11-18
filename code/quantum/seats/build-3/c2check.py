"""build-3 control C2: resonance null at kappa* = [2.810552549 +/- 3.83e-10] (verdict_M.json): the matching
determinant (D_HKA and D_new, engine transport + centre stage, same instrument as N6) must enclose 0 -> uniqueness refused.
Positive side: D_HKA(2) excludes 0.  -> c2.json"""
import sys, json, time
sys.path.insert(0, '<LAB>/quantum/seats/build-3')
import qb2stages as ST
from flint import arb, acb, ctx
from core import sym
old = ctx.prec
try:
    L = ST.engine(log=lambda *a: None); ctx.prec = 128
    ks = arb('2.810552549') + sym(arb('3.83e-10'))
    S = lambda a: [a.mid().str(15, radius=False), a.rad().str(3, radius=False)]
    out = {}
    for nm, k in (('kstar', ks), ('kappa2', arb(2))):
        t = time.time()
        dh = L.D(acb(k), 1, 0); dn = L.D(acb(k), 0, 1)
        out[nm] = dict(D_HKA=S(dh.real), D_new=S(dn.real), D_HKA_contains_0=bool(dh.real.contains(0)),
                       D_new_contains_0=bool(dn.real.contains(0)), sec=round(time.time() - t, 1))
    out['C2_refused'] = out['kstar']['D_HKA_contains_0'] and out['kstar']['D_new_contains_0']
    out['kappa2_unique'] = (not out['kappa2']['D_HKA_contains_0']) and (not out['kappa2']['D_new_contains_0'])
    json.dump(out, open('c2.json', 'w'), indent=1); print(json.dumps(out, indent=1))
finally:
    ctx.prec = old
