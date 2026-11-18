"""build-2 QB1b PQ7b: replay of D_new(2) with the sealed count instrument engine/cert/count_b4v3.py (install(),
imported, never edited) at its own matching point x_m = XBAR = -4.6.  Same calls as build-1 dnew.py.
usage (memcap 8): py -3.11 dnew2.py -> pq7b.json"""
import sys, os, time, json, hashlib
from flint import ctx as _ctx
_ctx.cap = 32
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import count_b4v3 as CV
import count_b4 as CB
import lin_stages_v2 as LS2
import sympy as sp
from flint import arb, acb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
S = lambda a: [a.mid().str(20, radius=False), a.rad().str(3, radius=False)]
T0 = time.time(); old = ctx.prec
out = dict(engine_sha256={f: hashlib.sha256(open('<LAB>/' + f, 'rb').read()).hexdigest() for f in
           ['engine/cert/count_b4v3.py', 'engine/cert/count_b4.py', 'engine/cert/count_b4v3_seal.py', 'engine/cert/lin_stages_v2.py',
            'engine/cert/lin_stages.py', 'engine/cert/lin_tube2_k1_3.pkl']})
try:
    CV.install(); assert ctx.cap == 32
    L = LS2.Linear2(sp.Rational(1, 3), log=print, nsub=CB.NSUB, tol=CB.TOL)
    out['linear'] = CV.INFO.get('linear'); out['prec'] = ctx.prec; out['x_m'] = L.xbar.str(20)
    d2 = L.D(acb(2), 0, 1); h2 = L.D(acb(2), 1, 0)
    out['D_new_2'] = dict(re=S(d2.real), im=S(d2.imag)); out['D_HKA_2'] = dict(re=S(h2.real), im=S(h2.imag))
    out['D_new_2_nonzero'] = bool(d2.real < 0 or d2.real > 0); out['D_new_2_sign'] = -1 if d2.real < 0 else 1
    out['abs_D_new_2_lower'] = abs(d2).lower().str(12, radius=False)
    b1 = json.load(open('<LAB>/quantum/seats/build-1/dnew_qb1.json'))['values']['2']['D_new']['re']
    out['build1_value'] = b1; out['overlaps_build1'] = bool(d2.real.overlaps(arb(b1[0]) + arb(0, arb(b1[1]).upper())))
finally:
    ctx.prec = old
out['time_s'] = round(time.time() - T0, 1)
json.dump(out, open(os.path.join(HERE, 'pq7b.json'), 'w'), indent=1); print(out)
