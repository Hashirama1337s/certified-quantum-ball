"""build-4 QB2b item 1: inward tube from the BG2 sonic box, by the pinned engine/cert/lin_tube_v2.build (imported,
NOT edited).  Two changes, both in this process only:
  * ENC[1/3]['V0'] = BG2 box from sideways/build-5/phase2/outward_C.pkl ('V0', width 7.58e-29, inside the E5 box);
  * the module's `open` is shadowed so its pickle lands in build-4/lin_tube_bg2_k1_3.pkl (never engine/cert).
W0, xc (centre side) stay the E5 boxes.  usage (memcap 8): py -3.11 tube_bg2.py"""
import sys, os, time, pickle, builtins
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import lin_tube_v2 as LT
HERE = '<LAB>/quantum/seats/build-4/'
OUT = HERE + 'lin_tube_bg2_k1_3.pkl'
bg2 = pickle.load(open('<LAB>/sideways/build-5/phase2/outward_C.pkl', 'rb'))
k = sp.Rational(1, 3)
old = LT.ENC[k]['V0']
LT.ENC[k]['V0'] = (bg2['V0'][0], bg2['V0'][1])
LT.ENC[k]['src'] = 'BG2 V0 box (sideways/build-5/phase2/outward_C.pkl), W0/xc E5 (runs6/target.log PE1)'


def _open(fn, *a, **kw):
    if isinstance(fn, str) and fn.replace(chr(92), '/').startswith('<LAB>/engine/cert/lin_tube2_'):
        fn = OUT
    return builtins.open(fn, *a, **kw)
LT.open = _open
# module constants of lin_tube_v2 are built at import precision (53 bits): rebuild the SAME values at prec 300
from flint import arb, ctx as _ctx
_old_prec = _ctx.prec; _ctx.prec = LT.PREC
LT.XBAR = arb('-4.6'); LT.X_SWITCH = arb('-0.05')
_ctx.prec = _old_prec
print('XBAR', LT.XBAR.str(30), 'X_SWITCH', LT.X_SWITCH.str(30), flush=True)
T0 = time.time()
print('V0 E5', old, '-> BG2', LT.ENC[k]['V0'], flush=True)
LT.build(k, log=lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True))
assert os.path.exists(OUT)
print('done', OUT)
