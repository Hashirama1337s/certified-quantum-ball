"""build-1 QB1 item 7: certified point values of the spherical characteristic function D_new(kappa) with the sealed
count instrument engine/cert/count_b4v3.py (BAR-LB-M HOLD b1b0ee1cde), imported and installed, never edited.

count_b4v3.install() substitutes its engine classes (cached Linear2 for k = 1/3 with sha256 pin, LinearTM3) exactly
as in the sealed count; L.D(kappa, a=0, b=1) is D_new (no gauge zero, PROOF-B4 Lemma 1.5) and L.D(kappa, 1, 0) is
D_HKA, both as acb balls.  The matching point of this instrument is fixed by its tube: x_m = XBAR = -4.6
(lin_tube2_k1_3.pkl).  The instrument has no x_m parameter; x_m = -2 is NOT available without editing it.
Float cross-check (NON-COUNTING): engine/cert/lin_float.py Problem(xm).D at x_m = -4.6 and -2.
usage (memcap 8): py -3.11 dnew.py  -> dnew_qb1.json
"""
import sys, os, time, json, hashlib
from flint import ctx as _ctx
_ctx.cap = 32                       # = the sealed count stack's cap (count_b4v3.py sets 32 before its imports)
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import count_b4v3 as CV
import count_b4 as CB
import lin_stages_v2 as LS2
import sympy as sp
from flint import arb, acb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))


def S(a):
    return [a.mid().str(20, radius=False), a.rad().str(3, radius=False)]


def main():
    T0 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    old = ctx.prec
    out = dict(instrument='engine/cert/count_b4v3.py (+count_b4.py, count_b4v3_seal.py, lin_stages_v2, lin_tm_v3), imported, install()',
               engine_sha256={f: hashlib.sha256(open('<LAB>/' + f, 'rb').read()).hexdigest()
                              for f in ['engine/cert/count_b4v3.py', 'engine/cert/count_b4.py', 'engine/cert/count_b4v3_seal.py',
                                        'engine/cert/lin_stages_v2.py', 'engine/cert/lin_stages.py', 'engine/cert/lin_tube2_k1_3.pkl']},
               script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
    try:
        CV.install()
        assert ctx.cap == 32
        L = LS2.Linear2(sp.Rational(1, 3), log=log, nsub=CB.NSUB, tol=CB.TOL)
        out['linear'] = CV.INFO.get('linear'); out['prec'] = ctx.prec
        out['x_m'] = L.xbar.str(20)
        vals = {}
        for kap in ('0', '1', '2'):
            t = time.time()
            dn = L.D(acb(arb(kap)), 0, 1)
            dh = L.D(acb(arb(kap)), 1, 0)
            vals[kap] = dict(D_new=dict(re=S(dn.real), im=S(dn.imag)), D_HKA=dict(re=S(dh.real), im=S(dh.imag)),
                             D_new_re_sign=(1 if dn.real > 0 else (-1 if dn.real < 0 else 0)), sec=round(time.time() - t, 1))
            log('kappa=%s: D_new = %s   D_HKA = %s  (%.0fs)' % (kap, dn.str(15), dh.str(15), time.time() - t))
        out['values'] = vals
        d2 = L.D(acb(2), 0, 1)
        out['D_new_2_nonzero_certified'] = bool(d2.real > 0 or d2.real < 0)
        out['D_new_2_abs_lower_bound'] = abs(d2).lower().str(10, radius=False) if hasattr(abs(d2), 'lower') else None
    finally:
        ctx.prec = old
    # float cross-check (NON-COUNTING)
    try:
        import lin_float as LF
        fl = {}
        for xm in (-4.6, -2.0):
            pr = LF.Problem(xm=xm)
            fl[str(xm)] = {k: complex(pr.D(float(k), 0, 1)[0]).real for k in (0.0, 1.0, 2.0)}
        out['float_lin_float_D_new_NONCOUNTING'] = fl
        log('float lin_float D_new:', fl)
    except Exception as e:
        out['float_lin_float_error'] = repr(e)[:200]
    out['time_s'] = round(time.time() - T0, 1)
    json.dump(out, open(os.path.join(HERE, 'dnew_qb1.json'), 'w'), indent=1)
    log('wrote dnew_qb1.json')


if __name__ == '__main__':
    main()
