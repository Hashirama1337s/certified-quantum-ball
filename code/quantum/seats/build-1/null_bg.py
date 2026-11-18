"""build-1 QB1 item 6 (null): an OFF-SOLUTION background, V0* = E5 box upper edge + 1e-15
(0.1124394013880983362 + 3.31e-20 + 1e-15, i.e. ~3e4 E5 radii outside the certified box), integrated the same way
as the certified records, so that the QB1 instrument can be run on it unchanged.

  in   inward leg = the control flow of engine/cert/lin_tube_v2.build (sonic manifold v32 with Dual V0, nu ladder,
       tau leg to x = -0.05, switch set, weighted-QR Lohner x-leg to -4.6, step records c / setrad / h), copied here
       because build() writes into engine/cert (it opens the pinned lin_tube2_k1_3.pkl for writing); no centre manifold
       (the shifted solution is not the regular-centre solution).  Same PREC 300, KS 80, KT 36, Q 0.25.
  out  outward leg = sideways/build-1/outward_rec.outward_rec (imported read-only, as build-5 b2_outward.py), BG2
       parameters PREC 256, K_t 50, X1 = 12.
Records are strings only.  usage (memcap 8): py -3.11 null_bg.py in|out  -> null_in.pkl / null_out.pkl
"""
import sys, os, time, pickle, hashlib
os.environ.setdefault('LB_NGRID', '24')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
mode = sys.argv[1]
if mode == 'out':
    sys.path.insert(0, '<LAB>/sideways/build-1')
    import outward_rec as OR          # sets ctx.prec 200, ctx.cap 32 at import (pinned module)
import driver_v32 as V32
import driver_v31 as V31
import sympy as sp
from flint import arb, ctx
import driver_v3 as V3
import driver_v2 as V2
import lohner as LO
import taylor as TY
import param as PM
from param import val
from core import sym, Dual
ctx.cap = 64
assert PM.arb_series is V31.arb_series_exact, 'R7 not installed'
assert LO.lohner_step is V3._step and V3.MODE == 'qrw', 'weighted-QR Lohner not installed'

V0STAR = arb('0.1124394013880983362') + arb('3.31e-20') + arb('1e-15')
V0RAD = arb('1e-30')
OM = '0.3256888894230684'
K13 = sp.Rational(1, 3)


def ser(a):
    return (a.mid().str(110, radius=False), a.rad().str(10, radius=False))


def inward(log):
    PREC, KS, KT, Q = 300, 80, 36, '0.25'
    NU_LADDER = ['0.06', '0.03', '0.02', '0.01', '0.004', '0.002', '0.001']
    XBAR = arb('-4.6'); X_SWITCH = arb('-0.05')
    ctx.prec = PREC
    T0 = time.time()
    V0 = V0STAR + sym(V0RAD)
    Pb = V2.Problem2(K13, prec=PREC, K_sonic=KS, K_centre=60, K_ty=KT)
    sm = None
    for nu in NU_LADDER:
        try:
            sm = V32.sonic_manifold_v32(Pb, Dual(V0, [arb(1)]), OM, nu, KS)
            if sm['eps'].is_finite() and sm['epsX'].is_finite(): break
            sm = None
        except (ArithmeticError, AssertionError) as ex:
            log('nu_s=%s: %s' % (nu, str(ex)[:90])); sm = None
    if sm is None: raise ArithmeticError('no nu_s contracts')
    nu_s = arb(nu)
    log('sonic manifold: nu_s=%s eps=%s epsX=%s (%.0fs)' % (nu, sm['eps'].str(3), sm['epsX'].str(3), time.time() - T0))
    th0 = -nu_s/20
    smc = sm['smc']
    c_thin = [PM.horner([val(a) for a in smc['co']['Y'][i]], th0) for i in range(3)] + [PM.horner([val(a) for a in smc['co']['X']], th0)]
    d_box = [PM.horner([a.d[0] for a in sm['co']['Y'][i]], th0) for i in range(3)] + [PM.horner([a.d[0] for a in sm['co']['X']], th0)]
    tails = [sm['eps']]*3 + [sm['epsX']]
    c_thin = [c_thin[i] + sym(tails[i]) for i in range(4)]
    S4 = V3.aligned_set_parts(c_thin, d_box, V0.rad())
    T4 = V2.tau_tape(K13)
    steps = []
    qa = arb(Q)

    def record(kind, Sset, F, ha, ncomp):
        c = [arb(t.mid()) for t in Sset.c]
        y1, Bx, Yc = TY.step_cauchy2(F, c, ha, KT, q=qa, conv=lambda cc: PM._c(cc, arb(0)))
        setrad = max(t.rad() for t in Sset.hull())
        steps.append(dict(kind=kind, h=ser(ha), c=[ser(t) for t in c], Rc=ser(abs(ha)/qa), setrad=ser(setrad)))
    xt = float(X_SWITCH.mid()); h = 0.02; sgn = None; n = 0
    while True:
        rate = float(T4.eval([arb(v.mid()) for v in S4.c], conv=lambda q_: PM._c(q_, arb(0)))[3].mid())
        xcur = float(S4.c[3].mid())
        if sgn is None: sgn = 1 if (xt - xcur)*rate > 0 else -1
        remaining = (xt - xcur)/rate*sgn
        last = remaining <= h
        ha = V2.land_tau(T4, S4.c, X_SWITCH, KT) if last else arb(repr(sgn*h))
        saved = (list(S4.c), [r_[:] for r_ in S4.B], list(S4.r))
        try:
            record('tau', S4, T4, ha, 4)
            LO.lohner_step(T4, S4, ha, KT, qa)
        except ArithmeticError:
            S4.c, S4.B, S4.r = saved
            if steps and steps[-1]['h'] == ser(ha): steps.pop()
            h /= 2
            if h < 1e-8: raise
            continue
        n += 1
        if last: break
        h *= 1.25
    log('tau leg: %d steps, width %s (%.0fs)' % (n, max(t.rad() for t in S4.hull()).str(3), time.time() - T0))
    S = V3.switch_set(Pb, S4, X_SWITCH)
    x = arb(X_SWITCH); h = 0.002; n = 0
    while float(x.mid()) > float(XBAR.mid()):
        rem = float((x - XBAR).mid()); hh = min(h, rem); last = hh >= rem - 1e-15
        ha = (XBAR - x) if last else arb(repr(-hh))
        saved = (list(S.c), [r_[:] for r_ in S.B], list(S.r))
        try:
            record('x', S, Pb.F, ha, 3)
            LO.lohner_step(Pb.F, S, ha, KT, qa)
        except ArithmeticError:
            S.c, S.B, S.r = saved
            if steps and steps[-1]['h'] == ser(ha): steps.pop()
            h = hh/2
            if h < 1e-8: raise
            continue
        x = XBAR if last else x + ha; n += 1; h = hh*1.25
        if n % 100 == 0: log('  x step %d x=%.3f width %s (%.0fs)' % (n, float(x.mid()), max(t.rad() for t in S.hull()).str(3), time.time() - T0))
    log('x leg: %d steps, end width %s (%.0fs)' % (n, max(t.rad() for t in S.hull()).str(3), time.time() - T0))
    out = dict(k='1/3', nu_s=nu, th0=ser(th0), XBAR=ser(XBAR), X_SWITCH=ser(X_SWITCH), KT=KT, Q=Q, prec=PREC,
               V0=ser(V0), source='build-1 QB1 null_bg.py inward (off-solution V0*)',
               sm=dict(Y=[[ser(val(a)) for a in sm['co']['Y'][i]] for i in range(3)], X=[ser(val(a)) for a in sm['co']['X']],
                       eps=ser(sm['eps']), epsX=ser(sm['epsX']), lam=ser(val(sm['lam']))),
               steps=steps, end=[ser(t) for t in S.hull()])
    return out


def outward(log):
    PREC, KT = 256, 50
    ctx.prec = PREC
    OR.KT = KT
    V0box = V0STAR + sym(V0RAD)
    Pb = V31.Problem2(K13, prec=PREC, K_sonic=OR.KS, K_centre=OR.KC, K_ty=KT)
    nu, _ = V2.choose_nu(Pb, V0box + sym(OR.NU_AT), OR.OM_T, OR.NU_LADDER, OR.KS, log=log)
    M = V32.MatcherV32(Pb, nu, OR.NUC, Ks=OR.KS, x_switch='-0.05', Kt=KT)
    r = OR.outward_rec(Pb, M, V0box, OR.OM_T, X1=12.0, log=log)
    log('outward: steps=%d D<0 %s V>0 %s lemma %s' % (r['steps'], r['ok_D'], r['ok_V'], r.get('lemma')))
    return dict(k='1/3', V0=(V0box.mid().str(60, radius=False), V0box.rad().str(10, radius=False)), nu_s=nu, X1=12.0, KT=KT,
                prec=PREC, steps=r['records'], sm=r['sm'], ok_D=r['ok_D'], ok_V=r['ok_V'], end=[OR.ser(t) for t in r['y1']],
                lemma=r.get('lemma'), source='build-1 QB1 null_bg.py outward (off-solution V0*), outward_rec.outward_rec')


if __name__ == '__main__':
    T00 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T00), *a, flush=True)
    old = ctx.prec
    try:
        out = inward(log) if mode == 'in' else outward(log)
    finally:
        ctx.prec = old
    out['instrument_sha256'] = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    out['time_s'] = round(time.time() - T00, 1)
    fn = os.path.join(HERE, 'null_%s.pkl' % mode)
    pickle.dump(out, open(fn, 'wb'))
    log('saved', fn)
