"""B4 build seat 2 (NEW file; lin_tube.py is left untouched): background tube for stage L, rebuilt
from the E5 certified enclosure (BAR-LB-E5 HOLD d0af75803f) with the repaired sonic manifold.

Changes against lin_tube.py (everything else is the same control flow, imported unchanged):
  T1  enclosure = runs6/target.log (k=1/3, PE1 K boxes at r=1e-15) and runs6/warmup.log (k=1/2, PE4),
      not runs3/RESULTS-E3.txt (E3 is VOID).
  T2  R7: driver_v32 is imported FIRST, which imports driver_v31, which rebinds param.arb_series to the
      untruncated series (driver_v31.py:53-59).  Asserted below.  lin_tube.py imported only driver_v3 /
      driver_v2, so its param.tail_bound inputs were truncated at ctx.cap = 10 (the R7 defect).
  T3  R8: the sonic manifold over the V0 box comes from driver_v32.sonic_manifold_v32 (mean-value head,
      driver_v32.py:104-127), not driver_v2.sonic_manifold_v2 (box head, the R8 defect).
usage: lin_tube_v2.py k  -> lin_tube2_k<k>.pkl
"""
import sys, os, time, pickle
os.environ.setdefault('LB_NGRID', '24')
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import driver_v32 as V32        # FIRST: installs R7 (via driver_v31) and the weighted-QR Lohner step
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
assert PM.arb_series is V31.arb_series_exact, 'R7 not installed'
assert LO.lohner_step is V3._step and V3.MODE == 'qrw', 'weighted-QR Lohner not installed'

PREC = 300
XBAR = arb('-4.6')
X_SWITCH = arb('-0.05')
NU_LADDER = ['0.06', '0.03', '0.02', '0.01', '0.004', '0.002', '0.001']
KS = 80; KT = 36; Q = '0.25'

# certified E5 enclosures (K boxes at the first contracting rung r = 1e-15)
ENC = {
    sp.Rational(1, 3): dict(V0=('0.1124394013880983362', '3.31e-20'), W0=('8.90132327775564352', '7.92e-18'),
                            xc=('0.212365646758945385', '7.95e-19'), om='0.3256888894230684',
                            src='runs6/target.log PE1 r=1e-15'),
    sp.Rational(1, 2): dict(V0=('-0.565685424949238020', '5.73e-19'), W0=('0.2962962962962963', '2.32e-17'),
                            xc=('0.0616119383369453', '4.08e-17'), om='0.7058823529411765',
                            src='runs6/warmup.log PE4 r=1e-15'),
}

def box(p):
    return arb(p[0]) + sym(arb(p[1]))

def ser(a):
    return (a.mid().str(110, radius=False), a.rad().str(10, radius=False))

def build(k, log=print):
    ctx.prec = PREC
    T0 = time.time()
    e = ENC[k]
    V0, W0, xc = box(e['V0']), box(e['W0']), box(e['xc'])
    Pb = V2.Problem2(k, prec=PREC, K_sonic=KS, K_centre=60, K_ty=KT)
    sm = None
    for nu in NU_LADDER:
        try:
            sm = V32.sonic_manifold_v32(Pb, Dual(V0, [arb(1)]), e['om'], nu, KS)
            if sm['eps'].is_finite() and sm['epsX'].is_finite(): break
            log('nu_s=%s: eps not finite' % nu); sm = None
        except (ArithmeticError, AssertionError) as ex:
            log('nu_s=%s: %s' % (nu, str(ex)[:90])); sm = None
    if sm is None: raise ArithmeticError('no nu_s contracts')
    nu_s = arb(nu)
    log('sonic manifold (R7+R8): nu_s=%s eps=%s epsX=%s (%.0fs)' % (nu, sm['eps'].str(3), sm['epsX'].str(3), time.time() - T0))
    th0 = -nu_s/20
    smc = sm['smc']
    c_thin = [PM.horner([val(a) for a in smc['co']['Y'][i]], th0) for i in range(3)] + [PM.horner([val(a) for a in smc['co']['X']], th0)]
    d_box = [PM.horner([a.d[0] for a in sm['co']['Y'][i]], th0) for i in range(3)] + [PM.horner([a.d[0] for a in sm['co']['X']], th0)]
    tails = [sm['eps']]*3 + [sm['epsX']]
    c_thin = [c_thin[i] + sym(tails[i]) for i in range(4)]
    S4 = V3.aligned_set_parts(c_thin, d_box, V0.rad())
    T4 = V2.tau_tape(k)
    steps = []
    qa = arb(Q)
    def record(kind, Sset, F, ha, ncomp):
        c = [arb(t.mid()) for t in Sset.c]
        y1, Bx, Yc = TY.step_cauchy2(F, c, ha, KT, q=qa, conv=lambda cc: PM._c(cc, arb(0)))
        nrv = [b.rad() for b in Bx]
        setrad = max(t.rad() for t in Sset.hull())
        steps.append(dict(kind=kind, h=ser(ha), c=[ser(t) for t in c], Yc=[[ser(a) for a in Yc[i][:KT + 1]] for i in range(ncomp)],
                          nrv=[ser(t) for t in nrv], Rc=ser(abs(ha)/qa), setrad=ser(setrad)))
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
    cm = Pb.centre_manifold(Dual(W0, [arb(1)]), arb('0.03'))
    log('centre manifold: eps=%s epsX=%s' % (cm['eps'].str(3), cm['epsX'].str(3)))
    out = dict(k=str(k), nu_s=nu, th0=ser(th0), XBAR=ser(XBAR), X_SWITCH=ser(X_SWITCH), KT=KT, Q=Q, prec=PREC,
               enc_src=e['src'], repairs='T1 E5 enclosure, T2 R7, T3 R8',
               sm=dict(Y=[[ser(val(a)) for a in sm['co']['Y'][i]] for i in range(3)], X=[ser(val(a)) for a in sm['co']['X']],
                       eps=ser(sm['eps']), epsX=ser(sm['epsX']), lam=ser(val(sm['lam'])), y0=[ser(val(t)) for t in sm['y0']]),
               cm=dict(Y=[[ser(val(a)) for a in cm['co']['Y'][i]] for i in range(3)], X=[ser(val(a)) for a in cm['co']['X']],
                       eps=ser(cm['eps']), epsX=ser(cm['epsX']), lam=ser(val(cm['lam'])), nu='0.03'),
               xc=ser(xc), W0=ser(W0), V0=ser(V0), steps=steps, end=[ser(t) for t in S.hull()])
    fn = '<LAB>/engine/cert/lin_tube2_k%s.pkl' % str(k).replace('/', '_')
    pickle.dump(out, open(fn, 'wb'))
    log('saved %s: %d steps (%.0fs)' % (fn, len(steps), time.time() - T0))
    return out

if __name__ == '__main__':
    T00 = time.time()
    build(sp.Rational(sys.argv[1]), log=lambda *a: print('[%7.1fs]' % (time.time() - T00), *a, flush=True))
