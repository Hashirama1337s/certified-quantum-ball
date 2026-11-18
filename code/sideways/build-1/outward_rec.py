"""outward_rec.py - build seat build-1 (B7), 2026-11-17.

Certified E5 outward leg (sonic point, theta > 0 side, to t = 0) re-run with STEP RECORDS.
Control flow = engine/cert/driver_v2.py outward_leg_v2 (l.262-317) as called by the sealed instrument
engine/cert/run_v32.py (target mode, l.96-101): same Problem2 / MatcherV32 parameters (PREC 200, K_s 16,
K_c 60, K_t 36, nu_c 0.03, nu_s ladder 0.004/0.002/0.001 at box radius 1e-9), same tau leg to x = +0.05,
same x-time Lohner integrate_rec to X1 = 12, same predicate checks (D_sonic < 0, V > 0) with bisection,
same far-side invariant-region lemma.  The only additions are records:

  for every accepted step: kind ('tau' or 'x'), h, the step-start enclosure hull (all certified solutions
  start the step inside it), ball Taylor jets Yc[i][0..KT] of ALL solutions starting in the hull
  (taylor.step_cauchy2 from the hull), the a priori box Bx and the Cauchy radius Rc = |h|/q.
  For t in [0, h] (complex |t| <= Rc):  y_i(t) in sum_n Yc[i][n] t^n  +/- rad(Bx_i) (|t|/Rc)^(KT+1)/(1-|t|/Rc).

V0 box = runs6/target.log PE1 K box at r = 1e-15 ([0.1124394013880983362 +/- 3.31e-20]), which contains the
sealed run's K ∩ X.  Records are strings only (arb.mid().str / arb.rad().str), never flint objects.
Failure-mode-3 control: the t = 0 m/r must overlap runs6/target.log PE3 [0.0595 +/- 8.24e-5].
usage (memcap 8):  py -3.11 outward_rec.py [X1]      -> outward_rec.pkl (X1 default 12 = sealed)
"""
import sys, os, time, pickle, hashlib
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
for f in ['LB_SMOKE', 'LB_EXPANDED', 'LB_MVSTART', 'LB_LADDER', 'LB_V3_LADDER', 'LB_NULADDER', 'LB_NU_AT', 'LB_X1',
          'LB_LOHNER_MODE', 'LB_V3ALIGN', 'LB_PREC', 'LB_KS', 'LB_KC', 'LB_KT', 'LB_NUC']:
    if os.environ.get(f) is not None:
        raise SystemExit('outward_rec: environment flag %s must be unset (as run_v32)' % f)
if os.environ.get('LB_NGRID', '24') != '24':
    raise SystemExit('outward_rec: LB_NGRID must be 24')
os.environ['LB_NGRID'] = '24'
import sympy as sp
from flint import arb, ctx
ctx.prec = 200
ctx.cap = 32
import driver_v31 as V31
import driver_v32 as V32
import driver_v3 as V3
import driver_v2 as V2
import driver2 as D2
import driver as DR
import lohner as LO
import param as PM
import taylor as TY
from param import val
from core import sym, sym_const, Tape
assert LO.lohner_step is V3._step and V3.MODE == 'qrw', 'weighted-QR Lohner not installed'
assert PM.arb_series is V31.arb_series_exact, 'R7 not installed'

PREC, KS, KC, KT, NUC = 200, 16, 60, 36, '0.03'
NU_LADDER = ['0.004', '0.002', '0.001']
NU_AT = arb('1e-9')
V0_MID, V0_RAD = '0.1124394013880983362', '3.31e-20'
OM_T = '0.3256888894230684'
PE3 = ('0.0595', '8.24e-5')
HERE = os.path.dirname(os.path.abspath(__file__))


def ser(a):
    return (a.mid().str(80, radius=False), a.rad().str(10, radius=False))


def record(kind, field, hull, ha, steps, q):
    y1, Bx, Yc = TY.step_cauchy2(field, list(hull), ha, KT, q=q, conv=lambda c: PM._c(c, arb(0)))
    steps.append(dict(kind=kind, h=ser(ha), hull=[ser(t) for t in hull], Yc=[[ser(a) for a in Yc[i][:KT + 1]] for i in range(len(hull))],
                      Bx=[ser(t) for t in Bx], Rc=ser(abs(ha)/q)))


def tau_leg_rec(T, S, x_target, K, steps, q='0.25', h0=0.02, grow=1.25, log=None):
    """driver_v2.tau_leg (l.75-102) + records"""
    qa = arb(q); T0 = time.time(); rec = []; n = 0
    sgn = None; h = h0
    xt = float(x_target.mid())
    while True:
        xc = float(S.c[3].mid())
        rate = float(T.eval([arb(v.mid()) for v in S.c], conv=lambda q_: PM._c(q_, arb(0)))[3].mid())
        if sgn is None:
            sgn = 1 if (xt - xc)*rate > 0 else -1
        remaining = (xt - xc)/rate*sgn
        last = remaining <= h
        ha = V2.land_tau(T, S.c, x_target, K) if last else arb(repr(sgn*h))
        start = S.hull()
        try:
            Bx, J = LO.lohner_step(T, S, ha, K, qa)
        except ArithmeticError:
            h = h/2
            if h < 1e-8: raise ArithmeticError('tau-leg collapse at x=%s' % xc)
            continue
        record('tau', T, start, ha, steps, qa)
        rec.append((xc, ha, Bx, start))
        n += 1
        if last:
            break
        h = h*grow
        if log and n % 50 == 0:
            log('  tau step %d x=%.5f h=%.4g width=%s t=%.1fs' % (n, float(S.c[3].mid()), float(ha.mid()), max(v.rad() for v in S.hull()).str(3), time.time() - T0))
    return dict(steps=n, rec=rec, time=time.time() - T0)


def integrate_rec2(F, S, x0, x1, K, steps, q='0.25', h0=0.002, grow=1.25, log=None):
    """driver2.integrate_rec (l.26-50) + records"""
    T0 = time.time(); qa = arb(q)
    x = arb(x0); x1 = arb(x1)
    sgn = 1 if float(x1.mid()) > float(x.mid()) else -1
    h = h0; n = 0; rec = []
    while (float(x1.mid()) - float(x.mid()))*sgn > 0:
        rem = abs(float((x1 - x).mid()))
        hh = min(h, rem); last = hh >= rem - 1e-15
        ha = (x1 - x) if last else arb(str(sgn*hh))
        start = S.hull()
        saved = (list(S.c), [row[:] for row in S.B], list(S.r), [v[:] for v in S.vh], [t[:] for t in S.t])
        try:
            Bx, J = LO.lohner_step(F, S, ha, K, qa)
        except ArithmeticError:
            S.c, S.B, S.r, S.vh, S.t = saved
            h = hh/2
            if h < 1e-7: raise ArithmeticError('Lohner step collapse at x=%s' % x)
            continue
        record('x', F, start, ha, steps, qa)
        steps[-1]['x0'] = ser(x)
        rec.append((float(x.mid()), ha, Bx, start))
        x = x1 if last else x + ha
        n += 1; h = hh*grow
        if log and n % 50 == 0:
            log('  lohner step %d x=%.4f h=%.4g width=%s t=%.1fs' % (n, float(x.mid()), hh, max(v.rad() for v in S.hull()).str(3), time.time() - T0))
    return dict(steps=n, rec=rec, time=time.time() - T0, x=x)


def outward_rec(Pb, M, V0box, om_guess, x_switch_out='0.05', X1=12.0, Dw='1e-3', log=None):
    """driver_v2.outward_leg_v2 (l.262-317) + records"""
    T0 = time.time()
    k_ = sp.Symbol('k'); N_, om_, V_ = sp.symbols('N om V')
    from pipeline import load_expr, FIELDS
    kk = {k_: Pb.k}
    ro = Tape([load_expr(FIELDS['ro'])], [N_, om_, V_], kk)
    rV = Tape([load_expr(FIELDS['rV'])], [N_, om_, V_], kk)
    aT = Tape([load_expr(FIELDS['Acons']) - 2 - (1 - k_)*om_], [N_, om_, V_], kk)
    sm = V2.sonic_manifold_v2(Pb, V0box, om_guess, M.nu_s, M.Ks2)
    th = -M.th1
    st = [PM.horner([val(a) for a in sm['co']['Y'][i]], th) + sym(sm['eps']) for i in range(3)]
    st.append(PM.horner([val(a) for a in sm['co']['X']], th) + sym(sm['epsX']))
    S4 = LO.LohnerSet(st)
    xs = arb(x_switch_out)
    steps = []
    tl = tau_leg_rec(M.T, S4, xs, M.Kt, steps, log=log)
    Y4 = S4.hull()
    predD = lambda Bv: bool(Pb.Dson.eval(Bv[:3], conv=sym_const)[0] < 0)
    predV = lambda Bv: bool(Bv[2] > 0)
    okD = V2.check_pred_bisect(M.T, tl['rec'], predD, K=M.Kt)
    okV = V2.check_pred_bisect(M.T, tl['rec'], predV, K=M.Kt)
    delta = xs - Y4[3]
    B, fB, half = DR.transport_bound(Pb, Y4[:3], delta.abs_upper())
    y = [Y4[i] + delta*fB[i] for i in range(3)]
    S = LO.LohnerSet(y)
    log('tau leg %d steps; x-gap transport to x = %s (%.0fs)' % (tl['steps'], xs.str(5), time.time() - T0))
    info = integrate_rec2(Pb.F, S, xs, arb(X1), M.Kt, steps, log=log)
    okD &= V2.check_pred_bisect(Pb.F, info['rec'], predD, K=M.Kt)
    okV &= V2.check_pred_bisect(Pb.F, info['rec'], predV, K=M.Kt)
    Vpiece = PM.horner([val(a) for a in sm['co']['Y'][2]], th/2 + sym(th/2)) + sym(sm['eps'])
    okV &= bool(Vpiece > 0)
    from core import Dual
    tI = Dual(th/2 + sym(th/2), [arb(1)])
    epsd = sm['eps']/(arb(M.nu_s) - th)
    Yd = [PM.horner([val(a) for a in sm['co']['Y'][i]], tI) for i in range(3)]
    Yd = [Dual(y_.v + sym(sm['eps']), [y_.d[0] + sym(epsd)]) for y_ in Yd]
    Dd = Pb.Dson.eval(Yd, conv=lambda c: Dual(sym_const(c), [arb(0)]))[0]
    okD_piece = bool(Dd.d[0] < 0)
    okD &= okD_piece
    N1, om1, V1 = S.hull()
    Dw = arb(Dw)
    Q = [N1.upper()/2 + sym(N1.upper()/2 + arb(10)**(-40)), om1 + sym(Dw), V1 + sym(Dw)]
    a = aT.eval(Q, conv=sym_const)[0]
    res = dict(steps=tl['steps'] + info['steps'], tau_steps=tl['steps'], x_steps=info['steps'], ok_D=okD, ok_V=okV,
               ok_D_manifold_piece=okD_piece, a=a, y1=S.hull())
    if not (a < 0):
        res.update(lemma=False, time=time.time() - T0)
    else:
        mu = -arb(a.upper())
        co_ = ro.eval(Q, conv=sym_const)[0].abs_upper(); cV = rV.eval(Q, conv=sym_const)[0].abs_upper()
        DQ = Pb.Dson.eval(Q, conv=sym_const)[0]
        dom = co_*N1.upper()/mu; dV = cV*N1.upper()/mu
        lemma = bool(dom < Dw) and bool(dV < Dw) and bool(DQ < 0) and bool(Q[2] > 0)
        om_inf = om1 + sym(dom); V_inf = V1 + sym(dV)
        A_inf = DR.centred_eval(Pb.Acons, [arb(0), om_inf, V_inf])
        res.update(lemma=lemma, mu=mu, om_inf=om_inf, V_inf=V_inf, A_inf=A_inf, m_over_r=(1 - 1/A_inf)/2, Q=Q, DQ=DQ,
                   time=time.time() - T0)
    res['records'] = steps
    res['sm'] = dict(Y=[[ser(val(a_)) for a_ in sm['co']['Y'][i]] for i in range(3)], X=[ser(val(a_)) for a_ in sm['co']['X']],
                     eps=ser(sm['eps']), epsX=ser(sm['epsX']), th=ser(th), nu_s=str(M.nu_s.mid()))
    return res


if __name__ == '__main__':
    T00 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T00), *a, flush=True)
    X1 = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0
    out_name = 'outward_rec.pkl' if X1 == 12.0 else 'outward_rec_probe.pkl'
    k13 = sp.Rational(1, 3)
    V0s = arb(V0_MID)
    Pb = V31.Problem2(k13, prec=PREC, K_sonic=KS, K_centre=KC, K_ty=KT)
    nu, _ = V2.choose_nu(Pb, V0s + sym(NU_AT), OM_T, NU_LADDER, KS, log=log)
    log('sealed rule: nu_s = %s' % nu)
    M = V32.MatcherV32(Pb, nu, NUC, Ks=KS, x_switch='-0.05', Kt=KT)
    V0box = V0s + sym(arb(V0_RAD))
    r = outward_rec(Pb, M, V0box, OM_T, X1=X1, log=log)
    log('outward: steps=%d (tau %d, x %d) D_sonic<0: %s V>0: %s lemma: %s' % (r['steps'], r['tau_steps'], r['x_steps'], r['ok_D'], r['ok_V'], r.get('lemma')))
    out = dict(k='1/3', V0=(V0_MID, V0_RAD), nu_s=nu, X1=X1, KT=KT, Q='0.25', prec=PREC, steps=r['records'], sm=r['sm'],
               ok_D=r['ok_D'], ok_V=r['ok_V'], ok_D_manifold_piece=r['ok_D_manifold_piece'], end=[ser(t) for t in r['y1']],
               lemma=r.get('lemma'), a=ser(r['a']))
    if r.get('lemma'):
        mr = r['m_over_r']
        lo_, hi_ = arb(PE3[0]) - arb(PE3[1]), arb(PE3[0]) + arb(PE3[1])
        overlap = bool(mr.upper() >= lo_.lower()) and bool(mr.lower() <= hi_.upper())
        out.update(mu=ser(r['mu']), om_inf=ser(r['om_inf']), V_inf=ser(r['V_inf']), A_inf=ser(r['A_inf']), m_over_r=ser(mr),
                   Q=[ser(t) for t in r['Q']], DQ=ser(r['DQ']), PE3_overlap=overlap)
        log('m/r(t=0) =', mr.str(10), ' overlaps runs6 PE3 [0.0595 +/- 8.24e-5]:', overlap)
    out['instrument_sha256'] = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    out['time_s'] = time.time() - T00
    pickle.dump(out, open(os.path.join(HERE, out_name), 'wb'))
    log('saved %s (%d step records)' % (out_name, len(out['steps'])))
