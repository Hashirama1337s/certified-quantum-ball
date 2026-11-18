"""BAR-LB-E2 (v2) driver: new module on top of the pinned v1 code.

Changes against v1 (see PROPOSAL-V2.md):
 V2-1  small sonic manifold: K_s small (default 16), nu_s from a sealed ladder, defect computed
       exactly to N2 = 12 K (polynomial composition) + Cauchy tail on R2 = 2 nu.  Box-valued
       coefficients stay finite (v1 failed with nan at K_s = 60 over a box).
 V2-2  mean-value start set over the V0 box:  y_s(V0) in y_s(V0~) + y_s'(X)(V0 - V0~).
 V2-3  tau-time leg: from the manifold point to x_switch with the desingularised polynomial field
       dy/dtau = G(y), dx/dtau = det M(y) (no singularity at the sonic point), Lohner/QR in 4D with
       tangents; the last tau-step lands on x_switch by Newton on the centre jets; residual x-gap
       transported in x-time with a second-order remainder.
 V2-4  x-time Lohner x_switch -> xbar, matching as v1 (driver2.Matcher3 logic).
 V2-5  AMEND-LB-E-2 items: radius ladder, PE2 bisection (both legs), null1 om guess.
"""
import sys, os, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import mpmath
from flint import arb, arb_mat
from core import Dual, sym, amax, sym_const, Tape
import param as PM
from param import val
import lohner as LO
import taylor as TY
import driver as DR
import driver2 as D2

def tau_tape(k):
    import system as S0
    Ms, bs = S0.Msym, S0.bsym
    det = Ms[0, 0]*Ms[1, 1] - Ms[0, 1]*Ms[1, 0]
    X_ = sp.Symbol('X')
    return Tape([det*S0.gN, Ms[1, 1]*bs[0] - Ms[0, 1]*bs[1], Ms[0, 0]*bs[1] - Ms[1, 0]*bs[0], det],
                [S0.N, S0.om, S0.V, X_], {S0.k: sp.Rational(k)})

def sonic_manifold_v2(Pb, V0, om_guess, nu, K, tails=True, N2fac=12):
    y0 = Pb.sonic_point(V0, om_guess)
    lam, v, DG0, others = Pb.eigen(y0)
    Jg, _ = PM.jac(Pb.Ss, y0)
    X1 = sum((Jg[0][j]*v[j] for j in range(3)), 0*lam)/lam
    v = [x/X1 for x in v]
    co = PM.coefficients(Pb.Gs, Pb.Ss, y0, lam, v, K)
    out = dict(y0=y0, lam=lam, v=v, co=co, others=others)
    if tails:
        Yc = [[val(a) for a in co['Y'][i]] for i in range(3)]
        Xc = [val(a) for a in co['X']]
        eps, epsX, info = PM.tail_bound(Pb.Gs, Pb.Ss, Yc, Xc, val(lam), co['DG0'], K, arb(nu), N2=N2fac*K)
        out.update(eps=eps, epsX=epsX, tinfo=info)
    return out

def choose_nu(Pb, V0box, om_guess, ladder, K, log=None):
    """first nu in the sealed ladder whose tail bound contracts over the V0 box (finite, no nan)."""
    for nu in ladder:
        try:
            sm = sonic_manifold_v2(Pb, V0box, om_guess, nu, K)
            if sm['eps'].is_finite() and sm['epsX'].is_finite():
                if log: log('nu_s=%s contracts over the box: eps=%s epsX=%s L*=%s' % (nu, sm['eps'].str(3), sm['epsX'].str(3), sm['tinfo']['L'].str(3)))
                return nu, sm
        except (ArithmeticError, AssertionError) as e:
            if log: log('nu_s=%s: %s' % (nu, str(e)[:100]))
    raise ArithmeticError('no nu_s in the ladder contracts')

def land_tau(T, c4, x_target, K):
    """tau-step h (arb, exact) landing the thin centre on x = x_target (Newton on centre jets)."""
    Yc = TY.jets(T, [arb(v.mid()) for v in c4], K, lambda q: PM._c(q, arb(0)))
    cx = [float(a.mid()) for a in Yc[3][:K + 1]]
    g = float(Yc[3][1].mid())
    h = (float(x_target.mid()) - cx[0])/g
    for it in range(30):
        p = sum(cx[n]*h**n for n in range(K + 1)) - float(x_target.mid())
        dp = sum(n*cx[n]*h**(n - 1) for n in range(1, K + 1))
        h -= p/dp
    return arb(repr(h))

def tau_leg(T, S, x_target, K=30, q='0.25', h0=0.02, grow=1.25, log=None):
    """Lohner in tau on (N, om, V, x) until x reaches x_target (direction from sign of det M)."""
    qa = arb(q); T0 = time.time(); rec = []; n = 0
    sgn = None; h = h0
    xt = float(x_target.mid())
    while True:
        xc = float(S.c[3].mid())
        rate = float(T.eval([arb(v.mid()) for v in S.c], conv=lambda q_: PM._c(q_, arb(0)))[3].mid())
        if sgn is None:
            sgn = 1 if (xt - xc)*rate > 0 else -1
        remaining = (xt - xc)/rate*sgn       # tau needed (positive) at current rate
        last = remaining <= h
        ha = land_tau(T, S.c, x_target, K) if last else arb(repr(sgn*h))
        start = S.hull()
        try:
            Bx, J = LO.lohner_step(T, S, ha, K, qa)
        except ArithmeticError:
            h = h/2
            if h < 1e-8: raise ArithmeticError('tau-leg collapse at x=%s' % xc)
            continue
        rec.append((xc, ha, Bx, start))
        n += 1
        if last:
            break
        h = h*grow
        if log and n % 50 == 0:
            log('  tau step %d x=%.5f h=%.4g width=%s t=%.1fs' % (n, float(S.c[3].mid()), float(ha.mid()), max(v.rad() for v in S.hull()).str(3), time.time() - T0))
    return dict(steps=n, rec=rec, time=time.time() - T0)

class MatcherV2(D2.Matcher3):
    def __init__(self, Pb, nu_s, nu_c, Ks=16, x_switch='-0.05', Kt=36, **kw):
        super().__init__(Pb, nu_s, nu_c, Kt=Kt, **kw)
        self.Ks2 = Ks; self.x_switch = arb(x_switch)
        self.T = tau_tape(Pb.k)

    def sonic_start(self, V0, om_guess, mv=True):
        """start point (Dual in V0) on the sonic manifold at theta1, 4D incl. x; tails; mean-value over box."""
        Pb = self.Pb
        V0d = Dual(V0, [arb(1)])
        sm = sonic_manifold_v2(Pb, V0d, om_guess, self.nu_s, self.Ks2)
        Y = [PM.horner(sm['co']['Y'][i], self.th1 + 0*V0d) for i in range(3)]
        X = PM.horner(sm['co']['X'], self.th1 + 0*V0d)
        st = Y + [X]
        if mv and V0.rad() > 0:
            Vc = arb(V0.mid())
            smc = sonic_manifold_v2(Pb, Dual(Vc, [arb(1)]), om_guess, self.nu_s, self.Ks2, tails=False)
            stc = [PM.horner(smc['co']['Y'][i], self.th1 + 0*V0d) for i in range(3)] + [PM.horner(smc['co']['X'], self.th1 + 0*V0d)]
            dV = V0 - Vc
            st = [Dual(stc[i].v + st[i].d[0]*dV, [st[i].d[0]]) for i in range(4)]
        return st, sm

    def residual(self, V0, W0, xc, om_guess, log=None, mvstart=True):
        Pb = self.Pb; T0 = time.time(); n = 3
        cvd = lambda like: (lambda c: PM._c(c, like))
        st, sm = self.sonic_start(V0, om_guess, mv=mvstart)
        tail4 = [sm['eps']]*3 + [sm['epsX']]
        # tau leg (4D Lohner with tangents: d/dV0 and unit vectors)
        S4 = LO.LohnerSet([val(v) for v in st],
                          vecs=[[v.d[0] for v in st]] + [[arb(1) if i == j else arb(0) for i in range(4)] for j in range(4)])
        tl = tau_leg(self.T, S4, self.x_switch, K=self.Kt, log=log)
        Y4 = S4.hull(); V4 = S4.vecs()
        Phi_t = [[V4[1 + j][i] for j in range(4)] for i in range(4)]
        e_t = [sum((Phi_t[i][j].abs_upper()*tail4[j] for j in range(4)), arb(0)) for i in range(4)]
        # transport the x-gap to the fixed abscissa x0 = x_switch (x-time, first order + remainder)
        self.x0 = self.x_switch
        y4 = [Dual(Y4[i], [V4[0][i]]) for i in range(4)]
        ys = y4[:3]; xs = y4[3]
        delta = self.x0 - xs
        fs = Pb.F.eval(ys, conv=cvd(ys[0]))
        y0 = [ys[i] + delta*fs[i] for i in range(3)]
        dmax = val(delta).abs_upper() + e_t[3]
        _, fB, half_acc = DR.transport_bound(Pb, [val(v) for v in ys], dmax, amax(e_t[:3]))
        e0 = amax(e_t[:3]) + e_t[3]*DR.fnorm(fB) + dmax**2*half_acc
        # centre side
        W0d = Dual(W0, [arb(1)])
        cm = Pb.centre_manifold(W0d, self.nu_c)
        Zc = [DR.emb(PM.horner(cm['co']['Y'][i], self.thm + 0*W0d), n, 1) for i in range(3)]
        chi = DR.emb(PM.horner(cm['co']['X'], self.thm + 0*W0d), n, 1)
        xcd = Dual(xc, [arb(0), arb(0), arb(1)])
        xm = self.thm.log() + chi + xcd
        if not hasattr(self, 'xbar'):
            self.xbar = arb(val(xm).mid())
            if os.environ.get('LB_SMOKE'):
                self.xbar = self.x0 - arb(os.environ['LB_SMOKE'])
        S = LO.LohnerSet([val(v) for v in y0],
                         vecs=[[v.d[0] for v in y0]] + [[arb(1) if i == j else arb(0) for i in range(3)] for j in range(3)])
        info = D2.integrate_rec(Pb.F, S, self.x0, self.xbar, K=self.Kt, q=self.q, h0=self.h0, log=log)
        Y = S.hull(); V = S.vecs()
        Phi = [[V[1 + j][i] for j in range(3)] for i in range(3)]
        e_end = amax([sum((Phi[i][j].abs_upper() for j in range(3)), arb(0)) for i in range(3)])*e0
        y_end = [Dual(Y[i], [V[0][i], arb(0), arb(0)]) for i in range(3)]
        dm = xm - self.xbar
        if os.environ.get('LB_SMOKE'):
            dm = dm - arb(val(dm).mid())
        fe = Pb.F.eval(y_end, conv=cvd(y_end[0]))
        dmx = val(dm).abs_upper() + cm['epsX']
        Bt, fBt, half_acc2 = DR.transport_bound(Pb, Y, dmx, e_end)
        ym = [y_end[i] + dm*fe[i] for i in range(3)]
        N, om, Vv = ym
        Zs = [1/N, om*N*N, Vv*N]
        R = [Zs[i] - Zc[i] for i in range(3)]
        Nb, ob, Vb = Bt
        JZ = amax([(1/(Nb*Nb)).abs_upper(), (2*ob*Nb).abs_upper() + (Nb*Nb).abs_upper(), Vb.abs_upper() + Nb.abs_upper()])
        eend = e_end + dmx**2*half_acc2 + DR.fnorm(fBt)*cm['epsX']
        Erad = JZ*eend + cm['eps']
        diag = dict(steps=info['steps'] + tl['steps'], tau_steps=tl['steps'], x_steps=info['steps'],
                    time=time.time() - T0, rec=info['rec'], tau_rec=tl['rec'],
                    boxes=[(x, float(h.mid()), Bx) for (x, h, Bx, s_) in info['rec']],
                    e0=e0, e_end=e_end, Erad=Erad, sm=sm, cm=cm, width_end=max(v.rad() for v in Y),
                    width_start=max(val(v).rad() for v in st), width_switch=max(v.rad() for v in Y4), Y=Y)
        return R, Erad, diag

# ---------------- PE2 on the tau leg ----------------
def classify_tau(T, Bx):
    Bv = [val(t) for t in Bx]
    if Bv[2] > 0: return '+'
    if Bv[2] < 0: return '-'
    try:
        g = T.eval(Bv, conv=sym_const)
        if (g[3] > 0 or g[3] < 0) and (g[2] > 0 or g[2] < 0): return 'mono'
    except Exception:
        pass
    return 'amb'

def pe2_tau_bisect(T, rec, K=30, q='0.25', depth=12):
    qa = arb(q)
    def sub(start, h, d):
        S = LO.LohnerSet(start); out = []
        for half in range(2):
            st = S.hull()
            Bx, _ = LO.lohner_step(T, S, h/2, K, qa)
            c = classify_tau(T, Bx)
            out += sub(st, h/2, d + 1) if (c == 'amb' and d < depth) else [c]
        return out
    seq = []; nb = 0
    for (x, ha, Bx, start) in rec:
        c = classify_tau(T, Bx)
        if c == 'amb':
            nb += 1
            try: seq += sub(start, ha, 1)
            except ArithmeticError: seq.append('amb')
        else:
            seq.append(c)
    return seq, nb

def pe2_full(Pb, M, d, K=36):
    """V sign structure over: sonic piece (theta in [th1,0]), tau leg, x leg, centre piece."""
    seq_t, nb_t = pe2_tau_bisect(M.T, d['tau_rec'], K=M.Kt)
    zx = D2.pe2_bisect(Pb, d['rec'], K=K)
    sm, cm = d['sm'], d['cm']
    Vs = PM.horner([val(a) for a in sm['co']['Y'][2]], M.th1/2 + sym(M.th1/2)) + sym(sm['eps'])
    Uc = PM.horner([val(a) for a in cm['co']['Y'][2]], M.thm/2 + sym(M.thm/2)) + sym(cm['eps'])
    # order along decreasing x: sonic piece, tau leg, x leg, centre piece
    seq = (['+'] if Vs > 0 else ['-'] if Vs < 0 else ['amb']) + seq_t
    cert_t = [s for s in seq if s in ('+', '-')]
    namb = seq.count('amb') + zx['amb']
    signs = cert_t + ([zx['first']] if zx['first'] else []) + ([zx['last']] if zx['last'] else []) + (['-'] if Uc < 0 else ['amb'])
    # total sign changes: tau-leg changes + x-leg changes + junction changes
    ch_t = sum(1 for a, b in zip(cert_t, cert_t[1:]) if a != b)
    junction = 1 if (cert_t and zx['first'] and cert_t[-1] != zx['first']) else 0
    end = 0 if (zx['last'] == '-' and Uc < 0) else (1 if zx['last'] == '+' and Uc < 0 else 0)
    total = ch_t + junction + zx['sign_changes'] + end
    return dict(total_sign_changes=total, amb=namb, tau_changes=ch_t, x_changes=zx['sign_changes'],
                tau_bisected=nb_t, x_bisected=zx['bisected_steps'], holds=(namb == 0 and total == 1 and bool(Uc < 0)))

def check_pred_bisect(T, rec, pred, K=30, q='0.25', depth=12):
    """pred(box) must hold on every step's a priori box; failing steps are bisected (Lohner restarted
    from the recorded step-start enclosure) up to depth."""
    qa = arb(q)
    def ok_sub(start, h, d):
        S = LO.LohnerSet(start)
        for half in range(2):
            st = S.hull()
            Bx, _ = LO.lohner_step(T, S, h/2, K, qa)
            if not pred([val(t) for t in Bx]):
                if d >= depth or not ok_sub(st, h/2, d + 1):
                    return False
        return True
    for (x, ha, Bx, start) in rec:
        if not pred([val(t) for t in Bx]):
            try:
                if not ok_sub(start, ha, 1): return False
            except ArithmeticError:
                return False
    return True

# ---------------- outward leg v2 ----------------
def outward_leg_v2(Pb, M, V0box, om_guess, x_switch_out='0.05', X1=12.0, Dw='1e-3', log=None):
    """theta>0 manifold point -> tau leg to +x_switch -> x-time Lohner to X1 -> far-side lemma (as v1)."""
    T0 = time.time()
    k_ = sp.Symbol('k'); N_, om_, V_ = sp.symbols('N om V')
    from pipeline import load_expr, FIELDS
    kk = {k_: Pb.k}
    ro = Tape([load_expr(FIELDS['ro'])], [N_, om_, V_], kk)
    rV = Tape([load_expr(FIELDS['rV'])], [N_, om_, V_], kk)
    aT = Tape([load_expr(FIELDS['Acons']) - 2 - (1 - k_)*om_], [N_, om_, V_], kk)
    sm = sonic_manifold_v2(Pb, V0box, om_guess, M.nu_s, M.Ks2)
    th = -M.th1
    st = [PM.horner([val(a) for a in sm['co']['Y'][i]], th) + sym(sm['eps']) for i in range(3)]
    st.append(PM.horner([val(a) for a in sm['co']['X']], th) + sym(sm['epsX']))
    S4 = LO.LohnerSet(st)
    xs = arb(x_switch_out)
    tl = tau_leg(M.T, S4, xs, K=M.Kt, log=log)
    Y4 = S4.hull()
    predD = lambda Bv: bool(Pb.Dson.eval(Bv[:3], conv=sym_const)[0] < 0)
    predV = lambda Bv: bool(Bv[2] > 0)
    okD = check_pred_bisect(M.T, tl['rec'], predD, K=M.Kt)
    okV = check_pred_bisect(M.T, tl['rec'], predV, K=M.Kt)
    # exact x: transport the x-gap (first order + remainder) and continue in x-time
    delta = xs - Y4[3]
    B, fB, half = DR.transport_bound(Pb, Y4[:3], delta.abs_upper())
    y = [Y4[i] + delta*fB[i] for i in range(3)]
    S = LO.LohnerSet(y)
    info = D2.integrate_rec(Pb.F, S, xs, arb(X1), K=M.Kt, log=log)
    okD &= check_pred_bisect(Pb.F, info['rec'], predD, K=M.Kt)
    okV &= check_pred_bisect(Pb.F, info['rec'], predV, K=M.Kt)
    Vpiece = PM.horner([val(a) for a in sm['co']['Y'][2]], th/2 + sym(th/2)) + sym(sm['eps'])
    okV &= bool(Vpiece > 0)
    # sonic-manifold piece theta in (0, th]: D(Y(0)) = 0 and d/dtheta D(Y(theta)) < 0 on [0, th]
    # (Y = Y_K + E, |E| <= eps on |theta| <= nu, |E'| <= eps/(nu - th) by Cauchy)
    tI = Dual(th/2 + sym(th/2), [arb(1)])
    epsd = sm['eps']/(arb(M.nu_s) - th)
    Yd = [PM.horner([val(a) for a in sm['co']['Y'][i]], tI) for i in range(3)]
    Yd = [Dual(y.v + sym(sm['eps']), [y.d[0] + sym(epsd)]) for y in Yd]
    Dd = Pb.Dson.eval(Yd, conv=lambda c: Dual(sym_const(c), [arb(0)]))[0]
    okD_piece = bool(Dd.d[0] < 0)
    okD &= okD_piece
    N1, om1, V1 = S.hull()
    Dw = arb(Dw)
    Q = [N1.upper()/2 + sym(N1.upper()/2 + arb(10)**(-40)), om1 + sym(Dw), V1 + sym(Dw)]
    a = aT.eval(Q, conv=sym_const)[0]
    res = dict(steps=tl['steps'] + info['steps'], ok_D=okD, ok_V=okV, ok_D_manifold_piece=okD_piece, a=a, y1=S.hull())
    if not (a < 0):
        res.update(lemma=False, time=time.time() - T0); return res
    mu = -arb(a.upper())
    co_ = ro.eval(Q, conv=sym_const)[0].abs_upper(); cV = rV.eval(Q, conv=sym_const)[0].abs_upper()
    DQ = Pb.Dson.eval(Q, conv=sym_const)[0]
    dom = co_*N1.upper()/mu; dV = cV*N1.upper()/mu
    lemma = bool(dom < Dw) and bool(dV < Dw) and bool(DQ < 0) and bool(Q[2] > 0)
    om_inf = om1 + sym(dom); V_inf = V1 + sym(dV)
    A_inf = DR.centred_eval(Pb.Acons, [arb(0), om_inf, V_inf])
    res.update(lemma=lemma, mu=mu, om_inf=om_inf, V_inf=V_inf, A_inf=A_inf, m_over_r=(1 - 1/A_inf)/2, time=time.time() - T0)
    return res


# ---------------- V2-6: tighter interval inflation in the sonic point and eigenpair ----------------
# v1 inflated the om0 interval-Newton start by x400 and the eigen Krawczyk box by x1000 of the DG
# width, which turned a 1e-14 V0 box into ~1e-9 manifold widths. Problem2 starts small and doubles
# only on failure (the Krawczyk/interval-Newton tests themselves are unchanged).
from pipeline import Problem
from flint import ctx

class Problem2(Problem):
    def sonic_point(self, V0, om_guess):
        c = self.c
        N0 = (1 - V0*c)/(c - V0)
        V0v, N0v = val(V0), val(N0)
        om_t = arb(om_guess)
        Nm, Vm = arb(N0v.mid()), arb(V0v.mid())
        for it in range(12):
            Fm = self.LB.eval([Nm, om_t, Vm])[0]
            J, _ = PM.jac(self.LB, [Nm, om_t, Vm])
            om_t = arb((om_t - Fm/J[0][1]).mid())
        r = arb(2)**(-ctx.prec + 20) + 4*(N0v.rad() + V0v.rad())
        for it in range(60):
            Iom = om_t + sym(r)
            Fm = self.LB.eval([N0v, om_t, V0v])[0]
            J, _ = PM.jac(self.LB, [N0v, Iom, V0v])
            newI = om_t - Fm/J[0][1]
            if Iom.contains_interior(newI):
                om_v = newI; break
            r = r*2
        else:
            raise ArithmeticError('sonic om0 interval Newton failed')
        if isinstance(V0, Dual):
            Jp, _ = PM.jac(self.LB, [N0v, om_v, V0v])
            om0 = Dual(om_v, [-(Jp[0][0]*N0.d[0] + Jp[0][2])/Jp[0][1]])
        else:
            om0 = om_v
        return [N0, om0, V0]

    def eigen(self, y0, branch='fast'):
        DG, _ = PM.jac(self.Gs, y0)
        DGv = [[val(DG[i][j]) for j in range(3)] for i in range(3)]
        Jm = mpmath.matrix([[float(DGv[i][j].mid()) for j in range(3)] for i in range(3)])
        ev, EV = mpmath.eig(Jm)
        ev = [complex(e) for e in ev]
        order = sorted(range(3), key=lambda i: -abs(ev[i]))
        ia = order[0]; lam_t = ev[ia].real
        vt = [float(mpmath.re(EV[i, ia])) for i in range(3)]
        jn = max(range(3), key=lambda i: abs(vt[i])); vt = [x/vt[jn] for x in vt]
        free = [i for i in range(3) if i != jn]
        def Fv(u, DGm):
            v = [None]*3; v[jn] = 1 + 0*u[0]; v[free[0]], v[free[1]] = u[0], u[1]
            return [sum((DGm[i][j]*v[j] for j in range(3)), 0*u[0]) - u[2]*v[i] for i in range(3)], v
        def JF(u, DGm):
            v = [None]*3; v[jn] = 1 + 0*u[0]; v[free[0]], v[free[1]] = u[0], u[1]
            return [[DGm[i][free[0]] - (u[2] if i == free[0] else 0), DGm[i][free[1]] - (u[2] if i == free[1] else 0), -v[i]]
                    for i in range(3)]
        ut = [arb(vt[free[0]]), arb(vt[free[1]]), arb(lam_t)]
        for it in range(8):
            Fm, _ = Fv(ut, DGv); Jm_ = JF(ut, DGv)
            du = PM.solve3([[arb(x.mid()) for x in row] for row in Jm_], [arb(x.mid()) for x in Fm])
            ut = [arb((ut[i] - du[i]).mid()) for i in range(3)]
        Fm, _ = Fv(ut, DGv)
        Cm = arb_mat([[arb(x.mid()) for x in row] for row in JF(ut, DGv)]).inv()
        CF = [sum((Cm[i, j]*Fm[j] for j in range(3)), arb(0)) for i in range(3)]
        rad = arb(2)**(-ctx.prec + 30) + 8*amax([DGv[i][j].rad() for i in range(3) for j in range(3)])
        for _try in range(40):
            U = [x + sym(rad) for x in ut]
            JU = JF(U, DGv)
            Kr = []
            for i in range(3):
                s = ut[i] - CF[i]
                for j in range(3):
                    cij = (1 if i == j else 0) - sum((Cm[i, l]*JU[l][j] for l in range(3)), arb(0))
                    s += cij*(U[j] - ut[j])
                Kr.append(s)
            if all(U[i].contains_interior(Kr[i]) for i in range(3)):
                break
            rad = rad*2
        else:
            raise ArithmeticError('eigen Krawczyk failed')
        uv = Kr
        if isinstance(y0[0], Dual):
            dDG = [[DG[i][j].d[0] if isinstance(DG[i][j], Dual) else arb(0) for j in range(3)] for i in range(3)]
            _, v = Fv(uv, DGv)
            rhs = [-sum((dDG[i][j]*v[j] for j in range(3)), arb(0)) for i in range(3)]
            du = PM.solve3(JF(uv, DGv), rhs)
            u = [Dual(uv[i], [du[i]]) for i in range(3)]
        else:
            u = uv
        v = [None]*3; one = 1 + 0*u[0]
        v[jn] = one; v[free[0]], v[free[1]] = u[0], u[1]
        return u[2], v, DG, [ev[i] for i in order[1:]]
