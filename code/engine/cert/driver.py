"""B3 certified driver (Lohner version): matching residual, generalised Krawczyk,
outward leg + far-side lemma, PE2 checks.  Builds on pipeline.Problem (manifolds,
sonic point, eigenpair) and lohner.py (wrapping-controlled integration).

Generalised Krawczyk (existence): with R(x) = Rh(x, tau(x)), tau = tails/transport remainders,
   R(x) in Rh(x~,0) + J_x(X)(x - x~) + E,   E >= |dRh/dtau| |tau|  on X,
   K = x~ - C Rh(x~,0) - C E + (I - C J_x(X))(X - x~)  in int X  ==>  R has a zero in X.
"""
import sys, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, arb_mat, ctx
from core import Dual, sym, amax, sym_const, Tape
import param as PM
from param import val
import lohner as LO
from pipeline import Problem, load_expr, FIELDS, N_, om_, V_, k_

def emb(x, n, i):
    if isinstance(x, Dual):
        d = [arb(0)]*n; d[i] = x.d[0]; return Dual(x.v, d)
    return Dual(x, [arb(0)]*n)

def fnorm(v): return amax([val(t) for t in v])

def centred_eval(tape, box):
    """mean-value form of a scalar tape over a box: g(c) + grad g(box).(box - c)"""
    c = [arb(b.mid()) for b in box]
    gc = tape.eval(c, conv=sym_const)[0]
    J, _ = PM.jac(tape, box)
    return gc + sum((J[0][j]*(box[j] - c[j]) for j in range(3)), arb(0))

def transport_bound(Pb, yv, dmax, extra=arb(0)):
    """box B containing the solution within |dx| <= dmax of yv, and sup|y''|/2 on it"""
    f = Pb.F.eval(yv, conv=sym_const)
    pad = 2*dmax*fnorm(f) + extra + arb(10)**(-60) + amax([v.rad() for v in yv])*arb('1e-3')
    B = [v + sym(pad) for v in yv]
    fB = Pb.F.eval(B, conv=sym_const)
    Jb, _ = PM.jac(Pb.F, B)
    acc = amax([sum((Jb[i][j]*fB[j] for j in range(3)), arb(0)) for i in range(3)])
    # consistency of the box: yv + [-dmax,dmax] f(B) inside B
    for i in range(3):
        if not B[i].contains_interior(yv[i] + sym(dmax*fB[i].abs_upper())):
            raise ArithmeticError('transport box not self-consistent')
    return B, fB, acc/2

class Matcher2:
    def __init__(self, Pb, nu_s, nu_c, frac_s='-0.9', frac_c='0.9', Kt=36, q='0.25', h0=0.002):
        self.Pb = Pb; self.nu_s = arb(nu_s); self.nu_c = arb(nu_c)
        self.th1 = arb(frac_s)*self.nu_s; self.thm = arb(frac_c)*self.nu_c
        self.Kt = Kt; self.q = q; self.h0 = h0
        self.x0 = None

    def residual(self, V0, W0, xc, om_guess, log=None):
        Pb = self.Pb; T0 = time.time(); n = 3
        cvd = lambda like: (lambda c: PM._c(c, like))
        V0d = Dual(V0, [arb(1)])
        sm = Pb.sonic_manifold(V0d, om_guess, self.nu_s)
        Yk, Xk = sm['co']['Y'], sm['co']['X']
        ys = [PM.horner(Yk[i], self.th1 + 0*V0d) for i in range(3)]
        xs = PM.horner(Xk, self.th1 + 0*V0d)
        if self.x0 is None:
            self.x0 = arb(val(xs).mid())
        delta = self.x0 - xs
        fs = Pb.F.eval(ys, conv=cvd(ys[0]))
        y0 = [ys[i] + delta*fs[i] for i in range(3)]
        dmax = val(delta).abs_upper() + sm['epsX']
        _, fB, half_acc = transport_bound(Pb, [val(v) for v in ys], dmax, sm['eps'])
        e0 = sm['eps'] + sm['epsX']*fnorm(fB) + dmax**2*half_acc
        # centre
        W0d = Dual(W0, [arb(1)])
        cm = Pb.centre_manifold(W0d, self.nu_c)
        Zc = [emb(PM.horner(cm['co']['Y'][i], self.thm + 0*W0d), n, 1) for i in range(3)]
        chi = emb(PM.horner(cm['co']['X'], self.thm + 0*W0d), n, 1)
        xcd = Dual(xc, [arb(0), arb(0), arb(1)])
        xm = self.thm.log() + chi + xcd
        if not hasattr(self, 'xbar'):
            import os
            self.xbar = arb(val(xm).mid())
            if os.environ.get('LB_SMOKE'):          # code-path smoke test only (wrong geometry, meaningless R)
                self.xbar = self.x0 - arb(os.environ['LB_SMOKE'])
        # Lohner integration of the set y0 with tangents (d/dV0, e1, e2, e3)
        S = LO.LohnerSet([val(v) for v in y0],
                         vecs=[[v.d[0] for v in y0]] + [[arb(1) if i == j else arb(0) for i in range(3)] for j in range(3)])
        info = LO.integrate(Pb.F, S, self.x0, self.xbar, K=self.Kt, q=self.q, h0=self.h0, log=log)
        Y = S.hull(); V = S.vecs()
        Phi = [[V[1 + j][i] for j in range(3)] for i in range(3)]
        e_end = amax([sum((Phi[i][j].abs_upper() for j in range(3)), arb(0)) for i in range(3)])*e0
        y_end = [Dual(Y[i], [V[0][i], arb(0), arb(0)]) for i in range(3)]
        # transport xbar -> xm
        dm = xm - self.xbar
        import os
        if os.environ.get('LB_SMOKE'):
            dm = dm - arb(val(dm).mid())           # smoke test only: zero the (huge) end transport, keep derivatives
        fe = Pb.F.eval(y_end, conv=cvd(y_end[0]))
        dmx = val(dm).abs_upper() + cm['epsX']
        Bt, fBt, half_acc2 = transport_bound(Pb, Y, dmx, e_end)
        ym = [y_end[i] + dm*fe[i] for i in range(3)]
        N, om, Vv = ym
        Zs = [1/N, om*N*N, Vv*N]
        R = [Zs[i] - Zc[i] for i in range(3)]
        Nb, ob, Vb = Bt
        JZ = amax([(1/(Nb*Nb)).abs_upper(), (2*ob*Nb).abs_upper() + (Nb*Nb).abs_upper(), Vb.abs_upper() + Nb.abs_upper()])
        eend = e_end + dmx**2*half_acc2 + fnorm(fBt)*cm['epsX']
        Erad = JZ*eend + cm['eps']
        diag = dict(steps=info['steps'], time=time.time() - T0, boxes=info['boxes'], e0=e0, e_end=e_end, Erad=Erad,
                    sm=sm, cm=cm, width_end=max(v.rad() for v in Y), Y=Y)
        return R, Erad, diag

def krawczyk(M, Xt, rad, om_guess, log=print):
    Rp, Ep, dp = M.residual(Xt[0], Xt[1], Xt[2], om_guess, log=log)
    Rv = [val(r) for r in Rp]
    Jt = [[Rp[i].d[j] for j in range(3)] for i in range(3)]
    C = arb_mat([[arb(Jt[i][j].mid()) for j in range(3)] for i in range(3)]).inv()
    Xb = [Xt[i] + sym(rad) for i in range(3)]
    if log: log('point run done %.1fs, |R|=%s, width=%s' % (dp['time'], amax(Rv).str(3), dp['width_end'].str(3)))
    Rb, Eb, db = M.residual(Xb[0], Xb[1], Xb[2], om_guess, log=log)
    JX = [[Rb[i].d[j] for j in range(3)] for i in range(3)]
    E = [sym(Eb) for i in range(3)]
    Kx = []
    for i in range(3):
        s = Xt[i] - sum((C[i, j]*(Rv[j] + E[j]) for j in range(3)), arb(0))
        for j in range(3):
            cij = (1 if i == j else 0) - sum((C[i, l]*JX[l][j] for l in range(3)), arb(0))
            s += cij*(Xb[j] - Xt[j])
        Kx.append(s)
    inside = all(Xb[i].contains_interior(Kx[i]) for i in range(3))
    return inside, Kx, Xb, dict(R=Rv, Rpoint=Rp, Ebox=Eb, Epoint=Ep, point=dp, box=db, J=JX, C=C)

def newton_refine(M, Xt, om_guess, iters=1, log=print):
    """midpoint Newton using the certified point run (non-rigorous step, only to centre the box)."""
    for it in range(iters):
        Rp, Ep, dp = M.residual(Xt[0], Xt[1], Xt[2], om_guess, log=None)
        Jm = arb_mat([[arb(Rp[i].d[j].mid()) for j in range(3)] for i in range(3)])
        Rm = arb_mat([[arb(val(Rp[i]).mid())] for i in range(3)])
        dX = Jm.solve(Rm)
        Xt = [arb((Xt[i] - dX[i, 0]).mid()) for i in range(3)]
        if log: log('newton %d: |R| = %s  step = %s  (%.0fs)' % (it, amax([val(r) for r in Rp]).str(3),
                     amax([dX[i, 0] for i in range(3)]).str(3), dp['time']))
    return Xt

# ---------------- PE2 checks ----------------
def saddle_check(sm):
    """eigenvalues of DG0 at the sonic point: {0 (sonic line, exact), lam, lam_b}; lam_b = tr - lam."""
    DG = sm['co']['DG0']
    tr = sum((val(DG[i][i]) for i in range(3)), arb(0))
    lam = val(sm['lam'])
    lam_b = tr - lam
    return dict(lam=lam, lam_b=lam_b, saddle=bool(lam*lam_b < 0))

def count_V_zeros(Pb, boxes):
    """over Lohner step boxes (x, h, Bx): sign changes of V, with f_V monotonicity where V straddles 0."""
    signs = []; ok = True; amb = 0
    for (x, h, Bx) in boxes:
        Vb = val(Bx[2])
        if Vb > 0: signs.append(1)
        elif Vb < 0: signs.append(-1)
        else:
            fv = Pb.F.eval([val(t) for t in Bx], conv=sym_const)[2]
            if fv > 0 or fv < 0:
                signs.append(0); amb += 1
            else:
                ok = False; signs.append(None)
    # zeros = sign changes among certified signs (ambiguous steps are monotone => at most one zero each)
    cert = [s for s in signs if s in (1, -1)]
    changes = sum(1 for a, b in zip(cert, cert[1:]) if a != b)
    return dict(sign_changes=changes, ambiguous_steps=amb, ok=ok, first=cert[0] if cert else None, last=cert[-1] if cert else None)

def outward_leg(Pb, V0box, om_guess, nu_s, X1=12.0, Kt=30, q='0.25', Dw='1e-3', log=None):
    """Lohner from the sonic manifold (theta > 0) to X1, D_sonic<0 and V>0 on every box; far-side lemma."""
    T0 = time.time()
    kk = {k_: Pb.k}
    ro = Tape([load_expr(FIELDS['ro'])], [N_, om_, V_], kk)
    rV = Tape([load_expr(FIELDS['rV'])], [N_, om_, V_], kk)
    aT = Tape([load_expr(FIELDS['Acons']) - 2 - (1 - k_)*om_], [N_, om_, V_], kk)
    sm = Pb.sonic_manifold(V0box, om_guess, arb(nu_s))
    th = arb('0.9')*arb(nu_s)
    y = [PM.horner([val(a) for a in sm['co']['Y'][i]], th) + sym(sm['eps']) for i in range(3)]
    x = PM.horner([val(a) for a in sm['co']['X']], th)
    xs = arb(x.mid())
    dmax = (xs - x).abs_upper() + sm['epsX']
    B, fB, _ = transport_bound(Pb, y, dmax)
    y = [y[i] + (xs - x)*fB[i] + sym(sm['epsX']*fB[i].abs_upper()) for i in range(3)]
    S = LO.LohnerSet(y)
    info = LO.integrate(Pb.F, S, xs, arb(X1), K=Kt, q=q, log=log)
    okD = True; okV = True; maxD = None
    for (_, _, Bx) in info['boxes']:
        Bv = [val(t) for t in Bx]
        Dv = Pb.Dson.eval(Bv, conv=sym_const)[0]
        okD &= bool(Dv < 0); okV &= bool(Bv[2] > 0)
    # sonic-manifold piece theta in [0, th]: V(theta) = V_K + E > 0 ?
    Vpiece = PM.horner([val(a) for a in sm['co']['Y'][2]], arb(0) + th/2 + sym(th/2)) + sym(sm['eps'])
    okV &= bool(Vpiece > 0)
    N1, om1, V1 = S.hull()
    Dw = arb(Dw)
    Q = [N1.upper()/2 + sym(N1.upper()/2 + arb(10)**(-40)), om1 + sym(Dw), V1 + sym(Dw)]
    a = aT.eval(Q, conv=sym_const)[0]
    res = dict(steps=info['steps'], ok_D=okD, ok_V=okV, a=a, y1=S.hull(), time=None)
    if not (a < 0):
        res.update(lemma=False, time=time.time() - T0); return res
    mu = -arb(a.upper())
    co_ = ro.eval(Q, conv=sym_const)[0].abs_upper(); cV = rV.eval(Q, conv=sym_const)[0].abs_upper()
    DQ = Pb.Dson.eval(Q, conv=sym_const)[0]
    dom = co_*N1.upper()/mu; dV = cV*N1.upper()/mu
    lemma = bool(dom < Dw) and bool(dV < Dw) and bool(DQ < 0) and bool(Q[2] > 0)
    om_inf = om1 + sym(dom); V_inf = V1 + sym(dV)
    A_inf = centred_eval(Pb.Acons, [arb(0), om_inf, V_inf])
    res.update(lemma=lemma, mu=mu, om_inf=om_inf, V_inf=V_inf, A_inf=A_inf, m_over_r=(1 - 1/A_inf)/2,
               time=time.time() - T0)
    return res
