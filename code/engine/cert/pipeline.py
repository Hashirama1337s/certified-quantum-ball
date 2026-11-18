"""Certified existence pipeline (B3) - matching residual and Krawczyk driver.

Unknowns X = (V0, W0, xc):  V0 = sonic velocity (KHA gauge, sonic point at x = 0),
W0 = centre density parameter (W = om/P^2 at the centre), xc = centre x-shift
(ln p = -xc).  Residual in the centre chart Z = (P, W, U):
    R(X) = Z_sonic-side(x_m') - Z_centre(theta_m)
Sonic side: certified 1D manifold (param.py) at the sonic equilibrium of the
desingularised polynomial field, then certified Taylor steps (taylor.py).
Centre side: certified 1D unstable manifold of the regular-centre equilibrium.
Tails enter as an interval addend E (generalised Krawczyk, see DESIGN-B3 notes):
    K(X) = x~ - C R^(x~) - C E + (I - C J(X)) (X - x~)  inside int X  =>  a zero exists.
"""
import sys, pickle, time
sys.path.insert(0, '<LAB>/engine/cert')
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import mpmath
from flint import arb, acb, arb_mat, ctx
from core import Tape, Dual, sym, csym, amax, sym_const
import param as PM
import taylor as TY
from param import val

N_, om_, V_, k_, P_, W_, U_ = sp.symbols('N om V k P W U')
FIELDS = pickle.load(open('<LAB>/engine/cert/fields.pkl', 'rb'))
def _ld(s): return sp.sympify(eval('sp.' + s) if False else sp.parse_expr(s) if False else sp.sympify(sp.srepr(sp.S(0))) )

def load_expr(s):
    return eval(s, {k: getattr(sp, k) for k in dir(sp)})

def D(v, n, i=None):
    """lift arb v to Dual with n directions (unit in direction i)"""
    z = arb(0)
    return Dual(v, [(arb(1) if j == i else z) for j in range(n)])

class Problem:
    def __init__(self, k, prec=200, K_sonic=40, K_centre=40, K_ty=24, log=print):
        ctx.prec = prec
        self.k = sp.Rational(k); self.log = log
        kk = {k_: self.k}
        self.kf = sym_const(self.k)
        self.c = self.kf.sqrt()
        ps = [load_expr(s) for s in FIELDS['pol_s']]
        pc = [load_expr(s) for s in FIELDS['pol_c']]
        f = [load_expr(s) for s in FIELDS['f']]
        # sonic field in STRUCTURED form (identical polynomial to fields.pol_s, checked in check_struct.py);
        # far less cancellation in ball jets than the expanded polynomial.
        import system as S0
        Ms, bs = S0.Msym, S0.bsym
        det = Ms[0, 0]*Ms[1, 1] - Ms[0, 1]*Ms[1, 0]
        gst = [det*S0.gN, Ms[1, 1]*bs[0] - Ms[0, 1]*bs[1], Ms[0, 0]*bs[1] - Ms[1, 0]*bs[0]]
        if __import__('os').environ.get('LB_EXPANDED'):
            gst, det = ps[:3], ps[3]
        self.Gs = Tape(gst, [N_, om_, V_], kk)
        self.Ss = Tape([det], [N_, om_, V_], kk)
        self.Gc = Tape(pc[:3], [P_, W_, U_], kk)
        self.Dc = Tape([pc[3]], [P_, W_, U_], kk)
        self.F = Tape(f, [N_, om_, V_], kk)
        import system as S
        lb = S.Msym[1, 0]*S.bsym[0] - S.Msym[0, 0]*S.bsym[1]
        self.LB = Tape([lb], [N_, om_, V_], kk)
        self.Acons = Tape([load_expr(FIELDS['Acons'])], [N_, om_, V_], kk)
        self.Dson = Tape([load_expr(FIELDS['Dsonic'])], [N_, om_, V_], kk)
        self.Ks, self.Kc, self.Kt = K_sonic, K_centre, K_ty
        self.Ustar = -arb(2)/(3*(1 + self.kf))
        self.lam_c = -self.kf*(1 + self.kf)
        self.cv = lambda c, like: PM._c(c, like)

    # ---------------- sonic equilibrium and eigenpair
    def sonic_point(self, V0, om_guess):
        """V0 arb or Dual(1 dir). Returns y0 list (same ring)."""
        c = self.c
        N0 = (1 - V0*c)/(c - V0)
        V0v, N0v = val(V0), val(N0)
        # interval Newton in om on LB(N0, om, V0) = 0
        om_t = arb(om_guess)
        Nm, Vm = arb(N0v.mid()), arb(V0v.mid())
        for it in range(12):          # midpoint Newton (non-rigorous, just to centre the box)
            Fm = self.LB.eval([Nm, om_t, Vm])[0]
            J, _ = PM.jac(self.LB, [Nm, om_t, Vm])
            om_t = arb((om_t - Fm/J[0][1]).mid())
        r = arb(2)**(-ctx.prec + 20) + 4*(N0v.rad() + V0v.rad())*100
        for it in range(40):          # rigorous interval Newton
            Iom = om_t + sym(r)
            Fm = self.LB.eval([N0v, om_t, V0v])[0]
            J, _ = PM.jac(self.LB, [N0v, Iom, V0v])
            newI = om_t - Fm/J[0][1]
            if Iom.contains_interior(newI):
                om_v = newI
                break
            r = r*4
        else:
            raise ArithmeticError('sonic om0 interval Newton failed')
        if isinstance(V0, Dual):
            Jp, _ = PM.jac(self.LB, [N0v, om_v, V0v])
            dN = N0.d[0]
            dom = -(Jp[0][0]*dN + Jp[0][2])/Jp[0][1]
            om0 = Dual(om_v, [dom])
        else:
            om0 = om_v
        return [N0, om0, V0]

    def eigen(self, y0, branch='fast'):
        """eigenpair of DGs(y0) for the nonzero eigenvalue of largest modulus (B1 sec.2)."""
        DG, _ = PM.jac(self.Gs, y0)
        DGv = [[val(DG[i][j]) for j in range(3)] for i in range(3)]
        Jm = mpmath.matrix([[float(DGv[i][j].mid()) for j in range(3)] for i in range(3)])
        ev, EV = mpmath.eig(Jm)
        ev = [complex(e) for e in ev]
        order = sorted(range(3), key=lambda i: -abs(ev[i]))
        ia = order[0]
        lam_t = ev[ia].real
        vt = [float(mpmath.re(EV[i, ia])) for i in range(3)]
        jn = max(range(3), key=lambda i: abs(vt[i]))
        vt = [x/vt[jn] for x in vt]
        free = [i for i in range(3) if i != jn]
        # Krawczyk on u = (v_free0, v_free1, lam)
        def Fv(u, DGm):
            v = [None]*3; v[jn] = 1 + 0*u[0]
            v[free[0]], v[free[1]] = u[0], u[1]
            return [sum((DGm[i][j]*v[j] for j in range(3)), 0*u[0]) - u[2]*v[i] for i in range(3)], v
        def JF(u, DGm):
            v = [None]*3; v[jn] = 1 + 0*u[0]; v[free[0]], v[free[1]] = u[0], u[1]
            return [[DGm[i][free[0]] - (u[2] if i == free[0] else 0), DGm[i][free[1]] - (u[2] if i == free[1] else 0), -v[i]]
                    for i in range(3)]
        ut = [arb(vt[free[0]]), arb(vt[free[1]]), arb(lam_t)]
        for it in range(8):   # Newton refine on midpoints
            Fm, _ = Fv(ut, DGv)
            Jm_ = JF(ut, DGv)
            du = PM.solve3([[arb(x.mid()) for x in row] for row in Jm_], [arb(x.mid()) for x in Fm])
            ut = [arb((ut[i] - du[i]).mid()) for i in range(3)]
        rad = arb(2)**(-ctx.prec + 30) + 1000*amax([DGv[i][j].rad() for i in range(3) for j in range(3)])
        U = [x + sym(rad) for x in ut]
        Fm, _ = Fv(ut, DGv)
        JU = JF(U, DGv)
        Cm = arb_mat([[arb(x.mid()) for x in row] for row in JF(ut, DGv)]).inv()
        CF = [sum((Cm[i, j]*Fm[j] for j in range(3)), arb(0)) for i in range(3)]
        Kr = []
        for i in range(3):
            s = ut[i] - CF[i]
            for j in range(3):
                cij = (1 if i == j else 0) - sum((Cm[i, l]*JU[l][j] for l in range(3)), arb(0))
                s += cij*(U[j] - ut[j])
            Kr.append(s)
        if not all(U[i].contains_interior(Kr[i]) for i in range(3)):
            raise ArithmeticError('eigen Krawczyk failed')
        uv = Kr
        if isinstance(y0[0], Dual):
            # implicit derivative: JF du = -(dDG/dV0) v
            dDG = [[DG[i][j].d[0] if isinstance(DG[i][j], Dual) else arb(0) for j in range(3)] for i in range(3)]
            _, v = Fv(uv, DGv)
            rhs = [-sum((dDG[i][j]*v[j] for j in range(3)), arb(0)) for i in range(3)]
            du = PM.solve3(JF(uv, DGv), rhs)
            u = [Dual(uv[i], [du[i]]) for i in range(3)]
        else:
            u = uv
        v = [None]*3
        one = 1 + 0*u[0]
        v[jn] = one; v[free[0]], v[free[1]] = u[0], u[1]
        lam = u[2]
        others = [ev[i] for i in order[1:]]
        return lam, v, DG, others

    # ---------------- sonic manifold
    def sonic_manifold(self, V0, om_guess, nu, tails=True):
        y0 = self.sonic_point(V0, om_guess)
        lam, v, DG0, others = self.eigen(y0)
        # normalise so that X_1 = 1  (theta ~ x near the sonic point)
        Jg, _ = PM.jac(self.Ss, y0)
        X1 = sum((Jg[0][j]*v[j] for j in range(3)), 0*lam)/lam
        v = [x/X1 for x in v]
        co = PM.coefficients(self.Gs, self.Ss, y0, lam, v, self.Ks)
        out = dict(y0=y0, lam=lam, v=v, co=co, others=others)
        if tails:
            Yc = [[val(a) for a in co['Y'][i]] for i in range(3)]
            Xc = [val(a) for a in co['X']]
            eps, epsX, info = PM.tail_bound(self.Gs, self.Ss, Yc, Xc, val(lam), co['DG0'], self.Ks, nu, log=self.log)
            out.update(eps=eps, epsX=epsX, tinfo=info)
        return out

    # ---------------- centre manifold
    def centre_manifold(self, W0, nu, tails=True):
        z = 0*W0
        y0 = [z, W0, self.Ustar + z]
        lam = self.lam_c + z
        v = [1 + z, z, z]
        co = PM.coefficients(self.Gc, self.Dc, y0, lam, v, self.Kc)
        out = dict(y0=y0, lam=lam, co=co)
        if tails:
            Yc = [[val(a) for a in co['Y'][i]] for i in range(3)]
            Xc = [val(a) for a in co['X']]
            eps, epsX, info = PM.tail_bound(self.Gc, self.Dc, Yc, Xc, val(lam), co['DG0'], self.Kc, nu, log=self.log)
            out.update(eps=eps, epsX=epsX, tinfo=info)
        return out


# ======================================================================
# integration, residual, Krawczyk
# ======================================================================
def emb(x, n, i):
    """embed a 1-direction Dual (or arb) into n directions at slot i"""
    if isinstance(x, Dual):
        d = [arb(0)]*n; d[i] = x.d[0]
        return Dual(x.v, d)
    return Dual(x, [arb(0)]*n)

def fnorm_inf(vec):
    return amax([val(v) for v in vec])

class Matcher:
    def __init__(self, Pb, nu_s, nu_c, frac_s='-0.9', frac_c='0.9', Kt=30, step_frac='0.2', q='0.25', h0='0.002'):
        self.q = arb(q); self.h0 = h0
        self.Pb = Pb; self.nu_s = arb(nu_s); self.nu_c = arb(nu_c)
        self.th1 = arb(frac_s)*self.nu_s
        self.thm = arb(frac_c)*self.nu_c
        self.Kt = Kt; self.step_frac = float(step_frac)
        self.sched = None
        self.cv = lambda like: (lambda c: PM._c(c, like))

    def _plan(self, y0, x0, x1):
        """step schedule from midpoint jets (non-rigorous choice, rigorous use)"""
        Pb = self.Pb; hs = []
        y = [arb(val(v).mid()) for v in y0]; x = float(x0)
        while x > float(x1) + 1e-15:
            Y = TY.jets(Pb.F, y, self.Kt, self.cv(y[0]))
            rho = min(min(abs(float(Y[i][n].mid()))**(-1.0/n) if float(Y[i][n].mid()) != 0 else 1e9
                          for i in range(3)) for n in range(self.Kt//2, self.Kt + 1))
            h = min(self.step_frac*rho, x - float(x1))
            h = float(arb(h).mid().str(6, radius=False)) if h > 1e-12 else h
            if x - h < float(x1) + 1e-9: h = x - float(x1)
            hs.append(h)
            y = [TY.horner(Y[i][:self.Kt + 1], arb(-h)) for i in range(3)]
            y = [arb(v.mid()) for v in y]
            x -= h
        return hs

    def residual(self, V0, W0, xc, om_guess, derivs=True, log=None):
        """R^(X, tails=0) as Dual(3) if derivs, plus E radius and diagnostics.
        V0, W0, xc: arb (points or boxes)."""
        Pb = self.Pb; T0 = time.time()
        n = 3
        V0d = Dual(V0, [arb(1)]) if derivs else V0
        sm = Pb.sonic_manifold(V0d, om_guess, self.nu_s)
        Yk = sm['co']['Y']; Xk = sm['co']['X']
        ys = [PM.horner(Yk[i], self.th1 + 0*V0d) for i in range(3)]
        xs = PM.horner(Xk, self.th1 + 0*V0d)
        if self.sched is None:
            self.x0 = arb(val(xs).mid())
        x0 = self.x0
        # transport start point from xs to the fixed x0 (first order + remainder)
        delta = x0 - xs
        fs = Pb.F.eval(ys, conv=self.cv(ys[0]))
        y0 = [ys[i] + delta*fs[i] for i in range(3)]
        dmax = val(delta).abs_upper() + sm['epsX']
        ysv = [val(v) for v in ys]
        B = [v + sym(2*dmax*fnorm_inf(fs) + arb(10)**(-60)) for v in ysv]
        fB = Pb.F.eval(B)
        Jb, _ = PM.jac(Pb.F, B)
        acc = amax([sum((Jb[i][j]*fB[j] for j in range(3)), arb(0)) for i in range(3)])
        R2s = dmax**2/2*acc
        e0 = sm['eps'] + sm['epsX']*fnorm_inf(fB) + R2s
        # centre side
        W0d = Dual(W0, [arb(1)]) if derivs else W0
        cm = Pb.centre_manifold(W0d, self.nu_c)
        Zc = [PM.horner(cm['co']['Y'][i], self.thm + 0*W0d) for i in range(3)]
        chi = PM.horner(cm['co']['X'], self.thm + 0*W0d)
        if derivs:
            Zc = [emb(z, n, 1) for z in Zc]; chi = emb(chi, n, 1)
            xcd = Dual(xc, [arb(0), arb(0), arb(1)])
            y0 = [emb(v, n, 0) for v in y0]
        else:
            xcd = xc
        xm = self.thm.log() + chi + xcd
        if self.sched is None:
            self.xbar = arb(val(xm).mid())
            self.sched = 'adaptive'
        # certified Taylor steps x0 -> xbar (adaptive; a priori via centred complex Picard)
        y = y0; x = arb(x0); Gam = arb(1); err = e0; detinfo = []
        h = float(self.h0); nst = 0; xt = float(self.xbar)
        while val(x).mid() > self.xbar.mid():
            rem = float((x - self.xbar).mid())
            hh = min(h, rem)
            last = hh >= rem - 1e-15
            ha = (x - self.xbar) if last else arb(str(hh))
            try:
                y1, Bx, Yc = TY.step_cauchy2(Pb.F, y, -ha, self.Kt, q=self.q, conv=self.cv(y[0]))
            except ArithmeticError:
                h = hh/2
                if h < 1e-7: raise ArithmeticError('step collapse at x=%s' % x)
                continue
            Bv = [val(b_) for b_ in Bx]
            Lh = TY.lipschitz(Pb.F, [b_ + sym(err*Gam) for b_ in Bv], self.cv(arb(0)))
            Gam = Gam*(Lh*ha).exp()
            detinfo.append((float(x.mid()), hh, Bv))
            y = y1; x = self.xbar if last else x - ha; nst += 1; h = hh*1.25
            if log and nst % 50 == 0: log('  step %d x=%.4f t=%.1fs' % (nst, float(x.mid()), time.time() - T0))
        self.nsub = nst
        # transport from xbar to xm
        dm = xm - self.xbar
        fe = Pb.F.eval(y, conv=self.cv(y[0]))
        yv = [val(v) for v in y]
        dmx = val(dm).abs_upper() + cm['epsX']
        B = [v + sym(2*dmx*fnorm_inf(fe) + err*Gam + arb(10)**(-60)) for v in yv]
        fB = Pb.F.eval(B); Jb, _ = PM.jac(Pb.F, B)
        acc = amax([sum((Jb[i][j]*fB[j] for j in range(3)), arb(0)) for i in range(3)])
        R2e = dmx**2/2*acc
        ym = [y[i] + dm*fe[i] for i in range(3)]
        N, om, V = ym
        Zs = [1/N, om*N*N, V*N]
        R = [Zs[i] - Zc[i] for i in range(3)]
        # E radius (tails, transports), in Z coordinates
        Nv, omv, Vv = [val(b) for b in B]
        JZ = amax([1/(Nv*Nv), (2*omv*Nv).abs_upper() + Nv*Nv, Vv.abs_upper() + Nv.abs_upper()])
        eend = Gam*e0 + R2e + fnorm_inf(fB)*cm['epsX']
        Erad = JZ*eend + cm['eps']
        diag = dict(Gamma=Gam, e0=e0, Erad=Erad, eps_s=sm['eps'], eps_c=cm['eps'], lam=sm['lam'], others=sm['others'],
                    nsteps=self.nsub, time=time.time() - T0, boxes=detinfo, y0=sm['y0'], sm=sm, cm=cm)
        return R, Erad, diag

def krawczyk(M, Xt, rad, om_guess, log=print):
    """Xt: list of 3 arb points (V0, W0, xc); rad: box half-width (arb)."""
    Rp, Ep, dp = M.residual(Xt[0], Xt[1], Xt[2], om_guess, derivs=True)
    Rv = [val(r) for r in Rp]
    Jt = [[val(Rp[i]).__class__(0) + Rp[i].d[j] for j in range(3)] for i in range(3)]
    C = arb_mat([[arb(Jt[i][j].mid()) for j in range(3)] for i in range(3)]).inv()
    Xb = [Xt[i] + sym(rad) for i in range(3)]
    Rb, Eb, db = M.residual(Xb[0], Xb[1], Xb[2], om_guess, derivs=True)
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
    return inside, Kx, Xb, dict(R=Rv, E=Eb, point=dp, box=db, J=JX, C=C)
