"""B4 build seat (new file): stages S (sonic Frobenius) and L (linear transport) of the certified
characteristic function, plus the theta-parameterised centre stage, assembled into D(kappa).

  D(kappa) = det[ s_kappa(xbar), rho_1(xbar), rho_2(xbar), e_V ]        (PROOF-B4 sec. 1.5)

Common analytic tool (RegSing): t dY/dt = (M0(t) + kappa M1(t)) Y, M analytic on |t| <= nu_c,
Taylor coefficients M0_n, M1_n (balls) for n <= K and circle sups
   supC0 >= sup_{|t| = nu_c} ||M0(t) - M0_0||_inf,  supC1 >= sup_{|t| = nu_c} ||M1(t)||_inf.
Series y_n: (n - M_0) y_n = sum_{m=1}^n M_m y_{n-m}; rigorous tail exactly as lin_cert.CentreLinear
(majorant with ||M_m|| <= S/nu_c^m, ||(n - M_0)^{-1}|| <= 1/(n - ||M_0||) for n > ||M_0||).

Stage S (sonic): parameter theta of the certified B3-type sonic manifold (param.py), x = X(theta),
   theta dY/dtheta = theta X'(theta) Q^{-1}(P0 + kappa P1) Y = [theta x' Ptop ; F/(detMx/theta)] Y
   (lin_mc.py universal matrix; detMx vanishes exactly at theta = 0 on the sonic line, so detMx/theta
   is the shifted series).  M_0 has rank one (exponent rho(kappa), PROOF-B4 Fact 1.3); y_0 solves
   [w ; I(kappa) ; e_A ; e_N] y_0 = [0, 0, a, b] (w = the nonzero row direction of M_0).
Stage L (transport): on each stored tube step (lin_tube.py) the linear system in the step variable t
   (tau or x) dY/dt = A(t) Y, A = [x' Ptop ; (F0 + kappa F1)/detMx] with background Taylor coefficients
   inflated by the certified deviation; transition matrix by Taylor series + Cauchy remainder
   (|Phi(t)| <= exp(a |t|) on |t| <= rho from a = sup ||A|| over the complex a priori box); vector
   carried as a complex Lohner set (weighted QR) so kappa-balls do not wrap.
Centre: theta_c parameter of the certified centre manifold, theta dY/dtheta = (Dc/lam_c) Mc(Z, kappa) Y,
   rho_1, rho_2 from y_0 = e_Nbar, e_wbar; evaluated at theta_m with ln theta_m + chi(theta_m) + xc = xbar.
"""
import sys, os, math, time, pickle
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, acb, acb_mat, arb_mat, ctx
from core import Tape, sym, csym, sym_const
import system as S

def deser(p):
    return arb(p[0]) + sym(arb(p[1]))

def load_expr(s):
    return eval(s, {k: getattr(sp, k) for k in dir(sp)})

def horner(c, t):
    r = 0*t
    for a in reversed(c): r = r*t + a
    return r

def deriv(c):
    return [(n + 1)*c[n + 1] for n in range(len(c) - 1)]

def smul(a, b, K):
    return [sum((a[i]*b[n - i] for i in range(n + 1) if i < len(a) and n - i < len(b)), arb(0)) for n in range(K + 1)]

def sdiv(a, b, K):
    """series a/b (b_0 != 0)"""
    qq = []
    for n in range(K + 1):
        s = a[n] if n < len(a) else arb(0)
        for i in range(n):
            if n - i < len(b): s = s - qq[i]*b[n - i]
        qq.append(s/b[0])
    return qq

def rownorm(Mrows):
    m = arb(0)
    for row in Mrows:
        m = m.max(sum((x.abs_upper() for x in row), arb(0)))
    return m

class LinFields:
    def __init__(self, k):
        self.k = sp.Rational(k)
        d = pickle.load(open('<LAB>/engine/cert/lin_fields2.pkl', 'rb'))
        omp, Vp, xp = sp.symbols('omp Vp xp')
        kk = {S.k: self.k}
        outs = [load_expr(s) for s in d['Ptop']] + [load_expr(s) for s in d['F0']] + [load_expr(s) for s in d['F1']] + [load_expr(d['detMx'])]
        self.step = Tape(outs, [S.N, S.om, S.V, omp, Vp, xp], kk)
        Mc = [load_expr(s) for s in d['Mc']]
        self.cen = Tape([e.subs(S.kappa, 0) for e in Mc] + [sp.diff(e, S.kappa) for e in Mc], [S.P, S.W, S.U], kk)
        Ir = [load_expr(s) for s in d['Irow']]
        self.irow = Tape([e.subs(S.kappa, 0) for e in Ir], [S.N, S.om, S.V], kk)
        self.abar_x = Tape([load_expr(d['Abar_x'])], [S.N, S.om, S.V], kk)
        self.kappa_g = Tape([load_expr(d['kappa_g'])], [S.N, S.om, S.V], kk)

    def series_out(self, tape, ins, K):
        z = arb(0)
        tape.jet_init(K + 1, z, conv=lambda c: sym_const(c))
        for n in range(K + 1):
            for i in range(len(ins)):
                tape.C[i][n] = ins[i][n] if n < len(ins[i]) else z
            tape.jet_order(n)
        return [[tape.C[o][n] for n in range(K + 1)] for o in tape.outs]

# ----------------------------------------------------------------------------------- RegSing
class RegSing:
    def __init__(self, M0, M1, supC0, supC1, nu_c, rank1=False):
        self.rank1 = rank1
        self.M0, self.M1 = M0, M1          # lists over n of 4x4 lists of arb
        self.supC0, self.supC1, self.nu = supC0, supC1, arb(nu_c)
        self.K = len(M0) - 1

    def series(self, kap, y0):
        kap = acb(kap); K = self.K
        M = [acb_mat([[acb(self.M0[n][i][j]) + kap*acb(self.M1[n][i][j]) for j in range(4)] for i in range(4)]) for n in range(K + 1)]
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        y = [acb_mat([[c] for c in y0])]
        for n in range(1, K + 1):
            r = acb_mat(4, 1)
            for m in range(1, n + 1):
                r += M[m]*y[n - m]
            y.append((I4*n - M[0]).solve(r))
        return y, M

    def evaluate(self, kap, y0, t):
        y, M = self.series(kap, y0)
        K = self.K
        t = acb(t); at = abs(t).upper()
        nup = ((self.nu + at)/2).min(2*at) if at > 0 else self.nu/2
        Snorm = self.supC0 + arb(abs(acb(kap)).upper())*self.supC1
        rat = nup/self.nu
        mup = Snorm*rat/(1 - rat)
        M0n = rownorm([[M[0][i, j] for j in range(4)] for i in range(4)])
        n0 = K
        if self.rank1:
            # M_0 = u w^T, w.u = tr M_0 = rho(kappa):  (n - M_0)^{-1} = (I + M_0/(n - rho))/n
            trM = sum((M[0][i, i] for i in range(4)), acb(0))
            def cn(n):
                gap = abs(acb(n) - trM).lower()
                return (1 + M0n/gap)/n if gap > 0 else None
        else:
            def cn(n):
                return 1/(arb(n) - M0n) if arb(n) > M0n else None
        while True:
            c = cn(n0 + 1)
            if c is not None and c*mup <= 1: break
            n0 += 1
            if n0 > 100000: raise ArithmeticError('RegSing tail: n0 too large')
        if n0 > K: raise ArithmeticError('RegSing tail needs K >= %d (have %d)' % (n0, K))
        B = arb(0)
        for n in range(K + 1):
            B = B.max(max((abs(y[n][i, 0]).upper() for i in range(4)))*nup**n)
        qq = at/nup
        tail = B*qq**(K + 1)/(1 - qq)
        out = []
        for i in range(4):
            s = acb(0)
            for n in range(K, -1, -1): s = s*t + y[n][i, 0]
            out.append(s + csym(tail))
        return out, dict(tail=tail, n0=n0, mup=mup)

def circle_sup(evalM, M0_0, nu_c, narc=128):
    """sup over |t| = nu_c of ||M0(t) - M0_0|| and ||M1(t)|| (row sums); evalM(tbox) -> (M0 rows, M1 rows) acb"""
    s0 = arb(0); s1 = arb(0)
    pi = arb.pi()
    for a in range(narc):
        th = (arb(a) + arb('0.5') + sym(arb('0.5')))*2*pi/narc
        tb = acb(th.cos(), th.sin())*arb(nu_c)
        A0, A1 = evalM(tb)
        r0 = rownorm([[A0[i][j] - acb(M0_0[i][j]) for j in range(4)] for i in range(4)])
        r1 = rownorm(A1)
        if not (r0.is_finite() and r1.is_finite()):
            raise ArithmeticError('circle enclosure not finite')
        s0 = s0.max(r0); s1 = s1.max(r1)
    return s0, s1

# ----------------------------------------------------------------------------------- assembly
class Linear:
    def __init__(self, k, tube_file=None, prec=200, KS=78, KC=60, KL=60, rho_fac='0.6', ncell=4, tol=1e-30, log=print):
        self.ncell = ncell; self.log_tol = math.log(tol)
        ctx.prec = prec
        self.log = log
        T0 = time.time()
        self.k = sp.Rational(k)
        fn = tube_file or '<LAB>/engine/cert/lin_tube_k%s.pkl' % str(self.k).replace('/', '_')
        tb = pickle.load(open(fn, 'rb'))
        self.tb = tb
        self.LF = LinFields(self.k)
        self.KS, self.KC, self.KL = KS, KC, KL
        self.rho_fac = arb(rho_fac)
        self.xbar = deser(tb['XBAR'])
        self._sonic_setup()
        log('stage S ready (%.0fs)' % (time.time() - T0))
        self._centre_setup()
        log('centre ready (%.0fs)' % (time.time() - T0))
        self._tube_setup()
        log('stage L ready: %d steps (%.0fs)' % (len(self.steps), time.time() - T0))

    # ---------- stage S
    def _sonic_setup(self):
        tb, LF = self.tb, self.LF
        sm = tb['sm']
        Y = [[deser(p) for p in sm['Y'][i]] for i in range(3)]
        X = [deser(p) for p in sm['X']]
        eps, epsX = deser(sm['eps']), deser(sm['epsX'])
        nu_s = arb(tb['nu_s'])
        K = min(self.KS, len(X) - 2)
        dY = [deriv(Y[i]) for i in range(3)]; dX = deriv(X)
        outs = LF.series_out(LF.step, [Y[0], Y[1], Y[2], dY[1], dY[2], dX], K + 1)
        Ptop, F0, F1, dM = outs[0:8], outs[8:16], outs[16:24], outs[24]
        if not dM[0].contains(0): raise ArithmeticError('detMx(theta=0) does not contain 0')
        dsh = dM[1:] + [arb(0)]            # detMx/theta (exact zero constant term on the sonic line)
        xPt = [smul(dX, Ptop[e], K + 1) for e in range(8)]
        M0 = []; M1 = []
        F0q = [sdiv(F0[e], dsh, K) for e in range(8)]; F1q = [sdiv(F1[e], dsh, K) for e in range(8)]
        for n in range(K + 1):
            top = [[(xPt[4*i + j][n - 1] if n >= 1 else arb(0)) for j in range(4)] for i in range(2)]
            fl0 = [[F0q[4*i + j][n] for j in range(4)] for i in range(2)]
            fl1 = [[F1q[4*i + j][n] for j in range(4)] for i in range(2)]
            M0.append(top + fl0); M1.append([[arb(0)]*4, [arb(0)]*4] + fl1)
        nuc = arb('0.95')*nu_s
        def evalM(tb_):
            yv = [horner([acb(a) for a in Y[i]], tb_) + csym(eps) for i in range(3)]
            dyv = [horner([acb(a) for a in dY[i]], tb_) + csym(eps/(nu_s - nuc)) for i in range(3)]
            dxv = horner([acb(a) for a in dX], tb_) + csym(epsX/(nu_s - nuc))
            o = LF.step.eval([yv[0], yv[1], yv[2], dyv[1], dyv[2], dxv], conv=lambda c: acb(sym_const(c)))
            dmt = o[24]/tb_
            A0 = [[tb_*dxv*o[4*i + j] for j in range(4)] for i in range(2)] + [[o[8 + 4*i + j]/dmt for j in range(4)] for i in range(2)]
            A1 = [[acb(0)]*4, [acb(0)]*4] + [[o[16 + 4*i + j]/dmt for j in range(4)] for i in range(2)]
            return A0, A1
        s0, s1 = circle_sup(evalM, M0[0], nuc)
        # rank-one check of M_0 at kappa = 0 and of M1_0 (both: fluid rows only, rank <= 1 by adj(Mx(0)))
        self.sonic = RegSing(M0, M1, s0, s1, nuc, rank1=True)
        self.y0bg = [deser(p) for p in sm['y0']]
        self.th0 = deser(tb['th0'])
        self.son_info = dict(supC0=s0, supC1=s1, nu=nuc, K=K)
        # exact sonic data: rho(kappa) pieces, kappa_g, Abar_x(0)
        self.kappa_g = LF.kappa_g.eval(self.y0bg, conv=sym_const)[0]
        self.abar_x0 = LF.abar_x.eval(self.y0bg, conv=sym_const)[0]
        self.I0 = LF.irow.eval(self.y0bg, conv=sym_const)

    def sonic_start(self, kap, a=0, b=1):
        kap = acb(kap)
        M0 = [[acb(self.sonic.M0[0][i][j]) + kap*acb(self.sonic.M1[0][i][j]) for j in range(4)] for i in range(4)]
        w = max((M0[2], M0[3]), key=lambda r: float(sum(abs(x).mid() for x in r)))
        I0 = [acb(self.I0[0]) + kap] + [acb(x) for x in self.I0[1:]]
        Mat = acb_mat([w, I0, [1, 0, 0, 0], [0, 1, 0, 0]])
        y0 = Mat.solve(acb_mat([[0], [0], [a], [b]]))
        y0 = [y0[i, 0] for i in range(4)]
        Y, info = self.sonic.evaluate(kap, y0, acb(self.th0))
        return Y, info

    # ---------- centre
    def _centre_setup(self):
        tb, LF = self.tb, self.LF
        cm = tb['cm']
        Z = [[deser(p) for p in cm['Y'][i]] for i in range(3)]
        chi = [deser(p) for p in cm['X']]
        eps, epsX = deser(cm['eps']), deser(cm['epsX'])
        lam = deser(cm['lam']); nu = arb(cm['nu'])
        K = min(self.KC, len(chi) - 1)
        import pipeline as PL
        Dtape = Tape([load_expr(PL.FIELDS['pol_c'][3])], [S.P, S.W, S.U], {S.k: self.k})
        outs = LF.series_out(LF.cen, Z, K)
        dc = LF.series_out(Dtape, Z, K)[0]
        fac = [x/lam for x in dc]
        M0 = []; M1 = []
        C0 = [smul(fac, outs[e], K) for e in range(16)]; C1 = [smul(fac, outs[16 + e], K) for e in range(16)]
        for n in range(K + 1):
            M0.append([[C0[4*i + j][n] for j in range(4)] for i in range(4)])
            M1.append([[C1[4*i + j][n] for j in range(4)] for i in range(4)])
        nuc = arb('0.95')*nu
        def evalM(tb_):
            zv = [horner([acb(a) for a in Z[i]], tb_) + csym(eps) for i in range(3)]
            o = LF.cen.eval(zv, conv=lambda c: acb(sym_const(c)))
            f = Dtape.eval(zv, conv=lambda c: acb(sym_const(c)))[0]/lam
            return [[f*o[4*i + j] for j in range(4)] for i in range(4)], [[f*o[16 + 4*i + j] for j in range(4)] for i in range(4)]
        s0, s1 = circle_sup(evalM, M0[0], nuc)
        self.centre = RegSing(M0, M1, s0, s1, nuc)
        # theta_m: ln theta + chi(theta) + xc = xbar  (interval Newton)
        xc = deser(tb['xc'])
        th = (self.xbar - xc).exp()
        dchi = deriv(chi)
        for it in range(40):
            g = th.log() + horner(chi, th) + xc - self.xbar
            dg = 1/th + horner(dchi, th)
            th = arb((th - g/dg).mid())
        r = arb('1e-30') + epsX*10 + xc.rad()*10
        for it in range(60):
            Th = th + sym(r)
            g = th.log() + horner(chi, th) + sym(epsX) + xc - self.xbar
            dg = 1/Th + horner(dchi, Th)
            Nt = th - g/dg
            if Th.contains_interior(Nt): break
            r *= 4
        else:
            raise ArithmeticError('theta_m interval Newton failed')
        self.theta_m = Nt
        self.cen_info = dict(supC0=s0, supC1=s1, nu=nuc, K=K, theta_m=Nt)

    def centre_regular(self, kap):
        out = []
        for e in (1, 2):
            y0 = [acb(1) if i == e else acb(0) for i in range(4)]
            Y, info = self.centre.evaluate(kap, y0, acb(self.theta_m))
            out.append(Y)
        return out

    # ---------- stage L
    def _tube_setup(self):
        LF = self.LF; KL = self.KL
        import driver_v2 as V2, pipeline as PL
        Ftape = Tape([load_expr(s) for s in PL.FIELDS['f']], [S.N, S.om, S.V], {S.k: self.k})
        T4 = V2.tau_tape(self.k)
        self.steps = []
        for st in self.tb['steps']:
            kind = st['kind']; h = deser(st['h']); Rc = deser(st['Rc'])
            c = [deser(p) for p in st['c']]
            nrv = [deser(p) for p in st['nrv']]
            setrad = deser(st['setrad'])
            ncomp = 4 if kind == 'tau' else 3
            Yc = [[deser(p) for p in st['Yc'][i][:KL + 2]] for i in range(ncomp)]
            rho = self.rho_fac*Rc
            # all certified backgrounds start the step in c + [-setrad, setrad]: recompute the jets and
            # the complex a priori box (step_cauchy2 theorem) with this BALL initial value, so the
            # coefficients and the disc box enclose every background member (no Lipschitz/Gronwall needed)
            import taylor as TY
            field = T4 if kind == 'tau' else Ftape
            y0b = [c[i] + sym(setrad*arb('1.01')) for i in range(ncomp)]
            from param import _c as _cc
            y1, Bx, Yc2 = TY.step_cauchy2(field, y0b, h, KL + 1, q=abs(h)/Rc, conv=lambda q_: _cc(q_, arb(0)))
            Yi = [[Yc2[i][n] for n in range(KL + 2)] for i in range(ncomp)]
            totrad = [Bx[i].rad() for i in range(ncomp)]
            dY = [deriv(Yi[i]) for i in range(ncomp)]
            xp = dY[3] if kind == 'tau' else [arb(1)] + [arb(0)]*KL
            outs = LF.series_out(LF.step, [Yi[0], Yi[1], Yi[2], dY[1], dY[2], xp], KL)
            Ptop, F0, F1, dM = outs[0:8], outs[8:16], outs[16:24], outs[24]
            xPt = [smul(xp, Ptop[e], KL) for e in range(8)]
            F0q = [sdiv(F0[e], dM, KL) for e in range(8)]; F1q = [sdiv(F1[e], dM, KL) for e in range(8)]
            A0 = []; A1 = []
            for n in range(KL + 1):
                A0.append(acb_mat([[xPt[4*i + j][n] for j in range(4)] for i in range(2)] + [[F0q[4*i + j][n] for j in range(4)] for i in range(2)]))
                A1.append(acb_mat([[0]*4, [0]*4] + [[F1q[4*i + j][n] for j in range(4)] for i in range(2)]))
            # sup ||A|| over the complex disc |t| <= rho: cells, background = Taylor polynomial on the cell
            # + Cauchy tails from the a priori box (|y_n| <= totrad/Rc^n for n >= 1)
            m_ = len(Yi[0]); r_ = rho/Rc
            tailY = [totrad[i]*r_**m_/(1 - r_) for i in range(ncomp)]
            tailD = [totrad[i]/Rc*m_*r_**(m_ - 1)/(1 - r_)**2 for i in range(ncomp)]
            dYp = [[acb(a) for a in deriv(Yi[i])] for i in range(ncomp)]
            Yp = [[acb(a) for a in Yi[i]] for i in range(ncomp)]
            def disc_sup(ncell, rho):
                hc = rho/ncell
                a0 = arb(0); a1 = arb(0)
                for ia in range(ncell):
                    for ib in range(ncell):
                        tcell = acb(-rho + (2*ia + 1)*hc, -rho + (2*ib + 1)*hc) + csym(hc*arb('1.415'))
                        yb = [horner(Yp[i], tcell) + csym(tailY[i]) for i in range(ncomp)]
                        dyb = [horner(dYp[i], tcell) + csym(tailD[i]) for i in range(ncomp)]
                        xpb = dyb[3] if kind == 'tau' else acb(1)
                        o = LF.step.eval([yb[0], yb[1], yb[2], dyb[1], dyb[2], xpb], conv=lambda cc: acb(sym_const(cc)))
                        if not o[24].is_finite() or (o[24].real.contains(0) and o[24].imag.contains(0)):
                            return None
                        a0 = a0.max(rownorm([[xpb*o[4*i + j] for j in range(4)] for i in range(2)] + [[o[8 + 4*i + j]/o[24] for j in range(4)] for i in range(2)]))
                        a1 = a1.max(rownorm([[o[16 + 4*i + j]/o[24] for j in range(4)] for i in range(2)]))
                return a0, a1
            res = None
            # a complex zero of detMx may sit inside |t| <= rho (seen at k = 1/2): refine the cells,
            # then shrink rho for this step (rho must stay > 1.3|h| for the Cauchy remainder)
            for shrink in ('1', '0.75', '0.55', '0.4'):
                rho_try = rho*arb(shrink)
                if not rho_try > 1.3*abs(h): break
                r_ = rho_try/Rc
                tailY = [totrad[i]*r_**m_/(1 - r_) for i in range(ncomp)]
                tailD = [totrad[i]/Rc*m_*r_**(m_ - 1)/(1 - r_)**2 for i in range(ncomp)]
                for nc in (self.ncell, 2*self.ncell, 4*self.ncell):
                    res = disc_sup(nc, rho_try)
                    if res is not None: break
                if res is not None:
                    rho = rho_try; break
            if res is None:
                raise ArithmeticError('detMx cell box contains 0 on tube step %d (%s, h=%s) even after shrinking rho' % (len(self.steps), kind, h.str(4)))
            a0, a1 = res
            if not (a0.is_finite() and a1.is_finite()): raise ArithmeticError('disc sup not finite')
            self.steps.append(dict(h=h, rho=rho, A0=A0, A1=A1, a0=a0, a1=a1, kind=kind))

    def transport(self, kap, Y0, dbg=None):
        kap = acb(kap); KL = self.KL
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        c = acb_mat([[acb(x.real.mid(), x.imag.mid())] for x in Y0])
        Bm = I4
        r = acb_mat([[Y0[i] - c[i, 0]] for i in range(4)])
        akap = arb(abs(kap).upper())
        for st in self.steps:
            a = st['a0'] + akap*st['a1']
            qq = abs(st['h'])/st['rho']
            # adaptive order: smallest K with exp(a rho) q^(K+1)/(1-q) <= tol (cap KL)
            la = float((a*st['rho']).upper()); lq = math.log(float(qq.upper()))
            Ks = KL
            for Kt_ in range(12, KL + 1):
                if la + (Kt_ + 1)*lq - math.log(1 - float(qq.upper())) < self.log_tol:
                    Ks = Kt_; break
            A = [st['A0'][n] + kap*st['A1'][n] for n in range(Ks + 1)]
            Phi = [I4]
            for n in range(Ks):
                s = acb_mat(4, 4)
                for m in range(n + 1):
                    s += A[m]*Phi[n - m]
                Phi.append(s*(acb(1)/(n + 1)))
            h = acb(st['h'])
            P = acb_mat(4, 4)
            for n in range(Ks, -1, -1):
                P = P*h + Phi[n]
            rem = (a*st['rho']).exp()*qq**(Ks + 1)/(1 - qq)
            E = csym(rem)
            P = P + acb_mat([[E]*4 for _ in range(4)])
            # complex Lohner step (weighted QR)
            p = P*c
            M = P*Bm
            w = [float(abs(r[j, 0]).upper()) + 1e-300 for j in range(4)]
            cols = [[complex(float(M[i, j].real.mid()), float(M[i, j].imag.mid())) for i in range(4)] for j in range(4)]
            order = sorted(range(4), key=lambda j: -math.sqrt(sum(abs(z)**2 for z in cols[j]))*w[j])
            Qv = []
            for j in order:
                v = cols[j][:]
                for qv in Qv:
                    dd = sum(x*y.conjugate() for x, y in zip(v, qv)); v = [x - dd*y for x, y in zip(v, qv)]
                nv = math.sqrt(sum(abs(z)**2 for z in v))
                Qv.append([z/nv for z in v])
            Qm = acb_mat([[acb(Qv[j][i].real, Qv[j][i].imag) for j in range(4)] for i in range(4)])
            Qi = Qm.inv()
            cn = acb_mat([[acb(p[i, 0].real.mid(), p[i, 0].imag.mid())] for i in range(4)])
            r = (Qi*M)*r + Qi*(p - cn)
            c = cn; Bm = Qm
            if dbg is not None:
                dbg.append((st['kind'], float(st['a0'].mid()), float(st['a1'].mid()), float(st['rho'].mid()), float(abs(st['h']).mid()), float(rem.mid()),
                            max(float(abs(r[i, 0]).upper()) for i in range(4)), max(float(abs(c[i, 0]).upper()) for i in range(4))))
        Y = c + Bm*r
        return [Y[i, 0] for i in range(4)]

    # ---------- D
    def D(self, kap, a=0, b=1, detector=3):
        Ys, si = self.sonic_start(kap, a, b)
        Yx = self.transport(kap, Ys)
        rho1, rho2 = self.centre_regular(kap)
        col = [acb(1) if i == detector else acb(0) for i in range(4)]
        m = acb_mat([[Yx[i], rho1[i], rho2[i], col[i]] for i in range(4)])
        return m.det()
