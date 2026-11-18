"""B4 build seat (new file): kappa-Taylor models of D(kappa).

Why: kappa-BALL evaluation of D wraps (measured at kappa_g: ball radius 0.0025 -> sonic width 3.7e-2,
D width 1.8e9), so segments of the argument-principle contour cannot be certified with balls.
Instead, on a disc |d| <= r around kappa_c (d = kappa - kappa_c), every quantity is carried as a
Taylor model  sum_{j<=m} c_j d^j + disc(R)  valid for ALL |d| <= r:
  * stage S: sonic Frobenius recursion with M(d) = Mc + d M1 (affine), the inverse
    (n - M_0(d))^{-1} applied by the exact Neumann hierarchy in d, remainder bounded with the
    uniform rank-one bound c_n(disc) >= ||(n - M_0(d))^{-1}||; theta-tail uniform over the disc
    (it depends on |kappa| only);
  * centre: same with y_0 = e_Nbar, e_wbar;
  * stage L: polynomial part by the exact d-hierarchy  Y_j' = A0(kc) Y_j + A1 Y_{j-1}  (Taylor in t),
    remainder function R (all orders > m, t-truncation, and the omitted d^{m+1} forcing) carried as a
    complex Lohner set with kappa-ball transition matrices (wrapping-controlled);
  * D = det[s, rho1, rho2, e_V] in Taylor-model arithmetic.
Then D(kappa_c + d) in D_tm(d) for every |d| <= r, and segment tests inside the disc cost nothing.
"""
import math, time
from flint import arb, acb, acb_mat
from core import sym, csym
from lin_stages import rownorm

class TM:
    """scalar Taylor model: coefficients c (list of acb, len m+1), remainder radius R (arb), radius r"""
    __slots__ = ('c', 'R', 'r')
    def __init__(self, c, R, r):
        self.c = c; self.R = R; self.r = r
    def maxnorm(self):
        s = self.R; p = arb(1)
        for a in self.c:
            s += abs(a).upper()*p; p *= self.r
        return s
    def __add__(self, o):
        return TM([a + b for a, b in zip(self.c, o.c)], self.R + o.R, self.r)
    def __sub__(self, o):
        return TM([a - b for a, b in zip(self.c, o.c)], self.R + o.R, self.r)
    def __mul__(self, o):
        m = len(self.c) - 1; r = self.r
        full = [acb(0)]*(2*m + 1)
        for i, a in enumerate(self.c):
            for j, b in enumerate(o.c):
                full[i + j] += a*b
        over = arb(0)
        for n in range(m + 1, 2*m + 1):
            over += abs(full[n]).upper()*r**n
        pa = TM(self.c, arb(0), r).maxnorm(); pb = TM(o.c, arb(0), r).maxnorm()
        R = over + self.R*pb + o.R*pa + self.R*o.R
        return TM(full[:m + 1], R, r)
    def neg(self):
        return TM([-a for a in self.c], self.R, self.r)
    def eval(self, d):
        s = acb(0)
        for a in reversed(self.c): s = s*d + a
        return s + csym(self.R)

class VTM:
    """vector Taylor model: C = acb_mat 4 x (m+1) (columns = powers of d), R = list of 4 arb radii"""
    def __init__(self, C, R, r):
        self.C = C; self.R = R; self.r = r
    @staticmethod
    def const(v, m, r):
        C = acb_mat(4, m + 1)
        for i in range(4): C[i, 0] = acb(v[i])
        return VTM(C, [arb(0)]*4, r)
    def comp(self, i):
        m = self.C.ncols() - 1
        return TM([self.C[i, j] for j in range(m + 1)], self.R[i], self.r)
    def norms(self):
        return [self.comp(i).maxnorm() for i in range(4)]

_SH = {}
def shmat(m1):
    if m1 not in _SH:
        S = acb_mat(m1, m1)
        for j in range(m1 - 1): S[j, j + 1] = 1
        _SH[m1] = S
    return _SH[m1]

def shift_cols(C):
    """columns j -> j+1, drop the last; returns (shifted, dropped column)"""
    m1 = C.ncols()
    dropped = [C[i, m1 - 1] for i in range(C.nrows())]
    return C*shmat(m1), dropped

def affine_matvec(Mc, M1, V):
    """(Mc + d M1) V  with Mc, M1 acb_mat 4x4"""
    r = V.r; m = V.C.ncols() - 1
    S, dropped = shift_cols(V.C)
    C = Mc*V.C + M1*S
    dd = M1*acb_mat([[x] for x in dropped])
    R = []
    for i in range(4):
        Ri = abs(dd[i, 0]).upper()*r**(m + 1)
        for j in range(4):
            Ri += (abs(Mc[i, j]).upper() + r*abs(M1[i, j]).upper())*V.R[j]
        R.append(Ri)
    return VTM(C, R, r)

def neumann_solve(Ac, A1, V, cn_disc):
    """X(d) = (Ac + d A1)^{-1} V(d):  Ac X_j = V_j - A1 X_{j-1};  remainder via the uniform bound
    cn_disc >= sup_disc ||(Ac + d A1)^{-1}||_inf."""
    r = V.r; m = V.C.ncols() - 1
    X = acb_mat(4, m + 1)
    prev = None
    for j in range(m + 1):
        b = acb_mat([[V.C[i, j]] for i in range(4)])
        if prev is not None: b = b - A1*prev
        xj = Ac.solve(b)
        for i in range(4): X[i, j] = xj[i, 0]
        prev = xj
    tail = A1*prev
    tnorm = max(abs(tail[i, 0]).upper() for i in range(4))*r**(m + 1)
    Rn = max(V.R) + tnorm
    R = [cn_disc*Rn]*4
    return VTM(X, R, r)

def det3_tm(cols, rows):
    """determinant of the 3x3 minor (rows) of three VTM columns, in TM arithmetic"""
    a = [[cols[j].comp(i) for j in range(3)] for i in rows]
    t1 = a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
    t2 = a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0])
    t3 = a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0])
    return t1 - t2 + t3

class LinearTM:
    """kappa-Taylor-model evaluation of D on a disc, reusing a lin_stages.Linear instance."""
    def __init__(self, L, m=24, rem_method='ball'):
        self.L = L; self.m = m; self.rem_method = rem_method

    # ----- RegSing in Taylor-model form
    def regsing_tm(self, RS, kc, r, y0, t, rank1):
        m = self.m; K = RS.K
        kc = acb(kc); rr = arb(r)
        akmax = arb(abs(kc).upper()) + rr
        Mc = [acb_mat([[acb(RS.M0[n][i][j]) + kc*acb(RS.M1[n][i][j]) for j in range(4)] for i in range(4)]) for n in range(K + 1)]
        M1 = [acb_mat([[acb(RS.M1[n][i][j]) for j in range(4)] for i in range(4)]) for n in range(K + 1)]
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        # uniform bound for (n - M_0(d))^{-1} over the disc
        M0n = rownorm([[Mc[0][i, j] for j in range(4)] for i in range(4)]) + rr*rownorm([[M1[0][i, j] for j in range(4)] for i in range(4)])
        trc = sum((Mc[0][i, i] for i in range(4)), acb(0)); tr1 = sum((M1[0][i, i] for i in range(4)), acb(0))
        nM1 = rownorm([[M1[0][i, j] for j in range(4)] for i in range(4)])
        def cn_tail(n):
            if rank1:
                gap = abs(acb(n) - trc).lower() - (rr*abs(tr1)).upper()
                return (1 + M0n/gap)/n if gap > 0 else None
            return 1/(arb(n) - M0n) if arb(n) > M0n else None
        def cn(n):
            if n > K or rank1: return cn_tail(n)      # rank one: exact formula valid for every n >= 1
            Ainv = (I4*n - Mc[0]).inv()
            na = rownorm([[Ainv[i, j] for j in range(4)] for i in range(4)])
            qn = na*rr*nM1
            return na/(1 - qn) if qn < 1 else None
        ys = [y0]
        for n in range(1, K + 1):
            rhs = None
            for mm in range(1, n + 1):
                term = affine_matvec(Mc[mm], M1[mm], ys[n - mm])
                rhs = term if rhs is None else VTM(rhs.C + term.C, [a + b for a, b in zip(rhs.R, term.R)], rhs.r)
            c = cn(n)
            if c is None: raise ArithmeticError('TM RegSing: no inverse bound at n=%d' % n)
            ys.append(neumann_solve(I4*n - Mc[0], -M1[0], rhs, c))
        # theta-tail, uniform over the disc (as RegSing.evaluate with |kappa| <= |kc| + r)
        t = acb(t); at = abs(t).upper()
        nup = ((RS.nu + at)/2).min(2*at)
        Snorm = RS.supC0 + akmax*RS.supC1
        rat = nup/RS.nu; mup = Snorm*rat/(1 - rat)
        n0 = K
        while True:
            c = cn_tail(n0 + 1)
            if c is not None and c*mup <= 1: break
            n0 += 1
            if n0 > 100000: raise ArithmeticError('TM tail n0')
        if n0 > K: raise ArithmeticError('TM tail needs K >= %d' % n0)
        B = arb(0)
        for n in range(K + 1):
            B = B.max(max(ys[n].norms())*nup**n)
        qq = at/nup
        tail = B*qq**(K + 1)/(1 - qq)
        C = acb_mat(4, m + 1); R = [tail]*4
        for n in range(K, -1, -1):
            C = C*t + ys[n].C
        Rs = [arb(0)]*4
        for n in range(K + 1):
            for i in range(4): Rs[i] += ys[n].R[i]*at**n
        return VTM(C, [R[i] + Rs[i] for i in range(4)], rr)

    def sonic_tm(self, kc, r, a=0, b=1):
        L = self.L; m = self.m; rr = arb(r); kc = acb(kc)
        RS = L.sonic
        M0c = [[acb(RS.M0[0][i][j]) + kc*acb(RS.M1[0][i][j]) for j in range(4)] for i in range(4)]
        M01 = [[acb(RS.M1[0][i][j]) for j in range(4)] for i in range(4)]
        iw = max((2, 3), key=lambda i: float(sum(abs(x).mid() for x in M0c[i])))
        I0 = [acb(L.I0[0]) + kc] + [acb(x) for x in L.I0[1:]]
        Ac = acb_mat([M0c[iw], I0, [1, 0, 0, 0], [0, 1, 0, 0]])
        A1 = acb_mat([M01[iw], [1, 0, 0, 0], [0]*4, [0]*4])
        # pole-free Cramer form: y0 = adj(A(d)) (0, 0, a, b)^T, A(d) = Ac + d A1 (affine): exact polynomials.
        # This multiplies s_kappa (hence D) by det A(kappa), an affine function whose only zero is
        # kappa_pole = -0.099 < 0 (PROOF-B4 Fact 1.3), so zeros and windings in Re kappa >= 0 are unchanged.
        z = [acb(0)]*(m + 1)
        def aff(c0, c1):
            c = list(z); c[0] = acb(c0); c[1] = acb(c1); return TM(c, arb(0), rr)
        Ad = [[aff(Ac[i, j], A1[i, j]) for j in range(4)] for i in range(4)]
        def minor3(rows, cols):
            x = [[Ad[i][j] for j in cols] for i in rows]
            return (x[0][0]*(x[1][1]*x[2][2] - x[1][2]*x[2][1]) - x[0][1]*(x[1][0]*x[2][2] - x[1][2]*x[2][0])
                    + x[0][2]*(x[1][0]*x[2][1] - x[1][1]*x[2][0]))
        bvec = [0, 0, a, b]
        C = acb_mat(4, m + 1)
        for i in range(4):
            acc = None
            for j in range(4):
                if bvec[j] == 0: continue
                # adj(A)_{ij} = (-1)^{i+j} M_{ji}
                cof = minor3([rr_ for rr_ in range(4) if rr_ != j], [cc for cc in range(4) if cc != i])
                if (i + j) % 2: cof = cof.neg()
                term = TM([x*bvec[j] for x in cof.c], cof.R*abs(bvec[j]), rr)
                acc = term if acc is None else acc + term
            for jj in range(m + 1): C[i, jj] = acc.c[jj]
        y0 = VTM(C, [arb(0)]*4, rr)
        return self.regsing_tm(RS, kc, r, y0, L.th0, rank1=True)

    def centre_tm(self, kc, r):
        L = self.L; m = self.m; rr = arb(r)
        out = []
        for e in (1, 2):
            y0 = VTM.const([1 if i == e else 0 for i in range(4)], m, rr)
            out.append(self.regsing_tm(L.centre, kc, r, y0, L.theta_m, rank1=False))
        return out

    def transport_tm(self, kc, r, Y):
        """Y: VTM at the tube start.  Polynomial part by the exact d-hierarchy; remainder R as a complex
        Lohner set propagated with kappa-ball transition matrices."""
        L = self.L; m = self.m; KL = L.KL; rr = arb(r); kc = acb(kc)
        kball = kc + csym(rr)
        akmax = arb(abs(kc).upper()) + rr
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        C = Y.C
        # remainder set (centre 0): r vector of complex boxes in basis Bm
        Bm = I4; rv = acb_mat([[csym(Y.R[i])] for i in range(4)])
        for st in L.steps:
            a = st['a0'] + akmax*st['a1']
            qq = abs(st['h'])/st['rho']
            la = float((a*st['rho']).upper()); lq = math.log(float(qq.upper()))
            Ks = KL
            for Kt_ in range(12, KL + 1):
                if la + (Kt_ + 1)*lq - math.log(1 - float(qq.upper())) < L.log_tol:
                    Ks = Kt_; break
            A0c = [st['A0'][n] + kc*st['A1'][n] for n in range(Ks + 1)]
            A1 = st['A1']
            h = acb(st['h'])
            # hierarchy series in t
            W = [C]
            colm = []          # column m of W_n, for the omitted forcing
            SH = shmat(m + 1)
            for n in range(Ks):
                Ts = acb_mat(4, m + 1); Us = acb_mat(4, m + 1)
                for l in range(n + 1):
                    Ts += A0c[l]*W[n - l]; Us += A1[l]*W[n - l]
                W.append((Ts + Us*SH)*(acb(1)/(n + 1)))
            Cn = acb_mat(4, m + 1)
            for n in range(Ks, -1, -1):
                Cn = Cn*h + W[n]
            # forcing into R: t-truncation of the full function + omitted d^{m+1} terms over the step
            rem_t = (a*st['rho']).exp()*qq**(Ks + 1)/(1 - qq)
            Ynorm = max(VTM(C, [arb(0)]*4, rr).norms()) + max(abs(rv[i, 0]).upper() for i in range(4))*rownorm([[Bm[i, j] for j in range(4)] for i in range(4)])
            sup_m = arb(0)
            for n in range(Ks + 1):
                sup_m += max(abs(W[n][i, m]).upper() for i in range(4))*abs(st['h']).upper()**n
            omit = abs(st['h']).upper()*(a*abs(st['h'])).exp()*rr**(m + 1)*st['a1']*sup_m*(1 + rem_t)
            force = rem_t*Ynorm + omit
            # remainder transport: thin transition at kappa_c (Lohner, true growth) plus an additive
            # perturbation ||Phi(kappa) - Phi(kappa_c)|| |R|, bounded by the first three kappa-derivative
            # matrices (exact hierarchy) and a one-step Cauchy tail with G = exp(a(r') |h|) on |d| <= r'
            if self.rem_method == 'ball':
                # remainder transport with the kappa-ball transition matrix (Lohner, centre 0)
                Akb = [st['A0'][n] + kball*st['A1'][n] for n in range(Ks + 1)]
                Phi = [I4]
                for n in range(Ks):
                    s_ = acb_mat(4, 4)
                    for l in range(n + 1): s_ += Akb[l]*Phi[n - l]
                    Phi.append(s_*(acb(1)/(n + 1)))
                P = acb_mat(4, 4)
                for n in range(Ks, -1, -1): P = P*h + Phi[n]
                P = P + acb_mat([[csym(rem_t)]*4 for _ in range(4)])
            else:
                NLV = 7
                H = [[I4]] + [[acb_mat(4, 4)] for _ in range(NLV - 1)]
                for n in range(Ks):
                    for lv in range(NLV):
                        s_ = acb_mat(4, 4)
                        for l in range(n + 1):
                            s_ += A0c[l]*H[lv][n - l]
                            if lv: s_ += A1[l]*H[lv - 1][n - l]
                        H[lv].append(s_*(acb(1)/(n + 1)))
                Phs = []
                for lv in range(NLV):
                    Pm = acb_mat(4, 4)
                    for n in range(Ks, -1, -1): Pm = Pm*h + H[lv][n]
                    Phs.append(Pm)
                P = Phs[0] + acb_mat([[csym(rem_t)]*4 for _ in range(4)])
                rp = rr*8
                ap = st['a0'] + (arb(abs(kc).upper()) + rp)*st['a1']
                G = (ap*abs(st['h'])).exp()
                dPhi = arb(0)
                for lv in range(1, NLV):
                    dPhi += rownorm([[Phs[lv][i, j] for j in range(4)] for i in range(4)])*rr**lv
                dPhi += G*(rr/rp)**NLV/(1 - rr/rp) + 2*rem_t
                Rbox = Bm*rv
                Rmax = max(abs(Rbox[i, 0]).upper() for i in range(4))
                force = force + dPhi*Rmax
            M = P*Bm
            w = [float(abs(rv[j, 0]).upper()) + 1e-300 for j in range(4)]
            cols = [[complex(float(M[i, j].real.mid()), float(M[i, j].imag.mid())) for i in range(4)] for j in range(4)]
            order = sorted(range(4), key=lambda j: -math.sqrt(sum(abs(z)**2 for z in cols[j]))*w[j])
            Qv = []
            for j in order:
                v = cols[j][:]
                for qv in Qv:
                    dd = sum(x*y.conjugate() for x, y in zip(v, qv)); v = [x - dd*y for x, y in zip(v, qv)]
                nv = math.sqrt(sum(abs(z)**2 for z in v)); Qv.append([z/nv for z in v])
            Qm = acb_mat([[acb(Qv[j][i].real, Qv[j][i].imag) for j in range(4)] for i in range(4)])
            Qi = Qm.inv()
            # midpoint + remainder: the ball radii of the propagated polynomial part go into R (Lohner),
            # so the d-polynomial is carried thin and does not wrap
            radv = []
            for i in range(4):
                ri = arb(0); pw = arb(1)
                for j in range(m + 1):
                    ri += (Cn[i, j].real.rad() + Cn[i, j].imag.rad())*pw; pw *= rr
                radv.append(ri)
            Cn = acb_mat([[acb(Cn[i, j].real.mid(), Cn[i, j].imag.mid()) for j in range(m + 1)] for i in range(4)])
            fv = acb_mat([[csym(force + radv[i])] for i in range(4)])
            rv = (Qi*M)*rv + Qi*fv
            Bm = Qm
            C = Cn
        Rend = Bm*rv
        return VTM(C, [abs(Rend[i, 0]).upper() for i in range(4)], rr)

    def D_tm(self, kc, r, a=0, b=1, detector=3):
        s = self.sonic_tm(kc, r, a, b)
        s = self.transport_tm(kc, r, s)
        r1, r2 = self.centre_tm(kc, r)
        if detector == 3:
            return det3_tm([s, r1, r2], (0, 1, 2))
        if detector == 0:
            return det3_tm([s, r1, r2], (1, 2, 3)).neg()
        raise ValueError('detector')

def winding_tm(LT, f_tm, corners, r, log=None, max_depth=30):
    """certified winding along the polygon: each edge covered by discs of radius r (centres 1.6 r apart);
    inside a disc the TM gives D on every sub-segment ball (half-plane test as lin_cert.winding)."""
    pts = list(corners) + [corners[0]]
    total = arb(0); ndisc = 0; nseg = 0
    for e in range(len(pts) - 1):
        a, b = pts[e], pts[e + 1]
        L_ = float(abs(b - a).mid())
        nd = max(1, math.ceil(L_/(1.6*r)))
        for i in range(nd):
            pa = a + (b - a)*acb(i)/nd; pb = a + (b - a)*acb(i + 1)/nd
            kc = (pa + pb)/2
            tm = f_tm(kc, r); ndisc += 1
            stack = [(pa - kc, pb - kc, 0)]
            while stack:
                u, v, dep = stack.pop()
                Du = tm.eval(u); Dv = tm.eval(v)
                cmid = (u + v)/2; rad = abs(v - u).upper()/2
                E = tm.eval(cmid + csym(rad))
                if (E*Du.conjugate()).real > 0 and (E*Dv.conjugate()).real > 0:
                    total += (Dv/Du).arg(); nseg += 1
                else:
                    if dep >= max_depth: raise ArithmeticError('TM segment not certified near %s' % (kc + cmid))
                    mid = (u + v)/2
                    stack.append((mid, v, dep + 1)); stack.append((u, mid, dep + 1))
            if log: log('   disc %d/%d on edge %d: centre %s, D_tm remainder %s' % (i + 1, nd, e, kc.str(5), tm.R.str(3)))
    w = total/(2*arb.pi())
    n = round(float(w.mid()))
    if not (w - n).abs_upper() < arb('0.25'):
        raise ArithmeticError('winding not near an integer: %s' % w)
    return n, ndisc, nseg, w

def _detA(self, kc):
    """det A(kappa) of the sonic y0 system (point value), for comparing TM c0 with point D"""
    L = self.L; RS = L.sonic; kc = acb(kc)
    M0c = [[acb(RS.M0[0][i][j]) + kc*acb(RS.M1[0][i][j]) for j in range(4)] for i in range(4)]
    iw = max((2, 3), key=lambda i: float(sum(abs(x).mid() for x in M0c[i])))
    I0 = [acb(L.I0[0]) + kc] + [acb(x) for x in L.I0[1:]]
    return acb_mat([M0c[iw], I0, [1, 0, 0, 0], [0, 1, 0, 0]]).det()
LinearTM.detA = _detA
