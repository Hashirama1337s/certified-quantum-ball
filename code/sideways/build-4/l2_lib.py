"""l2_lib.py - build seat build-4 (B7 SB2 "The Light Cone", BAR-SF-L2 seal 9077a35520), 2026-11-17.

Certified characteristic function D_2(kappa) of the axial (odd-parity) l = 2 gravitational-wave sector on the E5
solution, as kappa-Taylor models for the shared count engine (gravity_engine/argument_count.py: tm_at, point).

Object.  y = (Q, s_u u), u = g h_s (ax_sys.py; s_u > 0 a fixed balancing constant), y' = A(x, kappa) y with
A = A0 + kappa A1 + kappa^2 A2 (x-derivative), along the certified E5 path
    centre chart (theta_m, x = XBAR) -> 720 tube steps -> tube sonic manifold [th0, 0] -> outward sonic manifold
    [0, th_out] =: matching point M  <-  63 + (j_F - 63) outward steps  <-  light-cone Frobenius point x_lc.
    D(kappa) = det[ Y_lc(M), Y_c(M) ].

Pieces (each an enclosure valid for the TRUE background; all background data read-only, records as strings):
  centre   theta dy/dtheta = f(theta) A(Z(theta)) y, f = theta x_theta = Dc(Z)/lam_c (pipeline pol_c[3], lin_stages
           convention); y = theta ytil, theta dytil/dtheta = (f A - I) ytil, start ytil_0 = (1, -2/s_u) EXACT
           (ax_sys S5: M0 = [[-3, -2], [-2, 0]] exactly; kernel of M0 - 1), series + majorant tail (RegSing pattern of
           engine/cert/lin_stages.py, extended to kappa^2 and to kappa-Taylor models).
  steps    stored step records; per step the background Taylor jets are recomputed with step_cauchy2 from the
           certified start set (tube: c +- 1.01 setrad, as lin_stages/c1_cert; outward: the recorded hull), the
           transition P_s(kappa) = sum_{n<=K_s} Phi_n(kappa) h^n is an EXACT polynomial in kappa (ball coefficients),
           remainder exp(a(|kappa|) rho) q^(K_s+1)/(1-q), a = sup_{|t|<=rho} ||A0|| + |kappa| ||A1|| + |kappa|^2 ||A2||
           (Cauchy, lin_stages.transport).  Backward through a step: adj(P_s) = det(P_s) P_s^{-1}; det(P_s) =
           exp(int tr A) is entire and zero-free in kappa, so it changes no zero (stated in the normalisation proof).
  sonic    the two certified sonic manifolds as steps in theta (series around theta = 0, eps tails on |theta| <= nu_s).
  light    fresh step_cauchy2 expansion of the background at x_lc (state from L0.lc_root), N_0 := 1 exactly;
  cone     t y' = R(t) y, R = t A, rank-one residue R0 (ax_sys S3: det R0 = 0, tr R0 = sigma = (1 - kappa)/lam_lc).
           Analytic-class solution, X1 counting level, normalised by the factor prod_{j in S} (j - sigma):
             z_0 = v0 = (-R0_12, R0_11) (R0_12 = A/lam_lc, kappa-free, nonzero);
             z_n = adj(n - R0) [sum_m R_m z_{n-m} prod_{j in S, n-m < j <= n} (j - sigma)] / (n (n - sigma)^[n not in S])
             yhat_n = z_n prod_{j in S, j > n} (j - sigma);  S = {1} on every axial region (sigma = 1 iff kappa = 1 - lam).
           The 'bounded' null (PL2-4) is the same recursion with S = {} (normalised by the value at the cone).
Every function of kappa is carried as a Taylor model in d = kappa - kc of degree J with a remainder radius valid on
|d| <= r (class TM).  Frames: per node a constant float matrix (eigenvectors of A(node, kc)), used as an EXACT matrix;
absolute values are taken after the change of variables (lin_tm_v3 idea); any invertible frame is rigorous.
"""
import os, sys, math, time, pickle, hashlib, json
HERE = os.path.dirname(os.path.abspath(__file__))
B1 = os.path.normpath(os.path.join(HERE, '..', 'build-1'))
sys.path.insert(0, HERE); sys.path.insert(0, B1)
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import numpy as np
import sympy as sp
from flint import arb, acb, arb_poly, acb_poly, acb_mat, arb_series, ctx
import l0_cert as L0                      # read-only: lc_root, deser (sets ctx.cap = 32 at import; overridden below)
from core import Tape, sym, csym, sym_const
import taylor as TY
import param as PM
import system as S
import pipeline as PL
import driver_v2 as V2
import ax_sys
ctx.cap = 64                              # module top, once (doctrine #115): series below need at most CAPNEED terms

TUBE = '<LAB>/engine/cert/lin_tube2_k1_3.pkl'
OUTW = os.path.join(B1, 'outward_rec.pkl')
K13 = sp.Rational(1, 3)
UNIT = arb(0, 1)
KSCALE = 12                # transition polynomials are stored in khat = kappa/12 (all regions have |kappa| < 12)


def capguard(n):
    if ctx.cap < n:
        raise ArithmeticError('ctx.cap = %d < %d: series would be silently truncated' % (ctx.cap, n))


def coef(s, i):
    """strict coefficient access for arb_series / polys"""
    c = s.coeffs()
    if isinstance(s, arb_series):
        capguard(i + 1)
    return c[i] if i < len(c) else arb(0)


def deser(p):
    return L0.deser(p)


def horner(c, t):
    s = 0*t
    for a in reversed(c):
        s = s*t + a
    return s


def deriv(c):
    return [(n + 1)*c[n + 1] for n in range(len(c) - 1)]


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def up(x):
    return arb(x).abs_upper() if not isinstance(x, acb) else x.abs_upper()


# =========================================================================== Taylor models in d = kappa - kc
class TM:
    """p(d) (acb_poly, degree <= J) + E, |E| <= R, valid for |d| <= r"""
    __slots__ = ('p', 'R')

    def __init__(self, p, R):
        self.p = p if isinstance(p, acb_poly) else acb_poly([p])
        self.R = arb(R)


class Ctx:
    """per-disc constants: centre kc, radius r, degree J"""
    def __init__(self, kc, r, J):
        self.kc = acb(kc); self.r = arb(r); self.J = J
        self.kabs = abs(self.kc).upper() + self.r.upper()
        self.dpoly = acb_poly([self.kc, 1])         # kappa = kc + d
        self.dhat = acb_poly([self.kc/KSCALE, acb(1)/KSCALE])   # khat = kappa/KSCALE

    def psum(self, p):
        s = arb(0); rp = arb(1)
        for c in p.coeffs():
            s += c.abs_upper()*rp; rp *= self.r
        return s

    def trunc(self, p):
        cs = p.coeffs()
        if len(cs) <= self.J + 1:
            return p, arb(0)
        tail = arb(0); rp = self.r**(self.J + 1)
        for c in cs[self.J + 1:]:
            tail += c.abs_upper()*rp; rp *= self.r
        return acb_poly(cs[:self.J + 1]), tail

    def const(self, c):
        return TM(acb_poly([acb(c)]), 0)

    def kappa(self):
        return TM(self.dpoly, 0)

    def add(self, a, b):
        return TM(a.p + b.p, a.R + b.R)

    def sub(self, a, b):
        return TM(a.p - b.p, a.R + b.R)

    def scale(self, c, a):
        c = acb(c)
        return TM(a.p*c, a.R*c.abs_upper())

    def mul(self, a, b):
        q, tail = self.trunc(a.p*b.p)
        R = tail + self.psum(a.p)*b.R + self.psum(b.p)*a.R + a.R*b.R
        return TM(q, R)

    def sup(self, a):
        return self.psum(a.p) + a.R

    def evalb(self, a, u):
        return a.p(acb(u)) + csym(a.R)

    def recip_affine(self, c0, c1):
        """1/(c0 + c1 d) for |d| <= r (c0, c1 acb); refuses if the disc reaches the pole"""
        c0 = acb(c0); c1 = acb(c1)
        lo = c0.abs_lower()
        if not lo > 0:
            raise ArithmeticError('recip_affine: |c0| not bounded below')
        w = c1.abs_upper()*self.r/lo
        if not w < 1:
            raise ArithmeticError('recip_affine: disc reaches the pole (pole of the normalisation)')
        z = -c1/c0
        cs = [1/c0]
        for j in range(self.J):
            cs.append(cs[-1]*z)
        R = (1/lo)*w**(self.J + 1)/(1 - w)
        return TM(acb_poly(cs), R)

    def poly_kappa(self, P):
        """exact polynomial in kappa (arb_poly/acb_poly) -> TM in d"""
        q = acb_poly([acb(c) for c in P.coeffs()])(self.dhat) if P.length() > 1 else acb_poly([acb(c) for c in P.coeffs()])
        q2, tail = self.trunc(q)
        return TM(q2, tail)


def matvec(cx, M, v):
    return [cx.add(cx.mul(M[i][0], v[0]), cx.mul(M[i][1], v[1])) for i in range(2)]


def det2(cx, a, b):
    return cx.sub(cx.mul(a[0], b[1]), cx.mul(a[1], b[0]))


def rownorm(M):
    m = arb(0)
    for row in M:
        m = m.max(sum((acb(x).abs_upper() for x in row), arb(0)))
    return m


# =========================================================================== symbolic -> tapes
class Sys:
    """kappa-split A matrix entries as tapes.  kind = 'axial' (inputs N, A, om, V) or 'toy' (input N)."""
    def __init__(self, kind='axial', su=1):
        self.kind = kind
        self.su = sp.Rational(su)
        order = [(i, j) for i in range(2) for j in range(2)]
        if kind == 'axial':
            o = ax_sys.axial_structured(su=self.su)
            self.check = o['check']
            if not all(o['check'].values()):
                raise AssertionError('structured forms: %s' % o['check'])
            ins = [S.N, S.A, S.om, S.V]
            ent, num = o['step'], o['anum']
        else:
            o = ax_sys.toy_symbolic()
            ins = [o['N']]
            scale = {(0, 0): 1, (0, 1): self.su, (1, 0): 1/self.su, (1, 1): 1}   # y = (f, f'/s_u)
            ent = {ij: [c*scale[ij] for c in v] for ij, v in o['A'].items()}
            num = {ij: [c*scale[ij] for c in v] for ij, v in o['Anum'].items()}
        self.ins = ins
        self.tA = Tape([ent[ij][p] for p in range(3) for ij in order], ins, {S.k: K13})
        self.tAnum = Tape([num[ij][p] for p in range(3) for ij in order], ins, {S.k: K13})
        self.tden = Tape([ins[0]**2 - 1], ins, {S.k: K13})
        if kind == 'axial':
            self.tAcons = Tape([S.Acons], [S.N, S.om, S.V], {S.k: K13})
            self.tcen = Tape([o['cen'][ij][p] for p in range(3) for ij in order], [S.P, S.W, S.U], {S.k: K13})
            self.tDc = Tape([PL.load_expr(PL.FIELDS['pol_c'][3])], [S.P, S.W, S.U], {S.k: K13})
            self.fx = Tape([PL.load_expr(s) for s in PL.FIELDS['f']], [S.N, S.om, S.V], {S.k: K13})
            self.ftau = V2.tau_tape(K13)
        else:
            u = o['u']
            self.check = o['check']
            cen = [sp.cancel(o['cen'][ij][p]*scale[ij]) for p in range(3) for ij in order]
            self.tcen = Tape(cen, [u], {})

    @staticmethod
    def mats(vals):
        """12 outputs -> [M0, M1, M2] 2x2 lists"""
        return [[[vals[4*p + 2*i + j] for j in range(2)] for i in range(2)] for p in range(3)]


def conv_acb(c):
    return acb(sym_const(c))


# =========================================================================== background series per piece
def series_out(tape, ins, K):
    """Taylor coefficients 0..K of the tape outputs for input coefficient lists (ball, enclosing true coefficients).
    Plain Python lists (Tape.jet_order), no flint series: ctx.cap does not apply; every list has exactly K + 1 terms."""
    z = arb(0)
    tape.jet_init(K + 1, z, conv=lambda c: sym_const(c))
    for n in range(K + 1):
        for i in range(len(ins)):
            tape.C[i][n] = ins[i][n] if n < len(ins[i]) else z
        tape.jet_order(n)
    return [[tape.C[o][n] for n in range(K + 1)] for o in tape.outs]


def smul(a, b, K):
    return [sum((a[i]*b[n - i] for i in range(n + 1) if i < len(a) and n - i < len(b)), arb(0)) for n in range(K + 1)]


class Piece:
    """one transport piece: y' = A_t(t, kappa) y on t in [0, h] (t the piece variable), A_t = x'(t) A(state(t)).
    Data: A[p][n] (2x2 arb, p = kappa power, n = t-order), sups s[p] on |t| <= rho, h, rho; P[i][j] arb_poly in kappa."""

    def __init__(self, name, h, rho, Acoef, sups, K, tol):
        self.name = name; self.h = arb(h); self.rho = arb(rho)
        self.A = Acoef; self.s = sups; self.K = K; self.tol = tol
        self.q = abs(self.h)/self.rho
        if not self.q < 1:
            raise ArithmeticError('%s: |h| >= rho' % name)

    def rem(self, kabs, K=None):
        K = self.K if K is None else K
        a = self.s[0] + arb(kabs)*self.s[1] + arb(kabs)**2*self.s[2]
        return (a*self.rho).exp()*self.q**(K + 1)/(1 - self.q)

    def build_poly(self):
        """P(kappa) = sum_{n<=K} Phi_n h^n, exact polynomial in kappa (arb_poly entries)"""
        K = self.K
        KS = arb(KSCALE)          # polynomial in khat = kappa/KSCALE (|khat| <= 1 on every region): no 10^k amplification
        Am = [[[arb_poly([self.A[0][n][i][j], self.A[1][n][i][j]*KS, self.A[2][n][i][j]*KS*KS]) for j in range(2)] for i in range(2)] for n in range(K)]
        one = arb_poly([1]); zero = arb_poly([])
        Phi = [[[one, zero], [zero, one]]]
        for n in range(K):
            s_ = [[zero, zero], [zero, zero]]
            for m in range(n + 1):
                Am_ = Am[m]; Ph = Phi[n - m]
                for i in range(2):
                    for j in range(2):
                        s_[i][j] = s_[i][j] + Am_[i][0]*Ph[0][j] + Am_[i][1]*Ph[1][j]
            f = arb(1)/(n + 1)
            Phi.append([[s_[i][j]*f for j in range(2)] for i in range(2)])
        P = [[zero, zero], [zero, zero]]
        hp = arb(1)
        for n in range(K + 1):
            for i in range(2):
                for j in range(2):
                    P[i][j] = P[i][j] + Phi[n][i][j]*hp
            hp *= self.h
        self.P = P
        self.A0mid = [[float(self.A[p][0][i][j].mid()) for j in range(2)] for i in range(2) for p in range(3)]
        return P

    def frozen(self, kc):
        """float A(t = 0, kc) for the frame"""
        kc = complex(kc)
        M = np.zeros((2, 2), complex)
        for p in range(3):
            for i in range(2):
                for j in range(2):
                    M[i, j] += float(self.A[p][0][i][j].mid())*kc**p
        return M


def choose_K(sups_fn, rho, h, kmax, tol, Kcap):
    q = abs(h)/rho
    a = sups_fn[0] + arb(kmax)*sups_fn[1] + arb(kmax)**2*sups_fn[2]
    lg = float((a*rho).upper()) - math.log(1 - float(q.upper())) - math.log(tol)
    lq = -math.log(float(q.upper()))
    K = int(math.ceil(lg/lq))
    return max(4, K), K <= Kcap


def polar_cells(rho, nrad, nang, r_in=0):
    """acb boxes covering the closed disc |t| <= rho (annular sectors, inner disc |t| <= rho/nrad as one box when
    r_in == 0); returns (boxes, rmax) with rmax >= sup |t| over every box (the radius the tails must be valid on)"""
    pi = arb.pi()
    boxes = []
    rho = arb(rho)
    for ir in range(nrad):
        a = rho*ir/nrad; b = rho*(ir + 1)/nrad
        if ir == 0 and r_in == 0:
            boxes.append(csym(b))
            continue
        rr = (a + b)/2 + sym((b - a)/2)
        na = max(4, int(nang*(ir + 1)/nrad))
        for j in range(na):
            ph = (arb(j) + arb('0.5') + sym(arb('0.5')))*2*pi/na
            boxes.append(acb(rr*ph.cos(), rr*ph.sin()))
    rmax = arb(0)
    for bx in boxes:
        rmax = rmax.max(bx.abs_upper())
    return boxes, rmax


def disc_sups(evalA, boxes):
    """sup over the boxes (covering the disc) of row norms of the three kappa parts; None if any box fails"""
    s = [arb(0), arb(0), arb(0)]
    for bx in boxes:
        M = evalA(bx)
        if M is None:
            return None
        for p in range(3):
            v = rownorm(M[p])
            if not v.is_finite():
                return None
            s[p] = s[p].max(v)
    return s


# =========================================================================== E5 pieces
class E5Path:
    """the certified path; builds all Piece objects once (kappa-independent precompute)."""

    def __init__(self, SY, kmax=11.0, tol=1e-18, Kcap=90, rho_fac='0.88', nrad=3, nang=16, jF_dist='0.06', nuF='0.034',
                 d0='0.012', log=print, steps_slice=None):
        self.SY = SY; self.log = log; self.kmax = kmax; self.tol = tol; self.Kcap = Kcap
        self.rho_fac = arb(rho_fac); self.nrad = nrad; self.nang = nang
        self.tb = pickle.load(open(TUBE, 'rb'))
        self.ow = pickle.load(open(OUTW, 'rb'))
        self.info = {}
        T0 = time.time()
        self.lc = L0.lc_root(self.ow, lambda *a: None)
        self.info['x_lc'] = self.lc['x_lc'].str(20); self.info['root_step'] = self.lc['step']
        log('x_lc %s (root step %d) %.1fs' % (self.lc['x_lc'].str(12), self.lc['step'], time.time() - T0))
        self._choose_jF(arb(jF_dist))
        self.out_pieces = self._outward_pieces(steps_slice)
        log('outward pieces: %d (%.1fs)' % (len(self.out_pieces), time.time() - T0))
        self.tube_pieces = self._tube_pieces(steps_slice)
        log('tube pieces: %d (%.1fs)' % (len(self.tube_pieces), time.time() - T0))
        self.sm_tube = self._sonic_piece(self.tb['sm'], deser(self.tb['th0']), arb(self.tb['nu_s']), 'tube sonic manifold [th0, 0]')
        # matching point M = start of the outward tau leg: x_M = X_out(th_out) +- epsX_out (the true x there; both
        # charts put the sonic point at x = 0).  The outward manifold's tail disc (nu_s = 0.004, th_out = 0.0036) is
        # too thin for a transport disc, so the centre side reaches x_M on the TUBE manifold (nu_s = 0.06, 81 terms,
        # an enclosure of the same true solution): theta* with X_tube(theta*) = x_M (interval Newton).
        smo = self.ow['sm']; smt = self.tb['sm']
        Xo = [deser(p) for p in smo['X']]; tho = deser(smo['th'])
        xM = horner(Xo, tho) + sym(deser(smo['epsX']))
        Xt = [deser(p) for p in smt['X']]; dXt = deriv(Xt); epsXt = deser(smt['epsX'])
        th = arb(xM.mid())
        for it in range(40):
            th = arb((th - (horner(Xt, th) - xM)/horner(dXt, th)).mid())
        r = arb('1e-25') + xM.rad()*10
        for it in range(80):
            Th = th + sym(r)
            Nt = th - (horner(Xt, th) + sym(epsXt) - xM)/(horner(dXt, Th) + sym(epsXt/(arb(self.tb['nu_s']) - arb('0.01'))))
            if Th.contains_interior(Nt):
                break
            r *= 4
        else:
            raise ArithmeticError('theta*: interval Newton failed')
        self.theta_star = Nt
        self.info['matching point'] = dict(x_M=xM.str(20), theta_star=Nt.str(20), th_out=tho.str(10))
        self.sm_out = self._sonic_piece(smt, Nt, arb(self.tb['nu_s']), 'tube sonic manifold [0, theta*] (x = x_M)')
        log('sonic pieces (%.1fs)' % (time.time() - T0))
        self.centre = self._centre_setup()
        log('centre ready (%.1fs)' % (time.time() - T0))
        self.d0 = arb(d0)
        self.light = self._light_setup(arb(nuF))
        log('light-cone Frobenius ready (%.1fs)' % (time.time() - T0))
        self.near_pieces = self._near_lc_pieces(self.d0)
        log('near-cone pieces: %d, K %s (%.1fs)' % (len(self.near_pieces), [p.K for p in self.near_pieces], time.time() - T0))
        npoly = 0
        for pc in self.all_pieces():
            pc.build_poly(); npoly += 1
        log('transition polynomials built for %d pieces (%.1fs)' % (npoly, time.time() - T0))
        self.info['setup_s'] = round(time.time() - T0, 1)

    def all_pieces(self):
        return self.near_pieces + self.out_pieces + self.tube_pieces + [self.sm_tube, self.sm_out]

    # ---------------- A on a state
    def A_state_series(self, Ns, oms, Vs, xp, K):
        """series of x'(t) A(N, A, om, V) for state series (lists of arb)"""
        As = series_out(self.SY.tAcons, [Ns, oms, Vs], K)[0]
        out = series_out(self.SY.tA, [Ns, As, oms, Vs], K)
        out = [smul(xp, o, K) for o in out]
        return Sys.mats([[o[n] for n in range(K + 1)] for o in out])   # [p][i][j][n]

    def A_state_box(self, N, om, V, xp):
        A = self.SY.tAcons.eval([N, om, V], conv=conv_acb)[0]
        o = self.SY.tA.eval([N, A, om, V], conv=conv_acb)
        if not all(v.is_finite() for v in o):
            return None
        return Sys.mats([xp*v for v in o])

    @staticmethod
    def reorder(Mser, K):
        """[p][i][j][n] -> [p][n][i][j]"""
        return [[[[Mser[p][i][j][n] for j in range(2)] for i in range(2)] for n in range(K + 1)] for p in range(3)]

    # ---------------- generic step piece from a step record
    def _step_piece(self, name, kind, y0b, h, q):
        """background jets from step_cauchy2 on the certified start set y0b (contains the true state);
        Yc encloses the Taylor coefficients of every solution from y0b, Bx the complex a priori box on |t| <= Rc."""
        field = self.SY.ftau if kind == 'tau' else self.SY.fx
        ncomp = 4 if kind == 'tau' else 3
        KBG = 40
        Rc = abs(h)/q
        while True:
            y1, Bx, Yc = TY.step_cauchy2(field, list(y0b), h, KBG + 1, q=q, conv=lambda c: PM._c(c, arb(0)))
            Yi = [[Yc[i][n] for n in range(KBG + 2)] for i in range(ncomp)]
            totrad = [Bx[i].rad() for i in range(ncomp)]
            m_ = KBG + 2
            Yp = [[acb(a) for a in Yi[i]] for i in range(ncomp)]
            dYp = [[acb(a) for a in deriv(Yi[i])] for i in range(ncomp)]
            s = None
            for shrink in ('1', '0.8', '0.62', '0.48'):
                rho = self.rho_fac*Rc*arb(shrink)
                if not rho > arb('1.4')*abs(h):
                    break
                boxes, rmax = polar_cells(rho, self.nrad, self.nang)
                if not rmax < arb('0.98')*Rc:
                    continue
                r_ = rmax/Rc
                tailY = [totrad[i]*r_**m_/(1 - r_) for i in range(ncomp)]
                tailD = [totrad[i]/Rc*m_*r_**(m_ - 1)/(1 - r_)**2 for i in range(ncomp)]

                def evalA(tc):
                    yb = [horner(Yp[i], tc) + csym(tailY[i]) for i in range(ncomp)]
                    xp = (horner(dYp[3], tc) + csym(tailD[3])) if kind == 'tau' else acb(1)
                    return self.A_state_box(yb[0], yb[1], yb[2], xp)
                s = disc_sups(evalA, boxes)
                if s is not None:
                    break
            if s is None:
                raise ArithmeticError('%s: disc sup failed (a singularity of A inside every tried disc)' % name)
            K, ok = choose_K(s, rho, h, self.kmax, self.tol, self.Kcap)
            if K > 46:
                # finer cells give a tighter (still rigorous) sup near the cone; keep whichever needs the smaller K
                boxes2, rmax2 = polar_cells(rho, 3*self.nrad, 3*self.nang)
                if rmax2 < arb('0.98')*Rc:
                    r_ = rmax2/Rc
                    tailY = [totrad[i]*r_**m_/(1 - r_) for i in range(ncomp)]
                    tailD = [totrad[i]/Rc*m_*r_**(m_ - 1)/(1 - r_)**2 for i in range(ncomp)]
                    s2 = disc_sups(evalA, boxes2)
                    if s2 is not None:
                        K2, ok2 = choose_K(s2, rho, h, self.kmax, self.tol, self.Kcap)
                        if K2 < K:
                            s, K, ok = s2, K2, ok2
            if K <= KBG:
                break
            if not ok:
                raise ArithmeticError('%s: needs K = %d > Kcap' % (name, K))
            KBG = K + 2
        xp = deriv(Yi[3]) if kind == 'tau' else [arb(1)] + [arb(0)]*(K + 1)
        Mser = self.A_state_series(Yi[0][:K + 1], Yi[1][:K + 1], Yi[2][:K + 1], xp, K)
        pc = Piece(name, h, rho, self.reorder(Mser, K), s, K, self.tol)
        pc.kind = kind
        pc.y1 = y1
        return pc

    def _choose_jF(self, dist):
        """outward step index j_F whose START lies at distance about `dist` before x_lc; t_F = -(T + sum h)"""
        steps = self.ow['steps']; jr = self.lc['step']; T = self.lc['t']
        acc = T; j = jr
        while j > 0 and float(acc.mid()) < float(dist.mid()):
            j -= 1
            acc = acc + deser(steps[j]['h'])
            if steps[j]['kind'] != 'x':
                raise ArithmeticError('j_F would enter the tau leg')
        self.jF = j; self.tF = -acc
        self.info['j_F'] = j; self.info['t_F'] = self.tF.str(15)

    def _outward_pieces(self, sl):
        steps = self.ow['steps']; out = []
        idx = list(range(self.jF))
        if sl is not None:
            idx = idx[-sl:]
        for j in idx:
            st = steps[j]
            hull = [deser(p) for p in st['hull']]
            h = deser(st['h']); q = abs(h)/deser(st['Rc'])
            pc = self._step_piece('outward %d (%s)' % (j, st['kind']), st['kind'], hull, h, q)
            pc.index = j
            out.append(pc)
        return out

    def _tube_pieces(self, sl):
        steps = self.tb['steps']; out = []
        Q = arb(self.tb['Q'])
        idx = list(range(len(steps)))
        if sl is not None:
            idx = idx[:sl]
        for js in idx:
            st = steps[js]
            h = deser(st['h']); c = [deser(p) for p in st['c']]; setrad = deser(st['setrad'])
            ncomp = 4 if st['kind'] == 'tau' else 3
            y0b = [c[i] + sym(setrad*arb('1.01')) for i in range(ncomp)]
            pc = self._step_piece('tube %d (%s)' % (js, st['kind']), st['kind'], y0b, h, Q)
            pc.index = js
            out.append(pc)
        # true x at the tube end (autonomous system): x(tau-leg end) + sum of x-step h  (tau-leg end x from the last
        # tau step's own enclosure, not the stored X_SWITCH target)
        if sl is None:
            jt = max(i for i, pc in enumerate(out) if pc.kind == 'tau')
            xe = out[jt].y1[3]
            for pc in out[jt + 1:]:
                xe = xe + pc.h
            self.x_end = xe
            self.info['tube_end_x'] = xe.str(20)
            self.info['tube_end_x overlaps XBAR'] = bool(xe.overlaps(deser(self.tb['XBAR'])))
        else:
            self.x_end = deser(self.tb['XBAR'])
        return out

    def _sonic_piece(self, sm, th_end, nu_s, name):
        Y = [[deser(p) for p in sm['Y'][i]] for i in range(3)]
        X = [deser(p) for p in sm['X']]
        eps, epsX = deser(sm['eps']), deser(sm['epsX'])
        Kser = min(len(X) - 2, len(Y[0]) - 1)
        rho = arb('0.8')*nu_s
        dX = deriv(X)
        Yp = [[acb(a) for a in Y[i]] for i in range(3)]; dXp = [acb(a) for a in dX]
        boxes, rmax = polar_cells(rho, self.nrad, self.nang)
        if not rmax < arb('0.95')*nu_s:
            raise ArithmeticError('%s: cells leave the tail disc' % name)

        def evalA(tc):
            # eps bounds the series remainder on |theta| <= nu_s; Cauchy: |E'| <= epsX/(nu_s - |theta|)
            yb = [horner(Yp[i], tc) + csym(eps) for i in range(3)]
            xp = horner(dXp, tc) + csym(epsX/(nu_s - rmax))
            return self.A_state_box(yb[0], yb[1], yb[2], xp)
        s = disc_sups(evalA, boxes)
        if s is None:
            raise ArithmeticError('%s: disc sup failed' % name)
        if not abs(th_end) < rho/arb('1.4'):
            raise ArithmeticError('%s: |theta_end| too close to rho' % name)
        K, ok = choose_K(s, rho, th_end, self.kmax, self.tol, self.Kcap)
        if K > Kser:
            raise ArithmeticError('%s: needs K = %d > stored series order %d' % (name, K, Kser))
        Mser = self.A_state_series(Y[0][:K + 1], Y[1][:K + 1], Y[2][:K + 1], dX[:K + 1], K)
        pc = Piece(name, th_end, rho, self.reorder(Mser, K), s, K, self.tol)
        pc.kind = 'sonic'
        return pc

    # ---------------- centre
    def _centre_setup(self):
        tb = self.tb; cm = tb['cm']; SY = self.SY
        # coefficients stored as the exact string pair ('0', '0') are exact zeros (P(0) = 0 by the chart definition;
        # odd/even parity zeros, PROOF-SF-C1 sec. 0); everything else gets the usual deserialisation pad
        Z = [[arb(0) if (p[0] == '0' and p[1] == '0') else deser(p) for p in cm['Y'][i]] for i in range(3)]
        chi = [deser(p) for p in cm['X']]
        eps, epsX = deser(cm['eps']), deser(cm['epsX'])
        lam = deser(cm['lam']); nu = arb(cm['nu'])
        K = min(80, len(chi) - 1, len(Z[0]) - 1)
        capguard(K + 1)
        outs = series_out(SY.tcen, Z, K)
        dc = series_out(SY.tDc, Z, K)[0]
        fac = [x/lam for x in dc]
        Ms = [smul(fac, o, K) for o in outs]
        M = Sys.mats([[o[n] for n in range(K + 1)] for o in Ms])
        Mn = self.reorder(M, K)
        for i in range(2):                       # y = theta ytil:  M - I
            Mn[0][0][i][i] = Mn[0][0][i][i] - 1
        # the kappa parts vanish EXACTLY at theta = 0 (ax_sys S5b; P_0 = 0 exactly in the stored chart): M0 kappa-free
        for p in (1, 2):
            if not all(Mn[p][0][i][j] == 0 for i in range(2) for j in range(2)):
                raise ArithmeticError('centre: kappa part of M0 not exactly zero')
        nuc = arb('0.8')*nu
        Zp = [[acb(a) for a in Z[i]] for i in range(3)]
        M00 = [[[acb(Mn[p][0][i][j]) for j in range(2)] for i in range(2)] for p in range(3)]

        def evalM(tc):
            zv = [horner(Zp[i], tc) + csym(eps) for i in range(3)]
            o = SY.tcen.eval(zv, conv=conv_acb)
            f = SY.tDc.eval(zv, conv=conv_acb)[0]/lam
            m = Sys.mats([f*v for v in o])
            m[0] = [[m[0][i][j] - (1 if i == j else 0) for j in range(2)] for i in range(2)]
            return [[[m[p][i][j] - M00[p][i][j] for j in range(2)] for i in range(2)] for p in range(3)]
        sups = circle_sups(evalM, nuc)
        boxes, rmax = polar_cells(nuc, 4, 24)
        if not rmax < nu:
            raise ArithmeticError('centre cells leave the tail disc')
        if disc_sups(evalM, boxes) is None:
            raise ArithmeticError('centre: M not certified analytic on |theta| <= nu_c')
        # theta_m: ln theta + chi(theta) + xc = x_end (interval Newton, lin_stages pattern)
        xbar = self.x_end; xc = deser(tb['xc'])
        th = (xbar - xc).exp(); dchi = deriv(chi)
        for it in range(40):
            g = th.log() + horner(chi, th) + xc - xbar
            th = arb((th - g/(1/th + horner(dchi, th))).mid())
        r = arb('1e-30') + epsX*10 + xc.rad()*10 + xbar.rad()*10
        for it in range(80):
            Th = th + sym(r)
            g = th.log() + horner(chi, th) + sym(epsX) + xc - xbar
            Nt = th - g/(1/Th + horner(dchi, Th))
            if Th.contains_interior(Nt):
                break
            r *= 4
        else:
            raise ArithmeticError('theta_m interval Newton failed')
        # control: f(0) = Dc(Z(0))/lam = 1 (theta x_theta -> 1 at the centre)
        f0 = fac[0]
        self.info['centre'] = dict(K=K, nu_c=nuc.str(5), theta_m=Nt.str(15), f0=f0.str(10),
                                   M0til=[[Mn[0][0][i][j].str(8) for j in range(2)] for i in range(2)],
                                   sups=[s_.str(5) for s_ in sups])
        if not f0.contains(1):
            raise ArithmeticError('centre: theta x_theta at theta = 0 does not contain 1')
        return dict(M=Mn, K=K, sups=sups, nu=nuc, theta=Nt,
                    y0=[arb(1), arb(-2)/arb(int(self.SY.su.p))*arb(int(self.SY.su.q))], rank1=False)

    # ---------------- light cone
    def state_at_lc(self):
        lc = self.lc; ow = self.ow
        st = ow['steps'][lc['step']]; KT = ow['KT']
        Yc = [[deser(p) for p in st['Yc'][i]][:KT + 1] for i in range(3)]
        Rc = deser(st['Rc'])
        T = lc['t']
        tails = [deser(st['Bx'][i]).rad()*(T.abs_upper()/Rc)**(KT + 1)/(1 - T.abs_upper()/Rc) for i in range(3)]
        y0 = [horner(Yc[i], T) + sym(tails[i]) for i in range(3)]
        if not y0[0].contains(1):
            raise ArithmeticError('state at x_lc: N does not contain 1')
        return y0

    def _near_lc_pieces(self, d0, grow='0.08', rfrac='0.45'):
        """graded x-steps from t = -d0 to t = t_F(record node j_F) (t = x - x_lc): step k starts at distance d_k,
        length grow*d_k, transport disc rho = rfrac*d_k (keeps N = 1 outside), background re-expanded by
        step_cauchy2 from the enclosure carried along (state at x_lc from the root step polynomial)."""
        y = self.state_at_lc()
        yb = TY.step_cauchy2(self.SY.fx, y, -d0, 40, q=d0/arb('0.035'), conv=lambda c: PM._c(c, arb(0)))[0]
        D = -self.tF
        pieces = []
        d = d0
        k = 0
        while True:
            hlen = d*arb(grow)
            last = not (d + hlen*arb('1.15') < D)
            if last:
                hlen = D - d
            h = -hlen
            Rc = d*arb(rfrac)/self.rho_fac
            pc = self._step_piece('near-cone %d (x, d = %s)' % (k, d.str(4)), 'x', yb, h, abs(h)/Rc)
            pieces.append(pc)
            yb = pc.y1
            d = d + hlen; k += 1
            if last:
                break
        return pieces

    def _light_setup(self, nuF):
        SY = self.SY
        y0 = self.state_at_lc()
        KF = 120
        qF = arb('0.25')
        RcF = arb('0.04')                        # largest a priori disc step_cauchy2 certifies here (probe p_lc.py)
        hF = -RcF*qF
        _, Bx, Ycf = TY.step_cauchy2(SY.fx, y0, hF, KF, q=qF, conv=lambda c: PM._c(c, arb(0)))
        Yi = [[Ycf[i][n] for n in range(KF + 1)] for i in range(3)]
        Yi[0][0] = arb(1)                        # the true solution has N = 1 exactly at x_lc
        As = series_out(SY.tAcons, Yi, KF)[0]
        num = series_out(SY.tAnum, [Yi[0], As, Yi[1], Yi[2]], KF)
        den = series_out(SY.tden, [Yi[0], As, Yi[1], Yi[2]], KF)[0]
        if not den[0].contains(0):
            raise ArithmeticError('N^2 - 1 at x_lc')
        dsh = den[1:] + [arb(0)]                 # (N^2 - 1)/t : exact shift (constant term is exactly 0)
        KR = KF - 1
        Rser = []
        for o in num:
            qq = []
            for n in range(KR + 1):
                s_ = o[n]
                for i in range(n):
                    s_ = s_ - qq[i]*dsh[n - i]
                qq.append(s_/dsh[0])
            Rser.append(qq)
        M = Sys.mats([[o[n] for n in range(KR + 1)] for o in Rser])
        Mn = self.reorder(M, KR)
        nuc = nuF
        totrad = [Bx[i].rad() for i in range(3)]
        boxes, rmax = polar_cells(nuc, 6, 32)
        rmax = rmax.max(nuc*arb('1.05'))         # circle boxes (narc 128) stay inside 1.05 nu
        if not rmax < arb('0.98')*RcF:
            raise ArithmeticError('light cone: cells leave the a priori disc')
        rq = rmax/RcF
        tailY = [totrad[i]*rq**(KF + 1)/(1 - rq) for i in range(3)]
        # (N - 1)/t = sum_{n>=1} N_n t^(n-1) + E1,  |E1| <= sum_{n>KF} |N_n| |t|^(n-1) <= (totrad/RcF) rq^KF/(1 - rq)
        tailN1 = totrad[0]/RcF*rq**KF/(1 - rq)
        Yp = [[acb(a) for a in Yi[i]] for i in range(3)]
        N1p = Yp[0][1:]
        R00 = [[[acb(Mn[p][0][i][j]) for j in range(2)] for i in range(2)] for p in range(3)]

        def evalM(tc):
            yb = [horner(Yp[i], tc) + csym(tailY[i]) for i in range(3)]
            A = SY.tAcons.eval(yb, conv=conv_acb)[0]
            o = SY.tAnum.eval([yb[0], A, yb[1], yb[2]], conv=conv_acb)
            d_ = (horner(N1p, tc) + csym(tailN1))*(yb[0] + 1)          # (N^2 - 1)/t, valid also at t = 0
            if not d_.is_finite() or (d_.real.contains(0) and d_.imag.contains(0)):
                return None
            m = Sys.mats([v/d_ for v in o])
            return [[[m[p][i][j] - R00[p][i][j] for j in range(2)] for i in range(2)] for p in range(3)]
        sups = circle_sups(evalM, nuc)
        if disc_sups(evalM, boxes) is None:
            raise ArithmeticError('light cone: R not certified analytic on |t| <= nu')
        lam = -Yi[0][1]                           # N_x(x_lc) = N a, lambda_lc = -N_x
        tL = -self.d0
        self.info['light'] = dict(KF=KF, nu=nuc.str(5), RcF=RcF.str(5), lam_lc=lam.str(15), t_L=tL.str(15),
                                  sups=[s_.str(5) for s_ in sups])
        if not abs(tL) < nuc:
            raise ArithmeticError('t_L outside the Frobenius circle')
        return dict(M=Mn, K=KR, sups=sups, nu=nuc, t=tL, lam=lam, rank1=True, v0='col2')


def circle_sups(evalM, nu, narc=128):
    """sup over |t| = nu of ||M_p(t) - M_p(0)|| for p = 0 and ||M_p(t)|| for p = 1, 2 (evalM returns the differences
    for p = 0 and the full kappa parts for p = 1, 2 -- callers subtract the t = 0 values for every p, which is
    an upper bound for the kappa parts too since ||M_p - M_p0|| <= ||M_p|| + ||M_p0|| is not needed: we use the
    difference form for all p, valid for the Cauchy estimate of the coefficients n >= 1)"""
    s = [arb(0), arb(0), arb(0)]
    pi = arb.pi()
    for a in range(narc):
        th = (arb(a) + arb('0.5') + sym(arb('0.5')))*2*pi/narc
        tb_ = acb(th.cos(), th.sin())*nu
        M = evalM(tb_)
        if M is None:
            raise ArithmeticError('circle enclosure failed')
        for p in range(3):
            v = rownorm(M[p])
            if not v.is_finite():
                raise ArithmeticError('circle enclosure not finite')
            s[p] = s[p].max(v)
    return s


# =========================================================================== regular-singular series as TMs
def regsing_tm(cx, RS, y0, t, S_set=None, sigma=None):
    """series solution of t y' = M(t, kappa) y at the point t (acb or arb) as a TM vector.
    RS: dict(M=[p][n][i][j], K, sups, nu, rank1).  Centre: y0 = exact start (list of arb), S_set None.
    Light cone: y0 ignored; v0 = (-R0_12, R0_11); normalisation prod_{j in S_set} (j - sigma)."""
    K = RS['K']; Mx = RS['M']
    kap = cx.kappa()
    k2 = cx.mul(kap, kap)

    def Mtm(n):
        out = [[None, None], [None, None]]
        for i in range(2):
            for j in range(2):
                a0 = cx.const(acb(Mx[0][n][i][j]))
                a1 = cx.scale(acb(Mx[1][n][i][j]), kap)
                a2 = cx.scale(acb(Mx[2][n][i][j]), k2)
                out[i][j] = cx.add(cx.add(a0, a1), a2)
        return out
    Ms = [Mtm(n) for n in range(K + 1)]
    R0 = Ms[0]
    I2 = [[cx.const(1), cx.const(0)], [cx.const(0), cx.const(1)]]
    if RS['rank1']:
        sig = cx.add(R0[0][0], R0[1][1])                       # sigma(kappa), affine
        # sigma is affine: read its two coefficients exactly from the TM (R = 0 up to rounding)
        s0 = acb(sig.p.coeffs()[0]) if sig.p.length() > 0 else acb(0)
        s1 = acb(sig.p.coeffs()[1]) if sig.p.length() > 1 else acb(0)
        if not sig.R == 0:
            raise ArithmeticError('sigma TM not exact')
        S_set = sorted(S_set or [])
        # rank-one residue: both adjugate columns lie in ker R0; the problem fixes which one (proved nonvanishing)
        if RS['v0'] == 'col2':
            y = [[cx.scale(-1, R0[0][1]), R0[0][0]]]            # v0 = (-R0_12, R0_11)
        else:
            y = [[R0[1][1], cx.scale(-1, R0[1][0])]]            # v0 = (R0_22, -R0_21)
        jm = lambda j: TM(acb_poly([acb(j) - s0, -s1]), 0)      # (j - sigma) as TM in d (exact affine)
    else:
        y = [[cx.const(y0[0]), cx.const(y0[1])]]
    for n in range(1, K + 1):
        r = [cx.const(0), cx.const(0)]
        for m in range(1, n + 1):
            v = y[n - m]
            if RS['rank1']:
                for j in S_set:
                    if n - m < j < n:
                        v = [cx.mul(jm(j), v[0]), cx.mul(jm(j), v[1])]
            r = [cx.add(r[i], c) for i, c in enumerate(matvec(cx, Ms[m], v))]
        nM = [[cx.sub(cx.scale(n, I2[i][j]), R0[i][j]) for j in range(2)] for i in range(2)]
        if RS['rank1']:
            adj = [[nM[1][1], cx.scale(-1, nM[0][1])], [cx.scale(-1, nM[1][0]), nM[0][0]]]
            v = matvec(cx, adj, r)
            v = [cx.scale(acb(1)/n, c) for c in v]
            if n not in S_set:
                rc = cx.recip_affine(acb(n) - s0, -s1)
                v = [cx.mul(rc, c) for c in v]
        else:
            # centre: kappa-free M0 (exact up to the ball of the stored data): plain inverse
            M0c = [[acb(n)*(1 if i == j else 0) - acb(Mx[0][0][i][j]) for j in range(2)] for i in range(2)]
            inv = acb_mat(M0c).inv()
            v = [cx.add(cx.scale(inv[i, 0], r[0]), cx.scale(inv[i, 1], r[1])) for i in range(2)]
        y.append(v)
    if RS['rank1']:
        for n in range(K + 1):
            for j in S_set:
                if j > n:
                    y[n] = [cx.mul(jm(j), y[n][0]), cx.mul(jm(j), y[n][1])]
    # majorant tail (lin_stages.RegSing.evaluate, with |kappa| <= kabs on the disc)
    t = acb(t); at = abs(t).upper(); nu = RS['nu']
    nup = ((nu + at)/2).min(2*at) if at > 0 else nu/2
    kab = arb(cx.kabs)
    Snorm = RS['sups'][0] + kab*RS['sups'][1] + kab**2*RS['sups'][2]
    rat = nup/nu
    mup = Snorm*rat/(1 - rat)
    R0n = rownorm([[cx.sup(R0[i][j]) for j in range(2)] for i in range(2)])
    if RS['rank1']:
        # (n - R0)^{-1} = (I + R0/(n - sigma))/n ;  |n - sigma| >= n - |s0| - |s1| (r + 0) lower bound
        def cn(n):
            gap = (acb(n) - s0).abs_lower() - s1.abs_upper()*cx.r
            return (1 + R0n/gap)/n if gap > 0 else None
    else:
        def cn(n):
            return 1/(arb(n) - R0n) if arb(n) > R0n else None
    n0 = K
    while True:
        c = cn(n0 + 1)
        if c is not None and c*mup <= 1:
            break
        n0 += 1
        if n0 > 100000:
            raise ArithmeticError('regsing tail: n0 too large')
    if n0 > K:
        raise ArithmeticError('regsing tail needs K >= %d (have %d)' % (n0, K))
    B = arb(0)
    for n in range(K + 1):
        B = B.max(max(cx.sup(y[n][i]) for i in range(2))*nup**n)
    qq = at/nup
    tail = B*qq**(K + 1)/(1 - qq)
    out = []
    for i in range(2):
        s_ = cx.const(0)
        for n in range(K, -1, -1):
            s_ = cx.add(cx.scale(t, s_), y[n][i])
        out.append(TM(s_.p, s_.R + tail))
    return out, dict(tail=tail.str(5), n0=n0, mup=mup.str(5))


# =========================================================================== transport chain
def frame_for(M):
    """float eigenvector frame of a 2x2 complex matrix; identity if ill-conditioned"""
    try:
        w, V = np.linalg.eig(M)
        if np.linalg.cond(V) > 1e6 or not np.all(np.isfinite(V)):
            return np.eye(2, dtype=complex)
        return V
    except Exception:
        return np.eye(2, dtype=complex)


def exact_mat(V):
    T = [[acb(float(V[i, j].real), float(V[i, j].imag)) for j in range(2)] for i in range(2)]
    Ti = acb_mat(T).inv()
    return T, [[Ti[i, j] for j in range(2)] for i in range(2)]


def piece_tm(cx, pc, backward):
    """2x2 TM of P (forward) or adj(P) (backward) for this disc"""
    E = pc.rem(cx.kabs)
    P = [[cx.poly_kappa(pc.P[i][j]) for j in range(2)] for i in range(2)]
    P = [[TM(P[i][j].p, P[i][j].R + E) for j in range(2)] for i in range(2)]
    if backward:
        P = [[P[1][1], cx.scale(-1, P[0][1])], [cx.scale(-1, P[1][0]), P[0][0]]]
    return P


def conj_tm(cx, Ti_next, P, T_cur):
    """Ti_next P T_cur (exact constant matrices)"""
    PT = [[cx.add(cx.scale(T_cur[0][j], P[i][0]), cx.scale(T_cur[1][j], P[i][1])) for j in range(2)] for i in range(2)]
    return [[cx.add(cx.scale(Ti_next[i][0], PT[0][j]), cx.scale(Ti_next[i][1], PT[1][j])) for j in range(2)] for i in range(2)]


def _cmid(z):
    z = acb(z)
    return complex(float(z.real.mid()), float(z.imag.mid()))


def transport(cx, Y, chain, frames=True, dbg=None):
    """Y (TM vector in plain variables) through chain = [(piece, backward)].
    Frames (any invertible constant matrix is rigorous; it only picks the coordinates of the enclosure):
    T_s = V_s diag(p_s), V_s float eigenvectors of A(node, kc), p_s unit floats chosen so that the diagonal of
    N_s = T_{s+1}^{-1} P_s T_s is (nearly) real positive.  Why: acb balls are rectangles and every complex rotation of
    a vector box inflates it by up to |cos| + |sin|; on Re kappa = 0 the modes rotate by the total phase
    |kappa| int |mu| (probe p_dbg3: radius x 1e7 at kappa = 10i).  The phase now sits in the exact frame (applied to
    the thin matrix N, whose radius is ~1e-15), and the propagated vector only meets near-real-positive diagonals."""
    kc = _cmid(cx.kc)
    if frames:
        Tf = frame_for(chain[0][0].frozen(kc)) if chain else np.eye(2, dtype=complex)
        T, Ti = exact_mat(Tf)
        Z = [cx.add(cx.scale(Ti[i][0], Y[0]), cx.scale(Ti[i][1], Y[1])) for i in range(2)]
    else:
        Z = Y
    for k_, (pc, bw) in enumerate(chain):
        P = piece_tm(cx, pc, bw)
        if frames:
            if k_ + 1 < len(chain):
                V2 = frame_for(chain[k_ + 1][0].frozen(kc))
                Pf = np.array([[_cmid(P[i][j].p.coeffs()[0]) if P[i][j].p.length() else 0j for j in range(2)] for i in range(2)])
                M0 = np.linalg.solve(V2, Pf @ Tf)
                if abs(M0[0, 1]) + abs(M0[1, 0]) > abs(M0[0, 0]) + abs(M0[1, 1]):
                    V2 = V2[:, ::-1].copy(); M0 = np.linalg.solve(V2, Pf @ Tf)
                ph = [M0[i, i]/abs(M0[i, i]) if (np.isfinite(M0[i, i]) and abs(M0[i, i]) > 0) else 1.0 for i in range(2)]
                Tf2 = V2 @ np.diag(ph)
            else:
                Tf2 = np.eye(2, dtype=complex)
            T2, Ti2 = exact_mat(Tf2)
            P = conj_tm(cx, Ti2, P, T)
            T, Ti, Tf = T2, Ti2, Tf2
        Z = matvec(cx, P, Z)
        if dbg is not None:
            dbg.append((pc.name, max(float(Z[i].R.upper()) for i in range(2)), max(float(cx.sup(Z[i]).upper()) for i in range(2)), max(float(acb(Z[i].p.coeffs()[0]).rad().upper()) if Z[i].p.length() else 0.0 for i in range(2)), max(float(abs(acb(Z[i].p.coeffs()[0])).mid()) if Z[i].p.length() else 0.0 for i in range(2))))
    return Z


class Axial:
    """D_2(kappa) = det[Y_lc(M), Y_c(M)] on the E5 path"""

    def __init__(self, path, S_set=(1,), frames=True):
        self.path = path; self.S_set = list(S_set); self.frames = frames
        # Matching point M = the Frobenius evaluation point t_L = -d0 (x = x_lc - d0).  The centre-regular solution is the
        # DOMINANT one outward (r^(l-1) vs r^(-l-2)) and is carried the whole way; the light-cone-analytic solution is
        # subdominant inward (the t^sigma branch grows away from the cone), so it is NOT transported (probe p_dbg2:
        # carrying it inward lost every digit at kappa = 10i).
        # centre side: theta_m (x = x_end) -> tube steps 719 .. 0 backward -> tube sm th0 -> 0 (backward) -> 0 -> theta*
        #   (= outward tau-leg start) -> outward steps 0 .. j_F - 1 forward -> near-cone pieces backward (last first)
        self.chain_c = ([(pc, True) for pc in reversed(path.tube_pieces)] + [(path.sm_tube, True), (path.sm_out, False)]
                        + [(pc, False) for pc in path.out_pieces] + [(pc, True) for pc in reversed(path.near_pieces)])

    def D_tm(self, kc, r, J, dbg=None, info=None):
        cx = Ctx(kc, r, J)
        Ylc, ilc = regsing_tm(cx, self.path.light, None, self.path.light['t'], S_set=self.S_set)
        Yc, ic = regsing_tm(cx, self.path.centre, self.path.centre['y0'], self.path.centre['theta'])
        Yc = [cx.scale(acb(self.path.centre['theta']), c) for c in Yc]       # y = theta ytil
        Zc = transport(cx, Yc, self.chain_c, self.frames, dbg)
        if info is not None:
            info.update(light=ilc, centre=ic)
        return det2(cx, Ylc, Zc), cx


# =========================================================================== flat toy (G App. D), same modules
class ToyPath:
    """centre u-series at u_c, transport x-steps (exact background N = e^{-x}), light-cone Frobenius at x = 0."""

    def __init__(self, SY, u_c='0.25', x_F='-0.3', nsteps=16, kmax=3.5, tol=1e-18, Kcap=80, nuF='0.9', nuC='0.75', log=print):
        self.SY = SY; self.kmax = kmax; self.tol = tol
        uc = arb(u_c); xF = arb(x_F); xc = uc.log()
        self.info = {}
        # transport pieces from x_c to x_F (forward in x), h equal
        H = (xF - xc)/nsteps
        self.pieces = []
        for s in range(nsteps):
            x0 = xc + H*s
            self.pieces.append(self._piece(x0, H, Kcap, 'toy step %d' % s))
        # centre: u dy/du = M(u) y, M - 2I, start (1, 2/s_u... ) exact
        K = 90
        cu = [arb(0), arb(1)] + [arb(0)]*(K - 1)
        outs = series_out(SY.tcen, [cu], K)
        M = Sys.mats([[o[n] for n in range(K + 1)] for o in outs])
        Mn = E5Path.reorder(M, K)
        for i in range(2):
            Mn[0][0][i][i] = Mn[0][0][i][i] - 2
        for p in (1, 2):
            if not all(Mn[p][0][i][j] == 0 for i in range(2) for j in range(2)):
                raise ArithmeticError('toy centre: kappa part of M0 not exactly zero')
        nuc = arb(nuC)
        M00 = [[[acb(Mn[p][0][i][j]) for j in range(2)] for i in range(2)] for p in range(3)]

        def evalM(tc):
            o = SY.tcen.eval([tc], conv=conv_acb)
            if not all(v.is_finite() for v in o):
                return None
            m = Sys.mats(o)
            m[0] = [[m[0][i][j] - (2 if i == j else 0) for j in range(2)] for i in range(2)]
            return [[[m[p][i][j] - M00[p][i][j] for j in range(2)] for i in range(2)] for p in range(3)]
        su = SY.su
        if disc_sups(evalM, polar_cells(nuc, 4, 24)[0]) is None:
            raise ArithmeticError('toy centre: M not analytic on the disc')
        self.centre = dict(M=Mn, K=K, sups=circle_sups(evalM, nuc), nu=nuc, theta=uc,
                           y0=[arb(1), arb(2)*arb(int(su.q))/arb(int(su.p))], rank1=False)
        # light cone at x = 0: N = e^{-t}
        KF = 90
        Ns = [arb(-1)**n/arb(math.factorial(n)) for n in range(KF + 1)]
        num = series_out(SY.tAnum, [Ns], KF)
        den = series_out(SY.tden, [Ns], KF)[0]
        if not den[0].contains(0):
            raise ArithmeticError('toy N^2 - 1 at 0')
        dsh = den[1:] + [arb(0)]
        KR = KF - 1
        Rser = []
        for o in num:
            qq = []
            for n in range(KR + 1):
                s_ = o[n]
                for i in range(n):
                    s_ = s_ - qq[i]*dsh[n - i]
                qq.append(s_/dsh[0])
            Rser.append(qq)
        Mn = E5Path.reorder(Sys.mats([[o[n] for n in range(KR + 1)] for o in Rser]), KR)
        R00 = [[[acb(Mn[p][0][i][j]) for j in range(2)] for i in range(2)] for p in range(3)]
        nuf = arb(nuF)

        # (N^2 - 1)/t = (e^{-2t} - 1)/t = sum_{n>=0} (-2)^(n+1) t^n/(n+1)!  (exact coefficients), tail on |t| <= 1.1 nu
        KD = 80
        gco = [acb(arb(-2)**(n + 1)/arb(math.factorial(n + 1))) for n in range(KD + 1)]
        rr = nuf*arb('1.1')
        gtail = arb(2)**(KD + 2)*rr**(KD + 1)/arb(math.factorial(KD + 2))*2      # ratio of successive terms <= 1/2

        def evalR(tc):
            N = (-tc).exp()
            o = SY.tAnum.eval([N], conv=conv_acb)
            if not tc.abs_upper() <= rr:
                return None
            d_ = horner(gco, tc) + csym(gtail)
            if not d_.is_finite() or (d_.real.contains(0) and d_.imag.contains(0)):
                return None
            m = Sys.mats([v/d_ for v in o])
            return [[[m[p][i][j] - R00[p][i][j] for j in range(2)] for i in range(2)] for p in range(3)]
        if disc_sups(evalR, polar_cells(nuf, 6, 32)[0]) is None:
            raise ArithmeticError('toy light cone: R not analytic on the disc')
        self.light = dict(M=Mn, K=KR, sups=circle_sups(evalR, nuf), nu=nuf, t=xF, lam=arb(1), rank1=True, v0='col1')
        self.chain_c = [(pc, False) for pc in self.pieces]          # centre solution carried outward (dominant), as Axial

    def _piece(self, x0, h, Kcap, name):
        SY = self.SY
        rho = abs(h)*arb('3')
        n0 = (-x0).exp()

        def evalA(tc):
            N = n0*(-tc).exp()
            o = SY.tA.eval([N], conv=conv_acb)
            if not all(v.is_finite() for v in o):
                return None
            return Sys.mats(o)
        s = disc_sups(evalA, polar_cells(rho, 3, 16)[0])
        if s is None:
            raise ArithmeticError('toy disc sup failed')
        K, ok = choose_K(s, rho, h, self.kmax, self.tol, Kcap)
        if not ok:
            raise ArithmeticError('toy: K > Kcap')
        Ns = [n0*arb(-1)**n/arb(math.factorial(n)) for n in range(K + 1)]
        out = series_out(SY.tA, [Ns], K)
        pc = Piece(name, h, rho, E5Path.reorder(Sys.mats([[o[n] for n in range(K + 1)] for o in out]), K), s, K, self.tol)
        pc.kind = 'x'
        pc.build_poly()
        return pc


class Toy:
    def __init__(self, path, S_set=(1, 2), frames=True):
        self.path = path; self.S_set = list(S_set); self.frames = frames

    def D_tm(self, kc, r, J, dbg=None):
        cx = Ctx(kc, r, J)
        Ylc, _ = regsing_tm(cx, self.path.light, None, self.path.light['t'], S_set=self.S_set)
        Yc, _ = regsing_tm(cx, self.path.centre, self.path.centre['y0'], self.path.centre['theta'])
        u = acb(self.path.centre['theta'])
        Yc = [cx.scale(u*u, c) for c in Yc]
        Zc = transport(cx, Yc, self.path.chain_c, self.frames, dbg)
        return det2(cx, Ylc, Zc), cx


# =========================================================================== count-engine adapters
class DTM:
    """object for argument_count: .eval(u) encloses D(kc + u) on |u| <= r, .R remainder"""
    def __init__(self, tm, cx):
        self.tm = tm; self.cx = cx; self.R = tm.R

    def eval(self, u):
        """encloses D(kc + d) for every d in the box u with |d| <= r (the polynomial is evaluated on the whole box,
        a superset; segment_arg only needs the segment points, which lie in the disc by check_plan)"""
        return self.cx.evalb(self.tm, acb(u))


def make_callables(obj, J=24, J_point=0):
    def tm_at(kc, r):
        tm, cx = obj.D_tm(kc, arb(r), J)
        return DTM(tm, cx)

    def point(z):
        tm, cx = obj.D_tm(z, arb(0), J_point)
        return cx.evalb(tm, acb(0))
    return tm_at, point
