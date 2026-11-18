"""l2b_ext.py - build seat build-4, v2 extension for BAR-SF-L2b (seal fc89881caf), 2026-11-17.

New file; imports the PINNED v1 library l2_lib.py (40c831a3...) read-only and never edits it.

What v2 changes (DIAG-L2.md sec. 2, 4):
 T1  Tight background jets.  v1 computed the background Taylor coefficients by ball jets from the certified start BALL;
     dependency inflates the radius with the order (step 300: relative radius 2.6e-15 at n = 0, 5e-7 at n = 20, 1.6e-4 at
     n = 30), and the high orders dominate the kappa^k coefficients of the transition at |kappa| ~ 6.  v2: jets of the THIN
     centre y_c = mid(start ball) (rounding only) plus a rigorous deviation bound for the true solution:
       e(t) = y_true(t) - y_c(t),  |e(t)|_inf <= Delta e^{L |t|}  on |t| <= Rc   (Gronwall along rays; L = sup of the
       row-sum norm of the field Jacobian over the a priori complex box Bc, which contains both solutions for |t| <= Rc:
       step_cauchy2's a priori box for the start ball; mean-value form along the segment, Bc convex),
     so by Cauchy  |Y_n(true) - Y_n(thin)| <= Delta e^{L Rc} / Rc^n  (n >= 0), Delta = max radius of the start ball.
     Installed by rebinding taylor.step_cauchy2 while the v2 path is BUILT (restored in finally); the pinned l2_lib code
     calls it through the module attribute, so its own logic (cells, sups, K choice, Frobenius) is unchanged.
 T2  Near-cone background states read from the recorded outward steps (polynomial + Cauchy tail at the offset), not
     chained from step to step (v1 chained step_cauchy2 output balls, which grow).
 T3  Lohner-type carried centre solution: the transported vector is held as Y = c + B r (c thin midpoint vector of
     TMs-with-zero-remainder, B an orthonormal float frame used as an exact matrix, r a TM error vector) and every step
     re-orthonormalises B along the image of the dominant column (weighted QR, lin_stages.transport pattern), so the
     injected errors stay attached to their mode instead of filling the box.  The plain-frame transport of v1 remains
     available (lohner=False) for comparison.
 T4  Per-region matching distance: the light-cone solution may be carried outward through the first k near-cone pieces
     (forward) and the centre solution inward through the rest; the zero set of D does not depend on the matching point
     (PROOF-NORM-v2 N4: transport is an isomorphism, det P_s is zero-free).
 T5  Flat toy with the SINGULAR centre class r^{-l-1}: start (1, -3/s_u) (exponent -3 of M0 = [[0, 1], [6, -1]]),
     parity-aware centre recursion (odd orders are exactly zero: every coefficient of M at odd order is checked to be the
     exact zero; the resonance at order 5 = 2l + 1 is odd, so it never meets the recursion).  The light-cone side is the
     v1 module unchanged; for kappa in the PL2b-2 discs (Re kappa > 0) the system exponent rho = -kappa < 0 has no
     collision and the normaliser (1 - rho)(2 - rho) has no zero, so a zero of D is a parallelism of the two solutions.
"""
import os, sys, math, time, contextlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l2_lib as LL                       # pinned v1 (sets ctx.cap at its module top)
from l2_lib import TM, Ctx, matvec, det2, rownorm, polar_cells, disc_sups, circle_sups, Sys, E5Path, horner, deser
from flint import arb, acb, acb_mat, acb_poly, ctx
from core import Dual, sym, csym, sym_const
import taylor as TY

_ORIG_STEP = TY.step_cauchy2


def _conv(c):
    return acb(sym_const(c))


def jac_sup(tape, Bc):
    """sup over the complex box Bc of the row-sum norm of the field Jacobian (Dual evaluation)"""
    n = len(Bc)
    xs = [Dual(Bc[i], [acb(1) if j == i else acb(0) for j in range(n)]) for i in range(n)]
    out = tape.eval(xs, conv=_conv)
    L = arb(0)
    for o in out[:n]:
        row = sum((acb(d).abs_upper() for d in o.d), arb(0))
        if not row.is_finite():
            raise ArithmeticError('Jacobian bound not finite')
        L = L.max(row)
    return L


def tight_step_cauchy2(tape, y0, h, K, q=arb('0.25'), conv=None, maxit=60):
    """same contract as taylor.step_cauchy2 (y1, Bx, Yc): Yc enclose the Taylor coefficients of every solution from the
    ball y0, built as thin-centre jets + Delta e^{L Rc}/Rc^n (T1)."""
    y1b, Bx, _ = _ORIG_STEP(tape, list(y0), h, 2, q=q, conv=conv, maxit=maxit)
    q = arb(q)
    Rc = abs(h)/q
    yt = [arb(arb(v).mid()) for v in y0]
    Yt = TY.jets(tape, yt, K, conv or (lambda c: LL.PM._c(c, arb(0))))
    Delta = arb(0)
    for v in y0:
        Delta = Delta.max(arb(v).rad())
    Bc = [acb(arb(b.mid())) + csym(b.rad()) for b in Bx]
    if len(Bc) > len(Yt):
        Bc = Bc[:len(Yt)]
    L = jac_sup(tape, Bc)
    E = Delta*(L*Rc).exp()
    Yc = [[Yt[i][n] + sym(E/Rc**n) for n in range(len(Yt[i]))] for i in range(len(Yt))]
    geo = q**(K + 1)/(1 - q)
    y1 = [horner(Yc[i][:K + 1], h) + sym(Bx[i].rad()*geo) for i in range(len(Yc))]
    return y1, Bx, Yc


@contextlib.contextmanager
def tight_jets():
    TY.step_cauchy2 = tight_step_cauchy2
    try:
        yield
    finally:
        TY.step_cauchy2 = _ORIG_STEP


class E5PathB(E5Path):
    """v1 path with T1 (tight jets) and T2 (near-cone states from the records)"""

    def __init__(self, SY, tight=False, **kw):
        if tight:
            with tight_jets():
                super().__init__(SY, **kw)
        else:
            super().__init__(SY, **kw)
        self.info['v2'] = ('T1 tight jets ON' if tight else 'T1 OFF (v1 ball jets)') + '; T2 near-cone states from records'

    def state_at_dist(self, d):
        """background state at x = x_lc - d from the recorded outward step containing it (polynomial + Cauchy tail)"""
        steps = self.ow['steps']; jr = self.lc['step']; KT = self.ow['KT']
        acc = self.lc['t']                          # distance from the root step's start to x_lc
        j = jr
        while float((acc - d).mid()) < 0:
            j -= 1
            acc = acc + deser(steps[j]['h'])
        tau = acc - d                               # offset inside step j (from its start)
        st = steps[j]
        Yc = [[deser(p) for p in st['Yc'][i]][:KT + 1] for i in range(3)]
        Rc = deser(st['Rc'])
        if not (tau > -arb('1e-30')) or not (tau < deser(st['h'])*arb('1.000001')):
            raise ArithmeticError('state_at_dist: offset outside step %d' % j)
        r_ = abs(tau).upper()/Rc
        tails = [deser(st['Bx'][i]).rad()*r_**(KT + 1)/(1 - r_) for i in range(3)]
        return [horner(Yc[i], tau) + sym(tails[i]) for i in range(3)], j

    def _near_lc_pieces(self, d0, grow='0.08', rfrac='0.45'):
        D = -self.tF
        pieces = []
        d = d0; k = 0
        while True:
            hlen = d*arb(grow)
            last = not (d + hlen*arb('1.15') < D)
            if last:
                hlen = D - d
            h = -hlen
            Rc = d*arb(rfrac)/self.rho_fac
            yb, jrec = self.state_at_dist(d)
            pc = self._step_piece('near-cone %d (x, d = %s)' % (k, d.str(4)), 'x', yb, h, abs(h)/Rc)
            pc.d = d; pc.rec = jrec
            pieces.append(pc)
            d = d + hlen; k += 1
            if last:
                break
        return pieces


# =========================================================================== T3 Lohner-type carried vector
def _cm(z):
    return complex(float(acb(z).real.mid()), float(acb(z).imag.mid()))


def lohner_transport(cx, Y, chain, dbg=None, structured=False):
    """Y = c + B r through chain [(piece, backward)].  c: thin TM vector (remainder 0, coefficient balls kept thin by
    re-centring), B: float orthonormal 2x2 used as an exact acb matrix, r: TM error vector (holds every radius and every
    remainder).  Step: Y' = P Y = P c + (P B) r; new c = mid(P c); new B = QR-orthonormalised P B (dominant column first);
    r' = B'^{-1}(P B) r + B'^{-1}(P c - c').  All products are TM arithmetic (rigorous); B' is any invertible matrix."""
    def mid_tm(t):
        return TM(acb_poly([acb(acb(c).real.mid(), acb(c).imag.mid()) for c in t.p.coeffs()]) if t.p.length() else acb_poly([]), 0)

    def sub_tm(a, b):
        return cx.sub(a, b)
    c = [mid_tm(Y[i]) for i in range(2)]
    r = [sub_tm(Y[i], c[i]) for i in range(2)]               # error vector in the B = I frame
    Bf = np.eye(2, dtype=complex)
    B = [[acb(1), acb(0)], [acb(0), acb(1)]]
    for pc, bw in chain:
        P = LL.piece_tm(cx, pc, bw)
        Pc = matvec(cx, P, c)
        PB = [[cx.add(cx.scale(B[0][j], P[i][0]), cx.scale(B[1][j], P[i][1])) for j in range(2)] for i in range(2)]
        # float matrix of P B at d = 0 and the error weights -> new orthonormal frame (dominant direction first)
        Ma = [[acb(PB[i][j].p.coeffs()[0]) if PB[i][j].p.length() else acb(0) for j in range(2)] for i in range(2)]
        sc = max((abs(v).mid() for row in Ma for v in row), key=lambda a: float(a.log()) if a > 0 else -1e300)
        sc = sc if sc > 0 else arb(1)
        M = np.array([[_cm(Ma[i][j]/sc) for j in range(2)] for i in range(2)])     # scale-free float image (no overflow)
        w = [float(cx.sup(r[j]).upper()) + 1e-300 for j in range(2)]
        order = sorted(range(2), key=lambda j: -np.linalg.norm(M[:, j])*w[j])
        Q, _ = np.linalg.qr(M[:, order])
        if not np.all(np.isfinite(Q)) or abs(np.linalg.det(Q)) < 1e-12:
            Q = np.eye(2, dtype=complex)
        B2, Bi2 = LL.exact_mat(Q)
        c2 = [mid_tm(Pc[i]) for i in range(2)]
        dc = [sub_tm(Pc[i], c2[i]) for i in range(2)]
        Rm = LL.conj_tm(cx, Bi2, P, B)                       # B'^{-1} P B as ONE matrix: no rotation of the error box
        Rr = matvec(cx, Rm, r)
        dcb = [cx.add(cx.scale(Bi2[i][0], dc[0]), cx.scale(Bi2[i][1], dc[1])) for i in range(2)]
        r = [cx.add(Rr[i], dcb[i]) for i in range(2)]
        c, B, Bf = c2, B2, Q
        if dbg is not None:
            dbg.append((pc.name, max(float(cx.sup(r[i]).upper()) for i in range(2)), max(float(cx.sup(c[i]).upper()) for i in range(2))))
    if structured:
        return c, B, r                                         # Y = c + B r, kept apart for det_structured
    Br = [cx.add(cx.scale(B[i][0], r[0]), cx.scale(B[i][1], r[1])) for i in range(2)]
    return [cx.add(c[i], Br[i]) for i in range(2)]


def det_structured(cx, Zl, cBr):
    """det[Zl, c + B r] = det[Zl, c] + r_1 det[Zl, b_1] + r_2 det[Zl, b_2]: the same number, but the error coordinate
    along the dominant column b_1 (nearly parallel to the light-cone solution at the match) is multiplied by the small
    det[Zl, b_1] instead of being spread over both components by the box"""
    c, B, r = cBr
    d = det2(cx, Zl, c)
    for j in range(2):
        bj = [cx.const(B[0][j]), cx.const(B[1][j])]
        d = cx.add(d, cx.mul(r[j], det2(cx, Zl, bj)))
    return d


class AxialB:
    """D_2 with T3 (lohner) and T4 (match_k: number of near-cone pieces the light-cone solution is carried through)"""

    def __init__(self, path, S_set=(1,), lohner=True, match_k=0, centre='regular'):
        self.path = path; self.S_set = list(S_set); self.lohner = lohner; self.match_k = match_k; self.centre = centre
        near = path.near_pieces
        self.chain_l = [(pc, False) for pc in near[:match_k]]
        self.chain_c = ([(pc, True) for pc in reversed(path.tube_pieces)] + [(path.sm_tube, True), (path.sm_out, False)]
                        + [(pc, False) for pc in path.out_pieces] + [(pc, True) for pc in reversed(near[match_k:])])
        if centre == 'singular':
            RS = dict(path.centre); RS['shift'] = 1
            self.RSc = singular_centre(RS, exponent=-4, start=[arb(1), arb(1)/(2*int(path.SY.su))])   # ker(M0 + 4), M0 = [[-3, -2 su], [-2/su, 0]]

    def _tr(self, cx, Y, chain, dbg):
        if not chain:
            return Y
        return lohner_transport(cx, Y, chain, dbg) if self.lohner else LL.transport(cx, Y, chain, True, dbg)

    def D_tm(self, kc, r, J, dbg=None, info=None):
        if not (np.isfinite(_cm(kc).real) and np.isfinite(_cm(kc).imag)):
            raise ArithmeticError('non-finite kappa')
        cx = Ctx(kc, r, J)
        Ylc, ilc = LL.regsing_tm(cx, self.path.light, None, self.path.light['t'], S_set=self.S_set)
        th = acb(self.path.centre['theta'])
        if self.centre == 'regular':
            Yc, ic = LL.regsing_tm(cx, self.path.centre, self.path.centre['y0'], self.path.centre['theta'])
            Yc = [cx.scale(th, c) for c in Yc]
        else:
            Yc, ic = regsing_parity(cx, self.RSc, self.path.centre['theta'])
            Yc = [cx.scale(th**-4, c) for c in Yc]
        Zl = self._tr(cx, Ylc, self.chain_l, dbg)
        if info is not None:
            info.update(light=ilc, centre=ic)
        if self.lohner and self.chain_c:
            return det_structured(cx, Zl, lohner_transport(cx, Yc, self.chain_c, dbg, structured=True)), cx
        Zc = self._tr(cx, Yc, self.chain_c, dbg)
        return det2(cx, Zl, Zc), cx


# =========================================================================== T5 singular centre class
def singular_centre(RS, exponent, start):
    """RegSing data for y = theta^exponent ytil from a regular-class RS whose M already has -I (exponent 1) or -2I
    (toy, exponent 2) subtracted: shift the diagonal by (old exponent - new exponent)."""
    Mn = [[[[x for x in row] for row in mat] for mat in plist] for plist in RS['M']]     # structural copy (no pickling)
    old = RS['shift']
    for i in range(2):
        Mn[0][0][i][i] = Mn[0][0][i][i] + (old - exponent)
    d = dict(RS); d['M'] = Mn; d['y0'] = start; d['exponent'] = exponent
    return d


def regsing_parity(cx, RS, t):
    """centre series for the singular branch: every odd-order coefficient of M (all kappa powers) must be the exact
    zero; then y_odd = 0 solves the odd equations, and the even orders use (n - M0)^{-1} (n even, never the odd resonance).
    Tail: lin_stages majorant with c_n = 1/(n - ||M0||) for n > K (K even and >= ||M0|| + 2)."""
    K = RS['K']; Mx = RS['M']
    for p in range(3):
        for n in range(1, K + 1, 2):
            if not all(Mx[p][n][i][j] == 0 for i in range(2) for j in range(2)):
                raise ArithmeticError('parity: odd-order coefficient not exactly zero (p=%d, n=%d)' % (p, n))
    for p in (1, 2):
        if not all(Mx[p][0][i][j] == 0 for i in range(2) for j in range(2)):
            raise ArithmeticError('kappa part of M0 not exactly zero')
    kap = cx.kappa(); k2 = cx.mul(kap, kap)

    def Mtm(n):
        return [[cx.add(cx.add(cx.const(acb(Mx[0][n][i][j])), cx.scale(acb(Mx[1][n][i][j]), kap)),
                        cx.scale(acb(Mx[2][n][i][j]), k2)) for j in range(2)] for i in range(2)]
    Ms = [Mtm(n) for n in range(K + 1)]
    y = [[cx.const(RS['y0'][0]), cx.const(RS['y0'][1])]]
    for n in range(1, K + 1):
        if n % 2 == 1:
            y.append([cx.const(0), cx.const(0)]); continue
        rr = [cx.const(0), cx.const(0)]
        for m in range(2, n + 1, 2):
            rr = [cx.add(rr[i], c) for i, c in enumerate(matvec(cx, Ms[m], y[n - m]))]
        M0c = [[acb(n)*(1 if i == j else 0) - acb(Mx[0][0][i][j]) for j in range(2)] for i in range(2)]
        inv = acb_mat(M0c).inv()
        y.append([cx.add(cx.scale(inv[i, 0], rr[0]), cx.scale(inv[i, 1], rr[1])) for i in range(2)])
    # check the start vector: (0 - M0) y0 = 0 exactly (kernel)
    M00 = [[acb(Mx[0][0][i][j]) for j in range(2)] for i in range(2)]
    res = [sum((M00[i][j]*acb(RS['y0'][j]) for j in range(2)), acb(0)) for i in range(2)]
    if not all(v.contains(0) for v in res):
        raise ArithmeticError('singular start not in ker M0')
    t = acb(t); at = abs(t).upper(); nu = RS['nu']
    nup = ((nu + at)/2).min(2*at)
    kab = arb(cx.kabs)
    Snorm = RS['sups'][0] + kab*RS['sups'][1] + kab**2*RS['sups'][2]
    rat = nup/nu
    mup = Snorm*rat/(1 - rat)
    R0n = rownorm([[acb(Mx[0][0][i][j]) for j in range(2)] for i in range(2)])
    n0 = K
    while True:
        c = 1/(arb(n0 + 1) - R0n) if arb(n0 + 1) > R0n else None
        if c is not None and c*mup <= 1:
            break
        n0 += 1
        if n0 > 100000:
            raise ArithmeticError('tail: n0 too large')
    if n0 > K:
        raise ArithmeticError('singular centre tail needs K >= %d (have %d)' % (n0, K))
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
    return out, dict(tail=tail.str(5), n0=n0)


class ToySing:
    """flat toy with the singular centre class r^{-l-1} (exponent -3), centre carried outward, v1 light-cone module"""

    def __init__(self, tp, S_set=(1, 2)):
        self.tp = tp; self.S_set = list(S_set)
        RS = dict(tp.centre); RS['shift'] = 2
        self.RS = singular_centre(RS, exponent=-3, start=[arb(1), arb(-3)])

    def D_tm(self, kc, r, J, dbg=None):
        cx = Ctx(kc, r, J)
        Ylc, _ = LL.regsing_tm(cx, self.tp.light, None, self.tp.light['t'], S_set=self.S_set)
        Yc, _ = regsing_parity(cx, self.RS, self.tp.centre['theta'])
        u = acb(self.tp.centre['theta'])
        Yc = [cx.scale(u**-3, c) for c in Yc]
        Zc = LL.transport(cx, Yc, self.tp.chain_c, True, dbg)
        return det2(cx, Ylc, Zc), cx
