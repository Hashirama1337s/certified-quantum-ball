"""v3.1 repair driver (NEW module; every pinned file is imported, none is edited).

Repairs for REFEREE-E (seats/REFEREE-E.md) on top of driver_v3 (BAR-LB-E3, seal 8aec628fd2):
 4a  aligned_set_v31: first Lohner axis built rigorously.  r = Bi.(arb(v) dV) + Bi.b with Bi the
     ball inverse of the (exact, float-valued) basis B; the Gram-Schmidt leak (Bi v)_{1..3} dV is
     enclosed instead of dropped (driver_v3.py:122 assumed B[:,0].|v| = v exactly).
 4b  sonic tails sym(eps) (N, om, V) and sym(epsX) (x) are put into the thin start, so every Lohner
     set / a priori box of the tau leg and the x leg encloses the TRUE trajectory, and the mean-value
     segment between truncated and true start lies inside the set (e_t, e_end now rigorous).
 4c  V0 tangent at the tau->x switch is the derivative of the SAME map whose value switch_set
     encloses (exact crossing C(y,x) = phi_{x0-x}(y)):  dC = Dphi_delta dy - f(C) dx with
     Dphi_delta in I + delta Df(B) (I + sym(e^{L dmax} - 1))   (Gronwall + integral mean value).
     The end transport keeps the first-order surrogate R^ for BOTH value (point run) and derivative
     (box run); E bounds R - R^ (unchanged logic, see PROPOSAL-V31.md).
 4d  e^{L dmax} factor on the tails at the switch (e0) and, by the same argument, e^{L dmx} on
     e_end at the end transport (L = inf-norm of Df over the self-consistent transport box).
 2   asserts: S0 (companion at the sonic point) contains 0; Dc(Z0) - lam_c contains 0; exact sympy
     identities for both (exact_chart_checks).
 5   pe2_v31: det sign certified (one sign) on every inward tau box and x box (bisected if needed);
     x monotone on the inward manifold piece (X' > 0 incl. Cauchy tail derivative) and the centre
     piece (1 + theta chi' > 0); one f_V direction per run of 'mono' steps, consistent with the
     neighbouring certified signs; sign boxes inflated by max(e_t, e0, e_end) (margin; rigour
     already from 4b).
 R7  (new, found by the sealed-parameter tail control) param.tail_bound dropped the defect head:
     inputs were truncated to degree 9 by python-flint's default ctx.cap = 10 (see below).
 6   content_null: Krawczyk at a radius-r box centred off x~ must report K not inside X.
The sonic side (everything that depends on V0 only) is cached per exact V0 ball, so the W0 / xc
content nulls reuse the target's sonic legs bit-for-bit.
"""
import sys, os, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import driver_v3 as V3          # installs the weighted-QR Lohner step (V3-1) at import
import driver_v2 as V2
import driver2 as D2
import driver as DR
import lohner as LO
from flint import arb, arb_mat
from core import Dual, sym, amax, sym_const
import param as PM
from param import val

Problem2 = V2.Problem2

# ------------------------------------------------------------------ R7 (found by the v3.1 tail control)
# param.tail_bound (param.py:91) builds  ins = [arb_series(Yc[i]) ...]  BEFORE Tape.series_eval raises
# ctx.cap to N2.  python-flint's default ctx.cap is 10, so every input polynomial was truncated to degree
# 9 (series precision 10): the defect head n = K+1 .. N2-1 (param.py:105-109) was summed over an EMPTY
# range and the low-order check (param.py:94-101) covered n <= 9 only.  eps/epsX were then the
# Cauchy tail beyond N2 alone (sonic, sealed K=16: claimed 3.7e-59; actual |Y16 - Y24| = 4.3e-35 at
# 0.95 nu, matching mpmath at 90 digits).  Repair: the inputs are exact polynomials, so give them a
# series precision that never truncates; products are still truncated at ctx.cap = N2 inside
# series_eval.  Installed by rebinding param.arb_series (same pattern as V3-1); pinned file unchanged.
from flint import arb_series as _flint_arb_series
SERIES_PREC_EXACT = 1 << 16

def arb_series_exact(coeffs, prec=None):
    return _flint_arb_series(coeffs, prec=SERIES_PREC_EXACT if prec is None else prec)

PM.arb_series = arb_series_exact


def mat_inf(J):
    m = arb(0)
    for row in J:
        m = m.max(sum((x.abs_upper() for x in row), arb(0)))
    return m


def growth(L, d):
    """upper bound of e^{L d} as an exact arb"""
    return arb((L*d).exp().upper())


def _units(n):
    return [[arb(1) if i == j else arb(0) for i in range(n)] for j in range(n)]


# ------------------------------------------------------------------ 4a + 4b
def aligned_set_v31(c_thin, d_box, rV, tails):
    """{ P_K(V0) + tail : |V0 - V0~| <= rV },  P_K(V0) in c_thin + d (V0 - V0~),  |tail_i| <= tails[i].
    Basis B = [v/|v|, completion] (float Gram-Schmidt: B is just an exact matrix); r enclosed with the
    ball inverse Bi:  r in Bi.v . dV + Bi.((y_thin - c) + (d - v) dV),  v = float mid(d) (exact)."""
    n = len(c_thin)
    y_thin = [c_thin[i] + sym(tails[i]) for i in range(n)]
    c = [arb(t.mid()) for t in y_thin]
    S = LO.LohnerSet(c)
    v = [float(t.mid()) for t in d_box]
    nv = sum(x*x for x in v)**0.5
    if nv == 0 or not (rV > 0):
        S.r = [y_thin[i] - c[i] for i in range(n)]
        return S, arb(0)
    B = V3._gs_basis([[x/nv for x in v]] + [[1.0 if i == j else 0.0 for i in range(n)] for j in range(n)], n)
    Bi = V3._inv(B)
    dV = sym(rV)
    va = [arb(x) for x in v]
    w = LO.matvec(Bi, va)                         # ~ (|v|, O(1e-16|v|), ...): enclosed, not assumed
    b = [(y_thin[i] - c[i]) + (d_box[i] - va[i])*dV for i in range(n)]
    rb = LO.matvec(Bi, b)
    S.B = B
    S.r = [w[i]*dV + rb[i] for i in range(n)]
    leak = amax(w[1:])*arb(rV)                    # what driver_v3.py:122 dropped
    return S, leak


# ------------------------------------------------------------------ 4c + 4d at the switch
def switch_tangent(Pb, Y4, t4, x0, dmax, extra):
    """rigorous enclosure of d/dV0 of the exact crossing C(y, x) = phi_{x0 - x}(y), for every point of
    the hull Y4 (4D, incl. x) with V0-tangent t4 (4D).  Returns tan, L, B, fB, half."""
    B, fB, half = DR.transport_bound(Pb, Y4[:3], dmax, extra)
    J, _ = PM.jac(Pb.F, B)
    L = mat_inf(J)
    g = growth(L, dmax)
    W = [[(arb(1) if i == j else arb(0)) + sym(g - 1) for j in range(3)] for i in range(3)]
    JW = LO.matmul(J, W)
    delta = x0 - Y4[3]
    Dphi = [[(arb(1) if i == j else arb(0)) + delta*JW[i][j] for j in range(3)] for i in range(3)]
    tan = [sum((Dphi[i][j]*t4[j] for j in range(3)), arb(0)) - fB[i]*t4[3] for i in range(3)]
    return tan, L, g, B, fB, half


# ------------------------------------------------------------------ item 2 (exact identities)
def exact_chart_checks(k):
    """sympy: det M vanishes identically on the sonic line N0 = (1 - V c)/(c - V) (=> S0 = 0), and
    the centre companion Dc(P=0, W, U*) == lam_c = -k(1+k) identically in W."""
    import sympy as sp
    import system as S0
    from pipeline import load_expr, FIELDS, P_, W_, U_, k_
    k = sp.Rational(k)
    Ms = S0.Msym
    det = (Ms[0, 0]*Ms[1, 1] - Ms[0, 1]*Ms[1, 0]).subs(S0.k, k)
    c = sp.sqrt(k)
    N0 = (1 - S0.V*c)/(c - S0.V)
    e1 = sp.simplify(sp.factor(sp.together(det.subs(S0.N, N0))))
    pc3 = load_expr(FIELDS['pol_c'][3])
    e2 = sp.simplify(sp.expand(pc3.subs(k_, k).subs({P_: 0, U_: -sp.Integer(2)/(3*(1 + k))}) + k*(1 + k)))
    return dict(sonic_S0_identically_zero=(e1 == 0), centre_Dc_equals_lam_c=(e2 == 0))


# ------------------------------------------------------------------ matcher
class MatcherV31(V3.MatcherV3):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._cache = {}

    def _key(self, V0, om_guess):
        return (str(V0.mid().man_exp()), str(V0.rad().man_exp()), str(om_guess), str(self.xbar.mid().man_exp()))

    def sonic_side(self, V0, om_guess, log=None):
        Pb = self.Pb; T0 = time.time()
        c_thin, d_box, st, sm = self.start_parts(V0, om_guess)
        S0 = val(sm['co']['S0'])
        assert S0.contains(0), 'item 2: companion S0 at the sonic point does not contain 0: %s' % S0
        tail4 = [sm['eps']]*3 + [sm['epsX']]
        S4, leak = aligned_set_v31(c_thin, d_box, V0.rad(), tail4)          # 4a + 4b
        vecs4 = [list(d_box)] + _units(4)
        S4.vh = [[arb(x.mid()) for x in v] for v in vecs4]
        Bi4 = V3._inv(S4.B)
        S4.t = [LO.matvec(Bi4, [v[i] - S4.vh[kk][i] for i in range(4)]) for kk, v in enumerate(vecs4)]
        width_start = max(t.rad() for t in S4.hull())
        tl = V2.tau_leg(self.T, S4, self.x_switch, K=self.Kt, log=log)
        V4 = S4.vecs(); Y4 = S4.hull()
        Phi_t = [[V4[1 + j][i] for j in range(4)] for i in range(4)]
        e_t = [sum((Phi_t[i][j].abs_upper()*tail4[j] for j in range(4)), arb(0)) for i in range(4)]
        x0 = self.x_switch; self.x0 = x0
        dmax = (x0 - Y4[3]).abs_upper() + e_t[3]
        tan, L0, g0, B0, fB, half_acc = switch_tangent(Pb, Y4, V4[0], x0, dmax, amax(e_t[:3]))   # 4c
        e0 = g0*amax(e_t[:3]) + e_t[3]*DR.fnorm(fB) + dmax**2*half_acc                              # 4d
        # diagnostic only: the v3 surrogate tangent (not used)
        ys = [Dual(Y4[i], [V4[0][i]]) for i in range(4)]
        fs = Pb.F.eval(ys[:3], conv=lambda c: PM._c(c, ys[0]))
        tan_sur = [(ys[i] + (x0 - ys[3])*fs[i]).d[0] for i in range(3)]
        S = V3.switch_set(Pb, S4, x0)
        vecs = [tan] + _units(3)
        S.vh = [[arb(x.mid()) for x in v] for v in vecs]
        Bi3 = V3._inv(S.B)
        S.t = [LO.matvec(Bi3, [v[i] - S.vh[kk][i] for i in range(3)]) for kk, v in enumerate(vecs)]
        width_switch = max(t.rad() for t in S.hull())
        info = D2.integrate_rec(Pb.F, S, x0, self.xbar, K=self.Kt, q=self.q, h0=self.h0, log=log)
        Y = S.hull(); V = S.vecs()
        Phi = [[V[1 + j][i] for j in range(3)] for i in range(3)]
        e_end = amax([sum((Phi[i][j].abs_upper() for j in range(3)), arb(0)) for i in range(3)])*e0
        return dict(Y=Y, V0tan=V[0], e_end=e_end, e0=e0, e_t=e_t, et_max=amax(e_t), L0=L0, g0=g0,
                    leak4a=leak, sm=sm, rec=info['rec'], tau_rec=tl['rec'], tau_steps=tl['steps'],
                    x_steps=info['steps'], width_start=width_start, width_switch=width_switch,
                    tan=tan, tan_sur=tan_sur, time=time.time() - T0)

    def residual(self, V0, W0, xc, om_guess, log=None, mvstart=True):
        Pb = self.Pb; T0 = time.time(); n = 3
        cvd = lambda like: (lambda c: PM._c(c, like))
        W0d = Dual(W0, [arb(1)])
        cm = Pb.centre_manifold(W0d, self.nu_c)
        Sc = val(cm['co']['S0'])
        assert (Sc - Pb.lam_c).contains(0), 'item 2: Dc(Z0) does not contain lam_c: %s vs %s' % (Sc, Pb.lam_c)
        Zc = [DR.emb(PM.horner(cm['co']['Y'][i], self.thm + 0*W0d), n, 1) for i in range(3)]
        chi = DR.emb(PM.horner(cm['co']['X'], self.thm + 0*W0d), n, 1)
        xcd = Dual(xc, [arb(0), arb(0), arb(1)])
        xm = self.thm.log() + chi + xcd
        if not hasattr(self, 'xbar'):
            self.xbar = arb(val(xm).mid())
            if os.environ.get('LB_SMOKE'):
                self.xbar = self.x_switch - arb(os.environ['LB_SMOKE'])
        key = self._key(V0, om_guess)
        sd = self._cache.get(key)
        cached = sd is not None
        if not cached:
            sd = self.sonic_side(V0, om_guess, log=log)
            self._cache[key] = sd
        Y, Vt, e_end = sd['Y'], sd['V0tan'], sd['e_end']
        y_end = [Dual(Y[i], [Vt[i], arb(0), arb(0)]) for i in range(3)]
        dm = xm - self.xbar
        if os.environ.get('LB_SMOKE'):
            dm = dm - arb(val(dm).mid())
        fe = Pb.F.eval(y_end, conv=cvd(y_end[0]))
        dmx = val(dm).abs_upper() + cm['epsX']
        Bt, fBt, half_acc2 = DR.transport_bound(Pb, Y, dmx, e_end)
        Jt, _ = PM.jac(Pb.F, Bt)
        Lt = mat_inf(Jt); gt = growth(Lt, dmx)
        ym = [y_end[i] + dm*fe[i] for i in range(3)]
        N, om, Vv = ym
        Zs = [1/N, om*N*N, Vv*N]
        R = [Zs[i] - Zc[i] for i in range(3)]
        Nb, ob, Vb = Bt
        JZ = amax([(1/(Nb*Nb)).abs_upper(), (2*ob*Nb).abs_upper() + (Nb*Nb).abs_upper(), Vb.abs_upper() + Nb.abs_upper()])
        eend = e_end*gt + dmx**2*half_acc2 + DR.fnorm(fBt)*cm['epsX']                               # 4d (end)
        Erad = JZ*eend + cm['eps']
        diag = dict(sd)
        diag.update(steps=sd['tau_steps'] + sd['x_steps'], time=time.time() - T0, cached=cached,
                    boxes=[(x, float(h.mid()), Bx) for (x, h, Bx, s_) in sd['rec']],
                    Erad=Erad, cm=cm, Bt=Bt, gt=gt, width_end=max(v.rad() for v in Y),
                    tan_start=None, tan_end=max(abs(float(t.mid())) for t in Vt))
        return R, Erad, diag


# ------------------------------------------------------------------ item 5 (PE2)
def _sgn(a):
    try:
        if a > 0: return '+'
        if a < 0: return '-'
    except Exception:
        pass
    return None


def cls_tau(T, Bx, margin):
    """(V class, det sign) on a 4D tau-step box.  'mono' carries dV/ds, s = -x (traversal)."""
    Bv = [val(t) + sym(margin) for t in Bx]
    try:
        g = T.eval(Bv, conv=sym_const)
    except Exception:
        return 'amb', None
    ds = _sgn(g[3])
    vs = _sgn(Bv[2])
    if vs:
        return vs, ds
    g2 = _sgn(g[2])
    if ds and g2:
        dVdx = '+' if g2 == ds else '-'
        return ('mono', '-' if dVdx == '+' else '+'), ds
    return 'amb', ds


def cls_x(Pb, Bx, margin):
    Bv = [val(t) + sym(margin) for t in Bx[:3]]
    try:
        ds = _sgn(Pb.Ss.eval(Bv, conv=sym_const)[0])
    except Exception:
        ds = None
    vs = _sgn(Bv[2])
    if vs:
        return vs, ds
    try:
        f = _sgn(Pb.F.eval(Bv, conv=sym_const)[2])
    except Exception:
        f = None
    if f and ds:
        return ('mono', '-' if f == '+' else '+'), ds
    return 'amb', ds


def bisect_rec(F, rec, clsfn, K, q='0.25', depth=12):
    qa = arb(q)
    def bad(c, ds): return c == 'amb' or ds is None
    def sub(start, h, d):
        S = LO.LohnerSet(start); res = []
        for half in range(2):
            st = S.hull()
            Bx, _ = LO.lohner_step(F, S, h/2, K, qa)
            c, ds = clsfn(Bx)
            res += sub(st, h/2, d + 1) if (bad(c, ds) and d < depth) else [(c, ds)]
        return res
    out = []; nb = 0
    for (x, ha, Bx, start) in rec:
        c, ds = clsfn(Bx)
        if bad(c, ds):
            nb += 1
            try:
                out += sub(start, ha, 1)
            except ArithmeticError:
                out.append(('amb', None))
        else:
            out.append((c, ds))
    return out, nb


def mono_runs(classes):
    """every maximal run of ('mono', d) needs one d, certified neighbours, and d matching a sign change."""
    ok = True; runs = 0; i = 0
    while i < len(classes):
        if isinstance(classes[i], tuple):
            j = i; dirs = set()
            while j < len(classes) and isinstance(classes[j], tuple):
                dirs.add(classes[j][1]); j += 1
            runs += 1
            left = classes[i - 1] if i > 0 else None
            right = classes[j] if j < len(classes) else None
            if len(dirs) != 1 or left not in ('+', '-') or right not in ('+', '-'):
                ok = False
            elif left != right and dirs != {'-' if left == '+' else '+'}:
                ok = False
            i = j
        else:
            i += 1
    return ok, runs


def pe2_v31(Pb, M, d, depth=12):
    margin = d['et_max'].max(d['e0']).max(d['e_end'])
    sm, cm = d['sm'], d['cm']
    # inward sonic-manifold piece theta in [th1, 0] (th1 < 0)
    I1 = M.th1/2 + sym(M.th1/2)
    Vs = PM.horner([val(a) for a in sm['co']['Y'][2]], I1) + sym(sm['eps'])
    Xd = PM.horner([val(a) for a in sm['co']['X']], Dual(I1, [arb(1)]))
    Xp = Xd.d[0] + sym(sm['epsX']/(M.nu_s - M.th1.abs_upper()))          # Cauchy on |theta| <= nu_s
    lam = val(sm['lam'])
    mono_s = bool(Xp > 0)
    det_s = None
    if mono_s and _sgn(lam):
        det_s = '-' if _sgn(lam) == '+' else '+'                              # det = lam theta X', theta < 0
    # centre piece theta in [0, thm]
    Im = M.thm/2 + sym(M.thm/2)
    Uc = PM.horner([val(a) for a in cm['co']['Y'][2]], Im) + sym(cm['eps'])
    chi = PM.horner([val(a) for a in cm['co']['X']], Dual(Im, [arb(1)]))
    dxc = 1 + Im*chi.d[0] + sym(M.thm*cm['epsX']/(M.nu_c - M.thm))
    mono_c = bool(dxc > 0)
    tau_seq, nb_t = bisect_rec(M.T, d['tau_rec'], lambda Bx: cls_tau(M.T, Bx, margin), M.Kt, depth=depth)
    x_seq, nb_x = bisect_rec(Pb.F, d['rec'], lambda Bx: cls_x(Pb, Bx, margin), M.Kt, depth=depth)
    jt = cls_x(Pb, d['Bt'], arb(0))
    vs_c = _sgn(Vs) or 'amb'; uc_c = _sgn(Uc) or 'amb'
    classes = [vs_c] + [c for c, _ in tau_seq] + [c for c, _ in x_seq] + [jt[0]] + [uc_c]
    dets = [ds for _, ds in tau_seq] + [ds for _, ds in x_seq] + [jt[1]]
    det_set = set(dets)
    det_ok = (None not in det_set) and len(det_set) == 1 and det_s in det_set
    mono_ok, nruns = mono_runs(classes)
    namb = classes.count('amb')
    cert = [c for c in classes if c in ('+', '-')]
    changes = sum(1 for a, b in zip(cert, cert[1:]) if a != b)
    nmono = sum(1 for c in classes if isinstance(c, tuple))
    holds = (namb == 0 and det_ok and mono_ok and mono_s and mono_c and changes == 1 and uc_c == '-')
    return dict(holds=holds, sign_changes=changes, amb=namb, mono_steps=nmono, mono_runs=nruns, mono_ok=mono_ok,
                det_sign=sorted(det_set, key=str), det_ok=det_ok, det_manifold=det_s,
                x_monotone_manifold=mono_s, x_monotone_centre=mono_c,
                Xprime_manifold=Xp, dxdtheta_centre_factor=dxc, V_manifold=vs_c, U_centre=uc_c,
                tau_boxes=len(tau_seq), x_boxes=len(x_seq), tau_bisected=nb_t, x_bisected=nb_x, margin=margin)


# ------------------------------------------------------------------ item 6 (content null)
def content_null(M, Xt, axis, om_guess, off='1e-13', rad='1e-15', log=print):
    """box of radius rad centred off along `axis` from x~: Krawczyk must report K not inside X.
    An exception is NOT a pass (returned as void)."""
    Xs = list(Xt); Xs[axis] = Xt[axis] + arb(off)
    try:
        Rp, Ep, dp = M.residual(Xs[0], Xs[1], Xs[2], om_guess, log=None)
        Rv = [val(r) for r in Rp]
        C = arb_mat([[arb(Rp[i].d[j].mid()) for j in range(3)] for i in range(3)]).inv()
        r = arb(rad)
        Xb = [Xs[i] + sym(r) for i in range(3)]
        Rb, Eb, db = M.residual(Xb[0], Xb[1], Xb[2], om_guess, log=None)
        ok, Kx, Xb = D2.kr_test(Xs, Rv, C, Rb, Eb, r)
    except Exception as e:                       # noqa: BLE001  (void, never a pass)
        return dict(axis=axis, void=True, passed=False, error=repr(e)[:200])
    inside = [bool(Xb[i].contains_interior(Kx[i])) for i in range(3)]
    miss = max(float(((Kx[i].mid() - Xb[i].mid())/r).abs_upper()) for i in range(3))
    excess = max(float(Kx[i].rad()/Xb[i].rad()) for i in range(3))
    return dict(axis=axis, void=False, contracts=ok, passed=(not ok), inside=inside, miss_over_r=miss,
                excess=excess, Ebox=Eb, cached_point=dp['cached'], cached_box=db['cached'],
                K=[k.str(20) for k in Kx], X=[x.str(20) for x in Xb])
