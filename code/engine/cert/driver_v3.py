"""Attempt-3 driver (new module).  Root cause of the v2 box-run blow-up (diag_lohner*.py):
lohner.qr_float orders the QR columns by ||A_j|| only, ignoring the widths of r.  Lohner's
method needs the columns of A diag(|r|) ordered by norm; otherwise the first basis vector is not
aligned with the long axis of the set and wrapping accumulates.  Measured on x in [-1.5,-3.4],
k=1/3, box 1e-12:  growth 1.3e5 (v2) -> 7.8e2 (weighted QR) / 6.0e2 (parallelepiped);
true sigma_max of the flow = 564.

V3-1: every Lohner step (tau leg, x leg, PE2/outward bisection) uses lohner2.lohner_step2(mode='qrw').
      Installed by rebinding lohner.lohner_step at import time (the pinned files are not edited;
      all v2 callers resolve LO.lohner_step at call time).
Everything else is driver_v2 unchanged.
"""
import sys, os
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import lohner as LO
import lohner2 as L2

MODE = os.environ.get('LB_LOHNER_MODE', 'qrw')
_orig = LO.lohner_step
def _step(F, S, h, K, q):
    return L2.lohner_step2(F, S, h, K, q, mode=MODE)
LO.lohner_step = _step

import driver_v2 as V2          # noqa: E402  (uses LO.lohner_step via module attribute at call time)
import driver2 as D2            # noqa: E402
Problem2, MatcherV2 = V2.Problem2, V2.MatcherV2


# ---------------------------------------------------------------------------------------------
# V3-2: aligned start set and structure-preserving tau->x switch
#   start:  y(V0) in c + v (V0 - V0~) + b,  v = mid dy/dV0, b = small box (ball widths, derivative
#           spread); Lohner basis B = [v/|v|, completion], r = (|v| [-rV, rV], 0, 0, 0) + B^{-1} b.
#   switch: T(y, x) = y + (x0 - x) f(y);  T(c + B r) in T(c) + DT(hull) B r + R2,
#           DT = [I + (x0 - x) Df(y), -f(y)] over the hull; new 3D basis from weighted QR of mid(DT B).
# ---------------------------------------------------------------------------------------------
import math, time
from flint import arb, arb_mat
from core import Dual, sym, amax, sym_const
import param as PM
from param import val
import driver as DR

def _gs_basis(cols_float, n):
    Q = []
    for v in cols_float:
        for qv in Q:
            d = sum(a*b for a, b in zip(v, qv)); v = [a - d*b for a, b in zip(v, qv)]
        nv = math.sqrt(sum(x*x for x in v))
        if nv > 1e-8:
            Q.append([x/nv for x in v])
        if len(Q) == n:
            break
    return [[arb(Q[j][i]) for j in range(n)] for i in range(n)]

def _inv(B):
    n = len(B)
    Bi = arb_mat(B).inv()
    return [[Bi[i, j] for j in range(n)] for i in range(n)]

def aligned_set(st, rV):
    """LohnerSet from Dual start values st (derivative in V0 over the box), first axis along dy/dV0."""
    n = len(st)
    c = [arb(val(s).mid()) for s in st]
    v = [float(s.d[0].mid()) for s in st]
    nv = math.sqrt(sum(x*x for x in v))
    if nv == 0 or not (rV > 0):
        return LO.LohnerSet([val(s) for s in st])
    B = _gs_basis([[x/nv for x in v]] + [[1.0 if i == j else 0.0 for i in range(n)] for j in range(n)], n)
    Bi = _inv(B)
    dV = sym(rV)
    b = [val(st[i]) - c[i] - arb(v[i])*dV for i in range(n)]   # st.v encloses stc + d (V0 - V0~)
    rb = LO.matvec(Bi, b)
    S = LO.LohnerSet(c)
    S.B = B
    S.r = [rb[0] + arb(nv)*dV] + rb[1:]
    return S

def switch_set(Pb, S4, x0):
    """structure-preserving projection of the 4D tau-leg set onto x = x0 (3D Lohner set)."""
    Y4 = S4.hull()
    yb = Y4[:3]; delta = x0 - Y4[3]
    fB = Pb.F.eval(yb, conv=sym_const)
    J, _ = PM.jac(Pb.F, yb)
    DT = [[(1 if i == j else 0) + delta*J[i][j] for j in range(3)] + [-fB[i]] for i in range(3)]
    M = [[sum((DT[i][l]*S4.B[l][j] for l in range(4)), arb(0)) for j in range(4)] for i in range(3)]
    c4 = S4.c
    fc = Pb.F.eval(c4[:3], conv=sym_const)
    Tc = [c4[i] + (x0 - c4[3])*fc[i] for i in range(3)]
    dmax = delta.abs_upper()
    _, _, half = DR.transport_bound(Pb, yb, dmax)
    R2 = dmax**2*half
    w = [float(t.abs_upper()) + 1e-300 for t in S4.r]
    cols = [[float(M[i][j].mid()) for i in range(3)] for j in range(4)]
    order = sorted(range(4), key=lambda j: -math.sqrt(sum(x*x for x in cols[j]))*w[j])
    Q = _gs_basis([cols[j] for j in order] + [[1.0 if i == j else 0.0 for i in range(3)] for j in range(3)], 3)
    Qi = _inv(Q)
    c3 = [arb(t.mid()) for t in Tc]
    Mr = [sum((M[i][j]*S4.r[j] for j in range(4)), arb(0)) for i in range(3)]
    tot = [Mr[i] + (Tc[i] - c3[i]) + sym(R2) for i in range(3)]
    S3 = LO.LohnerSet(c3)
    S3.B = Q
    S3.r = LO.matvec(Qi, tot)
    return S3

def aligned_set_parts(c_thin, d_box, rV):
    """set { c + d (V0 - V0~) : |V0 - V0~| <= rV } enclosed WITHOUT the dependency loss of v1-v3a:
    y in c + v dV + (c_thin - c) + (d - v) dV,  v = mid(d); basis first axis along v."""
    n = len(c_thin)
    c = [arb(t.mid()) for t in c_thin]
    v = [float(t.mid()) for t in d_box]
    nv = math.sqrt(sum(x*x for x in v))
    S = LO.LohnerSet(c)
    if nv == 0 or not (rV > 0):
        S.r = [c_thin[i] - c[i] for i in range(n)]
        return S
    B = _gs_basis([[x/nv for x in v]] + [[1.0 if i == j else 0.0 for i in range(n)] for j in range(n)], n)
    Bi = _inv(B)
    dV = sym(rV)
    b = [(c_thin[i] - c[i]) + (d_box[i] - arb(v[i]))*dV for i in range(n)]
    rb = LO.matvec(Bi, b)
    S.B = B
    S.r = [rb[0] + arb(nv)*dV] + rb[1:]
    return S

class MatcherV3(V2.MatcherV2):
    def start_parts(self, V0, om_guess):
        """thin start at the box centre + derivative enclosure over the box (4D incl. x) + tails over box."""
        Pb = self.Pb
        V0d = Dual(V0, [arb(1)])
        sm = V2.sonic_manifold_v2(Pb, V0d, om_guess, self.nu_s, self.Ks2)
        stb = [PM.horner(sm['co']['Y'][i], self.th1 + 0*V0d) for i in range(3)] + [PM.horner(sm['co']['X'], self.th1 + 0*V0d)]
        Vc = arb(V0.mid())
        if V0.rad() > 0:
            smc = V2.sonic_manifold_v2(Pb, Dual(Vc, [arb(1)]), om_guess, self.nu_s, self.Ks2, tails=False)
            stc = [PM.horner(smc['co']['Y'][i], self.th1 + 0*V0d) for i in range(3)] + [PM.horner(smc['co']['X'], self.th1 + 0*V0d)]
        else:
            stc = stb
        c_thin = [val(t) for t in stc]
        d_box = [t.d[0] for t in stb]
        st = [Dual(c_thin[i] + d_box[i]*sym(V0.rad()), [d_box[i]]) for i in range(4)]   # for bookkeeping only
        return c_thin, d_box, st, sm

    def residual(self, V0, W0, xc, om_guess, log=None, mvstart=True):
        Pb = self.Pb; T0 = time.time(); n = 3
        cvd = lambda like: (lambda c: PM._c(c, like))
        c_thin, d_box, st, sm = self.start_parts(V0, om_guess)
        tail4 = [sm['eps']]*3 + [sm['epsX']]
        S4 = aligned_set_parts(c_thin, d_box, V0.rad())
        vecs4 = [[s.d[0] for s in st]] + [[arb(1) if i == j else arb(0) for i in range(4)] for j in range(4)]
        S4.vh = [[arb(x.mid()) for x in v] for v in vecs4]
        Bi4 = _inv(S4.B)
        S4.t = [LO.matvec(Bi4, [v[i] - S4.vh[kk][i] for i in range(4)]) for kk, v in enumerate(vecs4)]
        tl = V2.tau_leg(self.T, S4, self.x_switch, K=self.Kt, log=log)
        V4 = S4.vecs()
        Phi_t = [[V4[1 + j][i] for j in range(4)] for i in range(4)]
        e_t = [sum((Phi_t[i][j].abs_upper()*tail4[j] for j in range(4)), arb(0)) for i in range(4)]
        self.x0 = self.x_switch
        Y4 = S4.hull()
        y4 = [Dual(Y4[i], [V4[0][i]]) for i in range(4)]
        ys = y4[:3]; xs = y4[3]
        delta = self.x0 - xs
        fs = Pb.F.eval(ys, conv=cvd(ys[0]))
        y0 = [ys[i] + delta*fs[i] for i in range(3)]
        dmax = val(delta).abs_upper() + e_t[3]
        _, fB, half_acc = DR.transport_bound(Pb, [val(v) for v in ys], dmax, amax(e_t[:3]))
        e0 = amax(e_t[:3]) + e_t[3]*DR.fnorm(fB) + dmax**2*half_acc
        S = switch_set(Pb, S4, self.x0)
        vecs = [[v.d[0] for v in y0]] + [[arb(1) if i == j else arb(0) for i in range(3)] for j in range(3)]
        S.vh = [[arb(x.mid()) for x in v] for v in vecs]
        Bi3 = _inv(S.B)
        S.t = [LO.matvec(Bi3, [v[i] - S.vh[kk][i] for i in range(3)]) for kk, v in enumerate(vecs)]
        width_switch = max(t.rad() for t in S.hull())
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
                    width_start=max(val(v).rad() for v in st), width_switch=width_switch, Y=Y,
                    tan_start=max(abs(float(s.d[0].mid())) for s in st),
                    tan_switch=max(abs(float(t.mid())) for t in V4[0]),
                    tan_end=max(abs(float(t.mid())) for t in V[0]))
        return R, Erad, diag
