"""v3.2 driver (repair R8; NEW module; every pinned file and driver_v31 are imported, none is edited).

R8  (lead diagnosis runs4/diag_e4w_1e-15.log, verified in runs5/smoke/tails_r8.log):  param.tail_bound
    sums the defect head  sum_{K<n<N2} |[G(P_K)]_n| nu^n  with abs_upper of coefficients evaluated over
    the V0 BOX.  Their width is first order in the box radius r (dependency of ball arithmetic through the
    sonic point, eigenpair and coefficient recursion), so a quantity that is a smooth function of V0 was
    booked as tail ERROR: eps ~ r (k=1/2: 4.3e-22 at r=1e-15, 6.4e-13 at r=1e-9) instead of ~ point eps.
    Repair (mean-value form, the same form the start set already uses for P_K itself):
        |c_n(V0)| <= |c_n(V0~)| + r |c_n'(box)|      for every V0 in [V0~ - r, V0~ + r],
    c_n(V0~) from the point coefficients (V0~ = mid, exact), c_n'(box) by Dual (forward-mode) series
    arithmetic over the box (rigorous enclosure of dc_n/dV0 at every point of the box; its own width is
    O(r), so the r-term is O(r) |true derivative| + O(r^2)).  Per coefficient the bound used is
    min(box ball, mean-value bound): both are rigorous, so eps never exceeds the v3.1 value.
    Everything else in tail_bound (low-order check on the box series, Cauchy tail beyond N2 over the box,
    eigen-coordinate operator bound and Lipschitz constant over the box) is the pinned code, unchanged:
    the head coefficients reach it through a proxy tape (TMHeadTape) whose series_eval replaces ONLY the
    coefficients n > K by the ball [-h_n, h_n].
    The same proxy applies to the companion S (epsX).  The centre manifold over the W0 box is left as is
    (measured: its eps stays < 1e-50 up to r = 1e-9; see PROPOSAL-V32.md).
"""
import sys, os, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import driver_v31 as V31          # installs R7 (exact series inputs) and V3-1 (weighted QR) at import
import driver_v3 as V3
import driver_v2 as V2
import param as PM
from param import val
from core import Dual, sym, sym_const
from flint import arb, arb_series, ctx

assert PM.arb_series is V31.arb_series_exact, 'R7 not installed'
SER = V31.SERIES_PREC_EXACT


def _ser(coeffs):
    return arb_series(list(coeffs), prec=SER)      # never truncates at construction (lab doctrine #115)


def series_eval_dual(tape, ins, cap):
    """Tape.series_eval over Dual(arb_series, [arb_series]) inputs (forward-mode derivative of every
    output series).  ctx.cap is raised to cap BEFORE any series is built (doctrine #115)."""
    old = ctx.cap; ctx.cap = cap
    try:
        vals = [None]*len(tape.nodes)
        for i, (op, a, b) in enumerate(tape.nodes):
            if op == 'var': vals[i] = ins[a]
            elif op == 'const': vals[i] = Dual(_ser([sym_const(a)]), [_ser([arb(0)])])
            elif op == 'add': vals[i] = vals[a] + vals[b]
            elif op == 'mul': vals[i] = vals[a]*vals[b]
            elif op == 'div': vals[i] = vals[a]/vals[b]
        return [vals[j] for j in tape.outs]
    finally:
        ctx.cap = old


def _coeffs(s, N):
    cf = list(s.coeffs())
    return cf + [arb(0)]*(N - len(cf))


class TMHeadTape:
    """proxy for a pinned Tape inside param.tail_bound: eval / jac / outs / nodes delegate unchanged;
    series_eval returns the pinned box series with the coefficients n > K replaced by [-h_n, h_n],
    h_n = min(|box c_n|, |point c_n| + r |c_n'(box)|)."""
    def __init__(self, tape, K, r, pt_ins, dual_ins):
        self._t = tape; self.K = K; self.r = arb(r)
        self.pt_ins = pt_ins; self.dual_ins = dual_ins
        self.report = []

    def __getattr__(self, name):
        return getattr(self._t, name)

    def series_eval(self, ins, cap):
        box = self._t.series_eval(ins, cap)
        pt = self._t.series_eval(self.pt_ins, cap)
        du = series_eval_dual(self._t, self.dual_ins, cap)
        out = []
        for o in range(len(box)):
            cb = _coeffs(box[o], cap); cp = _coeffs(pt[o], cap); cd = _coeffs(du[o].d[0], cap)
            new = list(cb[:self.K + 1])
            nb = arb(0); nt = arb(0)
            for n in range(self.K + 1, cap):
                hb = cb[n].abs_upper()
                ht = cp[n].abs_upper() + self.r*cd[n].abs_upper()
                h = hb.min(ht)
                nb = nb.max(hb); nt = nt.max(ht)
                new.append(sym(h))
            self.report.append(dict(out=o, max_box=nb, max_tm=nt))
            old = ctx.cap; ctx.cap = cap
            try:
                out.append(_ser(new))
            finally:
                ctx.cap = old
        return out


def _dual_series(co_list):
    """[Dual(c_n, [dc_n])]_n -> Dual(series of values, [series of derivatives])"""
    v = [val(a) for a in co_list]
    d = [(a.d[0] if isinstance(a, Dual) else arb(0)) for a in co_list]
    return Dual(_ser(v), [_ser(d)])


def sonic_manifold_v32(Pb, V0, om_guess, nu, K, N2fac=12, tails=True):
    """sonic manifold over the V0 box (Dual in V0), with the R8 mean-value head.
    Returns the box dict of driver_v2.sonic_manifold_v2 (co, lam, v, ...) plus eps, epsX, tinfo, r8."""
    V0v = val(V0)
    V0d = V0 if isinstance(V0, Dual) else Dual(V0, [arb(1)])
    smb = V2.sonic_manifold_v2(Pb, V0d, om_guess, nu, K, tails=False)
    if not tails:
        return smb
    r = V0v.rad()
    Yc = [[val(a) for a in smb['co']['Y'][i]] for i in range(3)]
    Xc = [val(a) for a in smb['co']['X']]
    if r == 0:
        eps, epsX, info = PM.tail_bound(Pb.Gs, Pb.Ss, Yc, Xc, val(smb['lam']), smb['co']['DG0'], K, arb(nu), N2=N2fac*K)
        smb.update(eps=eps, epsX=epsX, tinfo=info, r8=None)
        return smb
    Vc = arb(V0v.mid())
    smc = V2.sonic_manifold_v2(Pb, Dual(Vc, [arb(1)]), om_guess, nu, K, tails=False)
    pt_ins = [_ser([val(a) for a in smc['co']['Y'][i]]) for i in range(3)]
    du_ins = [_dual_series(smb['co']['Y'][i]) for i in range(3)]
    Gp = TMHeadTape(Pb.Gs, K, r, pt_ins, du_ins)
    Sp = TMHeadTape(Pb.Ss, K, r, pt_ins, du_ins)
    eps, epsX, info = PM.tail_bound(Gp, Sp, Yc, Xc, val(smb['lam']), smb['co']['DG0'], K, arb(nu), N2=N2fac*K)
    smb.update(eps=eps, epsX=epsX, tinfo=info, r8=dict(G=Gp.report, S=Sp.report), smc=smc)
    return smb


# ------------------------------------------------------------------ matcher: R8 in the start set
class MatcherV32(V31.MatcherV31):
    """MatcherV31 with the sonic manifold of sonic_manifold_v32 (R8 head).  Everything downstream
    (4a-4d, item 2 asserts, the cache, residual, PE2 data) is inherited unchanged."""
    def start_parts(self, V0, om_guess):
        Pb = self.Pb
        V0d = Dual(V0, [arb(1)])
        sm = sonic_manifold_v32(Pb, V0d, om_guess, self.nu_s, self.Ks2)
        stb = [PM.horner(sm['co']['Y'][i], self.th1 + 0*V0d) for i in range(3)] + [PM.horner(sm['co']['X'], self.th1 + 0*V0d)]
        if V0.rad() > 0:
            smc = sm['smc']                                     # point manifold at V0~ (tails not needed)
            stc = [PM.horner(smc['co']['Y'][i], self.th1 + 0*V0d) for i in range(3)] + [PM.horner(smc['co']['X'], self.th1 + 0*V0d)]
        else:
            stc = stb
        c_thin = [val(t) for t in stc]
        d_box = [t.d[0] for t in stb]
        st = [Dual(c_thin[i] + d_box[i]*sym(V0.rad()), [d_box[i]]) for i in range(4)]
        return c_thin, d_box, st, sm


def kr_split(C, Rb, Eb, rad):
    """decomposition of rad(K)/rad(X) (diagnostic): error part sum_j |C_ij| E / r and
    derivative part sum_j |I - C J(X)|_ij, per row."""
    out = []
    for i in range(3):
        e = sum((C[i, j].abs_upper() for j in range(3)), arb(0))*Eb/arb(rad)
        jp = arb(0)
        for j in range(3):
            cij = (1 if i == j else 0) - sum((C[i, l]*Rb[l].d[j] for l in range(3)), arb(0))
            jp += cij.abs_upper()
        out.append((float(e.mid()), float(jp.mid())))
    return out
