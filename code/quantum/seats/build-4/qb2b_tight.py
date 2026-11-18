"""build-4: in-process tightening of the engine's homogeneous stage (engine files untouched).
T1  taylor.step_cauchy2 -> same (y1, Bx) a priori box, but the Taylor coefficients Yc2 are replaced by the
    mean-value enclosure  y_n(mid y0) + sum_j [d y_n/d y0_j](ball) [-r_j, r_j]  (point jets + Dual jets): no ball
    recursion wrapping.  Both are enclosures of the true coefficients for every start point in the ball.
T2  lin_stages.LinFields.series_out (the step / centre tapes on series inputs) -> mean-value form
    out_n in out_n(mid inputs) + [-R_n, R_n],  R_n = sum_k sum_{j<=n} |(d_k out)_{n-j}|(balls) r_{k,j},
    d_k out evaluated on the ball inputs as Dual series (seed at order 0): removes the dependency amplification.
"""
import sys
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
from flint import arb
from core import Dual, sym, sym_const
import taylor as TY
import lin_stages as LS
from param import _c as _cc

_orig_sc2 = TY.step_cauchy2
_orig_so = LS.LinFields.series_out


def step_cauchy2_tight(tape, y0, h, K, q=arb('0.25'), conv=None, maxit=60):
    y1, Bx, Yc2 = _orig_sc2(tape, y0, h, K, q=q, conv=conv, maxit=maxit)
    if not isinstance(y0[0], arb) or all(v.rad().is_zero() for v in y0):
        return y1, Bx, Yc2
    d = len(y0); n = len(Yc2[0]) - 1
    c = [arb(v.mid()) for v in y0]
    Yc = TY.jets(tape, c, n, lambda q_: _cc(q_, arb(0)))
    yd = [Dual(y0[i], [arb(1) if j == i else arb(0) for j in range(d)]) for i in range(d)]
    Yd = TY.jets(tape, yd, n, lambda q_: _cc(q_, yd[0]))
    rads = [v.rad() for v in y0]
    out = []
    for i in range(d):
        row = [y0[i]]
        for m in range(1, n + 1):
            v = Yc[i][m]
            for j in range(d):
                v = v + Yd[i][m].d[j]*sym(rads[j])
            assert v.overlaps(Yc2[i][m]), 'tight jet does not overlap the ball recursion'
            row.append(v)
        out.append(row)
    return y1, Bx, out


def series_out_mv(self, tape, ins, K):
    rad_any = any((not x.rad().is_zero()) for s in ins for x in s if isinstance(x, arb))
    if not rad_any:
        return _orig_so(self, tape, ins, K)
    z = arb(0); nk = len(ins); m = K + 1
    mid = [[arb(x.mid()) for x in s] for s in ins]
    rad = [[x.rad() for x in s] for s in ins]
    om = _orig_so(self, tape, mid, K)
    zd = Dual(z, [z]*nk)
    tape.jet_init(m, zd, conv=lambda c: Dual(sym_const(c), [z]*nk))
    for n in range(m):
        for i in range(nk):
            if n < len(ins[i]):
                tape.C[i][n] = Dual(ins[i][n], [arb(1) if (j == i and n == 0) else z for j in range(nk)])
            else:
                tape.C[i][n] = zd
        tape.jet_order(n)
    dd = [[tape.C[o][n] for n in range(m)] for o in tape.outs]
    out = []
    for o in range(len(dd)):
        row = []
        for n in range(m):
            R = arb(0)
            for k in range(nk):
                for j in range(min(n + 1, len(ins[k]))):
                    if not rad[k][j].is_zero():
                        R += dd[o][n - j].d[k].abs_upper()*rad[k][j]
            row.append(om[o][n] + sym(R))
        out.append(row)
    return out


TY.step_cauchy2 = step_cauchy2_tight
LS.LinFields.series_out = series_out_mv
