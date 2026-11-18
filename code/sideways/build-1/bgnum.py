"""Numeric background jets at a point (for identity checks only, never for certified numbers).
Any point (N, om, V) of the CSS ODE (A from the algebraic constraint) is a local exact solution; its x-jets
are produced by the engine's own background field (engine/cert/fields.pkl 'f', 'Acons', via core.Tape)."""
import sys
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import pickle, math
import sympy as sp
from flint import arb, arb_series, ctx
ctx.cap = 32                     # never the default 10 (E3 trap)
from core import Tape, sym_const
import lin_core as LC

FIELDS = pickle.load(open('<LAB>/engine/cert/fields.pkl', 'rb'))
def load_expr(s):
    return eval(s, {k_: getattr(sp, k_) for k_ in dir(sp)})
Ns, oms, Vs, ks = sp.symbols('N om V k')


def jets_at(Nv, omv, Vv, kv=sp.Rational(1, 3), order=5, prec=330):
    """returns dict jet-symbol -> sympy Float (dps ~ 90) for N, A, om, V up to `order` derivatives,
    plus sA, gV."""
    old = ctx.prec
    try:
        ctx.prec = prec
        F = Tape([load_expr(s) for s in FIELDS['f']], [Ns, oms, Vs], {ks: kv})
        Ac = Tape([load_expr(FIELDS['Acons'])], [Ns, oms, Vs], {ks: kv})
        y = [arb_series([arb(Nv)]), arb_series([arb(omv)]), arb_series([arb(Vv)])]
        # Picard on series: y' = f(y)
        for n in range(order + 1):
            fs = F.series_eval(y, order + 2)
            y = [arb_series([arb(v)] + [fs[i].coeffs()[j]/(j + 1) if j < len(fs[i].coeffs()) else arb(0)
                                         for j in range(order + 1)]) for i, v in enumerate((Nv, omv, Vv))]
        As = Ac.series_eval(y, order + 2)[0]
        ser = {'N': y[0], 'om': y[1], 'V': y[2], 'A': As}
        out = {}
        for f in LC.FN:
            cf = ser[f].coeffs()
            for i in range(5):
                c = cf[i] if i < len(cf) else arb(0)
                out[LC.BJ[f][i]] = sp.Float((c*math.factorial(i)).mid().str(95, radius=False), 95)
        A0v = ser['A'].coeffs()[0]
        out[LC.sA] = sp.Float(A0v.sqrt().mid().str(95, radius=False), 95)
        out[LC.gV] = sp.Float((1 - arb(Vv)**2).sqrt().mid().str(95, radius=False), 95)
        out[LC.k] = kv
        return out
    finally:
        ctx.prec = old
