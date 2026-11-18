"""build-4 helpers (no side effects on import)."""
import sys
sys.path.insert(0, '<LAB>/quantum/seats/build-3')
from flint import arb
from core import sym
import qb2lib as Q
de = lambda p: arb(p[0]) + sym(arb(p[1]))


def theta_sonic(tb, x):
    """theta on the sonic manifold with X(theta) = x (interval Newton, contraction checked); same as stages_b.py"""
    sm = tb['sm']; X = [de(p) for p in sm['X']]; epsX = de(sm['epsX']); nu = arb(tb['nu_s'])
    dX = Q.s_deriv(X); x = arb(x)
    th = arb(x.mid())
    for _ in range(60):
        th = arb((th - (Q.horner(X, th) - x)/Q.horner(dX, th)).mid())
    r = arb('1e-40')
    for _ in range(100):
        Th = th + sym(r)
        g = Q.horner(X, th) + sym(epsX) - x
        Nt = th - g/(Q.horner(dX, Th) + sym(epsX/(nu - Th.abs_upper())))
        if Th.contains_interior(Nt):
            return Nt
        r *= 2
    raise ArithmeticError('theta_sonic failed')


def horner_d(c, t, eps, nu):
    """p(t) + E(t) with |E| <= eps on |t| <= nu: value, first and second derivative (Cauchy on radius nu - |t|)"""
    d1 = Q.s_deriv(c); d2 = Q.s_deriv(d1)
    r = nu - t.abs_upper()
    return (Q.horner(c, t) + sym(eps), Q.horner(d1, t) + sym(eps/r), Q.horner(d2, t) + sym(2*eps/r**2))
