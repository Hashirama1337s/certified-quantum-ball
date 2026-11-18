"""referee-1 checks of C1 (own code ref_geom.py). Each test prints PASS/FAIL with the exact residual."""
import sys, sympy as sp
from ref_geom import Geom

which = sys.argv[1:] or ['hay', 'rn', 'pp', 'gen', 'flip']
xi = sp.Symbol('xi')
th, ph = sp.symbols('theta phi')


def zero(M):
    return all(sp.simplify(sp.expand_trig(sp.simplify(e)).rewrite(sp.exp)) == 0 for e in M)


if 'hay' in which:
    # Zahn eq. (23) reading: e^{-2u}(-2 dtau^2 + 2 du^2 + dOmega^2); target eq. (33) l.488-496 at xi = 0
    tau, u = sp.symbols('tau u')
    g = sp.diag(-2*sp.exp(-2*u), 2*sp.exp(-2*u), sp.exp(-2*u), sp.exp(-2*u)*sp.sin(th)**2)
    G = Geom(g, [tau, u, th, ph], simp=sp.simplify); G.tensors()
    print('hay R =', sp.simplify(G.R))
    V0 = G.V(0).applyfunc(sp.simplify)
    # Zahn's (6): V_{mu nu}(u,y) = e^{2u} V(y); evaluate at u = 0
    vals = [sp.simplify(V0[i, i].subs(u, 0)) for i in range(3)]
    print('hay V(xi=0) tt,uu,thth at u=0:', vals, 'expected 1/80, 17/240, -7/480',
          'PASS' if vals == [sp.Rational(1, 80), sp.Rational(17, 240), sp.Rational(-7, 480)] else 'FAIL')
    tr = sp.simplify(sum(G.gi[i, i]*V0[i, i] for i in range(4)))
    print('hay trace V(0) =', tr)
    Vb = G.V(sp.Rational(1, 6)).applyfunc(sp.simplify)
    print('hay V(1/6) - B/60 zero:', zero(Vb - G.B/60))
    print('hay T2 - T1/3 + 4B zero:', zero(G.T2 - G.T1/3 + 4*G.B))
    print('hay dV/dxi nonzero (R != 0):', [sp.simplify(sp.diff(G.V(xi)[i, i], xi)) for i in range(3)])

if 'rn' in which:
    t, r = sp.symbols('t r', positive=True)
    M, Q = sp.Integer(1), sp.Rational(1, 2)
    f = 1 - 2*M/r + Q**2/r**2
    G = Geom(sp.diag(-f, 1/f, r**2, r**2*sp.sin(th)**2), [t, r, th, ph]); G.tensors()
    print('rn R =', sp.simplify(G.R))
    Vx = G.V(xi).applyfunc(sp.simplify)
    print('rn dV/dxi == 0:', zero(Vx.diff(xi)))
    print('rn V - B/60 == 0:', zero(Vx - G.B/60))
    print('rn V_tt, V_rr, V_thth:', [sp.factor(Vx[i, i]) for i in range(3)])
    print('rn design C3 expected: (4r^2-8r+1)^2/(640 r^10), -1/(120 r^6), (4r^2-8r+1)/(240 r^6)')
    exp_ = [(4*r**2 - 8*r + 1)**2/(640*r**10), -1/(120*r**6), (4*r**2 - 8*r + 1)/(240*r**6)]
    print('rn match design C3:', [sp.simplify(Vx[i, i] - exp_[i]) == 0 for i in range(3)])
    Vr3 = [Vx[i, i].subs(r, 3) for i in range(3)]
    print('rn V at r=3:', [sp.N(v, 25) for v in Vr3])
    print('rn trace:', sp.simplify(sum(G.gi[i, i]*Vx[i, i] for i in range(4))))

if 'pp' in which:
    # POST HOC: pp-wave ds^2 = -2 du dv + H du^2 + dx^2 + dy^2, H = x^4 + x y^3 (R = 0, not Ricci flat, not static)
    uu, vv, x, y = sp.symbols('u v x y')
    H = x**4 + x*y**3 + sp.sin(uu)*x**2*y**2
    g = sp.Matrix([[H, -1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    G = Geom(g, [uu, vv, x, y]); G.tensors()
    print('pp R =', sp.simplify(G.R), ' Ric_uu =', sp.factor(G.Ric[0, 0]))
    Vx = G.V(xi).applyfunc(sp.simplify)
    L2 = sp.diff(H, x, 4) + 2*sp.diff(H, x, 2, y, 2) + sp.diff(H, y, 4)
    print('pp dV/dxi == 0:', zero(Vx.diff(xi)), ' V - B/60 == 0:', zero(Vx - G.B/60))
    print('pp V_uu =', sp.simplify(Vx[0, 0]), '; PAPER prediction -(1/240) Lap^2 H =', sp.simplify(-L2/240),
          'PASS' if sp.simplify(Vx[0, 0] + L2/240) == 0 else 'FAIL')
    others = [sp.simplify(Vx[i, j]) for i in range(4) for j in range(4) if (i, j) != (0, 0)]
    print('pp other components all zero:', all(o == 0 for o in others))

if 'gen' in which:
    # generic identity T2 - T1/3 = -4B on a non-symmetric diagonal metric depending on 2 coordinates (R != 0)
    t, x, y, z = sp.symbols('t x y z')
    g = sp.diag(-(1 + t**2*x + x**3), 1 + t*x**2, (1 + x*t)**2 + t**3, 2 + x - t*x**2)
    G = Geom(g, [t, x, y, z], simp=lambda e: e); G.tensors()
    pt = {t: sp.Rational(1, 3), x: sp.Rational(-1, 5)}
    res = [sp.nsimplify(sp.simplify((G.T2[i, j] - G.T1[i, j]/3 + 4*G.B[i, j]).subs(pt))) for i in range(4) for j in range(4)]
    print('gen T2 - T1/3 + 4B at rational point: all zero =', all(v == 0 for v in res), res[:6])
    print('gen R at point =', sp.simplify(G.R.subs(pt)))
    Bt = sp.simplify(sum(G.gi[i, i]*G.B[i, i] for i in range(4)).subs(pt))
    print('gen trace B at point =', Bt)

if 'flip' in which:
    # convention null: flip the sign of the R^{cd}R_{acbd} term in T2 -> Hayward must FAIL
    tau, u = sp.symbols('tau u')
    g = sp.diag(-2*sp.exp(-2*u), 2*sp.exp(-2*u), sp.exp(-2*u), sp.exp(-2*u)*sp.sin(th)**2)
    G = Geom(g, [tau, u, th, ph], simp=sp.simplify); G.tensors()
    n = 4
    RR = sp.Matrix(n, n, lambda a, b: sp.simplify(sum(G.Ricup[c, e]*G.Rl[(a, c, b, e)] for c in range(n) for e in range(n))))
    T2f = G.T2 + 8*RR
    Vf = (sp.Rational(1, 720) - sp.Rational(1, 288))*G.T1 - T2f/240
    vals = [sp.simplify(Vf[i, i].subs(u, 0)) for i in range(3)]
    print('flip-null Hayward:', vals, 'REFUSES' if vals != [sp.Rational(1, 80), sp.Rational(17, 240), sp.Rational(-7, 480)] else 'DOES NOT REFUSE')
