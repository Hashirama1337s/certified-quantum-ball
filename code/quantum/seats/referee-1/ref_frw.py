import sympy as sp
from ref_geom import Geom
tau = sp.Symbol('tau', positive=True); x, y, z, xi = sp.symbols('x y z xi')
a2 = tau**sp.Rational(8, 9)
G = Geom(sp.diag(-1, a2, a2, a2), [tau, x, y, z], simp=sp.simplify); G.tensors()
V = G.V(xi).applyfunc(sp.simplify)
print('FRW k=1/2: R =', sp.simplify(G.R), ' Bach all zero:', all(sp.simplify(b) == 0 for b in G.B))
print('V_tautau =', sp.factor(V[0, 0]), ' expected 2(6xi-1)^2/(243 tau^4):', sp.simplify(V[0, 0] - 2*(6*xi - 1)**2/(243*tau**4)) == 0)
print('V_xx =', sp.factor(V[1, 1]), ' expected 4(6xi-1)^2/(243 tau^(28/9)):', sp.simplify(V[1, 1] - 4*(6*xi - 1)**2/(243*tau**sp.Rational(28, 9))) == 0)
a2 = tau**sp.Rational(1, 1)   # k = 1/3 radiation, a ~ tau^(1/2)
G = Geom(sp.diag(-1, a2, a2, a2), [tau, x, y, z], simp=sp.simplify); G.tensors()
print('FRW k=1/3: R =', sp.simplify(G.R), ' V(xi) all zero:', all(sp.simplify(v) == 0 for v in G.V(xi)))
