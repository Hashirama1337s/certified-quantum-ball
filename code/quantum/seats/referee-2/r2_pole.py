"""referee-2 (exact sympy, POST HOC): at the sonic point, the identity row's fluid part annihilates r = ker M_x.
Hence Delta2(kappa) = const * (l E_u r - kappa l M_s r), i.e. kappa_pole = d exactly."""
import sympy as sp
N, om, V = sp.symbols('N om V', positive=True)
g = sp.Rational(4, 3); q = 1 - V**2
Mx = sp.Matrix([[1 + N*V, g*(N + V)/q], [(g - 1)*(N + V), g*(1 + N*V)/q]])
r = sp.Matrix([g*(N + V)/q, -(1 + N*V)])            # (Mx r)[0] == 0 identically
cf = sp.Matrix([[2*om*(1 + (g - 1)*V**2 + g*N*V)/q, 2*g*om*(N*(1 + V**2) + 2*V)/q**2]])
D = 3*N**2*V**2 - N**2 + 4*N*V - V**2 + 3
print('Mx r row0 == 0:', sp.simplify((Mx*r)[0]) == 0, '; Mx r row1 / D:', sp.factor((Mx*r)[1]))
print('det Mx =', sp.factor(Mx.det()))
print('c_f . r =', sp.factor((cf*r)[0]), ' ; ratio to D:', sp.factor((cf*r)[0]/D))
