import sympy as sp
from ref_geom import Geom
th, ph, tau, u = sp.symbols('theta phi tau u')
g = sp.diag(-2*sp.exp(-2*u), 2*sp.exp(-2*u), sp.exp(-2*u), sp.exp(-2*u)*sp.sin(th)**2)
G = Geom(g, [tau, u, th, ph], simp=sp.simplify); G.tensors()
print('B diag', [sp.simplify(G.B[i,i]) for i in range(4)])
print('T1 diag', [sp.simplify(G.T1[i,i]) for i in range(4)])
print('T2 diag', [sp.simplify(G.T2[i,i]) for i in range(4)])
print('resid', [sp.simplify(G.T2[i,i]-G.T1[i,i]/3+4*G.B[i,i]) for i in range(4)])
print('offdiag B', [sp.simplify(G.B[i,j]) for i in range(4) for j in range(4) if i!=j])
print('V0', [sp.simplify(G.V(0)[i,i]) for i in range(4)])
print('trace B', sp.simplify(sum(G.gi[i,i]*G.B[i,i] for i in range(4))))
