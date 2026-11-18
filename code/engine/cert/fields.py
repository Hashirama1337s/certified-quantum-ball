"""Certified pipeline, step 0: symbolic vector fields (exact, sympy) built from the
sealed derivation (engine/eqs_gen.py via engine/system.py).

  sonic field  (desingularised, polynomial):  y=(N,om,V)
      dy/dtau = Gs(y) = sv(y) * ( detM*N*g , adj(M) b )      sv = positive clearing factor
      dx/dtau = Ss(y) = sv(y) * detM
  centre field (blown up, polynomial):  Z=(P,W,U)
      dZ/dtau = Gc(Z) = Dc(Z) * Gcent(Z),   dx/dtau = Dc(Z)
  regular field (rational, for Taylor steps): y' = f(y) = (N g, M^-1 b)
  far-side reduced field: fN = N*(A-2-(1-k)om), f_om/N, f_V/N  (check divisibility)

Writes cert/fields.pkl (sympy srepr strings)."""
import sys, pickle
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import system as S

N, om, V, k, P, W, U = S.N, S.om, S.V, S.k, S.P, S.W, S.U

def clear(vec, posfactors):
    """multiply a rational vector by the lcm of denominators; report denominators"""
    vec = [sp.cancel(sp.together(e)) for e in vec]
    dens = [sp.denom(e) for e in vec]
    L = sp.lcm_list(dens)
    L = sp.factor(L)
    out = [sp.expand(sp.cancel(e*L)) for e in vec]
    return out, L

M, b = S.Msym, S.bsym
adj = M.adjugate()
w = adj*b
detM = M.det()
gN = S.gN
raw_s = [detM*gN, w[0], w[1], detM]
pol_s, Ls = clear(raw_s, None)
print('sonic clearing factor:', Ls)
for i, e in enumerate(pol_s):
    print('  comp', i, 'total degree', sp.Poly(e, N, om, V).total_degree(), 'terms', len(sp.Poly(e, N, om, V).terms()))

f = S.rhs_regular()
Gcb = S.rhs_centre_blown()
pol_c, Lc = clear(Gcb + [sp.Integer(1)], None)
print('centre clearing factor:', Lc)
for i, e in enumerate(pol_c):
    print('  comp', i, 'total degree', sp.Poly(e, P, W, U).total_degree())
print('Dc at P=0:', sp.factor(Lc.subs(P, 0)))

# far side: divisibility of f_om, f_V by N
fo = sp.cancel(sp.together(f[1])); fV = sp.cancel(sp.together(f[2]))
ro = sp.cancel(fo/N); rV = sp.cancel(fV/N)
print('f_om/N denominator has N factor?', sp.denom(ro).subs(N, 0) == 0, ' f_V/N:', sp.denom(rV).subs(N, 0) == 0)
print('f_om/N at N=0:', sp.simplify(ro.subs(N, 0)), ' f_V/N at N=0:', sp.simplify(rV.subs(N, 0)))

data = dict(pol_s=[sp.srepr(e) for e in pol_s], Ls=sp.srepr(Ls),
            pol_c=[sp.srepr(e) for e in pol_c], Lc=sp.srepr(Lc),
            f=[sp.srepr(sp.cancel(sp.together(e))) for e in f],
            ro=sp.srepr(ro), rV=sp.srepr(rV),
            Dsonic=sp.srepr(S.D_sonic), Acons=sp.srepr(S.Acons), detM=sp.srepr(sp.factor(detM)))
pickle.dump(data, open('<LAB>/engine/cert/fields.pkl', 'wb'))
print('saved')
