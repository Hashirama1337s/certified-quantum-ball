"""B4 build seat: general-k symbolic ingredients for stages S and L (new file).
HKA pencil (PROOF-B4 sec. 1.1, gamma = 1 + k), unknowns Y = (Abar_p, Nbar_p, wbar_p, V_p).

Universal step matrix, for any independent variable t along the background (t = x, tau, or the
sonic-manifold parameter theta), with x' = dx/dt and primes = d/dt:
     dY/dt = x' Q^{-1}(P0 + kappa P1)(y, y'/x') Y
           = [ x' Ptop ; (F0 + kappa F1)/detMx ] Y          (rows 0,1 ; rows 2,3)
     F0 = adj(Mx) (x' Pa + Pb om' + Pc V'),   F1 = adj(Mx)(-Ms) x'   (fluid columns 2,3)
where P_fluid = Pa + Pb omx + Pc Vx is affine in the background derivatives (checked).  No division
by x' appears, and near the sonic point theta*A = [theta x' Ptop ; F/(detMx/theta)] is analytic.

Centre (blown-up chart P, W, U): Mc(Z, kappa) = Q^{-1}(P0 + kappa P1) with omx, Vx from the field,
N = 1/P, om = W P^2, V = U P; checked against the proof seat's lin_fields.pkl at k = 1/3.
Writes lin_fields2.pkl.
"""
import sys, pickle
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import system as S

kk = S.k                     # same symbol as engine Acons / rhs_regular
N, om, V, kap = S.N, S.om, S.V, S.kappa
omp, Vp, xp = sp.symbols('omp Vp xp')
omx, Vx = sp.symbols('omx Vx')
g = 1 + kk
q = 1 - V**2
A = S.Acons            # in N, om, V, k

def hka(Ax_omx, Ax_Vx):
    wbx = Ax_omx/om
    Ms = sp.Matrix([[1, g*V/q], [(g - 1)*V, g/q]])
    Mx = sp.Matrix([[1 + N*V, g*(N + V)/q], [(g - 1)*(N + V), g*(1 + N*V)/q]])
    E1 = -(g + 2)/2*A*N*V
    E2 = (6 - 3*g)/2*N*V - (2 + g)/2*A*N*V + (2 - g)*om*N*V - N*V*wbx - g*N*Ax_Vx/q
    E3 = (2 - g)*om*N*V
    E4 = (6 - 3*g)/2*N - (2 + g)/2*A*N + (2 - g)*om*N - N*wbx - g*(1 + 2*N*V + V**2)*Ax_Vx/q**2
    F1 = (2 - 3*g)/2*A*N
    F2 = (2 - g)*(g - 1)*om*N + (7*g - 6)/2*N + (2 - 3*g)/2*A*N - (g - 1)*N*wbx - g*N*V*Ax_Vx/q
    F3 = (2 - g)*(g - 1)*om*N
    F4 = -(g - 1)*wbx - g*(N + 2*V + N*V**2)*Ax_Vx/q**2
    P0 = sp.Matrix([[-A, 0, 2*(1 + (g - 1)*V**2)*om/q, 4*g*om*V/q**2], [A, 0, (g - 2)*om, 0],
                    [E1, E2, E3, E4], [F1, F2, F3, F4]])
    P1 = sp.zeros(4, 4); P1[2:, 2:] = -Ms
    Irow = [kap - A, 2*g*N*V*om/q, 2*om*(1 + (g - 1)*V**2 + g*N*V)/q, 2*g*om*(N*(1 + V**2) + 2*V)/q**2]
    return Ms, Mx, P0, P1, Irow

Ms, Mx, P0, P1, Irow = hka(omx, Vx)
Ptop = P0[:2, :]
Pf = P0[2:, :]
Pa = Pf.subs({omx: 0, Vx: 0}); Pb = Pf.diff(omx); Pc = Pf.diff(Vx)
assert (Pf - Pa - Pb*omx - Pc*Vx).applyfunc(sp.expand) == sp.zeros(2, 4), 'P_fluid not affine'
assert Pb.diff(omx) == sp.zeros(2, 4) and Pc.diff(Vx) == sp.zeros(2, 4)
adjMx = Mx.adjugate(); detMx = sp.factor(Mx.det())
F0 = adjMx*(xp*Pa + Pb*omp + Pc*Vp)
F1 = sp.zeros(2, 4); F1[:, 2:] = adjMx*(-Ms)*xp
# check: [xp Ptop ; (F0 + kap F1)/detMx] == xp Q^-1 (P0 + kap P1) with omx = omp/xp, Vx = Vp/xp
Qm = sp.eye(4); Qm[2:, 2:] = Mx
lhs = xp*Qm.inv()*(P0 + kap*P1)
lhs = lhs.subs({omx: omp/xp, Vx: Vp/xp})
rhs = sp.Matrix.vstack(xp*Ptop, (F0 + kap*F1)/detMx)
ok = all(sp.cancel(sp.together(e)) == 0 for e in (lhs - rhs))
print('universal step matrix identity:', ok)
# centre
f = S.rhs_regular()
Mhka = Qm.inv()*(P0 + kap*P1)
Mhka = Mhka.subs({omx: f[1], Vx: f[2]})
cen = {N: 1/S.P, om: S.W*S.P**2, V: S.U*S.P}
Mc = Mhka.subs(cen).applyfunc(lambda e: sp.cancel(sp.together(e)))
d = pickle.load(open('<LAB>/engine/cert/lin_fields.pkl', 'rb'))
Mc_old = [sp.sympify(s) for s in d['Mc']]
chk = all(sp.cancel(sp.together(Mc[i] .subs(kk, sp.Rational(1, 3)) - Mc_old[i].subs(S.k, sp.Rational(1, 3)))) == 0 for i in range(16))
print('centre Mc(k) at k=1/3 == proof-seat lin_fields.pkl Mc:', chk)
Ic = [sp.cancel(sp.together(e.subs(cen))) for e in Irow]
# exact sonic quantities
Abar_x = 1 - A + 2*om*(1 + kk*V**2)/(1 - V**2)
kappa_g = 2 - A + (1 - kk)*om
pickle.dump(dict(Ptop=[sp.srepr(e) for e in Ptop], F0=[sp.srepr(e) for e in F0], F1=[sp.srepr(e) for e in F1],
                 detMx=sp.srepr(detMx), Irow=[sp.srepr(e) for e in Irow], Mc=[sp.srepr(e) for e in Mc],
                 Ic=[sp.srepr(e) for e in Ic], Abar_x=sp.srepr(Abar_x), kappa_g=sp.srepr(kappa_g)),
            open('<LAB>/engine/cert/lin_fields2.pkl', 'wb'))
print('detMx =', detMx)
print('saved lin_fields2.pkl')
