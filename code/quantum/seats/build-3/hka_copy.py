"""build-3: verbatim copy of hka() from engine/cert/lin_mc.py (sha256 of lin_mc.py recorded in controls_qb2.json).
lin_mc.py itself cannot be imported: at import it re-runs its checks and WRITES engine/cert/lin_fields2.pkl."""
import sys
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
