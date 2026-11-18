"""B4 (proof seat), new file.  Symbolic facts about the linearised Evans-Coleman problem.

Linear system in HKA bar variables  Y = (Abar_p, Nbar_p, wbar_p, V_p)  (Hara-Koike-Adachi
gr-qc/9607010, theory/src_tex/rflanl.tex:4991-5090 coefficients, 5132-5215 eigenmode ODE and
identity), written with the background in ENGINE variables (N, om, V), A = Acons (Maison (8)),
background derivatives from engine/system.py (sealed derivation).

Checks printed (each must print 0 / True):
  C1  background: d/dx Acons along the flow == A*(1 - A + 2 om (1+kV^2)/(1-V^2))   (KHA 18a)
  C2  HKA eigenmode system == linearisation of the sealed derive.py equations
      (eqs_gen P/Q rows 0,1,3,4) under Np = N Nbar_p, Ap = A Abar_p, op = om wbar_p, Vp = V_p
  C3  HKA identity row  ||  derive.py momentum row (row 2)
  C4  constraint propagation:  d/dx I(Y) = c(x) I(Y)  along solutions; prints c(x)
  C5  eigenvalues of Ms^-1 Mx are lambda_pm = 1 + N (V pm c)/(1 pm c V)
  C6  centre limit of the system matrix (HKA PPasmp1) in the blown-up chart
Writes lin_fields.pkl (srepr of the 4x4 pencil  Q(x) Y' = (P0(x) + kappa P1(x)) Y  and identity row).
"""
import sys, pickle
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import system as S

k = sp.Rational(1, 3); g = k + 1
N, om, V, kap = S.N, S.om, S.V, S.kappa
A = sp.cancel(S.Acons.subs(S.k, k))
f = [sp.cancel(sp.together(e.subs(S.k, k))) for e in S.rhs_regular()]
Nx, omx, Vx = f
Abar_x = 1 - A + 2*om*(1 + k*V**2)/(1 - V**2)

def zero(e):
    return sp.cancel(sp.together(e)) == 0

# C1
dA = sum(sp.diff(A, s)*d for s, d in zip((N, om, V), (Nx, omx, Vx)))
print('C1 background A_x consistency:', zero(dA - A*Abar_x))
Ax = A*Abar_x

# ---- HKA matrices (bar variables), gamma = 4/3
q = 1 - V**2
wbx = omx/om
Ms = sp.Matrix([[1, g*V/q], [(g - 1)*V, g/q]])
Mx = sp.Matrix([[1 + N*V, g*(N + V)/q], [(g - 1)*(N + V), g*(1 + N*V)/q]])
E1 = -(g + 2)/2*A*N*V
E2 = (6 - 3*g)/2*N*V - (2 + g)/2*A*N*V + (2 - g)*om*N*V - N*V*wbx - g*N*Vx/q
E3 = (2 - g)*om*N*V
E4 = (6 - 3*g)/2*N - (2 + g)/2*A*N + (2 - g)*om*N - N*wbx - g*(1 + 2*N*V + V**2)*Vx/q**2
F1 = (2 - 3*g)/2*A*N
F2 = (2 - g)*(g - 1)*om*N + (7*g - 6)/2*N + (2 - 3*g)/2*A*N - (g - 1)*N*wbx - g*N*V*Vx/q
F3 = (2 - g)*(g - 1)*om*N
F4 = -(g - 1)*wbx - g*(N + 2*V + N*V**2)*Vx/q**2
G1, G3, G4 = -A, 2*(1 + (g - 1)*V**2)*om/q, 4*g*om*V/q**2
H1, H3 = A, (g - 2)*om
Qm = sp.diag(1, 1, 1, 1); Qm[2:, 2:] = Mx
P0 = sp.Matrix([[G1, 0, G3, G4], [H1, 0, H3, 0], [E1, E2, E3, E4], [F1, F2, F3, F4]])
P1 = sp.zeros(4, 4); P1[2:, 2:] = -Ms
Irow = sp.Matrix([[kap - A, 2*g*N*V*om/q, 2*om*(1 + (g - 1)*V**2 + g*N*V)/q,
                   2*g*om*(N*(1 + V**2) + 2*V)/q**2]])

# ---- C2/C3: derive.py linearisation
sub = {S.Nb: N, S.Ab: A, S.ob: om, S.Vb: V, S.Nbx: Nx, S.Abx: Ax, S.obx: omx, S.Vbx: Vx, S.k: k}
def row(kind, ie):
    return sp.Matrix([[S.E['%s%d_%d' % (kind, ie, j)].subs(sub) for j in range(4)]])
# derive order hp = (Np, Ap, op, Vp); HKA order Y = (Ab, Nb, wb, V)
# hp = T Y  with  T mapping: Np = N*Y1, Ap = A*Y0, op = om*Y2, Vp = Y3
T = sp.zeros(4, 4); T[0, 1] = N; T[1, 0] = A; T[2, 2] = om; T[3, 3] = 1
Tx = sp.zeros(4, 4); Tx[0, 1] = Nx; Tx[1, 0] = Ax; Tx[2, 2] = omx
rowsP = sp.Matrix.vstack(*[row('P', i) for i in (0, 1, 3, 4)])
rowsQ = sp.Matrix.vstack(*[row('Q', i) for i in (0, 1, 3, 4)])
Lq = rowsQ*T                 # coefficient of Y'
Lp = rowsP*T + rowsQ*Tx      # coefficient of Y
Mder = -(Lq.inv())*Lp
Mhka = Qm.inv()*(P0 + kap*P1)
diff = (Mder - Mhka).applyfunc(lambda e: sp.cancel(sp.together(e)))
print('C2a HKA system == derive.py linearisation entrywise:', all(e == 0 for e in diff))
c2 = True
for i in range(4):
    r = [diff[i, j] for j in range(4)]
    if all(e == 0 for e in r):
        print('    row', i, ': identical'); continue
    rr = [sp.cancel(sp.together(r[j]/Irow[j])) for j in range(4)]
    okr = all(zero(rr[j] - rr[0]) for j in range(4)); c2 = c2 and okr
    print('    row', i, ': difference = multiple of identity row I:', okr, ' factor', sp.factor(rr[0]))
print('C2b systems agree modulo the identity I (equal on I = 0):', c2)
mom = row('P', 2)*T
assert all(zero(e) for e in row('Q', 2))
rat = [sp.cancel(sp.together(mom[j]/Irow[j])) for j in range(4)]
print('C3 momentum row / HKA identity row ratios equal:', all(zero(r - rat[0]) for r in rat), ' ratio =', sp.factor(rat[0]))

# ---- C4 constraint propagation
bgd = {N: Nx, om: omx, V: Vx}
dI = sp.Matrix([[sum(sp.diff(Irow[j], s)*d for s, d in bgd.items()) for j in range(4)]])
Mprop = dI + Irow*Mhka
cc = sp.cancel(sp.together(Mprop[0]/Irow[0]))
ok = all(zero(Mprop[j] - cc*Irow[j]) for j in range(4))
print('C4 dI/dx = c(x) I :', ok)
print('    c(x) =', sp.factor(sp.simplify(cc)))

# ---- C5 characteristic speeds
c = sp.sqrt(k)
lp = 1 + N*(V + c)/(1 + c*V); lm = 1 + N*(V - c)/(1 - c*V)
L = Ms.inv()*Mx
print('C5 trace, det of Ms^-1 Mx match lambda_pm:',
      zero(sp.simplify(L.trace() - lp - lm)), zero(sp.simplify(L.det() - lp*lm)))

# ---- C6 centre: blown-up chart P=1/N, om = W P^2, V = U P, independent u = e^x
Pb, W, U = S.P, S.W, S.U
cen = {N: 1/Pb, om: W*Pb**2, V: U*Pb}
Mc = Mhka.subs(cen).applyfunc(lambda e: sp.cancel(sp.together(e)))
U0 = -sp.Rational(2, 3)/(1 + k)
M00 = Mc.subs(Pb, 0).applyfunc(sp.simplify).subs(U, U0).applyfunc(sp.simplify)
print('C6 centre limit (P -> 0, U = -2/(3(1+k))):', M00.tolist(), ' eig:', M00.eigenvals())
# order of P in denominators at P=0 (regularity of the centre expansion)
bad = [(i, j) for i in range(4) for j in range(4) if sp.denom(Mc[i, j]).subs(Pb, 0) == 0]
print('C6 entries singular at P=0 (must be empty):', bad)
Ic = Irow.subs(cen).applyfunc(lambda e: sp.cancel(sp.together(e)))
print('    identity row at centre:', [sp.simplify(e.subs(Pb, 0)) for e in Ic])

pickle.dump(dict(Q=[sp.srepr(e) for e in Qm], P0=[sp.srepr(e) for e in P0], P1=[sp.srepr(e) for e in P1],
                 I=[sp.srepr(e) for e in Irow], cprop=sp.srepr(cc), Mc=[sp.srepr(e) for e in Mc],
                 Ic=[sp.srepr(e) for e in Ic]),
            open('<LAB>/engine/cert/lin_fields.pkl', 'wb'))
print('saved lin_fields.pkl')
