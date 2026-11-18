"""Load the generated equations (eqs_gen.py) as sympy expressions and build
the background ODE in the two charts used by B1:
  sonic chart  y = (N, om, V),   A = Acons(N,om,V)   (Maison eq. (8) / our A_s = 0)
  centre chart Y = (P, om, V),   P = 1/N, independent variable u = e^x
"""
import re, sympy as sp

SRC = '<LAB>/engine/eqs_gen.py'
N, A, om, V, k, P = sp.symbols('N A om V k P')
Nb, Ab, ob, Vb, Nbx, Abx, obx, Vbx, kappa = sp.symbols('Nb Ab ob Vb Nbx Abx obx Vbx kappa')

def load():
    txt = open(SRC).read()
    out = {}
    for m in re.finditer(r'def (\w+)\(([^)]*)\):\n    return (.*)\n', txt):
        name, args, body = m.group(1), m.group(2), m.group(3)
        out[name] = sp.sympify(body, locals={'N': N, 'A': A, 'om': om, 'V': V, 'k': k,
                                              'Nb': Nb, 'Ab': Ab, 'ob': ob, 'Vb': Vb, 'Nbx': Nbx,
                                              'Abx': Abx, 'obx': obx, 'Vbx': Vbx, 'kappa': kappa})
    return out

E = load()
Acons = E['Acons']
Msym = sp.Matrix([[E['M00'], E['M01']], [E['M10'], E['M11']]]).subs(A, Acons)
bsym = sp.Matrix([E['b0'], E['b1']]).subs(A, Acons)
gN = (N*(A + k*om - om - 2)).subs(A, Acons)          # N_x   (KHA (18b))
det = sp.factor(sp.simplify(Msym.det()))
D_sonic = k*(N + V)**2 - (1 + N*V)**2                 # det M ~ D_sonic

def rhs_regular():
    """explicit y' = f(y) away from the sonic point"""
    adj = Msym.adjugate()
    w = adj*bsym
    dt = Msym.det()
    return [gN, w[0]/dt, w[1]/dt]

def rhs_centre():
    """u dY/du = F(Y), Y = (P, om, V); P = 1/N"""
    f = rhs_regular()
    FP = sp.cancel(sp.together((-P**2*f[0]).subs(N, 1/P)))
    Fo = sp.cancel(sp.together(f[1].subs(N, 1/P)))
    FV = sp.cancel(sp.together(f[2].subs(N, 1/P)))
    return [FP, Fo, FV]

W, U = sp.symbols('W U')

def rhs_centre_blown():
    """u dZ/du = G(Z), Z = (P, W, U), om = W P^2, V = U P (blow-up of the centre)"""
    FP, Fo, FV = rhs_centre()
    sub = {om: W*P**2, V: U*P}
    FP = FP.subs(sub)
    GP = sp.cancel(sp.together(FP))
    GW = sp.cancel(sp.together(Fo.subs(sub)/P**2 - 2*W*FP/P))
    GU = sp.cancel(sp.together(FV.subs(sub)/P - U*FP/P))
    return [GP, GW, GU]

if __name__ == '__main__':
    F = rhs_centre_blown()
    Y = [P, W, U]
    G0 = [sp.simplify(Fi.subs(P, 0)) for Fi in F]
    print('G at P=0:', G0)
    Ueq = sp.solve(G0[2], U)
    print('U equilibria at P=0:', Ueq)
    for ue in Ueq:
        J = sp.Matrix([[sp.simplify(sp.diff(Fi, y).subs({P: 0}).subs(U, ue)) for y in Y] for Fi in F])
        print('U=', ue, ' J =', J, ' eig:', J.eigenvals())
    import sys; sys.exit()
    F = rhs_centre()
    Y = [P, om, V]
    J = sp.Matrix([[sp.simplify(sp.diff(Fi, y).subs({P: 0, om: 0, V: 0})) for y in Y] for Fi in F])
    print('centre Jacobian J0 =', J)
    print('eigen:', J.eigenvects())
    print('F at 0:', [sp.simplify(Fi.subs({P: 0, om: 0, V: 0})) for Fi in F])
    print('det M factor:', det)
