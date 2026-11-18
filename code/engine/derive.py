"""B1 derivation seat: Einstein-Euler, spherical, p = k rho, in KHA variables.

Independent symbolic derivation (sympy) from the metric, NOT transcribed from
the (garbled) PDF text.  Cross-checks against the readable parts of
Koike-Hara-Adachi gr-qc/9503007 eq. (18) and Maison gr-qc/9504008 eq. (8)
are printed at the end.

Coordinates: KHA eq. (16) ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2,
s = -ln(-t), x = ln(-r/t)  (KHA sec. 2),  N = alpha a^-1 e^-x, A = a^2,
omega = 4 pi r^2 A rho (KHA sec. 3), V = fluid 3-velocity w.r.t. Eulerian
observers, u^t = Gamma/alpha, u^r = Gamma V/a, Gamma = (1-V^2)^-1/2.

Output: engine/eqs_gen.py with plain-python expressions (only + - * / and
integer powers) usable with mpmath scalars, numpy, or our jet class.
"""
import sympy as sp

s, x, th, ph = sp.symbols('s x theta phi', real=True)
k = sp.symbols('k', positive=True)
N = sp.Function('N')(s, x)
A = sp.Function('A')(s, x)
om = sp.Function('om')(s, x)
V = sp.Function('V')(s, x)
X = [s, x, th, ph]

r2 = sp.exp(2*(x - s))                     # r^2
# ds^2 = r^2 { A[(1-N^2) ds^2 - 2 ds dx + dx^2] + dOmega^2 }
g = sp.zeros(4, 4)
g[0, 0] = r2*A*(1 - N**2)
g[0, 1] = g[1, 0] = -r2*A
g[1, 1] = r2*A
g[2, 2] = r2
g[3, 3] = r2*sp.sin(th)**2
gi = sp.simplify(g.inv())

def christoffel(g, gi):
    G = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                G[a][b][c] = sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b])
                                                         - sp.diff(g[b, c], X[d])) for d in range(4))/2)
    return G

Gam = christoffel(g, gi)

def ricci():
    R = sp.zeros(4, 4)
    for b in range(4):
        for c in range(4):
            R[b, c] = sp.simplify(sum(sp.diff(Gam[a][b][c], X[a]) for a in range(4))
                                  - sum(sp.diff(Gam[a][b][a], X[c]) for a in range(4))
                                  + sum(Gam[a][a][d]*Gam[d][b][c] for a in range(4) for d in range(4))
                                  - sum(Gam[a][c][d]*Gam[d][b][a] for a in range(4) for d in range(4)))
    return R

Ric = ricci()
Rs = sp.simplify(sum(gi[a, b]*Ric[a, b] for a in range(4) for b in range(4)))
Ein = sp.simplify(Ric - Rs*g/2)

# fluid
Gm = 1/sp.sqrt(1 - V**2)
es = sp.exp(s - x)
uup = [es*Gm/(N*sp.sqrt(A)), es*Gm/sp.sqrt(A)*(1/N + V), 0, 0]
# check normalisation
norm = sp.simplify(sum(g[a, b]*uup[a]*uup[b] for a in range(4) for b in range(4)))
print('u.u =', norm)
rho = om*sp.exp(-2*(x - s))/(4*sp.pi*A)
p = k*rho
Tup = sp.zeros(4, 4)
for a in range(4):
    for b in range(4):
        Tup[a, b] = (rho + p)*uup[a]*uup[b] + p*gi[a, b]
Tdn = sp.simplify(g*Tup*g)

E = sp.simplify(Ein - 8*sp.pi*Tdn)
divT = []
for nu in range(2):
    e = sum(sp.diff(Tup[m, nu], X[m]) for m in range(4))
    e += sum(Gam[m][m][l]*Tup[l, nu] for m in range(4) for l in range(4))
    e += sum(Gam[nu][m][l]*Tup[m, l] for m in range(4) for l in range(4))
    divT.append(sp.simplify(e))

# ---- solve Einstein constraints for A_x, N_x, A_s
Ax, Nx, As, Ns, omx, oms, Vx, Vs = sp.symbols('A_x N_x A_s N_s om_x om_s V_x V_s')
rep = {sp.Derivative(A, x): Ax, sp.Derivative(N, x): Nx, sp.Derivative(A, s): As,
       sp.Derivative(N, s): Ns, sp.Derivative(om, x): omx, sp.Derivative(om, s): oms,
       sp.Derivative(V, x): Vx, sp.Derivative(V, s): Vs}
Nn, An, on, Vn = sp.symbols('N A om V')
rep2 = {N: Nn, A: An, om: on, V: Vn}

def plain(e):
    return sp.simplify(e.subs(rep).subs(rep2))

Eeq = {(i, j): plain(E[i, j]) for i in range(4) for j in range(i, 4)}
for key, val in Eeq.items():
    if val != 0:
        fs = val.free_symbols & {Ax, Nx, As, Ns}
        print('E', key, 'contains', fs)
# pick: combos. Solve E_ss, E_sx, E_xx (first order) for Ax, Nx, As
sol = sp.solve([Eeq[(0, 0)], Eeq[(0, 1)], Eeq[(1, 1)]], [Ax, Nx, As], dict=True)
assert len(sol) == 1, sol
sol = {kk: sp.factor(sp.simplify(v)) for kk, v in sol[0].items()}
print('\nA_x/A =', sp.simplify(sol[Ax]/An))
print('N_x/N =', sp.simplify(sol[Nx]/Nn))
print('A_s   =', sol[As])

# fluid equations: substitute A_x, N_x, A_s
F = [plain(d) for d in divT]
F = [sp.simplify(f.subs(sol)) for f in F]
print('fluid eqs contain:', [f.free_symbols & {Ns, Nx, Ax, As} for f in F])
# solve the two fluid eqs for (om_s, V_s) as functions of x-derivatives
fs = sp.solve(F, [oms, Vs], dict=True)[0]
oms_e = sp.simplify(fs[oms]); Vs_e = sp.simplify(fs[Vs])

# ---- self-similar: drop s-derivatives. Fluid: M (om_x, V_x)^T = b
def strip_scale(f):
    """divide out the overall exp(c(s-x)) factor (scale invariance); assert it is overall"""
    f0 = f.subs({s: 0, x: 0})
    ratio = sp.simplify(f/f0)
    assert ratio.free_symbols <= {s, x}, ratio
    return sp.simplify(f0)

F = [strip_scale(f) for f in F]
F0 = [sp.expand(sp.numer(sp.together(f.subs({oms: 0, Vs: 0})))) for f in F]
M = sp.Matrix([[sp.diff(f, omx), sp.diff(f, Vx)] for f in F0])
b = sp.Matrix([-(f.subs({omx: 0, Vx: 0})) for f in F0])
assert sp.simplify(M*sp.Matrix([omx, Vx]) - b - sp.Matrix(F0)) == sp.zeros(2, 1)
M = sp.simplify(M); b = sp.simplify(b)
detM = sp.factor(M.det())
print('\ndet M =', detM)

# algebraic constraint from A_s = 0
cons = [c for c in sp.solve(sp.numer(sp.together(sol[As])), An) if c != 0]
print('A from A_s=0:', cons)

# ---- perturbations: h = hss(x) + eps e^{kap s} hp(x)
kap, eps = sp.symbols('kappa epsilon')
Nb, Ab, ob, Vb = sp.symbols('Nb Ab ob Vb')          # background values
Nbx, Abx, obx, Vbx = sp.symbols('Nbx Abx obx Vbx')  # background x-derivs
Np, Ap, op, Vp = sp.symbols('Np Ap op Vp')
Npx, Apx, opx, Vpx = sp.symbols('Npx Apx opx Vpx')
sub = {Nn: Nb + eps*Np, An: Ab + eps*Ap, on: ob + eps*op, Vn: Vb + eps*Vp,
       Nx: Nbx + eps*Npx, Ax: Abx + eps*Apx, omx: obx + eps*opx, Vx: Vbx + eps*Vpx,
       Ns: eps*kap*Np, As: eps*kap*Ap, oms: eps*kap*op, Vs: eps*kap*Vp}
# Full PDE set: Hamiltonian (Ax eq), slicing (Nx eq), momentum (As eq), 2 fluid
pde = [Ax - sol[Ax], Nx - sol[Nx], As - sol[As]] + [strip_scale(plain(d)) for d in divT]
lin = []
for e in pde:
    e1 = sp.diff(e.subs(sub, simultaneous=True), eps).subs(eps, 0)
    lin.append(sp.simplify(e1))
print('\nlinearised eqs built')

def code(e):
    return sp.pycode(e).replace('math.', '')

with open('<LAB>/engine/eqs_gen.py', 'w') as fh:
    fh.write('# GENERATED by derive.py - do not edit\n')
    fh.write('# self-similar background: M(N,A,om,V) [om_x,V_x] = b ; A_x, N_x algebraic\n')
    fh.write('def Ax(N,A,om,V,k):\n    return %s\n' % code(sol[Ax]))
    fh.write('def Nx(N,A,om,V,k):\n    return %s\n' % code(sol[Nx]))
    fh.write('def As_(N,A,om,V,k):\n    return %s\n' % code(sol[As]))
    for i in range(2):
        for j in range(2):
            fh.write('def M%d%d(N,A,om,V,k):\n    return %s\n' % (i, j, code(M[i, j])))
        fh.write('def b%d(N,A,om,V,k):\n    return %s\n' % (i, code(b[i])))
    fh.write('def detM(N,A,om,V,k):\n    return %s\n' % code(detM))
    fh.write('def Acons(N,om,V,k):\n    return %s\n' % code(cons[0]))
    # linear perturbation: each eq is linear in (Np,Ap,op,Vp,Npx,Apx,opx,Vpx)
    hp = [Np, Ap, op, Vp]; hpx = [Npx, Apx, opx, Vpx]
    args = 'Nb,Ab,ob,Vb,Nbx,Abx,obx,Vbx,kappa,k'
    for ie, e in enumerate(lin):
        e = sp.expand(e)
        for jv, v in enumerate(hp):
            fh.write('def P%d_%d(%s):\n    return %s\n' % (ie, jv, args, code(sp.simplify(e.coeff(v)))))
        for jv, v in enumerate(hpx):
            fh.write('def Q%d_%d(%s):\n    return %s\n' % (ie, jv, args, code(sp.simplify(e.coeff(v)))))
        rest = sp.simplify(e - sum(e.coeff(v)*v for v in hp + hpx))
        assert rest == 0, rest
print('wrote eqs_gen.py')

# ---- cross-checks vs sources
k3 = sp.Rational(1, 3)
print('\nCHECK KHA (18a) A_x/A - [1 - A + 2 om (1 + V^2/3)/(1-V^2)] at k=1/3:',
      sp.simplify((sol[Ax]/An).subs(k, k3) - (1 - An + 2*on*(1 + Vn**2/3)/(1 - Vn**2))))
print('CHECK KHA (18b) N_x/N - [-2 + A - 2 om/3] at k=1/3:',
      sp.simplify((sol[Nx]/Nn).subs(k, k3) - (-2 + An - 2*on/3)))
# Maison (8): A_M^2 = 1 - 2 omt (1 + (1+k) v/B + k v^2)/(1-v^2), A_M^2 = 1/A, omt = om/A, 1/B = N
print('CHECK Maison (8) with A_M^2=1/A, omt=om/A, v/B = N V:',
      sp.simplify(1/cons[0] - (1 - 2*(on/cons[0])*(1 + (1 + k)*Nn*Vn + k*Vn**2)/(1 - Vn**2))))
print('sonic factor check: detM / [k(N+V)^2-(1+NV)^2] =',
      sp.factor(sp.simplify(detM/(k*(Nn + Vn)**2 - (1 + Nn*Vn)**2))))
print('oms =', oms_e)
print('Vs =', Vs_e)
