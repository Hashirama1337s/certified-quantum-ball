"""referee-2 own derivation (exact sympy). Usage: r2_forced.py <k>   (k = 1/3 is the EC case; others POST HOC)

Forced KHA problem G - 8 pi T - Sigma = 0, div T = 0, perturbation in HKA BAR variables
  A = Ab (1 + eps e^{kap s} a), N = Nb (1 + eps e^{kap s} n), om = ob (1 + eps e^{kap s} w), V = Vb + eps e^{kap s} v,
  Sigma_{mu nu} = eps e^{kap s} sigma_{mu nu}(x)  (ss, sx, xx, thth).
Checks:
  K1  constraint row normalised so its kappa-coefficient is 1:  first entry kappa + c0;  is c0 == -A ?
  K2  propagation along the 4 evolution rows: J' = c J + Comp;  c == -A ?
  K3  Comp == lam * D_t,  D_t = e^{2(x-s)} e^{-kap s} [(div Sigma)_s + (div Sigma)_x];  print lam
  K4  (div Sigma)_x and _s versus the theta-theta residual: Bianchi with only E_thth != 0
  K5  (k = 1/3) HKA pencil source:  Q^{-1} S = F sigma_I + s4 sigma,  M4 - Mhka = F I;  print S rows 2,3
"""
import sys, time
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import system as SY
t0 = time.time()
kv = sp.Rational(sys.argv[1]) if len(sys.argv) > 1 else sp.Rational(1, 3)
print('k =', kv)
s, x, th, ph = sp.symbols('s x theta phi', real=True)
X = [s, x, th, ph]
Af, Nf, of, Vf = [sp.Function(nm)(s, x) for nm in ('A', 'N', 'om', 'V')]
r2 = sp.exp(2*(x - s))

def metric(A, N):
    g = sp.zeros(4, 4)
    g[0, 0] = r2*A*(1 - N**2); g[0, 1] = g[1, 0] = -r2*A; g[1, 1] = r2*A
    g[2, 2] = r2; g[3, 3] = r2*sp.sin(th)**2
    return g

def chris(g, gi):
    return [[[sp.together(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                              for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]

g = metric(Af, Nf); gi = g.inv().applyfunc(sp.simplify)
Gm = chris(g, gi)
def ricci(Gm):
    R = sp.zeros(4, 4)
    for b in range(4):
        for c in range(b, 4):
            e = sum(sp.diff(Gm[a][b][c], X[a]) - sp.diff(Gm[a][b][a], X[c]) for a in range(4))
            e += sum(Gm[a][a][d]*Gm[d][b][c] - Gm[a][c][d]*Gm[d][b][a] for a in range(4) for d in range(4))
            R[b, c] = R[c, b] = sp.simplify(e)
    return R
Ric = ricci(Gm)
Rs = sp.simplify(sum(gi[a, b]*Ric[a, b] for a in range(4) for b in range(4)))
G = Ric - Rs*g/2
# fluid: u normalised (checked in the referee report by hand), rho = om/(4 pi r^2 A), p = k rho
gam = 1/sp.sqrt(1 - Vf**2); es = sp.exp(s - x)
uup = [es*gam/(Nf*sp.sqrt(Af)), es*gam/sp.sqrt(Af)*(1/Nf + Vf), 0, 0]
nrm = sp.simplify(sum(g[a, b]*uup[a]*uup[b] for a in range(4) for b in range(4)))
assert nrm == -1, nrm
rho = of/(4*sp.pi*Af*r2); p = kv*rho
Tup = sp.Matrix(4, 4, lambda a, b: (rho + p)*uup[a]*uup[b] + p*gi[a, b])
Tdn = g*Tup*g
divT = []
for nu in range(2):
    e = sum(sp.diff(Tup[m, nu], X[m]) for m in range(4))
    e += sum(Gm[m][m][l]*Tup[l, nu] for m in range(4) for l in range(4))
    e += sum(Gm[nu][m][l]*Tup[m, l] for m in range(4) for l in range(4))
    divT.append(e)
print('geometry %.0fs' % (time.time() - t0), flush=True)

# plain symbols for first derivatives
Ax, Nx, As, Ns, ox, os_, Vx, Vs = sp.symbols('A_x N_x A_s N_s om_x om_s V_x V_s')
An, Nn, on, Vn = sp.symbols('A N om V')
rep1 = {sp.Derivative(Af, x): Ax, sp.Derivative(Nf, x): Nx, sp.Derivative(Af, s): As, sp.Derivative(Nf, s): Ns,
        sp.Derivative(of, x): ox, sp.Derivative(of, s): os_, sp.Derivative(Vf, x): Vx, sp.Derivative(Vf, s): Vs}
rep2 = {Af: An, Nf: Nn, of: on, Vf: Vn}
def plain(e):
    return sp.simplify(e.subs(rep1).subs(rep2))
S_ss, S_sx, S_xx = sp.symbols('S_ss S_sx S_xx')
Edn = G - 8*sp.pi*Tdn
e00 = plain(Edn[0, 0]) - S_ss; e01 = plain(Edn[0, 1]) - S_sx; e11 = plain(Edn[1, 1]) - S_xx
for nm, e in (('Ess', e00), ('Esx', e01), ('Exx', e11)):
    assert not (e.free_symbols & {s, x}), nm
sol = sp.solve([e00, e01, e11], [Ax, Nx, As], dict=True); assert len(sol) == 1; sol = sol[0]
fl = []
for d in divT:
    e = plain(d*r2**2)       # weight: div of T^{mu nu} ~ e^{-4(x-s)}
    e0 = e.subs({s: 0, x: 0}); assert sp.simplify(e - e0) == 0, 'fluid row not scale-free'
    fl.append(e0)
print('rows %.0fs' % (time.time() - t0), flush=True)

# linearise in BAR variables
eps, kap = sp.symbols('epsilon kappa')
Ab, Nb, ob, Vb, Abx, Nbx, obx, Vbx = sp.symbols('Ab Nb ob Vb Abx Nbx obx Vbx')
a, n, w, v, ax, nx, wx, vx = sp.symbols('a n w v a_x n_x w_x v_x')
sg = sp.symbols('sg_ss sg_sx sg_xx')
sub = {An: Ab*(1 + eps*a), Nn: Nb*(1 + eps*n), on: ob*(1 + eps*w), Vn: Vb + eps*v,
       Ax: Abx*(1 + eps*a) + eps*Ab*ax, Nx: Nbx*(1 + eps*n) + eps*Nb*nx,
       ox: obx*(1 + eps*w) + eps*ob*wx, Vx: Vbx + eps*vx,
       As: eps*kap*Ab*a, Ns: eps*kap*Nb*n, os_: eps*kap*ob*w, Vs: eps*kap*v,
       S_ss: eps*sg[0], S_sx: eps*sg[1], S_xx: eps*sg[2]}
rows = [Ax - sol[Ax], Nx - sol[Nx]] + fl + [As - sol[As]]   # Ham, slice, fl1, fl2, mom
lin = [sp.expand(sp.diff(e.subs(sub, simultaneous=True), eps).subs(eps, 0)) for e in rows]
Y = [a, n, w, v]; Yx = [ax, nx, wx, vx]
Pm = sp.Matrix(5, 4, lambda i, j: lin[i].coeff(Y[j]))
Qm = sp.Matrix(5, 4, lambda i, j: lin[i].coeff(Yx[j]))
Cm = sp.Matrix(5, 3, lambda i, j: lin[i].coeff(sg[j]))
for i in range(5):
    rest = sp.simplify(lin[i] - sum(Pm[i, j]*Y[j] + Qm[i, j]*Yx[j] for j in range(4)) - sum(Cm[i, j]*sg[j] for j in range(3)))
    assert rest == 0, (i, rest)
# on-shell background (engine sealed background ODE, general k)
Nsy, osy, Vsy = SY.N, SY.om, SY.V
Acons = sp.cancel(SY.Acons.subs(SY.k, kv))
f = [sp.cancel(sp.together(e.subs(SY.k, kv))) for e in SY.rhs_regular()]
def ddx(e):
    return sp.cancel(sp.together(sum(sp.diff(e, y)*d for y, d in zip((Nsy, osy, Vsy), f))))
dA = ddx(Acons)
bg = {Ab: Acons, Nb: Nsy, ob: osy, Vb: Vsy, Abx: dA, Nbx: f[0], obx: f[1], Vbx: f[2]}
def on(M):
    return M.subs(bg).applyfunc(lambda e: sp.cancel(sp.together(e)))
P = on(Pm); Q = on(Qm); C = on(Cm)
assert all(e == 0 for e in Q[4, :])
print('linearised + on-shell %.0fs' % (time.time() - t0), flush=True)

# K1: constraint row
Kp = P[4, :]; Kc = C[4, :]
nrmK = sp.cancel(sp.diff(Kp[0], kap))
assert all(sp.cancel(sp.diff(Kp[j], kap)) == 0 for j in range(1, 4))
Kn = (Kp/nrmK).applyfunc(sp.cancel); sI = (-Kc/nrmK).applyfunc(sp.cancel)   # J = Kn.Y - sI.sigma
c0 = sp.cancel(Kn[0] - kap)
print('K1 constraint first entry - kappa =', sp.factor(c0), ';  == -A :', sp.cancel(c0 + Acons) == 0)
# K2: propagation along the 4 evolution rows
ev = [0, 1, 2, 3]
Q4 = Q.extract(ev, list(range(4))); P4 = P.extract(ev, list(range(4))); C4 = C.extract(ev, list(range(3)))
Q4i = Q4.inv().applyfunc(lambda e: sp.cancel(sp.together(e)))
M4 = (-Q4i*P4).applyfunc(lambda e: sp.cancel(sp.together(e)))
s4 = (-Q4i*C4).applyfunc(lambda e: sp.cancel(sp.together(e)))
H = (Kn.applyfunc(ddx) + Kn*M4).applyfunc(lambda e: sp.cancel(sp.together(e)))
cc = sp.cancel(H[0]/Kn[0])
okprop = all(sp.cancel(H[j] - cc*Kn[j]) == 0 for j in range(4))
print('K2 dJ/dx = c J on the homogeneous system:', okprop, '; c =', sp.factor(cc), '; c == -A :', sp.cancel(cc + Acons) == 0)
sgx = sp.symbols('sgx_ss sgx_sx sgx_xx')
Comp = (Kn*s4*sp.Matrix(sg))[0] + cc*(sI*sp.Matrix(sg))[0] - (sI.applyfunc(ddx)*sp.Matrix(sg))[0] - (sI*sp.Matrix(sgx))[0]
Comp = sp.cancel(sp.together(Comp))
print('K2 kappa in Comp:', kap in Comp.free_symbols, flush=True)

# K3: covariant divergence of Sigma on the background
gb = metric(Af, Nf)   # background functions of x only after substitution below
sfun = [sp.Function(nm)(x) for nm in ('fss', 'fsx', 'fxx', 'fTT')]
ek = sp.exp(kap*s)
Sig = sp.zeros(4, 4)
Sig[0, 0] = ek*sfun[0]; Sig[0, 1] = Sig[1, 0] = ek*sfun[1]; Sig[1, 1] = ek*sfun[2]
Sig[2, 2] = ek*sfun[3]; Sig[3, 3] = ek*sfun[3]*sp.sin(th)**2
Smix = gi*Sig    # Sigma^mu_nu
div = []
for nu in range(4):
    e = sum(sp.diff(Smix[m, nu], X[m]) for m in range(4))
    e += sum(Gm[m][m][l]*Smix[l, nu] for m in range(4) for l in range(4))
    e -= sum(Gm[l][m][nu]*Smix[m, l] for m in range(4) for l in range(4))
    div.append(e)
sTT, sxTT = sp.symbols('sg_TT sgx_TT')
bgs = {sp.Derivative(Af, s): 0, sp.Derivative(Nf, s): 0, sp.Derivative(of, s): 0, sp.Derivative(Vf, s): 0}
fsub = {sp.Derivative(sfun[i], x): (sgx[i] if i < 3 else sxTT) for i in range(4)}
fval = {sfun[0]: sg[0], sfun[1]: sg[1], sfun[2]: sg[2], sfun[3]: sTT}
divp = []
for nu in range(2):
    e = div[nu].subs(bgs).subs(fsub).subs(fval)
    e = e.subs({sp.Derivative(Af, x): Abx, sp.Derivative(Nf, x): Nbx}).subs({Af: Ab, Nf: Nb})
    e = sp.simplify(e*sp.exp(2*(x - s))/ek)
    assert not (e.free_symbols & {s, x, th}), e.free_symbols
    divp.append(sp.cancel(sp.together(e.subs(bg))))
for nu in (2, 3):
    e = sp.simplify(div[nu].subs(bgs))
    print('    (div Sigma)_%d =' % nu, e)
Dt = sp.cancel(divp[0] + divp[1])
lam = sp.cancel(sp.together(Comp/Dt))
print('K3 Comp / D_t free of sigma:', not (lam.free_symbols & (set(sg) | set(sgx) | {sTT, sxTT})), '; lam =', sp.factor(lam))
print('K3 lam == -A :', sp.cancel(lam + Acons) == 0, '; lam == -A^2 :', sp.cancel(lam + Acons**2) == 0)
print('K3 sigma_TT in D_t:', sTT in Dt.free_symbols, '; in D_x:', sTT in divp[1].free_symbols, flush=True)

# K4: Bianchi with only E_thth != 0 (Eth = E_{theta theta}, E_phph = sin^2 Eth), background, lower-index divergence
Eth = sp.Function('Eth')(x)
Em = sp.zeros(4, 4); Em[2, 2] = Eth; Em[3, 3] = Eth*sp.sin(th)**2
Emix = gi*Em
dE = []
for nu in range(2):
    e = sum(sp.diff(Emix[m, nu], X[m]) for m in range(4))
    e += sum(Gm[m][m][l]*Emix[l, nu] for m in range(4) for l in range(4))
    e -= sum(Gm[l][m][nu]*Emix[m, l] for m in range(4) for l in range(4))
    dE.append(sp.simplify(e.subs(bgs)))
print('K4 div(E_thth only)_s =', dE[0], ';  _x =', dE[1], '; sum =', sp.simplify(dE[0] + dE[1]))

# K5: HKA pencil source (k = 1/3 only)
if kv == sp.Rational(1, 3):
    gg = kv + 1; q = 1 - Vsy**2; Nn_, om_, V_ = Nsy, osy, Vsy
    Ms = sp.Matrix([[1, gg*V_/q], [(gg - 1)*V_, gg/q]])
    Mx = sp.Matrix([[1 + Nn_*V_, gg*(Nn_ + V_)/q], [(gg - 1)*(Nn_ + V_), gg*(1 + Nn_*V_)/q]])
    Irow = sp.Matrix([[kap - Acons, 2*gg*Nn_*V_*om_/q, 2*om_*(1 + (gg - 1)*V_**2 + gg*Nn_*V_)/q,
                       2*gg*om_*(Nn_*(1 + V_**2) + 2*V_)/q**2]])
    print('K5 my constraint row == HKA identity row:', all(sp.cancel(Kn[j] - Irow[j]) == 0 for j in range(4)))
    import pickle
    d = pickle.load(open('<LAB>/engine/cert/lin_fields.pkl', 'rb'))
    loc = {'N': Nsy, 'om': osy, 'V': Vsy, 'kappa': SY.kappa}
    Qh = sp.Matrix(4, 4, [sp.sympify(sp.parse_expr(e, evaluate=True)) if False else eval(e, sp.__dict__.copy()) for e in d['Q']])
    P0h = sp.Matrix(4, 4, [eval(e, sp.__dict__.copy()) for e in d['P0']])
    P1h = sp.Matrix(4, 4, [eval(e, sp.__dict__.copy()) for e in d['P1']])
    rn = {sp.Symbol('kappa'): kap}
    Mh = (Qh.inv()*(P0h + kap*P1h)).subs(rn).applyfunc(lambda e: sp.cancel(sp.together(e)))
    Dm = (M4 - Mh).applyfunc(lambda e: sp.cancel(sp.together(e)))
    Fv = []
    for i in range(4):
        r = [Dm[i, j] for j in range(4)]
        if all(e == 0 for e in r):
            Fv.append(0); print('    row', i, 'M4 == Mhka'); continue
        fac = sp.cancel(r[0]/Kn[0])
        ok = all(sp.cancel(r[j] - fac*Kn[j]) == 0 for j in range(4))
        Fv.append(fac); print('    row', i, 'M4 - Mhka = F_i * I :', ok, ' F_i =', sp.factor(fac))
    Fv = sp.Matrix(Fv)
    QiS = (Fv*(sI*sp.Matrix(sg)) + s4*sp.Matrix(sg)).applyfunc(lambda e: sp.cancel(sp.together(e)))
    Sh = (Qh.subs(rn)*QiS).applyfunc(lambda e: sp.cancel(sp.together(e)))
    for i in range(4):
        print('K5 HKA source row %d:' % i, [sp.factor(sp.diff(Sh[i], sg[j])) for j in range(3)])
    print('K5 rows 2,3 identically zero:', Sh[2] == 0 and Sh[3] == 0)
    # Lemma (a) in HKA form: J' = -A J + I.Q^{-1}S - sI' - A sI ; compare with Comp
    CompH = sp.cancel(sp.together((Kn*QiS)[0] - Acons*(sI*sp.Matrix(sg))[0] - (sI.applyfunc(ddx)*sp.Matrix(sg))[0] - (sI*sp.Matrix(sgx))[0]))
    print('K5 (rate -A along HKA pencil) Comp via HKA pencil+full source == Comp via evolution rows:', sp.cancel(CompH - Comp) == 0)
    print('K5 sigma_I coefficients:', [sp.factor(e) for e in sI])
    print('K5 c_derive - (-A) == I.F :', sp.cancel(cc + Acons - (Kn*Fv)[0]) == 0)
    den = [sp.factor(sp.denom(sp.together(sp.diff(Sh[i], sg[j])))) for i in (2, 3) for j in range(3)]
    print('K5 denominators rows 2,3:', den)
print('total %.0fs' % (time.time() - t0))
