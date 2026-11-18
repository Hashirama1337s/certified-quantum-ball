"""ax_sys.py - build seat build-4 (B7 SB2 "The Light Cone", BAR-SF-L2 seal 9077a35520), 2026-11-17.

Symbolic layer of the l = 2 axial count.  Read-only inputs: build-1 eqs_l.pkl (own derivation, AX-Pi rows),
engine/eqs_gen.py via engine/system.py (background field), engine/cert/fields.pkl (x-field f, centre field pol_c).

First-order system used for the count (Q = kappa h_x - h_s', u = g h_s, g = (L - 2)/2 on shell, PROOF-SF-X1 PX2b-3):
    al_s Q' = -be_s Q - u                       (E_{s phi} row)
    u'      = -kappa (al_x Q' + be_x Q) - g Q   (E_{x phi} row and Q = kappa h_x - h_s')
i.e.  y' = (A0 + kappa A1 + kappa^2 A2) y,  y = (Q, u)  (x-derivative).  Only first background jets enter
(N_x = Nx(N, A, om, V), A_x = Ax(N, A, om, V), both algebraic and regular at the sonic point), so the system can be
evaluated on every chart (tau steps, sonic manifolds, x steps, centre chart) without second jets.

Checks (all exact rational-function identities, sympy cancel; the sqrt symbols sA, gV do not occur in these rows):
  S1  eliminating u reproduces the master equation  -g (a2 Q'' + b2 Q' + c2 Q)  with a2, b2, c2 built from the
      stored al/be rows with ga_s = ga_x = g (X1 on-shell identity);
  S2  stored eqs_l a2 (with ga_s -> g) equals al_s/g;
  S3  light-cone residue: det(t A)|_{t=0} = 0 and tr = (1 - kappa)/lambda_lc  (h_A-level exponent, X1 table);
  S4  Nx(N, Acons, om, V) == f_N (x-field), and grad(Acons) . f == Ax(N, Acons, om, V)  (jets are on-shell);
  S5  centre chart (N = 1/P, om = W P^2, V = U P): every entry of A analytic at P = 0, and
      M0 := (theta x_theta) A at theta = 0 is the exact matrix [[-3, -2], [-2, 0]] at U* = -1/2 for every W
      (theta x_theta -> 1 is checked numerically in l2_lib), exponents 1 = l - 1 and -4 = -l - 2 (C1 PC1);
  S6  flat toy (G App. D) first-order system (f, f') from the stored eqs_l toy rows, same kappa structure.
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
B1 = os.path.normpath(os.path.join(HERE, '..', 'build-1'))
sys.path.insert(0, B1); sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp

EQS = os.path.join(B1, 'eqs_l.pkl')
kap = sp.Symbol('kappa')
L_ = sp.Symbol('L')
Nj, Aj, omj, Vj = sp.symbols('N0 A0 om0 V0')
N1j, A1j = sp.symbols('N1 A1')
Q, Q1, Q2 = sp.symbols('Q Q1 Q2')


def ld(s):
    return eval(s, {k_: getattr(sp, k_) for k_ in dir(sp)})


def split_kappa(e):
    """e polynomial of degree <= 2 in kappa -> [c0, c1, c2] (exact)"""
    e = sp.cancel(e)
    num, den = sp.fraction(e)
    P = sp.Poly(sp.expand(num), kap)
    if P.degree() > 2:
        raise ValueError('kappa degree > 2')
    return [sp.cancel(P.coeff_monomial(kap**i)/den) for i in range(3)]


def axial_symbolic(Lval=6, checks=True):
    E = pickle.load(open(EQS, 'rb'))
    ax = {k_: ld(E['ax'][k_]) for k_ in ('al_s', 'be_s', 'al_x', 'be_x', 'a2', 'ga_s', 'ga_x')}
    g = (L_ - 2)/2
    al_s, be_s, al_x, be_x = ax['al_s'], ax['be_s'], ax['al_x'], ax['be_x']
    A = sp.Matrix([[-be_s/al_s, -1/al_s],
                   [kap*al_x*be_s/al_s - kap*be_x - g, kap*al_x/al_s]])
    out = dict(g=g)
    chk = {}
    if checks:
        import lin_core as LC
        assert LC.N0.name == 'N0' and LC.A1.name == 'A1'
        # S1: u = -al_s Q1 - be_s Q ; residual u' - (A21 Q + A22 u) == -g (a2 Q2 + b2 Q1 + c2 Q)
        Qs = LC.pert('Qm')
        u = -al_s*Qs[1] - be_s*Qs[0]
        res = LC.Dx(u) - (A[1, 0]*Qs[0] + A[1, 1]*u)
        a2 = al_s/g
        b2 = LC.Dx(al_s)/g + be_s/g - kap*al_x/g
        c2 = LC.Dx(be_s)/g - kap*be_x/g - 1
        master = a2*Qs[2] + b2*Qs[1] + c2*Qs[0]
        chk['S1 u-elimination == -g (a2 Q\'\' + b2 Q\' + c2 Q), ga_s = ga_x = g'] = sp.cancel(res + g*master) == 0
        # S2
        chk['S2 stored a2 == al_s/ga_s (so a2 = al_s/g on shell)'] = sp.cancel(ax['a2']*ax['ga_s'] - al_s) == 0
        # S3: light-cone residue (N0 = 1, al_s = (N^2-1)/(2 A N^2) ~ N_x t/A at x_lc)
        dal = LC.Dx(al_s)
        R0 = sp.Matrix([[-be_s/dal, -1/dal], [kap*al_x*be_s/dal, kap*al_x/dal]]).subs(Nj, 1)
        lam = -N1j
        chk['S3a det residue == 0 at N = 1'] = sp.cancel(R0.det()) == 0
        chk['S3b tr residue == (1 - kappa)/lambda_lc at N = 1 (h_A level)'] = sp.cancel(R0.trace() - (1 - kap)/lam) == 0
        chk['S3c residue entry 12 == -A/N_x = A/lambda_lc (kappa-free, nonzero: A > 0, N_x < 0)'] = sp.cancel(R0[0, 1] + Aj/N1j) == 0
    # background substitution: N_x, A_x algebraic (system.py: eqs_gen Nx, Ax), plain symbols N, A, om, V
    import system as S
    sub = {Nj: S.N, Aj: S.A, omj: S.om, Vj: S.V, N1j: S.E['Nx'], A1j: S.E['Ax'], L_: Lval}
    Ab = A.subs(sub).subs(S.k, sp.Rational(1, 3))
    ent = {}
    for i in range(2):
        for j in range(2):
            ent[(i, j)] = split_kappa(Ab[i, j])
    out['A'] = ent                                   # (i, j) -> [kappa^0, kappa^1, kappa^2] in (N, A, om, V)
    den = S.N**2 - 1
    out['Anum'] = {ij: [sp.cancel(c*den) for c in v] for ij, v in ent.items()}   # (N^2 - 1) A: analytic at N = 1
    if checks:
        for ij, v in out['Anum'].items():
            for c in v:
                d_ = sp.fraction(sp.cancel(c))[1]
                if sp.cancel(d_.subs(S.N, 1)) == 0:
                    raise AssertionError('Anum entry still singular at N = 1')
        # S4 jets on shell
        import pipeline as PL
        f = [PL.load_expr(s) for s in PL.FIELDS['f']]
        kk = {S.k: sp.Rational(1, 3)}
        Ac = S.Acons.subs(kk)
        chk['S4a Nx(N, Acons, om, V) == f_N'] = sp.cancel(S.E['Nx'].subs(S.A, Ac).subs(kk) - f[0].subs(kk)) == 0
        dA = sum(sp.diff(Ac, v)*fi.subs(kk) for v, fi in zip((S.N, S.om, S.V), f))
        chk['S4b grad(Acons).f == Ax(N, Acons, om, V)'] = sp.cancel(dA - S.E['Ax'].subs(S.A, Ac).subs(kk)) == 0
        # S5 centre chart
        P_, W_, U_ = S.P, S.W, S.U
        csub = {S.N: 1/P_, S.om: W_*P_**2, S.V: U_*P_}
        cen = {}
        for ij, v in ent.items():
            cen[ij] = [sp.cancel(c.subs(S.A, Ac).subs(csub)) for c in v]
        regular = all(sp.fraction(c)[1].subs(P_, 0) != 0 for v in cen.values() for c in v)
        chk['S5a centre entries analytic at P = 0 (denominators nonzero at P = 0)'] = bool(regular)
        M0 = sp.Matrix(2, 2, lambda i, j: sp.cancel(cen[(i, j)][0].subs(P_, 0).subs(U_, sp.Rational(-1, 2))))
        chk['S5b M0 at P = 0, U* = -1/2 == [[-3, -2], [-2, 0]] for every W (kappa parts vanish at P = 0)'] = \
            (M0 == sp.Matrix([[-3, -2], [-2, 0]])) and all(sp.cancel(cen[ij][p].subs(P_, 0)) == 0 for ij in cen for p in (1, 2))
        chk['S5c exponents 1 and -4 (eigenvalues of M0)'] = sorted(M0.eigenvals().keys()) == [-4, 1]
        chk['S5d regular start (1, -2) in ker(M0 - 1) exactly'] = (M0 - sp.eye(2))*sp.Matrix([1, -2]) == sp.zeros(2, 1)
        out['centre_check_M0'] = str(M0)
    out['check'] = chk
    return out


def axial_structured(Lval=6, su=1):
    """well-conditioned evaluation forms (no global expansion): kappa-split entries of the scaled A (y = (Q, u/su)),
    cancelled only in the JET symbols (N0, A0, om0, V0, N1, A1); background enters by xreplace (DAG sharing).
      step : inputs (N, A, om, V), N1 -> Nx = N (A + k om - om - 2), A1 -> Ax(A, om, V)   (eqs_gen)
      anum : (N0^2 - 1) * entry, cancelled in jet symbols (regular at N0 = 1)
      cen  : N0 -> 1/P, N1 -> a/P (a = N_x/N = A + (k-1) om - 2), cancelled in (P, A0, om0, V0, a, A1, kappa), then
             A0 -> Acen(P, W, U) = Acons(1/P, W P^2, U P) (cancelled, small), om0 -> W P^2, V0 -> U P."""
    import system as S
    E = pickle.load(open(EQS, 'rb'))
    ax = {k_: ld(E['ax'][k_]) for k_ in ('al_s', 'be_s', 'al_x', 'be_x')}
    g = sp.Rational(Lval - 2, 2)
    al_s, be_s, al_x, be_x = ax['al_s'], ax['be_s'], ax['al_x'], ax['be_x']
    su = sp.Rational(su)
    A = sp.Matrix([[-be_s/al_s, -su/al_s],
                   [(kap*al_x*be_s/al_s - kap*be_x - g)/su, kap*al_x/al_s]])
    kk = sp.Rational(1, 3)
    ent = {(i, j): [e.subs(sp.Symbol('k'), kk) for e in split_kappa(A[i, j])] for i in range(2) for j in range(2)}
    anum = {ij: [sp.cancel(c*(Nj**2 - 1)) for c in v] for ij, v in ent.items()}
    P_, a_ = sp.Symbol('P'), sp.Symbol('a')
    cen = {ij: [sp.cancel(c.subs({N1j: a_/P_}).subs({Nj: 1/P_})) for c in v] for ij, v in ent.items()}
    for ij, v in cen.items():
        for c in v:
            if sp.fraction(c)[1].subs(P_, 0) == 0:
                raise AssertionError('centre form singular at P = 0: %s' % (ij,))
    Nx = S.N*(S.A + S.k*S.om - S.om - 2)
    Ax = S.E['Ax']
    stepsub = {Nj: S.N, Aj: S.A, omj: S.om, Vj: S.V, N1j: Nx.subs(S.k, kk), A1j: Ax.subs(S.k, kk)}
    step = {ij: [c.xreplace(stepsub) for c in v] for ij, v in ent.items()}
    anum_s = {ij: [c.xreplace(stepsub) for c in v] for ij, v in anum.items()}
    Acen = sp.cancel(S.Acons.subs(S.k, kk).subs({S.N: 1/S.P, S.om: S.W*S.P**2, S.V: S.U*S.P}))
    if sp.fraction(Acen)[1].subs(S.P, 0) == 0:
        raise AssertionError('Acen singular at P = 0')
    omc, Vc = S.W*S.P**2, S.U*S.P
    ac = Acen + (kk - 1)*omc - 2
    A1c = Ax.subs(S.k, kk).xreplace({S.A: Acen, S.om: omc, S.V: Vc})
    censub = {P_: S.P, Aj: Acen, omj: omc, Vj: Vc, a_: ac, A1j: A1c}
    cen_s = {ij: [c.xreplace(censub) for c in v] for ij, v in cen.items()}
    M0 = sp.Matrix(2, 2, lambda i, j: sp.cancel(cen_s[(i, j)][0].subs(S.P, 0).subs(S.U, sp.Rational(-1, 2))))
    kfree = all(sp.cancel(cen_s[ij][p].subs(S.P, 0)) == 0 for ij in cen_s for p in (1, 2))
    chk = {'C1 structured centre M0 (P = 0, U* = -1/2), scaled: [[-3, -2 su], [-2/su, 0]] for every W':
           M0 == sp.Matrix([[-3, -2*su], [-2/su, 0]]),
           'C2 structured centre kappa parts vanish at P = 0': kfree,
           'C3 structured step form == ax_sys.axial_symbolic form (scaled), exact': all(
               sp.cancel(step[(i, j)][p] - axial_symbolic(Lval, checks=False)['A'][(i, j)][p]*([[1, su], [1/su, 1]][i][j])) == 0
               for i in range(2) for j in range(2) for p in range(3))}
    return dict(step=step, anum=anum_s, cen=cen_s, check=chk)


def toy_symbolic(Lval=6):
    """flat toy: f'' a + f' b + f c = 0 (eqs_l toy, n0 = e^{-x}); y = (f, f'):  A = [[0, 1], [-c/a, -b/a]]"""
    E = pickle.load(open(EQS, 'rb'))
    a, b, c = (ld(E['toy'][k_]) for k_ in ('a', 'b', 'c'))
    n0 = sp.Symbol('n0'); N = sp.Symbol('N')
    a, b, c = (sp.cancel(e.subs(n0, N).subs(L_, Lval)) for e in (a, b, c))
    A = sp.Matrix([[0, 1], [-c/a, -b/a]])
    ent = {(i, j): split_kappa(A[i, j]) for i in range(2) for j in range(2)}
    den = N**2 - 1
    anum = {ij: [sp.cancel(cc*den) for cc in v] for ij, v in ent.items()}
    u = sp.Symbol('u')
    cen = {ij: [sp.cancel(cc.subs(N, 1/u)) for cc in v] for ij, v in ent.items()}
    M0 = sp.Matrix(2, 2, lambda i, j: cen[(i, j)][0].subs(u, 0))
    chk = {'T1 toy M0 at u = 0 == [[0, 1], [L, -1]] (exponents l = 2, -l-1 = -3)': M0 == sp.Matrix([[0, 1], [Lval, -1]]),
           'T2 kappa parts vanish at u = 0': all(sp.cancel(cen[ij][p].subs(u, 0)) == 0 for ij in cen for p in (1, 2)),
           'T3 regular start (1, 2) in ker(M0 - 2)': (M0 - 2*sp.eye(2))*sp.Matrix([1, 2]) == sp.zeros(2, 1)}
    return dict(A=ent, Anum=anum, cen=cen, N=N, u=u, check=chk)


if __name__ == '__main__':
    import json, time
    T0 = time.time()
    o = axial_symbolic()
    for k_, v in o['check'].items():
        print(k_, v)
    for ij, v in o['A'].items():
        print(ij, [sp.count_ops(c) for c in v])
    t = toy_symbolic()
    print(t['check'])
    print('time %.1fs' % (time.time() - T0))
