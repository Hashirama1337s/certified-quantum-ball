"""Polar first-order systems (build-1).  From the derived linearised Einstein rows (lin_core), RW-type gauge:
  l >= 2 : h_Aa = 0, G = 0 (Regge-Wheeler); unknowns Hss, Hsx, Hxx, K; the traceless row is algebraic
           (Hss = -2 Hsx + (N^2 - 1) Hxx) and is used to eliminate Hss.
  l = 1  : h_Aa = 0 and K = 0 (our 'k = 0' gauge; the angular traceless harmonic vanishes for l = 1).
Fluid (dl = drho/rho, v = du_x/R, w = du_theta/R) enters the Einstein rows algebraically; the three
fluid-free combinations Z1 (ss,sx,xx), Z2 (s theta, x theta), Z3 (ss,sx,trace) are formed by exact
cross products (no division).  The fluid is then an algebraic function of the metric jets and the fluid
equations follow from the contracted Bianchi identity, so {Z1, Z2, Z3} is the complete mode system.
First-order form  C(x) eta' = (D0 + kappa D1) eta  (rows listed in build()).
"""
import sys
sys.path.insert(0, '<LAB>/sideways/build-1')
import sympy as sp
import lin_core as LC
from lin_core import R, E, Y0, Y1, sn


def einstein_rows(lsec):
    """lsec = 'l1' or 'lge2'. Returns dict of mode rows (sympy), perturbation jets."""
    B = LC.Lin()
    names = ['Hss', 'Hsx', 'Hxx', 'K', 'dl', 'v', 'w']
    P = {n: LC.pert(n) for n in names}
    h = sp.zeros(4, 4)
    h[0, 0] = R**2*E*P['Hss'][0]*Y0; h[0, 1] = h[1, 0] = R**2*E*P['Hsx'][0]*Y0; h[1, 1] = R**2*E*P['Hxx'][0]*Y0
    if lsec == 'lge2':
        h[2, 2] = R**2*E*P['K'][0]*Y0; h[3, 3] = R**2*E*P['K'][0]*Y0*sn**2
    du = [0, R*E*P['v'][0]*Y0, R*E*P['w'][0]*Y1, 0]
    du[0] = B.norm_du_s(h, du)
    drho = B.rho*E*P['dl'][0]*Y0
    comps = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2), (3, 3)]
    ein = B.einstein(h, drho, du, comps)
    rows = {}
    for key in [(0, 0), (0, 1), (1, 1)]:
        rows[key] = LC.finalize(ein[key]/Y0)
    for key in [(0, 2), (1, 2)]:
        rows[key] = LC.finalize(ein[key]/Y1)
    tt = LC.finalize(ein[(2, 2)]); pp = LC.finalize(ein[(3, 3)]/sn**2)
    tr = sp.cancel((tt + pp)/2); tl = sp.cancel((tt - pp)/2)
    rows['tr'] = sp.cancel(tr.subs({Y1: 0}).subs(Y0, 1))
    rows['tl'] = sp.cancel(sp.diff(tl, Y1)/(-sp.Rational(3, 4)))     # (cs, sn) = (3/5, 4/5) in finalize
    if lsec == 'l1':
        rows = {kk: v.subs(LC.L, 2) for kk, v in rows.items()}
    return rows, P


def cross3(ra, rb, rc, p, q):
    """combination of rows ra, rb, rc free of the two symbols p, q (exact, polynomial weights)"""
    pa, pb, pc = [sp.diff(r_, p) for r_ in (ra, rb, rc)]
    qa, qb, qc = [sp.diff(r_, q) for r_ in (ra, rb, rc)]
    return ra*(pb*qc - pc*qb) - rb*(pa*qc - pc*qa) + rc*(pa*qb - pb*qa)


def build(lsec):
    rows, P = einstein_rows(lsec)
    dl, v, w = P['dl'][0], P['v'][0], P['w'][0]
    Z1 = cross3(rows[(0, 0)], rows[(0, 1)], rows[(1, 1)], dl, v)
    Z2 = rows[(0, 2)]*sp.diff(rows[(1, 2)], w) - rows[(1, 2)]*sp.diff(rows[(0, 2)], w)
    Z3 = cross3(rows[(0, 0)], rows[(0, 1)], rows['tr'], dl, v)
    for z in (Z1, Z2, Z3):
        for f in (dl, v, w):
            assert sp.cancel(sp.diff(z, f)) == 0
    Hss, Hsx, Hxx, K = P['Hss'], P['Hsx'], P['Hxx'], P['K']
    eta = sp.symbols('e0:6')
    if lsec == 'lge2':
        # traceless: Hss = -2 Hsx + (N^2 - 1) Hxx  (checked against the derived row)
        tl = rows['tl']
        chk = sp.cancel(tl.subs({Hss[0]: -2*Hsx[0] + (LC.N0**2 - 1)*Hxx[0]}))
        assert chk == 0, chk
        HssE = [-2*Hsx[0] + (LC.N0**2 - 1)*Hxx[0]]
        for i in range(2):
            HssE.append(LC.Dx(HssE[-1]))
        sub = {Hss[i]: HssE[i] for i in range(3)}
        Z = [sp.expand(sp.numer(sp.together(z.subs(sub)))) for z in (Z1, Z2, Z3)]
        # eta = (K, K', Hsx, Hxx, Hxx') ; Hsx' from Z2 (first order)
        e = eta[:5]
        for z, nm in zip(Z, ('Z1', 'Z2', 'Z3')):
            print(nm, 'second derivatives:', [str(s_) for s_ in (Hsx[2], Hxx[2], K[2]) if sp.diff(z, s_) != 0],
                  'first:', [str(s_) for s_ in (Hsx[1], Hxx[1], K[1]) if sp.diff(z, s_) != 0])
        z2 = Z[1]
        assert all(sp.diff(z2, s_) == 0 for s_ in (Hsx[2], Hxx[2], K[2]))
        a2 = sp.diff(z2, Hsx[1])
        S = sp.cancel(-(z2 - a2*Hsx[1]).subs({K[1]: e[1], Hxx[1]: e[4], K[0]: e[0], Hsx[0]: e[2], Hxx[0]: e[3]})/a2)
        # S = Hsx' as a function of eta (algebraic); jets of S along the solution need eta' -> symbols
        ep = sp.symbols('ep0:5')
        def dtot(f):
            # total x-derivative of f(eta, background) along a solution: background via Dx, eta via ep
            out = LC.Dx(f)
            for i in range(5):
                out += sp.diff(f, e[i])*ep[i]
            return out
        Sp = dtot(S)
        mapZ = lambda z: z.subs({Hsx[2]: Sp}).subs({K[2]: ep[1], Hxx[2]: ep[4], Hsx[1]: S, Hxx[1]: e[4], K[1]: e[1],
                                                      K[0]: e[0], Hsx[0]: e[2], Hxx[0]: e[3]})
        rowsF = [ep[0] - e[1], ep[2] - S, ep[3] - e[4], mapZ(Z[0]), mapZ(Z[2])]
        n = 5
    else:
        Z = [sp.expand(z) for z in (Z1, Z2, Z3)]
        for z, nm in zip(Z, ('Z1', 'Z2', 'Z3')):
            print(nm, 'second derivatives:', [str(s_) for s_ in (Hss[2], Hsx[2], Hxx[2]) if sp.diff(z, s_) != 0],
                  'first:', [str(s_) for s_ in (Hss[1], Hsx[1], Hxx[1]) if sp.diff(z, s_) != 0])
        # eta = (Hss, Hss', Hsx, Hxx)
        e = eta[:4]; ep = sp.symbols('ep0:4')
        mp = {Hss[2]: ep[1], Hss[1]: e[1], Hsx[1]: ep[2], Hxx[1]: ep[3], Hss[0]: e[0], Hsx[0]: e[2], Hxx[0]: e[3]}
        for z in Z:
            assert all(sp.diff(z, s_) == 0 for s_ in (Hsx[2], Hxx[2]))
        rowsF = [ep[0] - e[1]] + [z.subs(mp) for z in Z]
        n = 4
    C = sp.Matrix([[sp.diff(r_, ep[j]) for j in range(n)] for r_ in rowsF])
    rest = [r_ - sum(C[i, j]*ep[j] for j in range(n)) for i, r_ in enumerate(rowsF)]
    D1 = sp.Matrix([[-sp.diff(sp.diff(r_, e[j]), LC.kap) for j in range(n)] for r_ in rest])
    Dfull = sp.Matrix([[-sp.diff(r_, e[j]) for j in range(n)] for r_ in rest])
    D0 = Dfull.subs(LC.kap, 0)
    # linearity checks
    import random
    rnd = random.Random(7)
    for i, r_ in enumerate(rest):
        lin = sum(sp.diff(r_, e[j])*e[j] for j in range(n))
        fs = sorted((r_ - lin).free_symbols, key=str)
        val_ = (r_ - lin).xreplace({s_: sp.Rational(rnd.randint(2, 40), rnd.randint(41, 90)) for s_ in fs})
        assert sp.nsimplify(val_) == 0, ('row not linear in eta', i)
    return dict(C=C, D0=D0, D1=D1, Dfull=Dfull, n=n, eta=e, Z=Z, P=P)

