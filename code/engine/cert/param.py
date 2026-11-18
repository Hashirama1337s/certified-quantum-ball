"""Certified parameterisation method for a 1D invariant manifold of a polynomial
vector field G at an equilibrium y0 (dim 3), plus a scalar companion S
(the x-coordinate: lam*theta*X'(theta) = S(Y(theta))).

Invariance equation:  G(Y(theta)) = lam * theta * Y'(theta),  Y(0)=y0, Y'(0)=v.
Coefficients: (n lam - DG0) Y_n = [G(Y)]_n |_{Y_n=0},  n >= 2   (ball / Dual arithmetic).
Tail (sup norm on |theta| <= nu, theta complex): Y = Y_K + E, E = O(theta^{K+1}),
  E = Lam^{-1} [ delta + H(Y_K+E) - H(Y_K) ],  delta = G(Y_K) - lam theta Y_K',
  H = G - DG0.  Lam^{-1} bound by the Schwarz lemma:
  for f = O(theta^{K+1}) and lam theta g' - mu g = f:  |g| <= sup|f| / |(K+1) lam - mu|
  (valid when (K+1) - mu/lam > 0), applied in approximate eigen-coordinates.
"""
import math, cmath
import mpmath
from flint import arb, acb, arb_series, arb_mat, acb_mat
from core import Dual, sym, csym, amax

def val(x):
    while isinstance(x, Dual): x = x.v
    return x

def solve3(A, b):
    """generic 3x3 Cramer (works for Dual / arb / acb)"""
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
                + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]))
    d = det3(A)
    out = []
    for j in range(3):
        m = [[(b[i] if c == j else A[i][c]) for c in range(3)] for i in range(3)]
        out.append(det3(m)/d)
    return out

def jac(tape, y):
    """DG(y) via one Dual evaluation with 3 directions (y entries may themselves be Dual)"""
    z = 0*y[0]
    one = z + 1
    ins = [Dual(y[i], [one if j == i else z for j in range(3)]) for i in range(3)]
    outs = tape.eval(ins, conv=lambda c: Dual(_c(c, y[0]), [z, z, z]))
    J = []
    for o in outs[:3]:
        J.append(o.d if isinstance(o, Dual) else [z, z, z])
    return J, [o.v if isinstance(o, Dual) else o for o in outs]

def _c(c, like):
    from core import sym_const
    v = sym_const(c)
    return v + 0*like

def coefficients(G, S, y0, lam, v, K, Xzero_order=True):
    """Y_n (n<=K) and companion X_n (n<=K) in the ring of y0 (arb or Dual)."""
    z = 0*y0[0]
    DG0, g0 = jac(G, y0)
    G.jet_init(K + 1, z, conv=lambda c: _c(c, y0[0]))
    S.jet_init(K + 1, z, conv=lambda c: _c(c, y0[0]))
    Y = [[z]*(K + 1) for _ in range(3)]
    X = [z]*(K + 1)
    def put(n):
        for i in range(3):
            G.C[i][n] = Y[i][n]; S.C[i][n] = Y[i][n]
        G.jet_order(n); S.jet_order(n)
    for i in range(3): Y[i][0] = y0[i]
    put(0)
    S0 = S.C[S.outs[0]][0]
    for i in range(3): Y[i][1] = v[i]
    put(1)
    X[1] = S.C[S.outs[0]][1]/lam
    for n in range(2, K + 1):
        put(n)                       # with Y_n = 0
        g = [G.C[G.outs[i]][n] for i in range(3)]
        A = [[(n*lam if i == j else z) - DG0[i][j] for j in range(3)] for i in range(3)]
        yn = solve3(A, g)
        for i in range(3): Y[i][n] = yn[i]
        put(n)
        X[n] = S.C[S.outs[0]][n]/(n*lam)
    return dict(Y=Y, X=X, DG0=DG0, G0=g0, S0=S0)

def horner(c, t):
    r = 0*t
    for a in reversed(c): r = r*t + a
    return r

def tail_bound(G, S, Yc, Xc, lam, DG0, K, nu, N2=None, R2=None, ngrid=None, narc=64, log=print):
    import os
    ngrid = ngrid or int(os.environ.get("LB_NGRID", "8"))
    """rigorous sup bounds eps (for Y) and epsX (companion) on |theta|<=nu. Inputs arb."""
    N2 = N2 or 3*K
    R2 = R2 or 2*nu
    nu = arb(nu); R2 = arb(R2)
    # 1. defect coefficients K<n<N2 via arb_series (polynomial composition is exact up to rounding)
    ins = [arb_series(Yc[i]) for i in range(3)]
    gser = G.series_eval(ins, N2)
    sser = S.series_eval(ins, N2)
    low_check = arb(0)
    for i in range(3):
        cf = gser[i].coeffs()
        for n in range(min(K + 1, len(cf))):
            d = cf[n] - (n*lam*Yc[i][n])
            if not d.contains(0):
                raise ArithmeticError('low-order defect not containing 0 at n=%d comp %d: %s' % (n, i, d))
            low_check = low_check.max(d.abs_upper())
    d_head = [arb(0)]*3; s_head = arb(0)
    for i in range(3):
        cf = gser[i].coeffs()
        for n in range(K + 1, min(N2, len(cf))):
            d_head[i] += cf[n].abs_upper()*nu**n
    cf = sser[0].coeffs()
    for n in range(K + 1, min(N2, len(cf))):
        s_head += cf[n].abs_upper()*nu**n
    # 2. Cauchy tails on |theta| = R2
    Mg = arb(0); Ms = arb(0)
    pi = arb.pi()
    for j in range(narc):
        phi = 2*pi*j/narc
        cen = acb(R2*phi.cos(), R2*phi.sin())
        box = cen + csym(R2*pi/narc*arb('1.01'))
        yb = [horner([acb(a) for a in Yc[i]], box) for i in range(3)]
        go = G.eval(yb, conv=lambda c: acb(_c(c, arb(0))))
        so = S.eval(yb, conv=lambda c: acb(_c(c, arb(0))))
        for o in go[:3]: Mg = Mg.max(o.abs_upper())
        Ms = Ms.max(so[0].abs_upper())
    q = nu/R2
    geo = q**N2/(1 - q)
    delta = amax([d_head[i] + Mg*geo for i in range(3)])
    sdef = s_head + Ms*geo
    # 3. operator bound in approximate eigen-coordinates
    Jm = mpmath.matrix([[float(val(DG0[i][j]).mid()) for j in range(3)] for i in range(3)])
    ev, EV = mpmath.eig(Jm)
    if any(abs(mpmath.im(e)) > 1e-12 for e in ev):
        raise ArithmeticError('complex eigenvalues at equilibrium: %s' % ev)
    St = arb_mat([[arb(float(mpmath.re(EV[i, j]))) for j in range(3)] for i in range(3)])
    Si = St.inv()
    DGm = arb_mat([[val(DG0[i][j]) for j in range(3)] for i in range(3)])
    T = Si*DGm*St
    Dd = [arb(T[i, i].mid()) for i in range(3)]
    Delta = arb(0)
    for i in range(3):
        row = arb(0)
        for j in range(3):
            row += (T[i, j] - (Dd[i] if i == j else 0)).abs_upper()
        Delta = Delta.max(row)
    c = arb(0)
    for i in range(3):
        den = (K + 1)*lam - Dd[i]
        if not (den/lam > 0):
            raise ArithmeticError('Schwarz condition fails for eigenvalue %s' % Dd[i])
        c = c.max(1/den.abs_lower())
    def mnorm(A):
        m = arb(0)
        for i in range(3):
            m = m.max(sum((A[i, j].abs_upper() for j in range(3)), arb(0)))
        return m
    if not (c*Delta < 1):
        raise ArithmeticError('c*Delta >= 1')
    cp = c/(1 - c*Delta)                       # |Lam~^{-1}| in eigen-coordinates g~ = S~^{-1} E
    nS, nSi = mnorm(St), mnorm(Si)
    rowS = [sum((St[i, j].abs_upper() for j in range(3)), arb(0)) for i in range(3)]
    Sl = [[St[i, j] for j in range(3)] for i in range(3)]
    Sil = [[Si[i, j] for j in range(3)] for i in range(3)]
    dtil = nSi*delta
    # 4. Lipschitz of S~^{-1}(DG - DG0)S~ and grad S over the complex neighbourhood
    def lips(epst):
        L = arb(0); LS = arb(0)
        h = nu/ngrid
        for a in range(ngrid):
            for b in range(ngrid):
                cx = -nu + (2*a + 1)*h; cy = -nu + (2*b + 1)*h
                cell = acb(cx, cy) + csym(h*arb('1.415'))
                yb = [horner([acb(t) for t in Yc[i]], cell) + csym(rowS[i]*epst) for i in range(3)]
                Jc, _ = jac(G, yb)
                Dm = [[Jc[i][j] - val(DG0[i][j]) for j in range(3)] for i in range(3)]
                DS = [[sum((Dm[i][l]*Sl[l][j] for l in range(3)), acb(0)) for j in range(3)] for i in range(3)]
                T2 = [[sum((Sil[i][l]*DS[l][j] for l in range(3)), acb(0)) for j in range(3)] for i in range(3)]
                for i in range(3):
                    L = L.max(sum((T2[i][j].abs_upper() for j in range(3)), arb(0)))
                Js, _ = jac(S, yb)
                LS = LS.max(sum((Js[0][j].abs_upper() for j in range(3)), arb(0)))
        return L, LS
    epst = 2*cp*dtil + arb(10)**(-80)
    for it in range(5):
        L, LS = lips(epst)
        if not (cp*L < 1):
            raise ArithmeticError('no contraction: cp=%s L*=%s c=%s normS=%s normSi=%s Delta=%s delta=%s' % (cp.str(4), L.str(4), c.str(4), nS.str(4), nSi.str(4), Delta.str(4), delta.str(4)))
        need = cp*dtil/(1 - cp*L)
        if need <= epst:
            break
        epst = 2*need
    else:
        raise ArithmeticError('tail iteration did not settle')
    assert cp*(dtil + L*epst) <= epst          # self-map of the ball |g~| <= epst
    eps = nS*epst                              # sup |E| (inf-norm, original coordinates)
    beta = cp
    epsX = (sdef + LS*eps)/((K + 1)*lam.abs_lower())
    info = dict(delta=delta, beta=beta, L=L, LS=LS, low_defect=low_check, Mg=Mg, eig=ev, c=c, Delta=Delta)
    return arb(eps.upper()), arb(epsX.upper()), info
