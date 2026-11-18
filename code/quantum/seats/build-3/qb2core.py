"""build-3 QB2 instrument, part 2: forcing columns in every chart, the augmented regular-singular series
(sonic, centre) and helpers.  Augmented state (NA = 12):
   0-3  Y = (Abar_p, Nbar_p, wbar_p, V_p)            HKA bar variables, rate kappa = 2
   4    J      scalar check variable for null N-a:   J' = -A J + (2 - A) S_0 z0     (Lemma (a), Consequence)
   5    Phi    C8 auxiliary, beta_nn^{C8} = A^{-1/2} Phi:   Phi' = -2 Phi + A^{1/2} g z8   (D_t = 0 <=> this)
   6 zR  real source (b = Bach of the background, beta = e^{-2x} b, four HKA rows)
   7 z0  row-0 part of the real source only (null N-a: row 0 x (1 + delta))
   8 zC  null N-c: b rebuilt with om -> om (1 + 1e-6) in A (Bach inputs only)
   9 z1  control C1: Q^{-1} S = phi' G_2,  phi = e^{2x}/(1+e^{2x}) - 1/2, G_2 = (Abar', Nbar' + 2, wbar', V')
   10 z8 control C8: explicit beta_nr, beta_rr (+ drives Phi, which enters Y through beta_nn = A^{-1/2} Phi)
   11 zE explicit beta_nr, beta_rr part of C8 alone (for C8' = (1+d) C8 - d E)
All forcing columns are in the x-form dY/dx = M Y + F; charts multiply by dx/dt (tube tau steps), theta x'
(sonic) or (1 + theta chi') (centre).
"""
import sys, math
sys.path.insert(0, '<LAB>/quantum/seats/build-3')
import qb2lib as Q
from qb2lib import s_mul, s_div, s_add, s_sub, s_scal, s_up, s_deriv, s_exp, s_sqrt, s_trunc, zlike
from flint import arb, acb, arb_mat, acb_mat, ctx
from core import sym, csym
Q.load_rows()

NA = 12
IJ, IP = 4, 5
ZC = dict(R=6, N0=7, C=8, C1=9, C8=10, E=11)
DELTA_NC = arb('1e-6')
EXTRA = list(range(4, NA))


def ser(x):
    m, e = x.mid().man_exp(); rm, re = x.rad().man_exp()
    return [str(int(m)), int(e), str(int(rm)), int(re)]


def des(t):
    return arb(arb(int(t[0]))*arb(2)**t[1], arb(int(t[2]))*arb(2)**t[3])


def one_like(z, n):
    return [z + 1] + [z]*(n - 1)


def apply_beta(R, sI, beta, n):
    """R[i][c] per-unit-beta row series (rows 0,1 final; rows 2,3 numerators), beta triple -> 4 rows, sigma_I"""
    z = zlike(R[0][0][0])
    rows = []
    for i in range(4):
        acc = [z]*n
        for c in range(3):
            if beta[c] is None:
                continue
            acc = s_add(acc, s_mul(R[i][c], beta[c], n))
        rows.append(acc)
    si = [z]*n
    for c in range(3):
        if beta[c] is None:
            continue
        si = s_add(si, s_mul(sI[c], beta[c], n))
    return rows, si


def generic_columns(Ns, oms, Vs, xp, ex, n_in, bach_override=None):
    """forcing columns (x-form) in the generic chart, series in the local parameter (length n = n_in - 4).
    Ns, oms, Vs: background series; xp = dx/dt; ex = e^{x} series.  Rows 2,3 of the beta columns are returned as
    NUMERATORS (adj(Mx).(S_2,S_3)); 'det' = det Mx; direct columns (C1) carry final rows 2,3."""
    n = n_in - 4
    z = zlike(Ns[0])
    if bach_override is None:
        bt = Q.bach_from_series(Ns, oms, Vs, xp, n_in)
        bw = Q.bach_from_series(Ns, oms, Vs, xp, n_in, om_fac=DELTA_NC)
    else:
        bt, bw = bach_override
    A = bt['A']; A1 = bt['A1']; l1 = bt['l1']
    N = s_trunc(Ns, n, z); V = s_trunc(Vs, n, z)
    e2x = s_div(one_like(z, n), s_mul(ex, ex, n), n)
    R, det, sI = Q.rows_coeff(N, V, A, n)
    out = dict(n=n, det=det, A=A, b=bt['b'], bw=bw['b'])
    bR = [s_mul(e2x, bt['b'][c], n) for c in range(3)]
    bC = [s_mul(e2x, bw['b'][c], n) for c in range(3)]
    out['R'] = ('beta',) + apply_beta(R, sI, bR, n)
    rows0 = [out['R'][1][0], [z]*n, [z]*n, [z]*n]
    out['N0'] = ('direct', rows0, [z]*n)
    out['C'] = ('beta',) + apply_beta(R, sI, bC, n)
    # C1: phi' G_2 (direct; sigma = 0)
    Dx = lambda s, m: s_div(s_deriv(s), xp, m)
    Nx = Dx(Ns, n); omx = Dx(oms, n); Vx = Dx(Vs, n)
    G2 = [s_div(A1, A, n), s_add(s_div(Nx, N, n), [z + 2] + [z]*(n - 1)), s_div(omx, s_trunc(oms, n, z), n), Vx]
    php = Q.phi_prime(ex, n)
    out['C1'] = ('direct', [s_mul(php, G2[i], n) for i in range(4)], [z]*n)
    # C8 explicit profiles
    bnr, brr = Q.c8_profiles(ex, n + 1)
    out['E'] = ('beta',) + apply_beta(R, sI, [None, bnr[:n], brr[:n]], n)
    out['C8'] = out['E']
    # coupling Phi -> Y : beta_nn = A^{-1/2} Phi
    rA = s_sqrt(A, n); irA = s_div(one_like(z, n), rA, n)
    out['P'] = ('beta',) + apply_beta(R, sI, [irA, None, None], n)
    # g = N beta_nr' + (2N + 2N' + N A'/A) beta_nr - (A'/(2A)) beta_rr,  N' = N l1
    dbnr = Dx(bnr, n)
    AA = s_div(A1, A, n)
    coef = s_mul(N, s_add(s_add([z + 2] + [z]*(n - 1), s_scal(l1, 2)), AA), n)
    g = s_sub(s_add(s_mul(N, dbnr, n), s_mul(coef, bnr[:n], n)), s_scal(s_mul(AA, brr[:n], n), arb(1)/2))
    out['Pdrv'] = s_mul(rA, g, n)
    out['Jdrv'] = s_mul(s_sub([z + 2] + [z]*(n - 1), A), out['R'][1][0], n)
    out['Jdiag'] = s_scal(A, -1)
    out['Pdiag'] = [z - 2] + [z]*(n - 1)
    # conservation D_t of the real source (kappa = 2; a3.py structure):
    # D_t = -bnn' + N bnr' - (2 + A'/(2A)) bnn + (2N + 2N' + N A'/A) bnr - (A'/(2A)) brr   on beta = e^{-2x} b
    m = n - 1
    Dt = s_add(s_add(s_scal(Dx(bR[0], m), -1), s_mul(N, Dx(bR[1], m), m)),
               s_sub(s_mul(coef, bR[1], m), s_add(s_mul(s_add([z + 2] + [z]*(m - 1), s_scal(AA, arb(1)/2)), bR[0], m),
                                                  s_scal(s_mul(AA, bR[2], m), arb(1)/2))))
    out['Dt'] = Dt
    out['trace'] = s_sub(s_sub(bt['b'][0], bt['b'][2]), s_scal(bt['b'][3], 2))
    return out


COLS = [('P', IP), ('R', 6), ('N0', 7), ('C', 8), ('C1', 9), ('C8', 10), ('E', 11)]


def x_form_matrix_entries(cols, n, rows23):
    """dict (i, j) -> series of the x-form augmented matrix entries (extras only).  rows23(num, det) -> final."""
    z = zlike(cols['A'][0])
    E = {}
    for name, j in COLS:
        kind, rows, si = cols[name]
        for i in range(4):
            if kind == 'beta' and i >= 2:
                E[(i, j)] = rows23(rows[i], cols['det'])
            else:
                E[(i, j)] = rows[i] if kind == 'direct' or i < 2 else None
    E[(IJ, IJ)] = cols['Jdiag']; E[(IJ, 7)] = cols['Jdrv']
    E[(IP, IP)] = cols['Pdiag']; E[(IP, 10)] = cols['Pdrv']
    return E


def sigma_I_of(cols, name):
    kind, rows, si = cols[name]
    return si


# ------------------------------------------------------------------ local series at a box (polynomial + tail)
def local_series(coefs, t0, L, tail_eps, rad_avail):
    """Taylor coefficients 0..L-1 about t0 of p(t) + E(t), p = sum coefs t^n, |E| <= tail_eps on |t| <= R
    with |t0| + rad_avail <= R: |E^{(j)}(t0)/j!| <= tail_eps / rad_avail^j (Cauchy)."""
    from math import comb
    out = []
    for j in range(L):
        acc = 0*t0
        for k in range(len(coefs) - 1, j - 1, -1):
            acc = acc*t0 + coefs[k]*comb(k, j)
        tb = tail_eps/rad_avail**j
        if isinstance(t0, acb):
            acc = acc + csym(tb)
        else:
            acc = acc + sym(tb)
        out.append(acc)
    return out


# ------------------------------------------------------------------ augmented regular-singular series
class RegSingAug:
    """t dY/dt = M(t) Y, M analytic on |t| <= nu, Taylor matrices M[n] (NA x NA arb_mat), Snorm >= sup_{|t|=nu}
    ||M(t) - M[0]||_inf.  cn(n) >= ||(n - M[0])^{-1}||_inf for n >= n_start.  Same majorant tail as
    engine lin_stages.RegSing.evaluate."""
    def __init__(self, M, Snorm, nu, cn):
        self.M = M; self.S = Snorm; self.nu = arb(nu); self.cn = cn; self.K = len(M) - 1

    def series(self, y0, fixed=None):
        K = self.K; n_ = self.M[0].nrows()
        I = arb_mat(n_, n_)
        for i in range(n_): I[i, i] = 1
        y = [y0]
        for n in range(1, K + 1):
            r = arb_mat(n_, 1)
            for m in range(1, n + 1):
                r += self.M[m]*y[n - m]
            y.append((I*n - self.M[0]).solve(r))
        return y

    def evaluate(self, y, t):
        K = self.K; n_ = y[0].nrows()
        at = abs(t).upper() if not isinstance(t, arb) else t.abs_upper()
        nup = ((self.nu + at)/2).min(2*at)
        rat = nup/self.nu
        mup = self.S*rat/(1 - rat)
        n0 = K
        while True:
            c = self.cn(n0 + 1)
            if c is not None and c*mup <= 1: break
            n0 += 1
            if n0 > 100000: raise ArithmeticError('tail n0')
        if n0 > K: raise ArithmeticError('RegSingAug tail needs K >= %d (have %d)' % (n0, K))
        B = arb(0)
        for n in range(K + 1):
            B = B.max(max(abs(y[n][i, 0]).upper() for i in range(n_))*nup**n)
        qq = at/nup
        tail = B*qq**(K + 1)/(1 - qq)
        out = []
        for i in range(n_):
            s = 0*t
            for n in range(K, -1, -1): s = s*t + y[n][i, 0]
            out.append(s + (sym(tail) if isinstance(t, arb) else csym(tail)))
        return out, dict(tail=tail, n0=n0, mup=mup)


def pow2_le(x):
    """2^-k with 2^k >= x (x > 0 float), exact arb"""
    k = max(0, math.ceil(math.log2(max(float(x), 1e-300))))
    return arb(2)**(-k)


def make_scaled(Mh, dev, baseS, order_rows=(IJ, IP)):
    """diagonal scaling D = diag(1,1,1,1, lam_4..lam_11) (powers of 2, exact).  Scaled M~ = D^-1 M D.
    dev[(i, j)] >= sup over the circle of |M_ij(t) - M_ij(0)| for the extras entries.
    Returns (Ms, Snorm, lam) with Snorm >= sup ||M~(t) - M~(0)||_inf (baseS covers the 4x4 block)."""
    lam = [arb(1)]*NA
    # rows J, P first (their own columns), then the z columns
    for j in (IP,):
        m = max([float(dev.get((i, j), arb(0)).upper()) for i in range(4)] + [0.0])
        lam[j] = pow2_le(8*m)
    for j in range(6, NA):
        m = max([float((dev.get((i, j), arb(0))/lam[i]).upper()) for i in range(NA) if (i, j) in dev] + [0.0])
        lam[j] = pow2_le(8*m)
    Ms = []
    for M in Mh:
        S_ = arb_mat(NA, NA)
        for i in range(NA):
            for j in range(NA):
                if not M[i, j].is_zero():
                    S_[i, j] = M[i, j]*lam[j]/lam[i]
        Ms.append(S_)
    rows = {}
    for (i, j), d in dev.items():
        rows[i] = rows.get(i, arb(0)) + d*lam[j]/lam[i]
    Snorm = baseS + max([rows.get(i, arb(0)) for i in range(4)])
    for i in range(4, NA):
        if i in rows:
            Snorm = Snorm.max(rows[i])
    return Ms, Snorm, lam


def rowsum_mat(M):
    m = arb(0)
    for i in range(M.nrows()):
        m = m.max(sum((M[i, j].abs_upper() for j in range(M.ncols())), arb(0)))
    return m
