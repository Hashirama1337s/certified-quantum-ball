"""build-3 QB2 instrument, part 1: series arithmetic, Tapes, the source b_ab and the forcing columns.

Everything is python-flint ball arithmetic (arb / acb).  No flint object is pickled.
Series = python list of balls (Taylor coefficients in the local parameter t), truncated to a stated length.

b_ab (node N0):  b = Bach_(ab)[gt] (bachsym.py, exact), as a Tape in (u = 1/N, A, A', A'', A''', A'''',
l', l'', l''', l''''),  l = ln N, primes = d/dx.  Along any parameter t with x' = dx/dt the x-derivatives are
series operations D = (1/x') d/dt (one order lost each), so b is obtained as a certified Taylor series in t
from the certified state series (N, om, V)(t) of the background record.  On complex/real boxes the x-jets come
from the background ODE y' = f(y) (engine/system.py rhs_regular, the record's own field).
"""
import sys, pickle, math
from flint import ctx
ctx.cap = 32                                     # = the sealed count stack's cap (count_b4v3.py)
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, acb
from core import Tape, sym_const, sym, csym
import system as SY

HERE = '<LAB>/quantum/seats/build-3/'
KV = sp.Rational(1, 3)


def load_expr(s):
    return eval(s, {k: getattr(sp, k) for k in dir(sp)})


# ------------------------------------------------------------------ series (generic ring)
def zlike(x):
    return 0*x


def s_trunc(a, n, z):
    a = list(a[:n]); return a + [z]*(n - len(a))


def s_mul(a, b, n):
    z = zlike(a[0])
    out = []
    for k in range(n):
        s = z
        for i in range(max(0, k - len(b) + 1), min(k, len(a) - 1) + 1):
            s = s + a[i]*b[k - i]
        out.append(s)
    return out


def s_div(a, b, n):
    z = zlike(a[0])
    q = []
    for k in range(n):
        s = a[k] if k < len(a) else z
        for i in range(max(0, k - len(b) + 1), k):
            s = s - q[i]*b[k - i]
        q.append(s/b[0])
    return q


def s_deriv(a):
    return [(k + 1)*a[k + 1] for k in range(len(a) - 1)]


def s_scal(a, c):
    return [c*x for x in a]


def s_add(a, b):
    n = min(len(a), len(b)); return [a[i] + b[i] for i in range(n)]


def s_sub(a, b):
    n = min(len(a), len(b)); return [a[i] - b[i] for i in range(n)]


def s_exp(a, n):
    """exp of a series: e' = a' e"""
    z = zlike(a[0]); e = [a[0].exp()]
    da = s_deriv(s_trunc(a, n + 1, z))
    for k in range(1, n):
        s = z
        for i in range(k):
            s = s + da[i]*e[k - 1 - i]
        e.append(s/k)
    return e


def s_sqrt(a, n):
    """sqrt of a series with a[0] > 0: s^2 = a"""
    z = zlike(a[0]); s = [a[0].sqrt()]
    for k in range(1, n):
        acc = a[k] if k < len(a) else z
        for i in range(1, k):
            acc = acc - s[i]*s[k - i]
        s.append(acc/(2*s[0]))
    return s


def s_up(a, m=1):
    """multiply by t^m"""
    z = zlike(a[0]); return [z]*m + list(a)


def s_shift_down(a, m=1):
    """divide by t^m (caller has checked/zeroed the first m coefficients)"""
    return list(a[m:])


def horner(c, t):
    r = 0*t
    for x in reversed(c):
        r = r*t + x
    return r


def series_out(tape, ins, n, z, conv):
    tape.jet_init(n, z, conv=conv)
    for k in range(n):
        for i in range(len(ins)):
            tape.C[i][k] = ins[i][k] if k < len(ins[i]) else z
        tape.jet_order(k)
    return [[tape.C[o][k] for k in range(n)] for o in tape.outs]


def conv_for(z):
    if isinstance(z, acb):
        return lambda c: acb(sym_const(c))
    return lambda c: sym_const(c)


# ------------------------------------------------------------------ tapes
Nn, omn, Vn = SY.N, SY.om, SY.V
A_s = sp.Symbol('A')
u_s = sp.Symbol('u'); q_s = sp.Symbol('q')
Asym = sp.symbols('A0:6'); lsym = sp.symbols('l0:6')
P_s, W_s, U_s = sp.symbols('P W U')

T_ACONS = Tape([sp.cancel(SY.Acons.subs(SY.k, KV))], [Nn, omn, Vn], {})
T_ACEN = Tape([sp.cancel(sp.together(SY.Acons.subs(SY.k, KV).subs({Nn: 1/P_s, omn: W_s*P_s**2, Vn: U_s*P_s})))],
              [P_s, W_s, U_s], {})
T_F = Tape([sp.cancel(sp.together(e.subs(SY.k, KV))) for e in SY.rhs_regular()], [Nn, omn, Vn], {})
T_BACH = Tape([load_expr(l) for l in open(HERE + 'bach_u.txt')], [u_s] + list(Asym[:5]) + list(lsym[1:5]), {})
T_ROWS = T_CROWS = None


def load_rows():
    """Tapes of the source rows (node N1 output, through qb2sym2.py -> qb2_exprs.pkl)"""
    global T_ROWS, T_CROWS, EX
    EX = pickle.load(open(HERE + 'qb2_exprs.pkl', 'rb'))
    _gr = [[load_expr(e) for e in r] for r in EX['gen_rows']]
    T_ROWS = Tape([e for r in _gr for e in r] + [load_expr(EX['gen_det'])] + [load_expr(e) for e in EX['gen_sI']],
                  [Nn, Vn, A_s], {})
    _cr = EX['cen_rows']
    T_CROWS = Tape([load_expr(_cr[(i, j)]) for i in range(4) for j in range(3)] + [load_expr(e) for e in EX['cen_sI']],
                   [u_s, q_s, A_s], {})


# ------------------------------------------------------------------ b from parameter series
def D_op(s, xp, n):
    """x-derivative of a series in t: (ds/dt)/x'(t), length n"""
    return s_div(s_deriv(s), xp, n)


def bach_from_series(Nser, omser, Vser, xp, n_in, om_fac=None, lN=None, user=None, mv=False):
    """Taylor series (length n_in - 4) of b = (b_nn, b_nr, b_rr, b_T) and of A, A', u, l' etc.
    Nser, omser, Vser: series of length >= n_in; xp: dx/dt series.  om_fac: N-c null, A built from om*(1+om_fac)
    (only the b-profile, i.e. the Bach inputs, use it).  lN: optional precomputed l' = (ln N)' series (centre chart);
    user: optional u = 1/N series (centre chart)."""
    z = zlike(Nser[0]); cv = conv_for(z)
    om_b = omser if om_fac is None else s_scal(omser, 1 + om_fac)
    A = series_out(T_ACONS, [Nser, om_b, Vser], n_in, z, cv)[0]
    Aj = [A]
    for j in range(4):
        Aj.append(D_op(Aj[-1], xp, n_in - 1 - j))
    if lN is None:
        l1 = s_div(D_op(Nser, xp, n_in - 1), Nser, n_in - 1)
    else:
        l1 = lN[:n_in - 1]
    lj = [l1]
    for j in range(3):
        lj.append(D_op(lj[-1], xp, n_in - 2 - j))
    m = n_in - 4
    u = s_div([z + 1], Nser, m) if user is None else user[:m]
    ins = [u] + [s_trunc(a, m, z) for a in Aj] + [s_trunc(a, m, z) for a in lj]
    b = bach_series_mv(ins, m) if mv else series_out(T_BACH, ins, m, z, cv)
    return dict(b=b, A=s_trunc(A, m, z), A1=s_trunc(Aj[1], m, z), l1=s_trunc(l1, m, z), n=m)


def rows_from_b(Nser, Vser, Aser, beta, n):
    """generic chart: (F0, F1, num2, num3, detMx, sigma_I) series for the beta-triple (length n)"""
    z = zlike(Nser[0]); cv = conv_for(z)
    R = series_out(T_ROWS, [Nser, Vser, Aser], n, z, cv)
    out = []
    for i in range(4):
        acc = [z]*n
        for c in range(3):
            acc = s_add(acc, s_mul(R[3*i + c], beta[c], n))
        out.append(acc)
    det = R[12]
    sI = [z]*n
    for c in range(3):
        sI = s_add(sI, s_mul(R[13 + c], beta[c], n))
    return out[0], out[1], out[2], out[3], det, sI


def rows_coeff(Nser, Vser, Aser, n):
    """generic chart: raw per-unit-beta row series R[i][c] (i = 0..3 numerators), det, sI[c]"""
    z = zlike(Nser[0]); cv = conv_for(z)
    R = series_out(T_ROWS, [Nser, Vser, Aser], n, z, cv)
    return [[R[3*i + c] for c in range(3)] for i in range(4)], R[12], [R[13 + c] for c in range(3)]


def ode_xseries(y, n):
    """local x-Taylor coefficients (length n) of (N, om, V) at the state y (ball/box) from y' = f(y)"""
    z = zlike(y[0]); cv = conv_for(z)
    T_F.jet_init(n, z, conv=cv)
    Y = [[z]*(n + 1) for _ in range(3)]
    for i in range(3):
        Y[i][0] = y[i]
    for k in range(n - 1):
        for i in range(3):
            T_F.C[i][k] = Y[i][k]
        T_F.jet_order(k)
        for i in range(3):
            Y[i][k + 1] = T_F.C[T_F.outs[i]][k]/(k + 1)
    return [Y[i][:n] for i in range(3)]


# explicit C8 profiles and phi, as functions of the series e^{x}
def c8_profiles(ex, n):
    """beta_nr = e^{3x}/(1+e^{2x})^2, beta_rr = e^{2x}/(1+e^{2x})^2 from the series of e^x"""
    e2 = s_mul(ex, ex, n); e3 = s_mul(e2, ex, n)
    one = [zlike(ex[0]) + 1] + [zlike(ex[0])]*(n - 1)
    d = s_add(one, e2); d2 = s_mul(d, d, n)
    return s_div(e3, d2, n), s_div(e2, d2, n)


def phi_prime(ex, n):
    """phi = e^{2x}/(1+e^{2x}) - 1/2,  phi' = 2 e^{2x}/(1+e^{2x})^2"""
    e2 = s_mul(ex, ex, n)
    one = [zlike(ex[0]) + 1] + [zlike(ex[0])]*(n - 1)
    d = s_add(one, e2)
    return s_scal(s_div(e2, s_mul(d, d, n), n), 2)


def phi_val(x):
    e2 = (2*x).exp()
    return e2/(1 + e2) - arb(1)/2


def bach_series_mv(ins, m):
    """mean-value form of the Bach tape on series inputs with ball coefficients (10 inputs):
    b_n in b_n(mid inputs) + [-R_n, R_n],  R_n = sum_k sum_{j<=n} |(d_k b)_{n-j}|(over the balls) r_{k,j},
    d_k b evaluated along the ball inputs as Dual series (seed at order 0)."""
    from core import Dual
    z = arb(0)
    nk = len(ins)
    mid = [[arb(c.mid()) for c in s[:m]] for s in ins]
    rad = [[c.rad() for c in s[:m]] for s in ins]
    bm = series_out(T_BACH, mid, m, z, conv_for(z))
    zd = Dual(z, [z]*nk)
    dins = []
    for k in range(nk):
        row = [Dual(ins[k][0], [arb(1) if j == k else z for j in range(nk)])]
        row += [Dual(ins[k][i], [z]*nk) for i in range(1, m)]
        dins.append(row)
    T_BACH.jet_init(m, zd, conv=lambda c: Dual(sym_const(c), [z]*nk))
    for n in range(m):
        for i in range(nk):
            T_BACH.C[i][n] = dins[i][n]
        T_BACH.jet_order(n)
    bd = [[T_BACH.C[o][n] for n in range(m)] for o in T_BACH.outs]
    out = []
    for c in range(4):
        row = []
        for n in range(m):
            R = arb(0)
            for k in range(nk):
                for j in range(n + 1):
                    if not rad[k][j].is_zero():
                        R += bd[c][n - j].d[k].abs_upper()*rad[k][j]
            row.append(bm[c][n] + sym(R))
        out.append(row)
    return out
