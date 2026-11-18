"""build-1 QB1: rigorous curvature on truncated multivariate Taylor jets in python-flint arb ball arithmetic.

Port of seats/smoke-1/jetcurv.py (float/mpmath) to arb.  Same formulas, same conventions (design-1 vcheck.py,
checked on Hayward and Reissner-Nordstrom): mostly plus, R^a_bcd = d_c Gam^a_db - d_d Gam^a_cb + Gam^a_ce Gam^e_db
- Gam^a_de Gam^e_cb, Ric_bd = R^a_bad.  Zahn (2509.10360v3) eq. 11-12 tensors T(1), T(2); V(xi) eq. 19;
Bach = -(T2 - T1/3)/4 (design-1 vcheck B); Route A (R = 0 on-shell form, DESIGN-QB1 sec. 1.2).

A Jet is a truncated Taylor polynomial in NV coordinate offsets (coordinates with index >= NV are Killing:
derivative 0), total degree <= DEG, with arb coefficients.  Every coefficient is a rigorous ball: the jet of the
product is the exact truncated product of the balls.  After k derivatives a jet is exact to degree DEG - k, so
Bach (4 metric derivatives) is exact at degree 0 for DEG = 4 and to degree 1 for DEG = 5.
setup(nv, deg) selects the space (module globals; one space per computation).
"""
from itertools import product
from flint import arb, ctx

ctx.cap = 64

NV = 3
DEG = 4
MONO = []
IDX = {}
DEGOF = []
TAB = []
ZERO = None


def setup(nv, deg):
    global NV, DEG, MONO, IDX, DEGOF, TAB, ZERO
    NV, DEG = nv, deg
    MONO = sorted([e for e in product(range(deg + 1), repeat=nv) if sum(e) <= deg], key=lambda e: (sum(e), e))
    IDX = {e: i for i, e in enumerate(MONO)}
    DEGOF = [sum(e) for e in MONO]
    TAB = []
    for e1 in MONO:
        row = []
        for e2 in MONO:
            s = tuple(a + b for a, b in zip(e1, e2))
            row.append(IDX.get(s, -1))
        TAB.append(row)
    ZERO = 0


setup(3, 4)


class Jet:
    __slots__ = ('c',)

    def __init__(self, c=None):
        self.c = c if c is not None else {}

    @staticmethod
    def const(v):
        v = arb(v)
        return Jet({0: v})

    @staticmethod
    def var(i):
        e = [0]*NV; e[i] = 1
        return Jet({IDX[tuple(e)]: arb(1)})

    def c0(self):
        return self.c.get(0, arb(0))

    def coef(self, e):
        return self.c.get(IDX[tuple(e)], arb(0))

    def __add__(self, o):
        if not isinstance(o, Jet):
            o = Jet.const(o)
        d = dict(self.c)
        for k, v in o.c.items():
            d[k] = d[k] + v if k in d else v
        return Jet(d)
    __radd__ = __add__

    def __neg__(self):
        return Jet({k: -v for k, v in self.c.items()})

    def __sub__(self, o):
        if not isinstance(o, Jet):
            o = Jet.const(o)
        d = dict(self.c)
        for k, v in o.c.items():
            d[k] = d[k] - v if k in d else -v
        return Jet(d)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Jet):
            o = arb(o)
            return Jet({k: v*o for k, v in self.c.items()})
        d = {}
        T = TAB
        for k1, v1 in self.c.items():
            row = T[k1]
            for k2, v2 in o.c.items():
                k = row[k2]
                if k < 0:
                    continue
                p = v1*v2
                if k in d:
                    d[k] += p
                else:
                    d[k] = p
        return Jet(d)
    __rmul__ = __mul__

    def nil(self):
        d = dict(self.c); d.pop(0, None)
        return Jet(d)

    def _series(self, coefs):
        """sum coefs[k] (self - c0)^k"""
        dl = self.nil(); out = Jet.const(coefs[0]); p = Jet.const(1)
        for k in range(1, DEG + 1):
            p = p*dl
            if not p.c:
                break
            out = out + p*coefs[k]
        return out

    def inv(self):
        a = self.c0()
        if a.contains(0):
            raise ArithmeticError('Jet.inv: constant term contains 0: %s' % a.str(5))
        ia = 1/a
        return self._series([(-1)**k*ia**(k + 1) for k in range(DEG + 1)])

    def __truediv__(self, o):
        if not isinstance(o, Jet):
            return self*(1/arb(o))
        return self*o.inv()

    def __rtruediv__(self, o):
        return self.inv()*o

    def sqrt(self):
        a = self.c0()
        if not a > 0:
            raise ArithmeticError('Jet.sqrt: constant term not certified positive')
        co = []; b = arb(1); sa = a.sqrt()
        for k in range(DEG + 1):
            co.append(b*sa/a**k)          # binom(1/2, k) a^(1/2 - k)
            b = b*(arb(1)/2 - k)/(k + 1)
        return self._series(co)

    def pow_half_neg(self):
        """self^(-1/2)"""
        a = self.c0()
        if not a > 0:
            raise ArithmeticError('Jet.pow_half_neg: constant term not certified positive')
        co = []; b = arb(1); s = 1/a.sqrt()
        for k in range(DEG + 1):
            co.append(b*s/a**k)           # binom(-1/2, k) a^(-1/2 - k)
            b = b*(-arb(1)/2 - k)/(k + 1)
        return self._series(co)

    def exp(self):
        a = self.c0().exp(); co = []; f = arb(1)
        for k in range(DEG + 1):
            co.append(a/f); f = f*(k + 1)
        return self._series(co)

    def log(self):
        a = self.c0()
        if not a > 0:
            raise ArithmeticError('Jet.log: constant term not certified positive')
        co = [a.log()]
        for k in range(1, DEG + 1):
            co.append(arb((-1)**(k + 1))/(k*a**k))
        return self._series(co)

    def sin(self):
        a = self.c0(); s, c = a.sin(), a.cos()
        ders = [s, c, -s, -c]; co = []; f = arb(1)
        for k in range(DEG + 1):
            co.append(ders[k % 4]/f); f = f*(k + 1)
        return self._series(co)

    def diff(self, i):
        if i >= NV:
            return Jet()
        d = {}
        for k, v in self.c.items():
            e = MONO[k]
            if e[i] > 0:
                ee = list(e); ee[i] -= 1
                d[IDX[tuple(ee)]] = v*e[i]
        return Jet(d)

    def iszero(self):
        return all(v.is_zero() for v in self.c.values())

    def maxrad(self):
        return max((v.rad() for v in self.c.values()), default=arb(0))


def horner(coeffs, dl):
    """sum_j coeffs[j] dl^j, dl a Jet (normally with zero constant term)"""
    out = Jet.const(0)
    for c in reversed(coeffs[:DEG + 1] if dl.c0().is_zero() else coeffs):
        out = out*dl + c
    return out


n4 = range(4)


def minv(g):
    """4x4 inverse of a jet matrix (Gauss-Jordan with pivoting on the largest |c0|)"""
    M = [[g[i][j] for j in n4] + [Jet.const(1) if i == j else Jet() for j in n4] for i in n4]
    for c in n4:
        piv = max(range(c, 4), key=lambda r: float(abs(M[r][c].c0()).mid()))
        M[c], M[piv] = M[piv], M[c]
        ip = M[c][c].inv()
        M[c] = [v*ip for v in M[c]]
        for r in n4:
            if r != c and not M[r][c].iszero():
                f = M[r][c]
                M[r] = [M[r][j] - f*M[c][j] for j in range(8)]
    return [[M[i][4 + j] for j in n4] for i in n4]


def _sum(it):
    s = Jet()
    for v in it:
        if v.c:
            s = s + v
    return s


def curv(g):
    gi = minv(g)
    dg = [[[g[a][b].diff(c) for c in n4] for b in n4] for a in n4]
    half = arb(1)/2
    Gam = [[[_sum(gi[a][d]*(dg[d][b][c] + dg[d][c][b] - dg[b][c][d]) for d in n4 if gi[a][d].c)*half
             for c in n4] for b in n4] for a in n4]
    Riem = [[[[Gam[a][b][d].diff(c) - Gam[a][b][c].diff(d)
               + _sum(Gam[a][c][e]*Gam[e][b][d] - Gam[a][d][e]*Gam[e][b][c] for e in n4)
               for d in n4] for c in n4] for b in n4] for a in n4]
    Ric = [[_sum(Riem[a][b][a][d] for a in n4) for d in n4] for b in n4]
    R = _sum(gi[a][b]*Ric[a][b] for a in n4 for b in n4)
    Rl = [[[[_sum(g[a][e]*Riem[e][b][c][d] for e in n4 if g[a][e].c) for d in n4] for c in n4] for b in n4] for a in n4]
    return dict(g=g, gi=gi, Gam=Gam, Ric=Ric, R=R, Rl=Rl)


def mm(A, B):
    return [[_sum(A[i][k]*B[k][j] for k in n4) for j in n4] for i in n4]


def hess(f, Gam):
    df = [f.diff(c) for c in n4]
    return [[df[a].diff(b) - _sum(Gam[c][a][b]*df[c] for c in n4) for b in n4] for a in n4]


def box2(T, gi, Gam):
    D = [[[T[a][b].diff(c) - _sum(Gam[e][c][a]*T[e][b] + Gam[e][c][b]*T[a][e] for e in n4)
           for c in n4] for b in n4] for a in n4]
    out = [[Jet() for _ in n4] for _ in n4]
    for a in n4:
        for b in n4:
            tot = Jet()
            for c in n4:
                for d in n4:
                    if not gi[c][d].c:
                        continue
                    v = D[a][b][c].diff(d) - _sum(Gam[e][d][a]*D[e][b][c] + Gam[e][d][b]*D[a][e][c]
                                                  + Gam[e][d][c]*D[a][b][e] for e in n4)
                    tot = tot + gi[c][d]*v
            out[a][b] = tot
    return out


def zahn(d, flip=False):
    """Zahn T(1), T(2).  flip=True flips the sign of the R^{rs}R_{m r n s} term (sign-flip null)."""
    g, gi, Gam, Ric, R, Rl = d['g'], d['gi'], d['Gam'], d['Ric'], d['R'], d['Rl']
    H = hess(R, Gam)
    boxR = _sum(gi[a][b]*H[a][b] for a in n4 for b in n4)
    boxRic = box2(Ric, gi, Gam)
    Ricup = mm(mm(gi, Ric), gi)
    Ric2 = _sum(Ric[a][b]*Ricup[a][b] for a in n4 for b in n4)
    RR = [[_sum(Ricup[r][s]*Rl[m][r][v][s] for r in n4 for s in n4) for v in n4] for m in n4]
    sg = -4 if not flip else 4
    T1 = [[4*H[a][b] - 4*g[a][b]*boxR - 4*R*Ric[a][b] + g[a][b]*R*R for b in n4] for a in n4]
    T2 = [[2*H[a][b] - g[a][b]*boxR - 2*boxRic[a][b] + g[a][b]*Ric2 + sg*RR[a][b] for b in n4] for a in n4]
    return T1, T2


def routeA(d, RicS):
    """V_A = box(RicS)/120 + RicS^{rs} C_{m r n s}/60 + g RicS^2/240 - (RicS RicS)_{mn}/60 (R = 0 form)"""
    g, gi, Gam, Rl = d['g'], d['gi'], d['Gam'], d['Rl']
    half = arb(1)/2
    C = [[[[Rl[a][b][c][e] - (g[a][c]*RicS[b][e] - g[a][e]*RicS[b][c] - g[b][c]*RicS[a][e] + g[b][e]*RicS[a][c])*half
            for e in n4] for c in n4] for b in n4] for a in n4]
    Ricup = mm(mm(gi, RicS), gi)
    Ric2 = _sum(RicS[a][b]*Ricup[a][b] for a in n4 for b in n4)
    boxRic = box2(RicS, gi, Gam)
    RC = [[_sum(Ricup[a][b]*C[m][a][v][b] for a in n4 for b in n4) for v in n4] for m in n4]
    RR = mm(mm(RicS, gi), RicS)
    return [[boxRic[m][v]/120 + RC[m][v]/60 + g[m][v]*Ric2/240 - RR[m][v]/60 for v in n4] for m in n4]


def Vxi(T1, T2, xi):
    xi = arb(xi)
    c1 = arb(1)/720 - (6*xi - 1)**2/288
    return [[T1[a][b]*c1 - T2[a][b]/240 for b in n4] for a in n4]


def bach(T1, T2):
    q = -arb(1)/4
    return [[(T2[a][b] - T1[a][b]/3)*q for b in n4] for a in n4]


def proj(T, u, v):
    """T_ab u^a v^b at the base point (u, v: lists of 4 arb)"""
    s = arb(0)
    for a in n4:
        if u[a] == 0:
            continue
        for b in n4:
            if v[b] == 0:
                continue
            s += T[a][b].c0()*u[a]*v[b]
    return s
