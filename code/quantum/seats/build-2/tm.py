"""build-2 QB1b: Taylor models in one cell variable eps (|eps| <= H) over python-flint arb, and their installation
into the pinned build-1 modules (ajet, aser, charts) WITHOUT editing those files.

A TM is  f(eps) in sum_{k<=P} c_k eps^k + [-r, r]  for every eps in [-H, H]   (c_k arb balls, r >= 0 an arb).
Semantics are those of ball arithmetic applied pointwise in eps; every operation returns a TM that contains the
exact result for every eps in the cell:
  * product: convolution, degrees > P bounded by sum |d_k| H^k, cross terms |P1| r2 + |P2| r1 + r1 r2;
  * analytic f(T) (1/x, sqrt, x^{-1/2} via 1/sqrt, log, exp, sin, cos): T = c0 + u with u's constant term exactly 0,
    f(T) = sum_{k<=P} f^(k)(c0)/k! u^k + M |u|^{P+1},  M >= sup_{xi in range(T)} |f^(P+1)(xi)|/(P+1)!
    (Lagrange remainder; c0 a ball: the Taylor coefficients are evaluated on the ball, which contains them at the
    true constant term).
Comparisons (> 0, contains) are decided on the certified range c0 + [-(sum_{k>=1}|c_k| H^k + r), +...].

install(): rebinds the module-global names that build-1 uses to build constants so that a TM passes through:
  ajet.arb, charts.arb -> ArbW (constructor that returns a TM unchanged; ArbW.pi = arb.pi);  aser.A -> passthrough.
The build-1 algorithms themselves are untouched; with arb inputs ArbW is exactly flint.arb.
"""
from math import factorial, comb
import flint
from flint import arb, ctx

ctx.cap = 64

P = 4
H = arb(0)


def setup(p, h):
    global P, H
    P = p; H = arb(h)


def _up(a):
    return arb(a).abs_upper()


class TM:
    __slots__ = ('c', 'r', '_b')

    def __init__(self, c, r=None):
        c = list(c[:P + 1]); c += [arb(0)]*(P + 1 - len(c))
        self.c = c; self.r = arb(0) if r is None else r; self._b = None

    @staticmethod
    def var(x0):
        """x0 + eps"""
        return TM([arb(x0), arb(1)])

    def B(self):
        """upper bound of |non-constant polynomial part| on the cell"""
        if self._b is None:
            s = arb(0); hk = arb(1)
            for k in range(1, P + 1):
                hk = hk*H
                if not self.c[k].is_zero():
                    s += self.c[k].abs_upper()*hk
            self._b = s
        return self._b

    def rng(self):
        return self.c[0] + arb(0, (self.B() + self.r).upper())

    def mag(self):
        return self.c[0].abs_upper() + self.B() + self.r

    def isconst(self):
        return self.r.is_zero() and all(v.is_zero() for v in self.c[1:])

    # ---- arithmetic
    def __add__(self, o):
        if isinstance(o, TM):
            return TM([a + b for a, b in zip(self.c, o.c)], self.r + o.r)
        c = list(self.c); c[0] = c[0] + o
        return TM(c, self.r)
    __radd__ = __add__

    def __neg__(self):
        return TM([-a for a in self.c], self.r)

    def __sub__(self, o):
        if isinstance(o, TM):
            return TM([a - b for a, b in zip(self.c, o.c)], self.r + o.r)
        c = list(self.c); c[0] = c[0] - o
        return TM(c, self.r)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, TM):
            o = arb(o)
            return TM([a*o for a in self.c], self.r*o.abs_upper())
        if o.isconst():
            return self*o.c[0]
        if self.isconst():
            return o*self.c[0]
        a, b = self.c, o.c
        d = [arb(0)]*(2*P + 1)
        for i in range(P + 1):
            if a[i].is_zero(): continue
            for j in range(P + 1):
                if b[j].is_zero(): continue
                d[i + j] += a[i]*b[j]
        high = arb(0); hk = H**(P + 1)
        for k in range(P + 1, 2*P + 1):
            if not d[k].is_zero():
                high += d[k].abs_upper()*hk
            hk = hk*H
        r = high + self.mag()*o.r + o.mag()*self.r - self.r*o.r     # |P1| r2 + |P2| r1 + r1 r2 with mag incl. r
        return TM(d[:P + 1], arb(r.upper()) if not r.is_zero() else arb(0))
    __rmul__ = __mul__

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = TM([arb(1)]); b = self
        while n:
            if n & 1: out = out*b
            n >>= 1
            if n: b = b*b
        return out

    def _compose(self, f, M):
        """sum_k f[k] u^k + M |u|^{P+1}, u = self - c0 (constant term exactly 0)"""
        u = TM([arb(0)] + self.c[1:], self.r)
        res = TM([f[P]])
        for k in range(P - 1, -1, -1):
            res = res*u + f[k]
        U = u.B() + u.r
        rem = arb(M).abs_upper()*U**(P + 1)
        return TM(res.c, res.r + rem)

    def inv(self):
        R = self.rng()
        if R.contains(0):
            raise ArithmeticError('TM.inv: range contains 0')
        c0 = self.c[0]
        f = [(-1)**k/c0**(k + 1) for k in range(P + 1)]
        M = 1/abs(R).lower()**(P + 2)
        return self._compose(f, M)

    def __truediv__(self, o):
        if isinstance(o, TM):
            return self*o.inv()
        return self*(1/arb(o))

    def __rtruediv__(self, o):
        return self.inv()*o

    def sqrt(self):
        R = self.rng()
        if not R > 0:
            raise ArithmeticError('TM.sqrt: range not positive')
        c0 = self.c[0]; a = arb(1)/2
        f = []; b = arb(1)
        for k in range(P + 1):
            f.append(b*c0.sqrt()/c0**k); b = b*(a - k)/(k + 1)
        # |binom(1/2, P+1)| x^{1/2 - (P+1)} is largest at the lower end of the range
        M = abs(b)*arb(R.lower()).sqrt()/arb(R.lower())**(P + 1)
        return self._compose(f, M)

    def log(self):
        R = self.rng()
        if not R > 0:
            raise ArithmeticError('TM.log: range not positive')
        c0 = self.c[0]
        f = [c0.log()] + [arb((-1)**(k + 1))/(k*c0**k) for k in range(1, P + 1)]
        M = 1/((P + 1)*arb(R.lower())**(P + 1))
        return self._compose(f, M)

    def exp(self):
        R = self.rng(); c0 = self.c[0]
        e = c0.exp()
        f = [e/factorial(k) for k in range(P + 1)]
        M = arb(R.upper()).exp()/factorial(P + 1)
        return self._compose(f, M)

    def sin(self):
        c0 = self.c[0]; s, c = c0.sin(), c0.cos(); d = [s, c, -s, -c]
        return self._compose([d[k % 4]/factorial(k) for k in range(P + 1)], arb(1)/factorial(P + 1))

    def cos(self):
        c0 = self.c[0]; s, c = c0.sin(), c0.cos(); d = [c, -s, -c, s]
        return self._compose([d[k % 4]/factorial(k) for k in range(P + 1)], arb(1)/factorial(P + 1))

    # ---- predicates / interop used by build-1 code
    def contains(self, z):
        return self.rng().contains(z)

    def overlaps(self, o):
        return self.rng().overlaps(o.rng() if isinstance(o, TM) else o)

    def __gt__(self, o):
        return self.rng() > (o.rng() if isinstance(o, TM) else o)

    def __lt__(self, o):
        return self.rng() < (o.rng() if isinstance(o, TM) else o)

    def is_zero(self):
        return self.r.is_zero() and all(v.is_zero() for v in self.c)

    def __eq__(self, o):
        if isinstance(o, (int, float)) and o == 0:
            return self.is_zero()
        return False

    __hash__ = None

    def __abs__(self):
        return abs(self.rng())

    def mid(self):
        return self.rng().mid()

    def rad(self):
        return self.rng().rad()

    def abs_upper(self):
        return self.rng().abs_upper()

    def str(self, *a, **k):
        return 'TM(' + self.rng().str(*a, **k) + ')'


class ArbW:
    """drop-in for the name `arb` inside the pinned build-1 modules: TMs pass through"""
    def __new__(cls, *a):
        if len(a) == 1 and isinstance(a[0], TM):
            return a[0]
        return arb(*a)
    pi = staticmethod(arb.pi)


_installed = False


def install():
    global _installed
    import ajet, aser, charts
    ajet.arb = ArbW
    charts.arb = ArbW
    aser.A = lambda v: v if isinstance(v, (arb, TM)) else arb(v)
    _installed = True


def rng(v):
    return v.rng() if isinstance(v, TM) else arb(v)
