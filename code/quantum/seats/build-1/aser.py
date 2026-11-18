"""build-1 QB1: univariate truncated power series over arb (exact ball algebra on coefficient lists), series
reversion and composition, and rigorous re-expansion of a 'curve piece' (polynomial + certified analytic tail)
at a parameter point or interval.

Ser(c, K): coefficients c[0..K] of f(h) = sum c_n h^n (h a local offset).  Every operation returns the exact
truncated result of ball arithmetic on the coefficients, so a true Taylor coefficient inside the input balls
stays inside the output balls (no truncation error: products/quotients/compositions of jets truncated at K are
exact to order K).
"""
from math import comb
from flint import arb, ctx

ctx.cap = 64


def A(v):
    return v if isinstance(v, arb) else arb(v)


class Ser:
    __slots__ = ('c', 'K')

    def __init__(self, c, K):
        c = [A(v) for v in c[:K + 1]]
        c += [arb(0)]*(K + 1 - len(c))
        self.c = c; self.K = K

    @staticmethod
    def const(v, K):
        return Ser([A(v)], K)

    @staticmethod
    def var(v0, K):
        return Ser([A(v0), arb(1)], K)

    def _co(self, o):
        return o if isinstance(o, Ser) else Ser.const(o, self.K)

    def __add__(self, o):
        o = self._co(o); return Ser([a + b for a, b in zip(self.c, o.c)], self.K)
    __radd__ = __add__

    def __sub__(self, o):
        o = self._co(o); return Ser([a - b for a, b in zip(self.c, o.c)], self.K)

    def __rsub__(self, o):
        return self._co(o) - self

    def __neg__(self):
        return Ser([-a for a in self.c], self.K)

    def __mul__(self, o):
        if not isinstance(o, Ser):
            o = A(o); return Ser([a*o for a in self.c], self.K)
        K = self.K; a, b = self.c, o.c
        return Ser([sum((a[i]*b[n - i] for i in range(n + 1)), arb(0)) for n in range(K + 1)], K)
    __rmul__ = __mul__

    def inv(self):
        a = self.c; K = self.K
        if a[0].contains(0):
            raise ArithmeticError('Ser.inv: c0 contains 0')
        out = [1/a[0]]
        for n in range(1, K + 1):
            out.append(-sum((a[i]*out[n - i] for i in range(1, n + 1)), arb(0))/a[0])
        return Ser(out, K)

    def __truediv__(self, o):
        if not isinstance(o, Ser):
            return self*(1/A(o))
        return self*o.inv()

    def __rtruediv__(self, o):
        return self.inv()*o

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        r = Ser.const(1, self.K); b = self
        while n:
            if n & 1: r = r*b
            n >>= 1
            if n: b = b*b
        return r

    def deriv(self):
        return Ser([self.c[n]*n for n in range(1, self.K + 1)], self.K)

    def compose_into(self, inner):
        """self(inner) where inner has zero constant term (inner a Ser)"""
        assert inner.c[0].is_zero(), 'compose: inner constant term must be exactly 0'
        out = Ser.const(0, inner.K)
        for a in reversed(self.c):
            out = out*inner + a
        return out

    def exp(self):
        # f = e^{c0} e^{g}, g = f - c0: solve f' = g' f
        K = self.K; g = self.c
        out = [g[0].exp()]
        for n in range(1, K + 1):
            out.append(sum((k*g[k]*out[n - k] for k in range(1, n + 1)), arb(0))/n)
        return Ser(out, K)

    def log(self):
        K = self.K; a = self.c
        if not a[0] > 0:
            raise ArithmeticError('Ser.log: c0 not positive')
        # (log f)' = f'/f
        d = (self.deriv()/self).c
        return Ser([a[0].log()] + [d[n - 1]/n for n in range(1, K + 1)], K)

    def sqrt(self):
        K = self.K; a = self.c
        if not a[0] > 0:
            raise ArithmeticError('Ser.sqrt: c0 not positive')
        s = [a[0].sqrt()]
        for n in range(1, K + 1):
            acc = a[n] - sum((s[i]*s[n - i] for i in range(1, n)), arb(0))
            s.append(acc/(2*s[0]))
        return Ser(s, K)

    def shift0(self):
        """series with the constant term removed"""
        return Ser([arb(0)] + self.c[1:], self.K)


def revert(f):
    """f(h) = a1 h + a2 h^2 + ... (exact zero constant, a1 != 0) -> g with f(g(d)) = d, to order K"""
    K = f.K
    assert f.c[0].is_zero()
    if f.c[1].contains(0):
        raise ArithmeticError('revert: linear coefficient contains 0')
    g = Ser([arb(0), 1/f.c[1]], K)
    # Newton-free order-by-order: g_n from f(g) = d
    for n in range(2, K + 1):
        e = f.compose_into(g)
        # coefficient n of f(g) must vanish; d/dg_n of coefficient n is a1
        g.c[n] = -e.c[n]/f.c[1]
    return g


# ----------------------------------------------------------------------------------- curve pieces
class Piece:
    """components y_i(p) = sum_n C[i][n] p^n + E_i(p) on the complex disc |p| <= R, with
       kind 'rec':  |E_i(p)| <= M_i (|p|/R)^(K+1)/(1 - |p|/R)     (Cauchy step record, R = Rc)
       kind 'sch':  |E_i(p)| <= M_i (|p|/R)^(K+1)                 (|E_i| <= eps_i on |p| <= R and E_i = O(p^(K+1)))
       C[i] has K+1 entries (degree K)."""

    def __init__(self, C, M, R, kind, K=None):
        self.C = C; self.M = [A(m) for m in M]; self.R = A(R); self.kind = kind
        self.K = K if K is not None else len(C[0]) - 1
        for c in C:
            assert len(c) == self.K + 1

    def tail_on(self, rho, i):
        """sup |E_i| on |p| <= rho (rho < R)"""
        if self.kind == 'sch' and rho is self.R:
            return self.M[i]
        q = A(rho)/self.R
        if self.kind == 'rec':
            if not q < 1: raise ArithmeticError('Piece: rho >= R')
            return self.M[i]*q**(self.K + 1)/(1 - q)
        if not q <= 1: raise ArithmeticError('Piece: rho > R')
        return self.M[i]*q**(self.K + 1)

    def jets(self, ps, J, comps=None):
        """Taylor coefficients (order 0..J) of each component about the parameter value ps (arb; point or
        interval).  For every p in ps the true coefficients lie in the returned balls."""
        ps = A(ps)
        pa = ps.abs_upper()
        if not pa < self.R:
            raise ArithmeticError('Piece.jets: |p| not inside the disc')
        # Tail coefficients: both kinds give |e_n| <= M/R^n for n > K (Cauchy on the full disc; for 'rec' this is
        # exactly how the record bound M (|p|/R)^(K+1)/(1-|p|/R) arises), e_n = 0 for n <= K.  Hence the j-th
        # Taylor coefficient of E at p* is bounded by (M/R^j) sum_{n>K} C(n,j) q^(n-j), q = |p*|/R, and the sum
        # is <= C(K+1,j) q^(K+1-j)/(1 - rq), rq = q (K+2)/(K+2-j) (term ratio bound), valid when rq < 1.
        q = pa/self.R
        comps = range(len(self.C)) if comps is None else comps
        exact0 = ps.is_zero()              # at p = 0 the tail E = O(p^(K+1)) has zero coefficients up to K
        out = []
        for i in comps:
            c = self.C[i]; K = self.K
            row = []
            for j in range(J + 1):
                # polynomial part: sum_n C(n,j) c_n ps^(n-j)  (Horner in ps)
                acc = arb(0)
                for n in range(K, j - 1, -1):
                    acc = acc*ps + c[n]*comb(n, j)
                if exact0 and j <= K:
                    row.append(acc)
                else:
                    # explicit terms n = K+1 .. n1 until the term ratio q (n+1)/(n+1-j) (decreasing in n)
                    # is <= 0.9, then a geometric bound for the rest
                    s = arb(0); n = K + 1
                    term = comb(n, j)*q**(n - j)
                    while True:
                        rq = q*(n + 1)/(n + 1 - j)
                        if rq.upper() <= arb('0.9') or n > K + 20000:
                            break
                        s += term; n += 1; term = term*rq
                    if not rq < 1:
                        raise ArithmeticError('Piece.jets: tail ratio not < 1')
                    s = s + term/(1 - rq)
                    tj = self.M[i]/self.R**j*s
                    row.append(acc + arb(0, tj.upper()))
            out.append(row)
        return out


def ball(m, r):
    return arb(m) + arb(0, arb(r).abs_upper())
