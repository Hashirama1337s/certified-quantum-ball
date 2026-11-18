"""Certified core: generic straight-line tapes over arb / acb / Dual, and helpers.
All arithmetic is python-flint ball arithmetic (rigorous enclosures); no floats."""
import sympy as sp
import flint
from flint import arb, acb, arb_series, ctx

def set_prec(bits):
    ctx.prec = bits

def sym_const(c):
    """exact sympy constant -> arb"""
    c = sp.nsimplify(c) if not isinstance(c, sp.Basic) else c
    if c.is_Rational:
        return arb(int(c.p))/arb(int(c.q))
    raise ValueError('non-rational constant %s' % c)

class Dual:
    """value + gradient (list) ; components any ring type (arb/acb)"""
    __slots__ = ('v', 'd')
    def __init__(self, v, d):
        self.v = v; self.d = d
    @staticmethod
    def lift(x, n):
        return x if isinstance(x, Dual) else Dual(x, [0*x]*n if not isinstance(x, (int,)) else [arb(0)]*n)
    def _o(self, o):
        return o if isinstance(o, Dual) else Dual(o, [0*self.v for _ in self.d])
    def __add__(self, o):
        o = self._o(o); return Dual(self.v + o.v, [a + b for a, b in zip(self.d, o.d)])
    __radd__ = __add__
    def __sub__(self, o):
        o = self._o(o); return Dual(self.v - o.v, [a - b for a, b in zip(self.d, o.d)])
    def __rsub__(self, o):
        return self._o(o) - self
    def __neg__(self):
        return Dual(-self.v, [-a for a in self.d])
    def __mul__(self, o):
        o = self._o(o); return Dual(self.v*o.v, [a*o.v + self.v*b for a, b in zip(self.d, o.d)])
    __rmul__ = __mul__
    def __truediv__(self, o):
        o = self._o(o); q = self.v/o.v
        return Dual(q, [(a - q*b)/o.v for a, b in zip(self.d, o.d)])
    def __rtruediv__(self, o):
        return self._o(o)/self
    def __repr__(self):
        return 'Dual(%s, %s)' % (self.v, self.d)

class Tape:
    """straight-line program from sympy expressions (hash-consed)."""
    def __init__(self, exprs, inputs, params=None):
        self.inputs = list(inputs); self.params = dict(params or {})
        self.nodes = []; self.memo = {}
        for i, v in enumerate(self.inputs):
            self._add(('var', i, None), ('var', v))
        self.outs = [self._b(sp.sympify(e)) for e in exprs]
    def _add(self, node, key):
        if key in self.memo: return self.memo[key]
        self.nodes.append(node); self.memo[key] = len(self.nodes) - 1
        return len(self.nodes) - 1
    def _c(self, c):
        return self._add(('const', c, None), ('const', str(c)))
    def _bin(self, op, a, b):
        key = (op,) + (tuple(sorted((a, b))) if op in ('add', 'mul') else (a, b))
        return self._add((op, a, b), key)
    def _b(self, e):
        if e in self.inputs: return self.memo[('var', e)]
        if e.is_Symbol:
            if e in self.params: return self._c(sp.Rational(self.params[e]))
            raise ValueError('free symbol %s' % e)
        if e.is_Rational: return self._c(sp.Rational(e))
        if e.is_Add:
            args = [self._b(a) for a in e.args]; cur = args[0]
            for a in args[1:]: cur = self._bin('add', cur, a)
            return cur
        if e.is_Mul:
            num, den = [], []
            for a in e.args:
                if a.is_Pow and a.exp.is_Integer and a.exp < 0: den.append(self._b(a.base**(-a.exp)))
                else: num.append(self._b(a))
            cur = num[0] if num else self._c(sp.Integer(1))
            for a in num[1:]: cur = self._bin('mul', cur, a)
            for d in den: cur = self._bin('div', cur, d)
            return cur
        if e.is_Pow and e.exp.is_Integer:
            n = int(e.exp)
            if n < 0: return self._bin('div', self._c(sp.Integer(1)), self._b(e.base**(-n)))
            b = self._b(e.base); res = None; sq = b
            while n:
                if n & 1: res = sq if res is None else self._bin('mul', res, sq)
                n >>= 1
                if n: sq = self._bin('mul', sq, sq)
            return res
        raise ValueError('unsupported %s' % e)

    # ---- pointwise evaluation (any ring type)
    def eval(self, xs, conv=lambda c: sym_const(c)):
        vals = [None]*len(self.nodes)
        for i, (op, a, b) in enumerate(self.nodes):
            if op == 'var': vals[i] = xs[a]
            elif op == 'const': vals[i] = conv(a)
            elif op == 'add': vals[i] = vals[a] + vals[b]
            elif op == 'mul': vals[i] = vals[a]*vals[b]
            elif op == 'div': vals[i] = vals[a]/vals[b]
        return [vals[j] for j in self.outs]

    # ---- Taylor jets order by order (any ring type)
    def jet_init(self, K, zero, conv=lambda c: sym_const(c)):
        self.K = K; self.zero = zero
        self.C = [[zero]*(K + 2) for _ in self.nodes]
        for i, (op, a, b) in enumerate(self.nodes):
            if op == 'const':
                self.C[i][0] = conv(a)
    def jet_order(self, n):
        C = self.C
        for i, (op, a, b) in enumerate(self.nodes):
            if op in ('var', 'const'): continue
            Ca, Cb = C[a], C[b]
            if op == 'add':
                C[i][n] = Ca[n] + Cb[n]
            elif op == 'mul':
                s = Ca[0]*Cb[n]
                for j in range(1, n + 1): s = s + Ca[j]*Cb[n - j]
                C[i][n] = s
            elif op == 'div':
                Ci = C[i]
                if n == 0: Ci[0] = Ca[0]/Cb[0]
                else:
                    s = Ca[n]
                    for j in range(n): s = s - Ci[j]*Cb[n - j]
                    Ci[n] = s/Cb[0]

    # ---- whole-series evaluation with arb_series (fast; for defects of polynomials)
    def series_eval(self, ins, cap):
        old = ctx.cap; ctx.cap = cap
        vals = [None]*len(self.nodes)
        for i, (op, a, b) in enumerate(self.nodes):
            if op == 'var': vals[i] = ins[a]
            elif op == 'const': vals[i] = arb_series([sym_const(a)])
            elif op == 'add': vals[i] = vals[a] + vals[b]
            elif op == 'mul': vals[i] = vals[a]*vals[b]
            elif op == 'div': vals[i] = vals[a]/vals[b]
        out = [vals[j] for j in self.outs]
        ctx.cap = old
        return out

def ub(x):
    """upper bound of |x| as arb (exact-ish)"""
    if isinstance(x, acb): return x.abs_upper()
    return x.abs_upper()

def amax(xs):
    m = arb(0)
    for x in xs:
        m = m.max(ub(x))
    return m

UNIT = arb(0, 1)          # the ball [-1, 1]

def sym(r):
    """the interval [-|r|_up, |r|_up] (rigorous, r any arb/number)"""
    return UNIT*arb(r).abs_upper()

def csym(r):
    """complex box containing the disc |z| <= r"""
    s = sym(r)
    return acb(s, s)
