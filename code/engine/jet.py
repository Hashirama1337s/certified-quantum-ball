"""Tiny Taylor-jet engine (non-rigorous, mpmath).  A Tape is a straight-line
program of binary ops built from sympy expressions with hash-consing; order(n)
fills the n-th Taylor coefficient of every node from the inputs' coefficients
0..n.  The same structure is what B3 will run in Arb ball arithmetic.
"""
import sympy as sp
from mpmath import mp, mpf, fdot

class Tape:
    def __init__(self, exprs, inputs, params=None, K=40):
        self.inputs = list(inputs)
        self.params = dict(params or {})
        self.nodes = []          # (op, a, b) ; op in var,const,add,sub,mul,div
        self.memo = {}
        self.K = K
        for i, v in enumerate(self.inputs):
            self._add(('var', i, None), key=('var', v))
        self.outs = [self._build(sp.sympify(e)) for e in exprs]

    def _add(self, node, key):
        if key in self.memo:
            return self.memo[key]
        self.nodes.append(node)
        idx = len(self.nodes) - 1
        self.memo[key] = idx
        return idx

    def _const(self, c):
        return self._add(('const', c, None), key=('const', str(c)))

    def _build(self, e):
        if e in self.inputs:
            return self.memo[('var', e)]
        if e.is_Symbol:
            if e in self.params:
                return self._const(self.params[e])
            raise ValueError('free symbol %s' % e)
        if e.is_Number:
            return self._const(sp.Rational(e) if e.is_Rational else sp.Float(e, 80))
        if e.is_NumberSymbol:
            return self._const(sp.Float(sp.N(e, 80), 80))
        if e.is_Add:
            args = [self._build(a) for a in e.args]
            cur = args[0]
            for a in args[1:]:
                cur = self._add(('add', cur, a), key=('add',) + tuple(sorted((cur, a))))
            return cur
        if e.is_Mul:
            num, den = [], []
            for a in e.args:
                if a.is_Pow and a.exp.is_Integer and a.exp < 0:
                    den.append(self._build(a.base**(-a.exp)))
                else:
                    num.append(self._build(a))
            cur = num[0] if num else self._const(1)
            for a in num[1:]:
                cur = self._add(('mul', cur, a), key=('mul',) + tuple(sorted((cur, a))))
            for d in den:
                cur = self._add(('div', cur, d), key=('div', cur, d))
            return cur
        if e.is_Pow and e.exp.is_Integer:
            n = int(e.exp)
            if n < 0:
                return self._add(('div', self._const(1), self._build(e.base**(-n))),
                                 key=('div', self._const(1), self._build(e.base**(-n))))
            b = self._build(e.base)
            # binary powering
            res = None; sq = b; m = n
            while m:
                if m & 1:
                    res = sq if res is None else self._add(('mul', res, sq), key=('mul',) + tuple(sorted((res, sq))))
                m >>= 1
                if m:
                    sq = self._add(('mul', sq, sq), key=('mul', sq, sq))
            return res
        raise ValueError('unsupported %s' % e)

    def reset(self, K=None):
        if K: self.K = K
        self.C = [[mpf(0)]*(self.K + 2) for _ in self.nodes]
        for i, (op, a, b) in enumerate(self.nodes):
            if op == 'const':
                self.C[i][0] = mpf(a.p)/mpf(a.q) if isinstance(a, sp.Rational) else mpf(str(a))

    def set_input(self, i, n, val):
        self.C[i][n] = val

    def order(self, n):
        C = self.C
        for i, (op, a, b) in enumerate(self.nodes):
            if op == 'var' or op == 'const':
                continue
            Ca, Cb = C[a], C[b]
            if op == 'add':
                C[i][n] = Ca[n] + Cb[n]
            elif op == 'sub':
                C[i][n] = Ca[n] - Cb[n]
            elif op == 'mul':
                C[i][n] = fdot(Ca[:n + 1], Cb[n::-1])
            elif op == 'div':
                Ci = C[i]
                if n == 0:
                    Ci[0] = Ca[0]/Cb[0]
                else:
                    Ci[n] = (Ca[n] - fdot(Ci[:n], Cb[n:0:-1]))/Cb[0]

    def out(self, j, n):
        return self.C[self.outs[j]][n]
