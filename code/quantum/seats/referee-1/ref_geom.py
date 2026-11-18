"""referee-1: independent curvature / Zahn-V / Bach code (sympy). Written from scratch for REFEREE-QB1.
Conventions (Wald = Zahn l.227): signature mostly plus,
  R^a_{bcd} = d_c Gamma^a_{db} - d_d Gamma^a_{cb} + Gamma^a_{ce} Gamma^e_{db} - Gamma^a_{de} Gamma^e_{cb},
  R_{bd} = R^a_{bad}  (positive on spheres).
Zahn eq. (11), (12), (19) (src/zahn.txt l.250, l.256, l.329-343):
  T1 = 4 D_a D_b R - 4 g box R - 4 R R_ab + g R^2
  T2 = 2 D_a D_b R - g box R - 2 box R_ab + g Ric.Ric - 4 R^{cd} R_{acbd}
  V(xi) = (1/720 - (6 xi - 1)^2/288) T1 - (1/240) T2
Bach: B_ab = D^c D^d C_{acbd} + (1/2) R^{cd} C_{acbd}.
"""
import itertools
import sympy as sp

S = sp.simplify


class Geom:
    def __init__(self, g, X, simp=sp.cancel):
        self.g = sp.Matrix(g); self.X = list(X); self.n = len(X); self.simp = simp
        self.gi = self.g.inv() if not self.g.is_diagonal() else sp.diag(*[1/self.g[i, i] for i in range(self.n)])
        self.gi = self.gi.applyfunc(simp)
        n = self.n; d = sp.diff
        # Gamma^a_{bc}
        self.G = [[[simp(sum(self.gi[a, e]*(d(self.g[e, b], X[c]) + d(self.g[e, c], X[b]) - d(self.g[b, c], X[e]))
                             for e in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
        G = self.G
        # Riemann R^a_{bcd}
        Rm = {}
        for a, b, c, dd in itertools.product(range(n), repeat=4):
            if c >= dd:
                continue
            v = d(G[a][dd][b], X[c]) - d(G[a][c][b], X[dd]) + sum(G[a][c][e]*G[e][dd][b] - G[a][dd][e]*G[e][c][b]
                                                                  for e in range(n))
            v = simp(v)
            Rm[(a, b, c, dd)] = v; Rm[(a, b, dd, c)] = -v
        for a, b, c in itertools.product(range(n), repeat=3):
            Rm[(a, b, c, c)] = 0
        self.Rup = Rm
        # all-lower R_{abcd}
        self.Rl = {k: simp(sum(self.g[k[0], e]*Rm[(e, k[1], k[2], k[3])] for e in range(n)))
                   for k in itertools.product(range(n), repeat=4)}
        self.Ric = sp.Matrix(n, n, lambda b, dd: simp(sum(Rm[(a, b, a, dd)] for a in range(n))))
        self.R = simp(sum(self.gi[a, b]*self.Ric[a, b] for a in range(n) for b in range(n)))
        self.Ricup = (self.gi*self.Ric*self.gi).applyfunc(simp)

    # covariant derivative of an all-lower tensor stored as dict {index tuple: expr}, rank r -> rank r+1 (new index first)
    def cov(self, T, r):
        n = self.n; out = {}
        for idx in itertools.product(range(n), repeat=r):
            for l in range(n):
                v = sp.diff(T[idx], self.X[l])
                for p in range(r):
                    v -= sum(self.G[e][l][idx[p]]*T[idx[:p] + (e,) + idx[p + 1:]] for e in range(n))
                out[(l,) + idx] = self.simp(v)
        return out

    def tensors(self):
        n = self.n; g = self.g; gi = self.gi; simp = self.simp
        Rs = {(): self.R}
        dR = self.cov(Rs, 0)                                    # D_a R
        ddR = self.cov(dR, 1)                                   # D_b D_a R  (symmetric)
        boxR = simp(sum(gi[a, b]*ddR[(a, b)] for a in range(n) for b in range(n)))
        Ric = {(a, b): self.Ric[a, b] for a in range(n) for b in range(n)}
        dRic = self.cov(Ric, 2); ddRic = self.cov(dRic, 3)
        boxRic = sp.Matrix(n, n, lambda a, b: simp(sum(gi[c, e]*ddRic[(c, e, a, b)] for c in range(n) for e in range(n))))
        RicRic = simp(sum(self.Ricup[a, b]*self.Ric[a, b] for a in range(n) for b in range(n)))
        RR = sp.Matrix(n, n, lambda a, b: simp(sum(self.Ricup[c, e]*self.Rl[(a, c, b, e)] for c in range(n) for e in range(n))))
        DDR = sp.Matrix(n, n, lambda a, b: ddR[(b, a)])
        R = self.R
        T1 = (4*DDR - 4*g*boxR - 4*R*self.Ric + g*R**2).applyfunc(simp)
        T2 = (2*DDR - g*boxR - 2*boxRic + g*RicRic - 4*RR).applyfunc(simp)
        # Weyl, n = 4
        C = {}
        for a, b, c, d in itertools.product(range(n), repeat=4):
            C[(a, b, c, d)] = simp(self.Rl[(a, b, c, d)]
                                   - sp.Rational(1, 2)*(g[a, c]*self.Ric[b, d] - g[a, d]*self.Ric[b, c]
                                                        - g[b, c]*self.Ric[a, d] + g[b, d]*self.Ric[a, c])
                                   + R/6*(g[a, c]*g[b, d] - g[a, d]*g[b, c]))
        dC = self.cov(C, 4)                                     # D_l C_{a c b d}
        # Dv_{a c b} = g^{ld} D_l C_{a c b d}
        Dv = {(a, c, b): simp(sum(gi[l, d]*dC[(l, a, c, b, d)] for l in range(n) for d in range(n)))
              for a, c, b in itertools.product(range(n), repeat=3)}
        dDv = self.cov(Dv, 3)                                   # D_k Dv_{a c b}
        B = sp.Matrix(n, n, lambda a, b: simp(sum(gi[k, c]*dDv[(k, a, c, b)] for k in range(n) for c in range(n))
                                              + sp.Rational(1, 2)*sum(self.Ricup[c, d]*C[(a, c, b, d)]
                                                                      for c in range(n) for d in range(n))))
        self.T1, self.T2, self.B, self.C = T1, T2, B, C
        self.boxR = boxR
        return T1, T2, B

    def V(self, xi):
        xi = sp.sympify(xi)
        return ((sp.Rational(1, 720) - (6*xi - 1)**2/sp.Integer(288))*self.T1 - self.T2/240)
