"""B4 build seat 2 (NEW file; lin_tm.py is imported, never edited): kappa-Taylor models of D, v2.

LinearTM2 = lin_tm.LinearTM with a new transport_tm:
  M1  C1 (derivative) transition ball.  The remainder set is still a complex Lohner set, but the ball
      transition matrix over the disc |d| <= r is no longer the Taylor recursion run on the kappa BALL
      (which wraps: its radius carries exp(h |kappa| ||A1||) from the absolute-value recursion).  It is
         P(kc + d) in P_0 + sum_{j=1..J} P_j d^j + tail,
      with the thin kappa-derivative matrices P_j = (1/j!) d^j P/dkappa^j at kc from the exact hierarchy
      H_j' = A(kc) H_j + A1 H_{j-1} (Taylor in t, same order Ks), and the tail by the Cauchy estimate on
      |d| <= rp:  |P_j| <= G'/rp^j,  G' = exp(a(rp)|h|) + rem_t(rp)  (a(rp) = a0 + (|kc|+rp) a1),
      tail <= G' (r/rp)^(J+1)/(1 - r/rp); rp picked from {4r,8r,16r,32r} by the smallest bound.
      The t-truncation of Phi (uniform over the disc) is added as before.
  M2  optimal Cauchy radius for every t-truncation remainder (lin_stages_v2.rem_opt).
  The polynomial part (exact d-hierarchy on the solution, lin_tm.py:248-262), the omitted d^(m+1)
  forcing, the rounding radii and the sonic/centre stages are unchanged.
rem_method='ball' keeps the v1 kappa-ball matrix (with M2) for comparison.
"""
import math
from flint import arb, acb, acb_mat
from core import csym
import lin_tm as LT
from lin_tm import VTM, TM, shmat
from lin_stages import rownorm
import lin_stages_v2 as LS2
from lin_stages_v2 import qr_frame


class LinearTM2(LT.LinearTM):
    def __init__(self, L, m=24, rem_method='c1', J=6):
        super().__init__(L, m=m, rem_method=rem_method)
        self.J = J
        self.stats = {}

    def transport_tm(self, kc, r, Y):
        L = self.L; m = self.m; KL = L.KL; rr = arb(r); kc = acb(kc)
        kball = kc + csym(rr)
        akc = arb(abs(kc).upper())
        akmax = akc + rr
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        C = Y.C
        Bm = I4; rv = acb_mat([[csym(Y.R[i])] for i in range(4)])
        J = self.J
        st_wrap = 0.0; st_force = 0.0; st_K = 0
        for st in L.steps:
            a = st['a0'] + akmax*st['a1']
            Ks, rem_t = L.step_order(st, a)
            st_K += Ks
            A0c = [st['A0'][n] + kc*st['A1'][n] for n in range(Ks + 1)]
            A1 = st['A1']
            h = acb(st['h']); ah = abs(st['h']).upper()
            W = [C]
            SH = shmat(m + 1)
            for n in range(Ks):
                Ts = acb_mat(4, m + 1); Us = acb_mat(4, m + 1)
                for l in range(n + 1):
                    Ts += A0c[l]*W[n - l]; Us += A1[l]*W[n - l]
                W.append((Ts + Us*SH)*(acb(1)/(n + 1)))
            Cn = acb_mat(4, m + 1)
            for n in range(Ks, -1, -1):
                Cn = Cn*h + W[n]
            Ynorm = max(VTM(C, [arb(0)]*4, rr).norms()) + max(abs(rv[i, 0]).upper() for i in range(4))*rownorm([[Bm[i, j] for j in range(4)] for i in range(4)])
            sup_m = arb(0)
            for n in range(Ks + 1):
                sup_m += max(abs(W[n][i, m]).upper() for i in range(4))*ah**n
            omit = ah*(a*ah).exp()*rr**(m + 1)*st['a1']*sup_m*(1 + rem_t)
            force = rem_t*Ynorm + omit
            if self.rem_method == 'ball':
                Akb = [st['A0'][n] + kball*st['A1'][n] for n in range(Ks + 1)]
                Phi = [I4]
                for n in range(Ks):
                    s_ = acb_mat(4, 4)
                    for l in range(n + 1): s_ += Akb[l]*Phi[n - l]
                    Phi.append(s_*(acb(1)/(n + 1)))
                P = acb_mat(4, 4)
                for n in range(Ks, -1, -1): P = P*h + Phi[n]
                P = P + acb_mat([[csym(rem_t)]*4 for _ in range(4)])
            else:
                H = [[I4]] + [[acb_mat(4, 4)] for _ in range(J)]
                for n in range(Ks):
                    for lv in range(J + 1):
                        s_ = acb_mat(4, 4)
                        for l in range(n + 1):
                            s_ += A0c[l]*H[lv][n - l]
                            if lv: s_ += A1[l]*H[lv - 1][n - l]
                        H[lv].append(s_*(acb(1)/(n + 1)))
                Ph = []
                for lv in range(J + 1):
                    Pm = acb_mat(4, 4)
                    for n in range(Ks, -1, -1): Pm = Pm*h + H[lv][n]
                    Ph.append(Pm)
                # Cauchy tail of the kappa-expansion of the (t-truncated) transition polynomial
                best = None
                for c_ in (4, 8, 16, 32):
                    rp = rr*c_
                    ap = st['a0'] + (akc + rp)*st['a1']
                    Gp = (ap*ah).exp() + LS2.rem_opt(ap, st['h'], st['rho'], Ks)
                    tl = Gp*(rr/rp)**(J + 1)/(1 - rr/rp)
                    if best is None or tl < best: best = tl
                tail = best + rem_t
                P = Ph[0]
                pw = arb(1)
                for lv in range(1, J + 1):
                    pw = pw*rr
                    P = P + acb_mat([[csym(abs(Ph[lv][i, j]).upper()*pw) for j in range(4)] for i in range(4)])
                P = P + acb_mat([[csym(tail)]*4 for _ in range(4)])
            st_wrap = max(st_wrap, max(float(P[i, j].real.rad()) + float(P[i, j].imag.rad()) for i in range(4) for j in range(4)))
            M = P*Bm
            Qm, Qi = qr_frame(M, rv)
            radv = []
            for i in range(4):
                ri = arb(0); pw = arb(1)
                for j in range(m + 1):
                    ri += (Cn[i, j].real.rad() + Cn[i, j].imag.rad())*pw; pw *= rr
                radv.append(ri)
            Cn = acb_mat([[acb(Cn[i, j].real.mid(), Cn[i, j].imag.mid()) for j in range(m + 1)] for i in range(4)])
            fv = acb_mat([[csym(force + radv[i])] for i in range(4)])
            st_force += float(force.upper()) + float(max(radv).upper())
            rv = (Qi*M)*rv + Qi*fv
            Bm = Qm
            C = Cn
        Rend = Bm*rv
        self.stats = dict(max_wrap=st_wrap, sum_force=st_force, mean_K=st_K/max(1, len(L.steps)))
        return VTM(C, [abs(Rend[i, 0]).upper() for i in range(4)], rr)
