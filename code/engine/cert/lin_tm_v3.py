"""B4 speedup seat (NEW file; lin_tm.py and lin_tm_v2.py are imported, never edited).

LinearTM3 = lin_tm_v2.LinearTM2 with the kappa-disc remainder carried in a CHARACTERISTIC frame.

Why.  In v2 the transition matrix over the disc |d| <= r is the ball  P(kc+d) in P_0 + sum_j |P_j| r^j + tail,
with the absolute values taken entrywise in the ORIGINAL variables (Abar, Nbar, omega, V).  The Lohner
remainder vector then grows every step by about  r |h| ||A1||_inf  times the frame conditioning, which
compounds to exp(c r Int ||A1||) over the 1440 steps: the measured exp(c r) law (PROPOSAL-B4-BUILD2 sec. 3).
The d-dependence of P is, to first order, d * Int Phi A1 Phi; in the eigenbasis of the large-kappa principal
part A1 (the characteristic variables w_+-, plus the two kappa-free directions, PROOF-B4-v4 sec. 2) this is
DIAGONAL, so its entrywise modulus is |h| max|mu| r instead of |h| ||A1|| r times frame factors.

What changes (everything else is v2: polynomial part by the exact d-hierarchy, C1 derivative matrices P_j
from the exact kappa-hierarchy, Cauchy tail, t-truncation, omitted d^(m+1) forcing, rounding radii):
  F1  per step s a constant invertible matrix T_s (float, used as an EXACT matrix) is chosen from the step
      data at the midpoint t = h/2:
        frame='A1':   T = [[I, 0], [-B^{-1} C, U]] with A1 = [[0, 0], [C, B]], U = eigenvectors of B
                      (columns 3, 4 = the characteristic directions mu_+-, columns 1, 2 = the kernel of A1);
        frame='full': eigenvectors of A0(h/2) + kc A1(h/2) (the full frozen matrix at the disc centre).
      Any invertible T is rigorous: it is only the basis in which the remainder set is written.
  F2  the remainder set is R_Y = T_s B rv (rv a complex box, B the Lohner frame, both in T_s coordinates).
      One step:  N(d) = T_{s+1}^{-1} P(kc+d) T_s  is enclosed as
        N_0 + sum_{j=1..J} csym(|T_{s+1}^{-1} P_j T_s| r^j) + T_{s+1}^{-1} csym(tail) T_s,
      i.e. the absolute values are taken AFTER the change of variables (this is the whole point);
      the forcing box f (Y coordinates) enters as T_{s+1}^{-1} f;  then the usual weighted-QR Lohner step.
      Start: rv = T_0^{-1} R_Y(start) (one wrap), B = I.  End: T_last := I, so R_Y(end) = B rv.
  The rigour argument is the v2 one applied to the linear map Z_s -> Z_{s+1}, Z_s = T_s^{-1} Y: the
  enclosure of N(d) holds for every |d| <= r because each P_j is thin and the conjugation is exact ball
  arithmetic; the forcing and the Lohner step are unchanged linear algebra.
"""
import math
import numpy as np
from flint import arb, acb, acb_mat
from core import csym
import lin_tm as LT
from lin_tm import VTM, shmat
from lin_stages import rownorm
import lin_stages_v2 as LS2
from lin_stages_v2 import qr_frame
import lin_tm_v2 as LT2


def _cm(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))


def _np(M):
    return np.array([[_cm(M[i, j]) for j in range(4)] for i in range(4)])


def _at(coef, t):
    s = np.zeros((4, 4), complex)
    for n in range(len(coef) - 1, -1, -1):
        s = s*t + coef[n]
    return s


def _exact(Tf):
    return acb_mat([[acb(float(Tf[i, j].real), float(Tf[i, j].imag)) for j in range(4)] for i in range(4)])


class LinearTM3(LT2.LinearTM2):
    def __init__(self, L, m=16, J=4, frame='A1', nfit=12):
        super().__init__(L, m=m, rem_method='c1', J=J)
        self.frame = frame
        self.nfit = nfit
        self._fA = None
        self._Tcache = {}

    def _float_steps(self):
        if self._fA is None:
            self._fA = [([_np(st['A0'][n]) for n in range(self.nfit + 1)], [_np(st['A1'][n]) for n in range(self.nfit + 1)],
                         _cm(acb(st['h']))) for st in self.L.steps]
        return self._fA

    def frames(self, kc):
        """list of (T exact acb_mat, T^{-1} ball acb_mat) per step (+ identity at the end)"""
        key = (self.frame, None if self.frame == 'A1' else (round(kc.real, 12), round(kc.imag, 12)))
        if key in self._Tcache:
            return self._Tcache[key]
        out = []
        for A0c, A1c, h in self._float_steps():
            A1m = _at(A1c, h/2)
            if self.frame == 'A1':
                B = A1m[2:4, 2:4]; C = A1m[2:4, 0:2]
                w, U = np.linalg.eig(B)
                U = U/np.linalg.norm(U, axis=0)
                Tf = np.zeros((4, 4), complex)
                Tf[0:2, 0:2] = np.eye(2); Tf[2:4, 0:2] = -np.linalg.solve(B, C); Tf[2:4, 2:4] = U
            elif self.frame == 'full':
                A = _at(A0c, h/2) + kc*A1m
                w, Tf = np.linalg.eig(A)
                Tf = Tf/np.linalg.norm(Tf, axis=0)
            elif self.frame == 'id':
                Tf = np.eye(4, dtype=complex)
            else:
                raise ValueError(self.frame)
            T = _exact(Tf)
            out.append((T, T.inv()))
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        out.append((I4, I4))
        self._Tcache[key] = out
        return out

    def transport_tm(self, kc, r, Y):
        L = self.L; m = self.m; rr = arb(r); kc = acb(kc)
        akc = arb(abs(kc).upper())
        akmax = akc + rr
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        FR = self.frames(_cm(kc))
        C = Y.C
        T0, T0i = FR[0]
        rv = T0i*acb_mat([[csym(Y.R[i])] for i in range(4)])
        Bm = I4
        J = self.J
        st_wrap = 0.0; st_force = 0.0; st_K = 0; st_cond = 0.0
        for s, st in enumerate(L.steps):
            Ts, _ = FR[s]; _, Tni = FR[s + 1]
            a = st['a0'] + akmax*st['a1']
            Ks, rem_t = L.step_order(st, a)
            st_K += Ks
            A0c = [st['A0'][n] + kc*st['A1'][n] for n in range(Ks + 1)]
            A1 = st['A1']
            h = acb(st['h']); ah = abs(st['h']).upper()
            W = [C]
            SH = shmat(m + 1)
            for n in range(Ks):
                Ts_ = acb_mat(4, m + 1); Us = acb_mat(4, m + 1)
                for l in range(n + 1):
                    Ts_ += A0c[l]*W[n - l]; Us += A1[l]*W[n - l]
                W.append((Ts_ + Us*SH)*(acb(1)/(n + 1)))
            Cn = acb_mat(4, m + 1)
            for n in range(Ks, -1, -1):
                Cn = Cn*h + W[n]
            # |R_Y| <= |T_s B rv| for the forcing scale (as v2's Ynorm, Y coordinates)
            RY = (Ts*Bm)*rv
            Ynorm = max(VTM(C, [arb(0)]*4, rr).norms()) + max(abs(RY[i, 0]).upper() for i in range(4))
            sup_m = arb(0)
            for n in range(Ks + 1):
                sup_m += max(abs(W[n][i, m]).upper() for i in range(4))*ah**n
            omit = ah*(a*ah).exp()*rr**(m + 1)*st['a1']*sup_m*(1 + rem_t)
            force = rem_t*Ynorm + omit
            # exact kappa-hierarchy of the transition polynomial (as v2)
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
            best = None
            for c_ in (4, 8, 16, 32):
                rp = rr*c_
                ap = st['a0'] + (akc + rp)*st['a1']
                Gp = (ap*ah).exp() + LS2.rem_opt(ap, st['h'], st['rho'], Ks)
                tl = Gp*(rr/rp)**(J + 1)/(1 - rr/rp)
                if best is None or tl < best: best = tl
            tail = best + rem_t
            # F2: conjugate the thin matrices first, take absolute values in the characteristic frame
            N = Tni*Ph[0]*Ts
            pw = arb(1)
            for lv in range(1, J + 1):
                pw = pw*rr
                Nj = Tni*Ph[lv]*Ts
                N = N + acb_mat([[csym(abs(Nj[i, j]).upper()*pw) for j in range(4)] for i in range(4)])
            N = N + Tni*acb_mat([[csym(tail)]*4 for _ in range(4)])*Ts
            st_wrap = max(st_wrap, max(float(N[i, j].real.rad()) + float(N[i, j].imag.rad()) for i in range(4) for j in range(4)))
            M = N*Bm
            Qm, Qi = qr_frame(M, rv)
            radv = []
            for i in range(4):
                ri = arb(0); pw = arb(1)
                for j in range(m + 1):
                    ri += (Cn[i, j].real.rad() + Cn[i, j].imag.rad())*pw; pw *= rr
                radv.append(ri)
            Cn = acb_mat([[acb(Cn[i, j].real.mid(), Cn[i, j].imag.mid()) for j in range(m + 1)] for i in range(4)])
            fv = Tni*acb_mat([[csym(force + radv[i])] for i in range(4)])
            st_force += float(force.upper()) + float(max(radv).upper())
            rv = (Qi*M)*rv + Qi*fv
            Bm = Qm
            C = Cn
        Rend = Bm*rv          # T_last = I: Y coordinates
        self.stats = dict(max_wrap=st_wrap, sum_force=st_force, mean_K=st_K/max(1, len(L.steps)))
        return VTM(C, [abs(Rend[i, 0]).upper() for i in range(4)], rr)
