"""B4 build seat 2 (NEW file; lin_stages.py is imported, never edited).

Linear2 = lin_stages.Linear with three changes:
  S1  R7/R8 routing: driver_v32 is imported first (it installs R7 through driver_v31; asserted), and the
      default tube is lin_tube2_k*.pkl (lin_tube_v2.py: E5 enclosure + sonic_manifold_v32).  Stages S, C
      and the RegSing machinery are reused unchanged (lin_stages.py:90-291).
  S2  step splitting (C3 repair).  lin_stages._tube_setup accepted any step whose detMx cell boxes avoid 0,
      however close a complex zero of detMx came: the k=1/2 tube then had steps with a0*rho of order 1e5,
      and the Cauchy remainder exp(a rho) q^(K+1) produced D = +-1e121700 (b4c3.txt, lin_c3_probe.py).
      Here a step is accepted only if a0*rho <= AR_MAX; otherwise it is split in two halves, each started
      from a rigorous ball (Taylor polynomial of the full-step jets at the split point + Cauchy tail
      totrad (|t|/Rc)^(K+2)/(1-|t|/Rc)), and processed recursively (depth <= DMAX).  Every substep starts
      from a ball that contains every certified background, as in the original (lin_stages.py:308-316).
  S3  optimal Cauchy radius in the transport remainder.  For |Phi(t)| <= exp(a|t|) on |t| <= rho, the bound
      exp(a rho') (|h|/rho')^(K+1)/(1-|h|/rho') holds for EVERY rho' in (|h|, rho] (the sup over the
      smaller disc is <= a).  rho' = (K+1)/a minimises it; at large |kappa| this lowers the t-truncation
      remainder by many orders without raising the order (it is (e a|h|/(K+1))^(K+1) when interior).
"""
import sys, os, math, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import driver_v32 as V32                    # FIRST: R7 via driver_v31
import driver_v31 as V31
import param as PM
assert PM.arb_series is V31.arb_series_exact, 'R7 not installed'
from flint import arb, acb, acb_mat
from core import Tape, sym, csym, sym_const
import system as S
import lin_stages as LS
from lin_stages import deser, deriv, smul, sdiv, rownorm, horner, load_expr

AR_MAX = 40          # accept a step only if a0*rho <= AR_MAX
DMAX = 7             # at most 2^7 substeps per tube step
REM_MAX = 1e-18      # and only if its kappa = 0 Cauchy remainder at order KL is <= REM_MAX (C3: q = |h|/rho
                     # reached 0.75 on k=1/2 tau steps whose rho was shrunk, rem 3e6; b4build2/c3diag2.log)
Q_SUB = arb('0.25')


def rem_opt(a, h, rho, K):
    """rigorous min over rho' in (|h|, rho] of exp(a rho') (|h|/rho')^(K+1)/(1-|h|/rho'), with rho'
    chosen in floating point (any choice in the interval is rigorous) and evaluated in arb."""
    ah = abs(h).upper() if hasattr(abs(h), 'upper') else abs(h)
    af = float(a.upper()); hf = float(ah); rf = float(rho.lower())
    rp = (K + 1)/af if af > 0 else rf
    rp = min(rf, max(rp, 1.05*hf))
    rp = arb(repr(rp))
    if not (rp <= rho): rp = rho
    q = ah/rp
    return (a*rp).exp()*q**(K + 1)/(1 - q)


def log_rem_opt(af, hf, rf, K):
    rp = min(rf, max((K + 1)/af if af > 0 else rf, 1.05*hf))
    q = hf/rp
    return af*rp + (K + 1)*math.log(q) - math.log(1 - q)


class Linear2(LS.Linear):
    def __init__(self, k, tube_file=None, ar_max=AR_MAX, dmax=DMAX, nsub=1, **kw):
        import sympy as sp
        self.ar_max = arb(ar_max); self.dmax = dmax; self.nsub = nsub
        self.log_rem_max = math.log(REM_MAX)
        self.split_log = []
        fn = tube_file or '<LAB>/engine/cert/lin_tube2_k%s.pkl' % str(sp.Rational(k)).replace('/', '_')
        super().__init__(k, tube_file=fn, **kw)

    # ---------- stage L with splitting (S2)
    def _tube_setup(self):
        LF = self.LF; KL = self.KL
        import driver_v2 as V2, pipeline as PL, taylor as TY
        from param import _c as _cc
        Ftape = Tape([load_expr(s) for s in PL.FIELDS['f']], [S.N, S.om, S.V], {S.k: self.k})
        T4 = V2.tau_tape(self.k)
        self.steps = []
        nsplit = 0

        def jets(field, y0b, h, Rc):
            y1, Bx, Yc2 = TY.step_cauchy2(field, y0b, h, KL + 1, q=abs(h)/Rc, conv=lambda q_: _cc(q_, arb(0)))
            return Bx, Yc2

        def process(kind, field, ncomp, y0b, h, depth, Rc=None):
            nonlocal nsplit
            if Rc is not None:
                try:
                    Bx, Yc2 = jets(field, y0b, h, Rc)
                except ArithmeticError:
                    Rc = None
            if Rc is None:
                Rc = abs(h)/Q_SUB
                Bx, Yc2 = jets(field, y0b, h, Rc)
            Yi = [[Yc2[i][n] for n in range(KL + 2)] for i in range(ncomp)]
            totrad = [Bx[i].rad() for i in range(ncomp)]
            dY = [deriv(Yi[i]) for i in range(ncomp)]
            xp = dY[3] if kind == 'tau' else [arb(1)] + [arb(0)]*KL
            outs = LF.series_out(LF.step, [Yi[0], Yi[1], Yi[2], dY[1], dY[2], xp], KL)
            Ptop, F0, F1, dM = outs[0:8], outs[8:16], outs[16:24], outs[24]
            ok_dm = not dM[0].contains(0)
            res = None; rho = self.rho_fac*Rc
            m_ = len(Yi[0])
            if ok_dm:
                dYp = [[acb(a) for a in deriv(Yi[i])] for i in range(ncomp)]
                Yp = [[acb(a) for a in Yi[i]] for i in range(ncomp)]
                for shrink in ('1', '0.75', '0.55'):
                    rho_try = rho*arb(shrink)
                    if not rho_try > 1.3*abs(h): break
                    r_ = rho_try/Rc
                    tailY = [totrad[i]*r_**m_/(1 - r_) for i in range(ncomp)]
                    tailD = [totrad[i]/Rc*m_*r_**(m_ - 1)/(1 - r_)**2 for i in range(ncomp)]
                    for nc in (self.ncell, 2*self.ncell):
                        res = self._disc_sup(kind, ncomp, Yp, dYp, tailY, tailD, nc, rho_try)
                        if res is not None: break
                    if res is not None and (res[0]*rho_try) <= self.ar_max and \
                            log_rem_opt(float(res[0].upper()), float(abs(h).upper()), float(rho_try.lower()), KL) <= self.log_rem_max:
                        rho = rho_try; break
                    res = None
            if res is None:
                if depth >= self.dmax:
                    raise ArithmeticError('step %d (%s, h=%s): no acceptable disc after %d splits' % (len(self.steps), kind, h.str(4), depth))
                nsplit += 1
                h2 = h/2
                tq = abs(h2)/Rc
                tail = [totrad[i]*tq**(KL + 2)/(1 - tq) for i in range(ncomp)]
                ymid = [horner(Yi[i], h2) + sym(tail[i]) for i in range(ncomp)]
                process(kind, field, ncomp, y0b, h2, depth + 1)
                process(kind, field, ncomp, ymid, h - h2, depth + 1)
                return
            a0, a1 = res
            A0 = []; A1 = []
            xPt = [smul(xp, Ptop[e], KL) for e in range(8)]
            F0q = [sdiv(F0[e], dM, KL) for e in range(8)]; F1q = [sdiv(F1[e], dM, KL) for e in range(8)]
            for n in range(KL + 1):
                A0.append(acb_mat([[xPt[4*i + j][n] for j in range(4)] for i in range(2)] + [[F0q[4*i + j][n] for j in range(4)] for i in range(2)]))
                A1.append(acb_mat([[0]*4, [0]*4] + [[F1q[4*i + j][n] for j in range(4)] for i in range(2)]))
            if not (a0.is_finite() and a1.is_finite()): raise ArithmeticError('disc sup not finite')
            self.steps.append(dict(h=h, rho=rho, A0=A0, A1=A1, a0=a0, a1=a1, kind=kind, depth=depth))

        for st in self.tb['steps']:
            kind = st['kind']; h = deser(st['h'])
            c = [deser(p) for p in st['c']]
            setrad = deser(st['setrad'])
            ncomp = 4 if kind == 'tau' else 3
            field = T4 if kind == 'tau' else Ftape
            y0b = [c[i] + sym(setrad*arb('1.01')) for i in range(ncomp)]
            if self.nsub <= 1:
                process(kind, field, ncomp, y0b, h, 0)
                continue
            # forced uniform split (speed: q = |h_sub|/rho drops by nsub at the same disc radius Rc)
            Rc = abs(h)/Q_SUB
            Bx, Yc2 = jets(field, y0b, h, Rc)
            Yi = [[Yc2[i][n] for n in range(KL + 2)] for i in range(ncomp)]
            totrad = [Bx[i].rad() for i in range(ncomp)]
            for s in range(self.nsub):
                ts = h*s/self.nsub; hs = h*(s + 1)/self.nsub - ts
                if s == 0:
                    ys = y0b
                else:
                    tq = abs(ts)/Rc
                    ys = [horner(Yi[i], ts) + sym(totrad[i]*tq**(KL + 2)/(1 - tq)) for i in range(ncomp)]
                process(kind, field, ncomp, ys, hs, 0, Rc=Rc)
        self.nsplit = nsplit

    def _disc_sup(self, kind, ncomp, Yp, dYp, tailY, tailD, ncell, rho):
        LF = self.LF
        hc = rho/ncell
        a0 = arb(0); a1 = arb(0)
        for ia in range(ncell):
            for ib in range(ncell):
                tcell = acb(-rho + (2*ia + 1)*hc, -rho + (2*ib + 1)*hc) + csym(hc*arb('1.415'))
                yb = [horner(Yp[i], tcell) + csym(tailY[i]) for i in range(ncomp)]
                dyb = [horner(dYp[i], tcell) + csym(tailD[i]) for i in range(ncomp)]
                xpb = dyb[3] if kind == 'tau' else acb(1)
                o = LF.step.eval([yb[0], yb[1], yb[2], dyb[1], dyb[2], xpb], conv=lambda cc: acb(sym_const(cc)))
                if not o[24].is_finite() or (o[24].real.contains(0) and o[24].imag.contains(0)):
                    return None
                a0 = a0.max(rownorm([[xpb*o[4*i + j] for j in range(4)] for i in range(2)] + [[o[8 + 4*i + j]/o[24] for j in range(4)] for i in range(2)]))
                a1 = a1.max(rownorm([[o[16 + 4*i + j]/o[24] for j in range(4)] for i in range(2)]))
        return a0, a1

    # ---------- transport with the optimal Cauchy radius (S3); otherwise lin_stages.Linear.transport
    def step_order(self, st, a):
        """smallest order K (12..KL) with the optimal-radius remainder below tol; returns (K, rem arb)"""
        af = float(a.upper()); hf = float(abs(st['h']).upper()); rf = float(st['rho'].lower())
        Ks = self.KL
        for K in range(12, self.KL + 1):
            if log_rem_opt(af, hf, rf, K) < self.log_tol:
                Ks = K; break
        return Ks, rem_opt(a, st['h'], st['rho'], Ks)

    def transport(self, kap, Y0, dbg=None):
        kap = acb(kap)
        I4 = acb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
        c = acb_mat([[acb(x.real.mid(), x.imag.mid())] for x in Y0])
        Bm = I4
        r = acb_mat([[Y0[i] - c[i, 0]] for i in range(4)])
        akap = arb(abs(kap).upper())
        for st in self.steps:
            a = st['a0'] + akap*st['a1']
            Ks, rem = self.step_order(st, a)
            A = [st['A0'][n] + kap*st['A1'][n] for n in range(Ks + 1)]
            Phi = [I4]
            for n in range(Ks):
                s = acb_mat(4, 4)
                for m in range(n + 1):
                    s += A[m]*Phi[n - m]
                Phi.append(s*(acb(1)/(n + 1)))
            h = acb(st['h'])
            P = acb_mat(4, 4)
            for n in range(Ks, -1, -1):
                P = P*h + Phi[n]
            P = P + acb_mat([[csym(rem)]*4 for _ in range(4)])
            p = P*c
            M = P*Bm
            Qm, Qi = qr_frame(M, r)
            cn = acb_mat([[acb(p[i, 0].real.mid(), p[i, 0].imag.mid())] for i in range(4)])
            r = (Qi*M)*r + Qi*(p - cn)
            c = cn; Bm = Qm
            if dbg is not None:
                dbg.append((st['kind'], float(st['a0'].mid()), float(st['a1'].mid()), float(st['rho'].mid()), float(abs(st['h']).mid()), float(rem.mid()), Ks))
        Y = c + Bm*r
        return [Y[i, 0] for i in range(4)]


def qr_frame(M, rv):
    """weighted-QR Lohner frame (float Gram-Schmidt; the frame is an exact matrix, its inverse a ball)"""
    w = [float(abs(rv[j, 0]).upper()) + 1e-300 for j in range(4)]
    cols = [[complex(float(M[i, j].real.mid()), float(M[i, j].imag.mid())) for i in range(4)] for j in range(4)]
    order = sorted(range(4), key=lambda j: -math.sqrt(sum(abs(z)**2 for z in cols[j]))*w[j])
    Qv = []
    for j in order:
        v = cols[j][:]
        for qv in Qv:
            dd = sum(x*y.conjugate() for x, y in zip(v, qv)); v = [x - dd*y for x, y in zip(v, qv)]
        nv = math.sqrt(sum(abs(z)**2 for z in v))
        Qv.append([z/nv for z in v])
    Qm = acb_mat([[acb(Qv[j][i].real, Qv[j][i].imag) for j in range(4)] for i in range(4)])
    return Qm, Qm.inv()
