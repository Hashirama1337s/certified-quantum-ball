"""build-1 QB1: certified Evans-Coleman k = 1/3 background -> profile series (F2, At, Wt, Vt) for charts.py.

Sources (all read-only; strings deserialised into rigorous balls, never pickled flint objects):
  TUBE  engine/cert/lin_tube2_k1_3.pkl  (E5 enclosure, BAR-LB-E5 HOLD d0af75803f): inward x-leg steps
        x in [-4.6, -0.05] (step-start ball c +- 1.01 setrad, Cauchy jets recomputed with taylor.step_cauchy2 as in
        sideways/build-1/s1_cert.tube_steps), the sonic manifold sm (K 80, |theta| <= nu_s = 0.06, eps 1e-63) and the
        blown-up centre manifold cm (P, W, U, X)(theta), K 60, |theta| <= 0.03, eps 3.3e-56, x = log theta + X + xc.
  OUTW  sideways/build-5/phase2/outward_C.pkl (BG2, BAR-SF-BG2 HOLD ed61d09268): outward x-steps x in [0.05, 12]
        (Yc, a priori box Bx, Cauchy radius Rc per step) and the sharp sonic manifold (K 16, nu_s 0.004).
The low-order Taylor coefficients of a parameterisation-method manifold at theta = 0 are the stored coefficient
balls themselves (the tail is O(theta^(K+1))); elsewhere tails enter through Cauchy estimates (aser.Piece).
"""
import sys, os, pickle
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, ctx
import driver_v32 as V32          # noqa: F401  (first, as lin_tube_v2 / s1_cert: installs R7; not used for tails here)
import taylor as TY
import param as PM
from core import Tape, sym_const
from pipeline import load_expr, FIELDS
from aser import Ser, revert, Piece

ctx.cap = 64

TUBE = '<LAB>/engine/cert/lin_tube2_k1_3.pkl'
OUTW = '<LAB>/sideways/build-5/phase2/outward_C.pkl'
K13 = sp.Rational(1, 3)
N_, om_, V_, k_ = sp.symbols('N om V k')
def KAPPA():
    return arb(1)/3            # at the current precision (never at import time)


def deser(p):
    """(mid, rad) strings -> rigorous ball, inflated for the decimal rounding of both strings (as s1_cert.deser)"""
    if p[0] == '0' and p[1] == '0':
        return arb(0)                    # stored exact zero (structural zeros of the manifolds, e.g. P(0), X(0))
    m = arb(p[0]); r = arb(p[1])
    return m + arb(0, (r*arb('1.00000001') + abs(m)*arb('1e-75') + arb('1e-300')).upper())


_T = {}


def tapes():
    if not _T:
        kk = {k_: K13}
        _T['F'] = Tape([load_expr(s) for s in FIELDS['f']], [N_, om_, V_], kk)
        _T['A'] = Tape([load_expr(FIELDS['Acons'])], [N_, om_, V_], kk)
    return _T


def cconv(c):
    return PM._c(c, arb(0))


# ------------------------------------------------------------------ (N, om, V) x-series -> profile series
def profiles_from_NomV(x0, cN, co, cV, K):
    """x0 arb (point or interval: every x in x0), cN/co/cV x-Taylor lists about x0 -> F2, At, Wt, Vt lists"""
    N = Ser(cN, K); om = Ser(co, K); V = Ser(cV, K)
    ex = Ser([arb(0), arb(1)], K).exp()*x0.exp()        # e^{x0 + h}
    e2 = ex*ex
    k = KAPPA()
    one = Ser.const(1, K)
    At = om*2*(one + N*V*(k + 1) + V*V*k)/((one - V*V)*e2)
    A = one + At*e2
    # F2 = A M^2 with M = N e^x.  Computing N^2 e^{2x} by series products cancels catastrophically where
    # N ~ e^{-x} (centre) or N ~ e^{-0.84x} (far side); instead use the field equation obeyed by the certified
    # solution, N_x = N (A + k om - om - 2) (KHA 18b, engine system.gN):  (ln M)_x = A - 1 - (1 - k) om.
    s = A - 1 - om*(1 - k)
    lnM = Ser([(N.c[0]).log() + x0] + [s.c[n]/(n + 1) for n in range(K)], K)
    M = lnM.exp()
    F2 = A*M*M
    F2_direct = A*N*N*e2
    chk = all((a - b).contains(0) for a, b in zip(F2.c, F2_direct.c))
    Wt = om/e2
    Vt = V/ex
    # independent check of the A identity: Acons(N, om, V) (engine field) - (1 + At e^{2x}) must contain 0
    Ac = tapes()['A'].eval([N, om, V], conv=lambda c: Ser.const(cconv(c), K))[0]
    chk = chk and all((a - b).contains(0) for a, b in zip(Ac.c, A.c))
    return [F2.c, At.c, Wt.c, Vt.c], chk


def xseries_from_param(pj, xj, K):
    """pj: list of p-Taylor lists of components at p*, xj: p-Taylor list of x(p) at p*  ->  x0, x-Taylor lists"""
    X = Ser(xj, K)
    g = revert(X.shift0())
    return xj[0], [Ser(c, K).compose_into(g).c for c in pj]


# ------------------------------------------------------------------ loaders
class Background:
    def __init__(self, load_tube=True, load_outw=True):
        self.tb = pickle.load(open(TUBE, 'rb')) if load_tube else None
        self.ow = pickle.load(open(OUTW, 'rb')) if load_outw else None
        self._tube_pieces = None

    # ---------- sonic manifolds
    def sm_piece(self, which):
        sm = (self.ow if which == 'bg2' else self.tb)['sm']
        Y = [[deser(p) for p in sm['Y'][i]] for i in range(3)]
        X = [deser(p) for p in sm['X']]
        eps = deser(sm['eps']).abs_upper(); epsX = deser(sm['epsX']).abs_upper()
        # tail bounds hold on |theta| <= nu_s (tube: '0.06'; BG2: [0.004 +/- 1.6e-80]); use a slightly smaller disc
        nu = arb('0.0599999999') if which == 'tube' else arb('0.0039999999')
        return Piece(Y + [X], [eps, eps, eps, epsX], nu, 'sch')

    def sonic_series(self, which, ps, K):
        P = self.sm_piece(which)
        jj = P.jets(ps, K)
        x0, (cN, co, cV) = xseries_from_param(jj[:3], jj[3], K)
        return x0, cN, co, cV

    # ---------- centre manifold
    def cm_piece(self):
        cm = self.tb['cm']
        Y = [[deser(p) for p in cm['Y'][i]] for i in range(3)]
        X = [deser(p) for p in cm['X']]
        eps = deser(cm['eps']).abs_upper(); epsX = deser(cm['epsX']).abs_upper()
        return Piece(Y + [X], [eps, eps, eps, epsX], arb(cm['nu']), 'sch'), deser(self.tb['xc'])

    def centre_profiles_theta(self, ths, K):
        """theta* > 0 (point or interval): p-jets -> profile x-series (F2, At, Wt, Vt) about x(theta*)"""
        P, xc = self.cm_piece()
        jj = P.jets(ths, K)
        Pj, Wj, Uj, Xj = (Ser(c, K) for c in jj)
        th = Ser.var(ths, K)
        Q = Pj/th                                        # P/theta
        E = (Xj + xc).exp()                              # z/theta
        k = KAPPA()
        one = Ser.const(1, K)
        PU = Pj*Uj
        At = Wj*(Q/E)*(Q/E)*2*(one + Uj*(k + 1) + PU*PU*k)/(one - PU*PU)
        z = th*E
        A = one + At*z*z
        F2 = A*(E/Q)*(E/Q)
        Wt = Wj*(Q/E)*(Q/E)
        Vt = Uj*Q/E
        xs = th.log() + Xj + xc
        x0, out = xseries_from_param([F2.c, At.c, Wt.c, Vt.c], xs.c, K)
        return x0, out

    def centre_profiles_w(self, K=8):
        """exact w = z^2 series at the centre (theta = 0) of F2, At, Wt, Vt from the stored coefficients
        (no tail: low-order Taylor coefficients of the manifold are exact).  Returns w-lists and odd-part check."""
        P, xc = self.cm_piece()
        Kz = 2*K + 2
        C = [c[:Kz + 2] for c in P.C]
        Pj = Ser(C[0], Kz + 1); Wj = Ser(C[1], Kz + 1); Uj = Ser(C[2], Kz + 1); Xj = Ser(C[3], Kz + 1)
        assert C[0][0].is_zero(), 'P(0) must be exactly 0'
        Q = Ser(C[0][1:], Kz + 1)                        # P/theta, exact coefficient shift
        E = (Xj + xc).exp()
        k = KAPPA()
        one = Ser.const(1, Kz + 1)
        PU = Pj*Uj
        At = Wj*(Q/E)*(Q/E)*2*(one + Uj*(k + 1) + PU*PU*k)/(one - PU*PU)
        zth = Ser([arb(0)] + E.c[:-1], Kz + 1)           # z(theta) = theta E(theta)
        F2 = (one + At*zth*zth)*(E/Q)*(E/Q)              # A (z/P)^2, A = 1 + At z^2
        Wt = Wj*(Q/E)*(Q/E)
        Vt = Uj*Q/E
        th_of_z = revert(zth)
        outs = []; oddmax = arb(0)
        for f in (F2, At, Wt, Vt):
            fz = f.compose_into(th_of_z)
            ev = [fz.c[2*n] for n in range(K + 1)]
            od = [fz.c[2*n + 1] for n in range(K)]
            for o in od:
                if not o.contains(0):
                    raise ArithmeticError('centre: odd z-coefficient excludes 0 (parity broken)')
                oddmax = max(oddmax, o.abs_upper())
            outs.append(ev)
        return outs, oddmax, xc

    # ---------- tube x-steps (E5)
    def tube_xsteps(self):
        if self._tube_pieces is None:
            tb = self.tb
            xs = deser(tb['X_SWITCH'])
            x = arb('-0.05')
            assert xs.contains(x)
            KT = tb['KT']
            lst = []
            for st in tb['steps']:
                if st['kind'] != 'x':
                    continue
                h = deser(st['h'])
                lst.append(dict(x0=x, h=h, st=st))
                x = x + h
            assert deser(tb['XBAR']).overlaps(x), 'tube x-leg does not end at XBAR'
            self._tube_pieces = lst
        return self._tube_pieces

    def tube_piece(self, rec):
        if 'piece' in rec:
            return rec['piece']
        tb = self.tb; st = rec['st']; KT = tb['KT']; Q = arb(tb['Q'])
        h = rec['h']
        c = [deser(p) for p in st['c']]
        setrad = deser(st['setrad'])
        y0b = [c[i] + arb(0, (setrad*arb('1.01')).upper()) for i in range(3)]
        y1, Bx, Yc = TY.step_cauchy2(tapes()['F'], y0b, h, KT, q=Q, conv=cconv)
        Rc = abs(h)/Q
        P = Piece([Yc[i][:KT + 1] for i in range(3)], [Bx[i].rad() for i in range(3)], Rc, 'rec', K=KT)
        rec['piece'] = P
        return P

    # ---------- outward x-steps (BG2)
    def outw_xsteps(self):
        out = []
        KT = self.ow['KT']
        for st in self.ow['steps']:
            if st['kind'] != 'x':
                continue
            out.append(dict(x0=deser(st['x0']), h=deser(st['h']), st=st))
        return out

    def outw_piece(self, rec):
        if 'piece' in rec:
            return rec['piece']
        st = rec['st']; KT = self.ow['KT']
        Yc = [[deser(p) for p in st['Yc'][i]][:KT + 1] for i in range(3)]
        M = [deser(p).rad() for p in st['Bx']]
        Rc = deser(st['Rc'])
        Rc = arb(Rc.lower())
        P = Piece(Yc, M, Rc, 'rec', K=KT)
        rec['piece'] = P
        return P

    # ---------- point / interval series from x-step records
    def step_series(self, which, x, K):
        """x (arb point or interval inside ONE step of the chosen leg) -> x-lists of N, om, V (method a)"""
        steps = self.tube_xsteps() if which == 'tube' else self.outw_xsteps()
        rec = self.find_step(steps, x)
        P = self.tube_piece(rec) if which == 'tube' else self.outw_piece(rec)
        t = x - rec['x0']
        jj = P.jets(t, K)
        return jj

    def point_series(self, which, x, K):
        """method b: certified state ball at x (record polynomial + order-0 tail) -> x-Taylor lists from the ODE
        y' = f(y) (engine field FIELDS['f'], regular for |x| >= 0.05).  Valid for every x in the ball x."""
        jj = self.step_series(which, x, 0)
        ybox = [jj[i][0] for i in range(3)]
        return ybox, self.ode_jets(ybox, K)

    @staticmethod
    def find_step(steps, x):
        for rec in steps:
            a = rec['x0']; b = rec['x0'] + rec['h']
            lo, hi = (a, b) if b.mid() > a.mid() else (b, a)
            # x may touch the step ends (|t| stays far inside the Cauchy disc R = 4|h|)
            if bool(x.lower() >= lo.lower()) and bool(x.upper() <= hi.upper()):
                return rec
        raise ArithmeticError('no single step contains %s' % x.str(10))

    # ---------- state box over an x-interval (several steps) -> ODE jets (method b)
    def state_box(self, which, xa, xb):
        steps = self.tube_xsteps() if which == 'tube' else self.outw_xsteps()
        box = None
        for rec in steps:
            a = rec['x0']; b = rec['x0'] + rec['h']
            lo, hi = (a, b) if b.mid() > a.mid() else (b, a)
            if hi.mid() < xa.mid() - 1e-30 or lo.mid() > xb.mid() + 1e-30:
                continue
            # overlap of [lo, hi] with [xa, xb]
            l2 = lo if lo.mid() > xa.mid() else xa
            h2 = hi if hi.mid() < xb.mid() else xb
            if h2.mid() < l2.mid():
                continue
            P = self.tube_piece(rec) if which == 'tube' else self.outw_piece(rec)
            ta = l2 - rec['x0']; tb_ = h2 - rec['x0']
            tl = arb(min(ta.lower(), tb_.lower())) if False else ta.union(tb_)
            jj = P.jets(tl, 0)
            v = [jj[i][0] for i in range(3)]
            box = v if box is None else [box[i].union(v[i]) for i in range(3)]
        if box is None:
            raise ArithmeticError('state_box: no step meets [%s, %s]' % (xa.str(8), xb.str(8)))
        # coverage check: the union of the used step intervals must cover [xa, xb]
        return box

    @staticmethod
    def ode_jets(ybox, K):
        Y = TY.jets(tapes()['F'], ybox, K, cconv)
        return [Y[i][:K + 1] for i in range(3)]
