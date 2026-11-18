"""AMEND-LB-E-2 driver extension (new module; imports the pinned driver/lohner/pipeline).

 (1) PE1 radius ladder: one point run at x~ (reused), then box runs for r in the sealed ladder,
     stop at the first Krawczyk contraction.
 (2) PE2: every ambiguous Lohner step is bisected in x (restarting Lohner steps from the step's
     recorded start enclosure), to depth 12. Classes per (sub)step box:
        '+' V>0 certified, '-' V<0 certified,
        'mono' V straddles 0 but f_V has a certified sign on the box (at most one zero, monotone),
        'amb' anything else  -> bisected; PE2 needs 0 'amb' left and exactly one certified sign change.
 (3) null1: sonic-om guess at the shifted V0 by mpmath Newton (engine/bg.py) before interval Newton.

Optional (OFF by default, NOT part of AMEND-LB-E-2): LB_MVSTART=1 uses a mean-value start set
   y_s(V0) in y_s(V0~) + y_s'(X).(V0 - V0~)  instead of evaluating the manifold over the box.
"""
import sys, os, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
from flint import arb, arb_mat
from core import Dual, sym, amax, sym_const
import param as PM
from param import val
import lohner as LO
import driver as DR

LADDER = ['1e-12', '1e-11', '1e-10', '1e-9', '1e-8', '1e-7']

def integrate_rec(F, S, x0, x1, K=36, q='0.25', h0=0.002, grow=1.25, log=None):
    """as lohner.integrate, but records (x, h_arb, Bx, start_hull) for every accepted step."""
    T0 = time.time(); q = arb(q)
    x = arb(x0); x1 = arb(x1)
    sgn = 1 if float(x1.mid()) > float(x.mid()) else -1
    h = h0; n = 0; rec = []
    while (float(x1.mid()) - float(x.mid()))*sgn > 0:
        rem = abs(float((x1 - x).mid()))
        hh = min(h, rem); last = hh >= rem - 1e-15
        ha = (x1 - x) if last else arb(str(sgn*hh))
        start = S.hull()
        saved = (list(S.c), [row[:] for row in S.B], list(S.r), [v[:] for v in S.vh], [t[:] for t in S.t])
        try:
            Bx, J = LO.lohner_step(F, S, ha, K, q)
        except ArithmeticError:
            S.c, S.B, S.r, S.vh, S.t = saved
            h = hh/2
            if h < 1e-7: raise ArithmeticError('Lohner step collapse at x=%s' % x)
            continue
        rec.append((float(x.mid()), ha, Bx, start))
        x = x1 if last else x + ha
        n += 1; h = hh*grow
        if log and n % 50 == 0:
            log('  lohner step %d x=%.4f h=%.4g width=%s t=%.1fs' % (n, float(x.mid()), hh, max(v.rad() for v in S.hull()).str(3), time.time() - T0))
    return dict(steps=n, rec=rec, time=time.time() - T0, x=x)

class Matcher3(DR.Matcher2):
    def residual(self, V0, W0, xc, om_guess, log=None, mvstart=None):
        Pb = self.Pb; T0 = time.time(); n = 3
        mvstart = (os.environ.get('LB_MVSTART') == '1') if mvstart is None else mvstart
        cvd = lambda like: (lambda c: PM._c(c, like))
        V0d = Dual(V0, [arb(1)])
        sm = Pb.sonic_manifold(V0d, om_guess, self.nu_s)
        Yk, Xk = sm['co']['Y'], sm['co']['X']
        ys = [PM.horner(Yk[i], self.th1 + 0*V0d) for i in range(3)]
        xs = PM.horner(Xk, self.th1 + 0*V0d)
        if mvstart and V0.rad() > 0:
            Vc = arb(V0.mid())
            smc = Pb.sonic_manifold(Dual(Vc, [arb(1)]), om_guess, self.nu_s, tails=False)
            ysc = [PM.horner(smc['co']['Y'][i], self.th1 + 0*V0d) for i in range(3)]
            xsc = PM.horner(smc['co']['X'], self.th1 + 0*V0d)
            dV = V0 - Vc
            ys = [Dual(ysc[i].v + ys[i].d[0]*dV, [ys[i].d[0]]) for i in range(3)]
            xs = Dual(xsc.v + xs.d[0]*dV, [xs.d[0]])
        if self.x0 is None:
            self.x0 = arb(val(xs).mid())
        delta = self.x0 - xs
        fs = Pb.F.eval(ys, conv=cvd(ys[0]))
        y0 = [ys[i] + delta*fs[i] for i in range(3)]
        dmax = val(delta).abs_upper() + sm['epsX']
        _, fB, half_acc = DR.transport_bound(Pb, [val(v) for v in ys], dmax, sm['eps'])
        e0 = sm['eps'] + sm['epsX']*DR.fnorm(fB) + dmax**2*half_acc
        W0d = Dual(W0, [arb(1)])
        cm = Pb.centre_manifold(W0d, self.nu_c)
        Zc = [DR.emb(PM.horner(cm['co']['Y'][i], self.thm + 0*W0d), n, 1) for i in range(3)]
        chi = DR.emb(PM.horner(cm['co']['X'], self.thm + 0*W0d), n, 1)
        xcd = Dual(xc, [arb(0), arb(0), arb(1)])
        xm = self.thm.log() + chi + xcd
        if not hasattr(self, 'xbar'):
            self.xbar = arb(val(xm).mid())
            if os.environ.get('LB_SMOKE'):
                self.xbar = self.x0 - arb(os.environ['LB_SMOKE'])
        S = LO.LohnerSet([val(v) for v in y0],
                         vecs=[[v.d[0] for v in y0]] + [[arb(1) if i == j else arb(0) for i in range(3)] for j in range(3)])
        info = integrate_rec(Pb.F, S, self.x0, self.xbar, K=self.Kt, q=self.q, h0=self.h0, log=log)
        Y = S.hull(); V = S.vecs()
        Phi = [[V[1 + j][i] for j in range(3)] for i in range(3)]
        e_end = amax([sum((Phi[i][j].abs_upper() for j in range(3)), arb(0)) for i in range(3)])*e0
        y_end = [Dual(Y[i], [V[0][i], arb(0), arb(0)]) for i in range(3)]
        dm = xm - self.xbar
        if os.environ.get('LB_SMOKE'):
            dm = dm - arb(val(dm).mid())
        fe = Pb.F.eval(y_end, conv=cvd(y_end[0]))
        dmx = val(dm).abs_upper() + cm['epsX']
        Bt, fBt, half_acc2 = DR.transport_bound(Pb, Y, dmx, e_end)
        ym = [y_end[i] + dm*fe[i] for i in range(3)]
        N, om, Vv = ym
        Zs = [1/N, om*N*N, Vv*N]
        R = [Zs[i] - Zc[i] for i in range(3)]
        Nb, ob, Vb = Bt
        JZ = amax([(1/(Nb*Nb)).abs_upper(), (2*ob*Nb).abs_upper() + (Nb*Nb).abs_upper(), Vb.abs_upper() + Nb.abs_upper()])
        eend = e_end + dmx**2*half_acc2 + DR.fnorm(fBt)*cm['epsX']
        Erad = JZ*eend + cm['eps']
        diag = dict(steps=info['steps'], time=time.time() - T0, rec=info['rec'],
                    boxes=[(x, float(h.mid()), Bx) for (x, h, Bx, st) in info['rec']],
                    e0=e0, e_end=e_end, Erad=Erad, sm=sm, cm=cm, width_end=max(v.rad() for v in Y),
                    width_start=max(val(v).rad() for v in y0), Y=Y)
        return R, Erad, diag

def kr_test(Xt, Rv, C, Rb, Eb, rad):
    Xb = [Xt[i] + sym(rad) for i in range(3)]
    JX = [[Rb[i].d[j] for j in range(3)] for i in range(3)]
    E = [sym(Eb) for i in range(3)]
    Kx = []
    for i in range(3):
        s = Xt[i] - sum((C[i, j]*(Rv[j] + E[j]) for j in range(3)), arb(0))
        for j in range(3):
            cij = (1 if i == j else 0) - sum((C[i, l]*JX[l][j] for l in range(3)), arb(0))
            s += cij*(Xb[j] - Xt[j])
        Kx.append(s)
    return all(Xb[i].contains_interior(Kx[i]) for i in range(3)), Kx, Xb

def krawczyk_ladder(M, Xt, om_guess, ladder=LADDER, log=print):
    """one point run at Xt (reused), then box runs up the ladder; stop at first contraction."""
    Rp, Ep, dp = M.residual(Xt[0], Xt[1], Xt[2], om_guess, log=log)
    Rv = [val(r) for r in Rp]
    C = arb_mat([[arb(Rp[i].d[j].mid()) for j in range(3)] for i in range(3)]).inv()
    if log: log('point run: %.1fs |R|=%s width_end=%s Epoint=%s' % (dp['time'], amax(Rv).str(3), dp['width_end'].str(3), Ep.str(3)))
    tried = []
    for r in ladder:
        rad = arb(r)
        Xb = [Xt[i] + sym(rad) for i in range(3)]
        try:
            Rb, Eb, db = M.residual(Xb[0], Xb[1], Xb[2], om_guess, log=log)
        except ArithmeticError as e:
            tried.append((r, 'evaluation failed: %s' % str(e)[:120]))
            if log: log('rung r=%s: box evaluation failed: %s' % (r, str(e)[:120]))
            continue
        ok, Kx, Xb = kr_test(Xt, Rv, C, Rb, Eb, rad)
        excess = max(float(Kx[i].rad()/Xb[i].rad()) for i in range(3))
        tried.append((r, ok, excess, float(Eb.mid())))
        if log: log('rung r=%s: contraction=%s  max rad(K)/rad(X)=%.3g  Ebox=%s  width_start=%s width_end=%s (%.0fs)'
                    % (r, ok, excess, Eb.str(3), db['width_start'].str(3), db['width_end'].str(3), db['time']))
        if ok:
            return dict(ok=True, r=r, Kx=Kx, Xb=Xb, point=dp, box=db, tried=tried, R=Rv)
    return dict(ok=False, tried=tried, point=dp, R=Rv)

def classify(Pb, Bx):
    Bv = [val(t) for t in Bx]
    if Bv[2] > 0: return '+'
    if Bv[2] < 0: return '-'
    try:
        fv = Pb.F.eval(Bv, conv=sym_const)[2]
        if fv > 0 or fv < 0: return 'mono'
    except Exception:
        pass
    return 'amb'

def pe2_bisect(Pb, rec, K=36, q='0.25', depth=12, log=None):
    """bisect ambiguous Lohner steps from their recorded start enclosures."""
    qa = arb(q)
    def sub(start, h, d):
        cls = []
        S = LO.LohnerSet(start)
        h2 = h/2
        out = []
        for half in range(2):
            st = S.hull()
            Bx, _ = LO.lohner_step(Pb.F, S, h2, K, qa)
            c = classify(Pb, Bx)
            if c == 'amb' and d < depth:
                out += sub(st, h2, d + 1)
            else:
                out.append(c)
        return out
    seq = []; nb = 0
    for (x, ha, Bx, start) in rec:
        c = classify(Pb, Bx)
        if c == 'amb':
            nb += 1
            try:
                seq += sub(start, ha, 1)
            except ArithmeticError:
                seq.append('amb')
        else:
            seq.append(c)
    cert = [s for s in seq if s in ('+', '-')]
    changes = sum(1 for a, b in zip(cert, cert[1:]) if a != b)
    namb = seq.count('amb'); nmono = seq.count('mono')
    return dict(sign_changes=changes, amb=namb, mono=nmono, bisected_steps=nb,
                first=cert[0] if cert else None, last=cert[-1] if cert else None,
                holds=(namb == 0 and changes == 1))

def null1_om_guess(V0shift, om_target, k):
    """mpmath Newton for the sonic om at the shifted V0 (engine/bg.py, independent float code)."""
    import sympy as sp
    from mpmath import mp, mpf
    import bg as BG
    B = BG.Background(sp.Rational(k), dps=30, K=10)
    N0, om0 = B.sonic_point(mpf(V0shift), om_guess=mpf(om_target))
    return mp.nstr(om0, 25)
