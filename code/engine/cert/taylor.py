"""Certified Taylor integrator (Moore / high-order-enclosure method) over arb balls,
optionally carrying first derivatives (Dual) = variational equation.

Theorem used (Corliss & Rihm 1996; Nedialkov-Jackson-Corliss 1999, "HOE"):
 if  sum_{i<k} t^i f^[i](Y0) + t^k f^[k](B)  is contained in B for all t in [0,h],
 then the solution from any y0 in Y0 exists on [0,h] and stays in B.
Remainder: y(h) = sum_{n<=K} f^[n](y0) h^n + f^[K+1](y(xi)) h^{K+1},  y(xi) in B.
f^[n] = n-th Taylor coefficient function (computed by jets). Applied to the
augmented system (y, dy/dparam) when inputs are Dual.
"""
from flint import arb
from core import Dual, sym
from param import val, jac

def _is_dual(x): return isinstance(x, Dual)

def jets(tape, y0, K, conv):
    z = 0*y0[0]
    tape.jet_init(K + 1, z, conv=conv)
    d = len(y0)
    Y = [[z]*(K + 2) for _ in range(d)]
    for i in range(d): Y[i][0] = y0[i]
    for n in range(K + 1):
        for i in range(d): tape.C[i][n] = Y[i][n]
        tape.jet_order(n)
        for i in range(d): Y[i][n + 1] = tape.C[tape.outs[i]][n]/(n + 1)
    return Y

def horner(c, t):
    r = 0*c[0]
    for a in reversed(c): r = r*t + a
    return r

def hull_contains(B, X):
    """X inside interior of B (componentwise, incl. Dual parts)"""
    if _is_dual(B):
        return hull_contains(B.v, X.v) and all(hull_contains(b, x) for b, x in zip(B.d, X.d))
    return B.contains_interior(X)

def inflate(X, fac=arb('0.5'), add=arb(10)**(-40)):
    if _is_dual(X):
        return Dual(inflate(X.v, fac, add), [inflate(a, fac, add) for a in X.d])
    return X + sym(X.rad()*fac + add)

def step(tape, y0, h, K, Ka=8, conv=None, maxit=12):
    from param import _c
    conv = conv or (lambda c: _c(c, y0[0]))
    d = len(y0)
    tI = h/2 + sym(h/2)                                      # [0,h] (either sign)
    Yc = jets(tape, y0, K, conv)
    # a priori enclosure, low order Ka
    head = [horner(Yc[i][:Ka], tI) for i in range(d)]
    B = [inflate(head[i], arb(1), arb(10)**(-30)) for i in range(d)]
    for it in range(maxit):
        Yb = jets(tape, B, Ka, conv)
        new = [head[i] + Yb[i][Ka]*tI**Ka for i in range(d)]
        if all(hull_contains(B[i], new[i]) for i in range(d)):
            break
        B = [inflate(_union(B[i], new[i]), arb(1), arb(10)**(-30)) for i in range(d)]
    else:
        raise ArithmeticError('a priori enclosure failed (h=%s)' % h)
    # remainder of the main expansion
    YB = jets(tape, B, K + 1, conv) if K + 1 > Ka else Yb
    y1 = [horner(Yc[i][:K + 1], h) + YB[i][K + 1]*h**(K + 1) for i in range(d)]
    return y1, B, Yc

def _union(a, b):
    if _is_dual(a):
        return Dual(a.v.union(b.v), [x.union(y) for x, y in zip(a.d, b.d)])
    return a.union(b)

def lipschitz(tape, B, conv):
    """|| Df(B) ||_inf over the a priori box (values only)"""
    Bv = [val(b) for b in B]
    J, _ = jac(tape, Bv)
    m = arb(0)
    for row in J:
        m = m.max(sum((x.abs_upper() for x in row), arb(0)))
    return m


# ----------------------------------------------------------------------
# Analytic (Cauchy-estimate) remainder: avoids interval jets over wide boxes.
# Complex a priori enclosure (Picard on the complex disc |t| <= Rc):
#   y0 + D_Rc * f(Bc)  inside int(Bc)   ==>  the solution is analytic on |t| < Rc
#   and stays in Bc (continuation along rays; f analytic on Bc).
# Then |y_n| <= M / Rc^n with M = Rc*|f(Bc)|, and for |h| = q Rc:
#   |sum_{n>K} y_n h^n| <= M q^{K+1}/(1-q).
# Same for the variational part Phi (Dual derivative): Phi' = Df(y) Phi.
# ----------------------------------------------------------------------
from flint import acb
from core import csym

def _to_acb(x):
    if isinstance(x, Dual):
        return Dual(_to_acb(x.v), [_to_acb(a) for a in x.d])
    return acb(x)

def _abs_up(x):
    return x.abs_upper()

def step_cauchy(tape, y0, h, K, q=arb('0.25'), conv=None, maxit=30):
    from param import _c
    conv = conv or (lambda c: _c(c, y0[0]))
    d = len(y0)
    Yc = jets(tape, y0, K, conv)
    Rc = abs(h)/q
    D = csym(Rc)
    dual = isinstance(y0[0], Dual)
    y0c = [_to_acb(v) for v in y0]
    f0 = tape.eval(y0c, conv=lambda c: _to_acb(conv(c)))
    def mag(x): return _abs_up(x.v) if isinstance(x, Dual) else _abs_up(x)
    # initial radii
    rv = [2*Rc*_abs_up(val(f0[i]) if not dual else f0[i].v) + arb(10)**(-50) for i in range(d)]
    rd = [[2*Rc*_abs_up(f0[i].d[j]) + arb(10)**(-50) for j in range(len(f0[i].d))] for i in range(d)] if dual else None
    for it in range(maxit):
        if dual:
            Bc = [Dual(y0c[i].v + csym(rv[i]), [y0c[i].d[j] + csym(rd[i][j]) for j in range(len(rd[i]))]) for i in range(d)]
        else:
            Bc = [y0c[i] + csym(rv[i]) for i in range(d)]
        fb = tape.eval(Bc, conv=lambda c: _to_acb(conv(c)))
        ok = True
        nrv = []; nrd = []
        for i in range(d):
            fv = fb[i].v if dual else fb[i]
            if not fv.is_finite():
                ok = False; nrv.append(rv[i]*2); 
                if dual: nrd.append([r*2 for r in rd[i]])
                continue
            m = Rc*_abs_up(fv)
            nrv.append(m)
            if not (m < rv[i]): ok = False
            if dual:
                row = []
                for j in range(len(rd[i])):
                    mj = Rc*_abs_up(fb[i].d[j])
                    row.append(mj)
                    if not (mj < rd[i][j]): ok = False
                nrd.append(row)
        if ok:
            break
        rv = [max_(2*nrv[i], rv[i]*arb('1.5')) for i in range(d)]
        if dual: rd = [[max_(2*nrd[i][j], rd[i][j]*arb('1.5')) for j in range(len(rd[i]))] for i in range(d)]
    else:
        raise ArithmeticError('complex a priori enclosure failed (h=%s)' % h)
    qq = abs(h)/Rc
    geo = qq**(K + 1)/(1 - qq)
    y1 = []
    for i in range(d):
        p = horner(Yc[i][:K + 1], h)
        if dual:
            p = Dual(p.v + sym(nrv[i]*geo), [p.d[j] + sym(nrd[i][j]*geo) for j in range(len(p.d))])
        else:
            p = p + sym(nrv[i]*geo)
        y1.append(p)
    # real box containing y(t), t in [0,h]
    Breal = [(val(y0[i]) + sym(nrv[i])) for i in range(d)]
    return y1, Breal, Yc

def max_(a, b):
    return a.max(b)


class QuotientField:
    """f_i = n_i / d with polynomial n_i, d (tape Pq: outputs n0,n1,n2,d).
    Box evaluation in centred (mean-value) form, which keeps d(B) away from 0
    near the sonic surface where naive interval evaluation fails:
        p(B) in p(c) + grad p(B) . (B - c)        for p in {n_i, d}
        f_i(B) = n_i(B)/d(B);  Df(B) = (grad n_i d - n_i grad d)/d^2   (all over B)."""
    def __init__(self, Pq):
        self.Pq = Pq
    def box(self, B, c, conv=None):
        from core import sym_const
        ka = lambda q: acb(sym_const(q))
        dual = isinstance(B[0], Dual)
        Bv = [b.v for b in B] if dual else B
        cv_ = [x.v for x in c] if dual else c
        pc = self.Pq.eval(cv_, conv=ka)
        J, pB = _jac4(self.Pq, Bv)
        dx = [Bv[j] - cv_[j] for j in range(3)]
        cen = [pc[i] + sum((J[i][j]*dx[j] for j in range(3)), acb(0)) for i in range(4)]
        dB = cen[3]
        out = []
        for i in range(3):
            fi = cen[i]/dB
            if dual:
                grads = [(J[i][l]*dB - cen[i]*J[3][l])/(dB*dB) for l in range(3)]
                dd = [sum((grads[l]*B[l].d[p] for l in range(3)), acb(0)) for p in range(len(B[0].d))]
                out.append(Dual(fi, dd))
            else:
                out.append(fi)
        return out, dB

def _jac4(tape, y):
    from core import sym_const
    z = acb(0); one = acb(1)
    ins = [Dual(y[i], [one if j == i else z for j in range(3)]) for i in range(3)]
    outs = tape.eval(ins, conv=lambda q: Dual(acb(sym_const(q)), [z, z, z]))
    return [o.d for o in outs], [o.v for o in outs]

def step_cauchy_q(tape, Q, y0, h, K, q=arb('0.25'), conv=None, maxit=40):
    """as step_cauchy, but the a priori Picard check uses the centred quotient box Q."""
    from param import _c
    conv = conv or (lambda c: _c(c, y0[0]))
    d = len(y0)
    Yc = jets(tape, y0, K, conv)
    Rc = abs(h)/q
    dual = isinstance(y0[0], Dual)
    y0c = [_to_acb(v) for v in y0]
    f0, _ = Q.box(y0c, y0c, conv)
    rv = [2*Rc*(f0[i].v if dual else f0[i]).abs_upper() + arb(10)**(-50) for i in range(d)]
    rd = [[2*Rc*f0[i].d[j].abs_upper() + arb(10)**(-50) for j in range(len(f0[i].d))] for i in range(d)] if dual else None
    for it in range(maxit):
        if dual:
            Bc = [Dual(y0c[i].v + csym(rv[i]), [y0c[i].d[j] + csym(rd[i][j]) for j in range(len(rd[i]))]) for i in range(d)]
        else:
            Bc = [y0c[i] + csym(rv[i]) for i in range(d)]
        fb, dB = Q.box(Bc, y0c, conv)
        ok = dB.is_finite() and not dB.contains(0)
        nrv = []; nrd = []
        for i in range(d):
            fv = fb[i].v if dual else fb[i]
            m = Rc*fv.abs_upper() if fv.is_finite() else rv[i]*4
            nrv.append(m)
            if not (m < rv[i]): ok = False
            if dual:
                row = []
                for j in range(len(rd[i])):
                    x = fb[i].d[j]
                    mj = Rc*x.abs_upper() if x.is_finite() else rd[i][j]*4
                    row.append(mj)
                    if not (mj < rd[i][j]): ok = False
                nrd.append(row)
        if ok: break
        rv = [max_(nrv[i]*arb('1.3'), rv[i]*arb('1.3')) for i in range(d)]
        if dual: rd = [[max_(nrd[i][j]*arb('1.3'), rd[i][j]*arb('1.3')) for j in range(len(rd[i]))] for i in range(d)]
    else:
        raise ArithmeticError('complex a priori (quotient) failed (h=%s)' % h)
    qq = abs(h)/Rc
    geo = qq**(K + 1)/(1 - qq)
    y1 = []
    for i in range(d):
        p = horner(Yc[i][:K + 1], h)
        if dual:
            p = Dual(p.v + sym(nrv[i]*geo), [p.d[j] + sym(nrd[i][j]*geo) for j in range(len(p.d))])
        else:
            p = p + sym(nrv[i]*geo)
        y1.append(p)
    Breal = [(val(y0[i]) + sym(nrv[i])) for i in range(d)]
    return y1, Breal, Yc


def jacd(tape, y):
    """Jacobian (all outputs x all inputs) over acb/arb inputs y (plain numbers)."""
    from core import sym_const
    d = len(y)
    z = 0*y[0]; one = z + 1
    mk = (lambda q: acb(sym_const(q))) if isinstance(z, acb) else (lambda q: sym_const(q))
    ins = [Dual(y[i], [one if j == i else z for j in range(d)]) for i in range(d)]
    outs = tape.eval(ins, conv=lambda q: Dual(mk(q), [z]*d))
    return [o.d for o in outs], [o.v for o in outs]

def centred_box(tape, Bc, c):
    """f(Bc) in f(c) + Df(Bc)(Bc - c); handles Dual(value, derivative) boxes:
    derivative part = Df(Bc) . dBc."""
    from core import sym_const
    dual = isinstance(Bc[0], Dual)
    Bv = [b.v for b in Bc] if dual else Bc
    cv_ = [x.v for x in c] if dual else c
    d = len(Bv)
    fc = tape.eval(cv_, conv=lambda q: acb(sym_const(q)))
    J, _ = jacd(tape, Bv)
    dx = [Bv[j] - cv_[j] for j in range(d)]
    out = []
    for i in range(len(J)):
        fi = fc[i] + sum((J[i][j]*dx[j] for j in range(d)), acb(0))
        if dual:
            np_ = len(Bc[0].d)
            out.append(Dual(fi, [sum((J[i][j]*Bc[j].d[p] for j in range(d)), acb(0)) for p in range(np_)]))
        else:
            out.append(fi)
    return out

def step_cauchy_c(tape, y0, h, K, q=arb('0.25'), conv=None, maxit=40):
    """Cauchy-remainder step with centred-form complex a priori enclosure."""
    from param import _c
    conv = conv or (lambda c: _c(c, y0[0]))
    d = len(y0)
    Yc = jets(tape, y0, K, conv)
    Rc = abs(h)/q
    dual = isinstance(y0[0], Dual)
    y0c = [_to_acb(v) for v in y0]
    f0 = centred_box(tape, y0c, y0c)
    gv = lambda x: x.v if dual else x
    rv = [2*Rc*gv(f0[i]).abs_upper() + arb(10)**(-50) for i in range(d)]
    rd = [[2*Rc*f0[i].d[j].abs_upper() + arb(10)**(-50) for j in range(len(f0[i].d))] for i in range(d)] if dual else None
    for it in range(maxit):
        if dual:
            Bc = [Dual(y0c[i].v + csym(rv[i]), [y0c[i].d[j] + csym(rd[i][j]) for j in range(len(rd[i]))]) for i in range(d)]
        else:
            Bc = [y0c[i] + csym(rv[i]) for i in range(d)]
        fb = centred_box(tape, Bc, y0c)
        ok = True; nrv = []; nrd = []
        for i in range(d):
            fv = gv(fb[i])
            m = Rc*fv.abs_upper() if fv.is_finite() else rv[i]*4
            nrv.append(m)
            if not (m < rv[i]): ok = False
            if dual:
                row = []
                for j in range(len(rd[i])):
                    x = fb[i].d[j]
                    mj = Rc*x.abs_upper() if x.is_finite() else rd[i][j]*4
                    row.append(mj)
                    if not (mj < rd[i][j]): ok = False
                nrd.append(row)
        if ok: break
        rv = [max_(nrv[i]*arb('1.25'), rv[i]*arb('1.25')) for i in range(d)]
        if dual: rd = [[max_(nrd[i][j]*arb('1.25'), rd[i][j]*arb('1.25')) for j in range(len(rd[i]))] for i in range(d)]
    else:
        raise ArithmeticError('complex a priori (centred) failed (h=%s)' % h)
    qq = abs(h)/Rc
    geo = qq**(K + 1)/(1 - qq)
    y1 = []
    for i in range(d):
        p = horner(Yc[i][:K + 1], h)
        if dual:
            p = Dual(p.v + sym(nrv[i]*geo), [p.d[j] + sym(nrd[i][j]*geo) for j in range(len(p.d))])
        else:
            p = p + sym(nrv[i]*geo)
        y1.append(p)
    Breal = [(val(y0[i]) + sym(nrv[i])) for i in range(d)]
    return y1, Breal, Yc


def hessian_box(tape, B):
    """second derivatives d2f_i/dy_j dy_l over the (acb) box B via hyper-duals."""
    from core import sym_const
    d = len(B); z = acb(0); one = acb(1)
    zi = Dual(z, [z]*d)
    ins = []
    for i in range(d):
        inner = Dual(B[i], [one if l == i else z for l in range(d)])
        ins.append(Dual(inner, [Dual(one if j == i else z, [z]*d) for j in range(d)]))
    outs = tape.eval(ins, conv=lambda q: Dual(Dual(acb(sym_const(q)), [z]*d), [zi]*d))
    H = [[[o.d[j].d[l] for l in range(d)] for j in range(d)] for o in outs]
    return H

def centred2_box(tape, Bc, c, Jc=None):
    """f(Bc) in f(c) + Df(c)(Bc-c) + 1/2 D2f(Bc)[Bc-c]^2 ; Dual part: Df(c)dB + D2f(Bc)[Bc-c, dB]."""
    from core import sym_const
    dual = isinstance(Bc[0], Dual)
    Bv = [b.v for b in Bc] if dual else Bc
    cv_ = [x.v for x in c] if dual else c
    d = len(Bv)
    fc = tape.eval(cv_, conv=lambda q: acb(sym_const(q)))
    if Jc is None:
        Jc, _ = jacd(tape, cv_)
    H = hessian_box(tape, Bv)
    dx = [Bv[j] - cv_[j] for j in range(d)]
    rx = [x.abs_upper() for x in dx]
    out = []
    for i in range(len(fc)):
        lin = fc[i] + sum((Jc[i][j]*dx[j] for j in range(d)), acb(0))
        quad = sum((H[i][j][l].abs_upper()*rx[j]*rx[l] for j in range(d) for l in range(d)), arb(0))/2
        fi = lin + csym(quad)
        if dual:
            np_ = len(Bc[0].d)
            dd = []
            for p in range(np_):
                lp = sum((Jc[i][j]*c[j].d[p] for j in range(d)), acb(0))
                # Phi in box Bc.d ; Df(Bc) Phi in Df(c)Phi + D2f(Bc)[dx, Phi] + Df(c)(Phi - Phi_c)
                rphi = [Bc[j].d[p] - c[j].d[p] for j in range(d)]
                lin2 = lp + sum((Jc[i][j]*rphi[j] for j in range(d)), acb(0))
                q2 = sum((H[i][j][l].abs_upper()*rx[j]*Bc[l].d[p].abs_upper() for j in range(d) for l in range(d)), arb(0))
                dd.append(lin2 + csym(q2))
            out.append(Dual(fi, dd))
        else:
            out.append(fi)
    return out

def step_cauchy2(tape, y0, h, K, q=arb('0.25'), conv=None, maxit=60):
    """Cauchy-remainder step; a priori Picard with second-order centred box evaluation."""
    from param import _c
    conv = conv or (lambda c: _c(c, y0[0]))
    d = len(y0)
    Yc = jets(tape, y0, K, conv)
    Rc = abs(h)/q
    dual = isinstance(y0[0], Dual)
    y0c = [_to_acb(v) for v in y0]
    cv_ = [x.v for x in y0c] if dual else y0c
    Jc, _ = jacd(tape, cv_)
    f0 = centred2_box(tape, y0c, y0c, Jc)
    gv = lambda x: x.v if dual else x
    rv = [Rc*gv(f0[i]).abs_upper()*arb('1.2') + arb(10)**(-50) for i in range(d)]
    rd = [[Rc*f0[i].d[j].abs_upper()*arb('1.2') + arb(10)**(-50) for j in range(len(f0[i].d))] for i in range(d)] if dual else None
    for it in range(maxit):
        if dual:
            Bc = [Dual(y0c[i].v + csym(rv[i]), [y0c[i].d[j] + csym(rd[i][j]) for j in range(len(rd[i]))]) for i in range(d)]
        else:
            Bc = [y0c[i] + csym(rv[i]) for i in range(d)]
        fb = centred2_box(tape, Bc, y0c, Jc)
        ok = True; nrv = []; nrd = []
        for i in range(d):
            fv = gv(fb[i])
            m = Rc*fv.abs_upper() if fv.is_finite() else rv[i]*4
            nrv.append(m)
            if not (m < rv[i]): ok = False
            if dual:
                row = []
                for j in range(len(rd[i])):
                    x = fb[i].d[j]
                    mj = Rc*x.abs_upper() if x.is_finite() else rd[i][j]*4
                    row.append(mj)
                    if not (mj < rd[i][j]): ok = False
                nrd.append(row)
        if ok: break
        rv = [max_(nrv[i]*arb('1.1'), rv[i]) for i in range(d)]
        if dual: rd = [[max_(nrd[i][j]*arb('1.1'), rd[i][j]) for j in range(len(rd[i]))] for i in range(d)]
    else:
        raise ArithmeticError('complex a priori (centred2) failed (h=%s)' % h)
    qq = abs(h)/Rc
    geo = qq**(K + 1)/(1 - qq)
    y1 = []
    for i in range(d):
        p = horner(Yc[i][:K + 1], h)
        if dual:
            p = Dual(p.v + sym(nrv[i]*geo), [p.d[j] + sym(nrd[i][j]*geo) for j in range(len(p.d))])
        else:
            p = p + sym(nrv[i]*geo)
        y1.append(p)
    Breal = [(val(y0[i]) + sym(nrv[i])) for i in range(d)]
    return y1, Breal, Yc
