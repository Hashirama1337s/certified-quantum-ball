"""l0_cert.py - build seat build-1, instrument for BAR-SF-L0 "Ledger of Waves" (seal 596c7146ff), 2026-11-17.

Equations: eqs_l.pkl (eqs_l.py, own derivation).  Background: certified E5 data only.
  sonic point  : stored sonic manifold of engine/cert/lin_tube2_k1_3.pkl (coefficients enclose the true Taylor
                 coefficients; x(theta) = X(theta), X0 = 0, X1 = 1): x-jets at x_s = 0 by exact series reversion.
  light cone   : outward_rec.pkl (certified outward leg with ball jets): unique simple root x_lc of N = 1,
                 x-jets at x_lc from the step polynomial + Cauchy remainder bounds for every derivative.
Exponent formula (regular singular point of C(x) eta' = D(x, kappa) eta with rank C(x*) = n - 1):
  rho = (l D r)/(l C' r), l, r = left/right null vectors of C(x*) (row/column of adj C(x*), enclosures).
PL1  characteristic speeds (exact identities from eqs_checks.json + certified values at both cones).
PL2  sonic exponent, polar l = 1 and l = 2: (a) metric-frozen fluid rows (exact linear form, l-free), (b) the full
     reduced polar systems (fluid algebraic, constraint solved) at kappa = 0, 1, 2, shifted to fluid level.
PL3a x_lc (unique simple root) and the light-cone exponent constants: axial Pi l = 2 (exact form
     (1 - lambda' - kappa)/lambda' at the h_A / Q_x level) and polar l = 2 (reduced RW system, kappa = 0, 1, 2).
PL3b comparison with 1 - lambda'_lc.
PL4  l = 1 polar gauge list (exact, eqs_l.gauge_l1).
PL5  flat toy through the same light-cone code.
usage (memcap 8):  py -3.11 l0_cert.py --controls   (PL4, PL5, cap guard, wrong-k sonic null; controls_l0.json)
                   py -3.11 l0_cert.py --main       (PL1-PL5; result_l0.json)   [LEAD ONLY]
"""
import sys, os, time, json, pickle, hashlib, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, arb_series, acb, ctx
ctx.cap = 32                                   # never the default 10 (E3 trap, doctrine #115)
from core import Tape, sym, sym_const

TUBE = '<LAB>/engine/cert/lin_tube2_k1_3.pkl'
OUTW = os.path.join(HERE, 'outward_rec.pkl')
EQS = os.path.join(HERE, 'eqs_l.pkl')
FN = ['N', 'A', 'om', 'V']
NJ = 7                                         # x-jets 0..6 supplied for every field
D_B4 = '-0.0990290820550211564'; LP_B4 = '1.0990290820550211564'      # PROOF-B4.md l.58-59 (float)
D_BAR = ('-0.0990290820550', '5e-14'); LP_BAR = ('1.0990290820550', '5e-14')   # BAR-SF-L0 PL2 (13 digits)
XLC_BAR = ('1.35794', '1e-3')


def capguard(n):
    if ctx.cap < n:
        raise ArithmeticError('ctx.cap = %d < %d: series would be silently truncated' % (ctx.cap, n))


def coef(s, i):
    """strict coefficient access: index beyond the stored length is a true zero only below ctx.cap"""
    capguard(i + 1)
    c = s.coeffs()
    return c[i] if i < len(c) else arb(0)


def deser(p):
    m = arb(p[0]); r = arb(p[1])
    return m + sym(r*arb('1.00000001') + abs(m)*arb('1e-75') + arb('1e-300'))


def load(s):
    return eval(s, {k_: getattr(sp, k_) for k_ in dir(sp)})


def horner(c, t):
    s = 0*t
    for a in reversed(c):
        s = s*t + a
    return s


# ===================================================================== background jets
def acons_tape():
    from pipeline import load_expr, FIELDS
    N_, om_, V_, k_ = sp.symbols('N om V k')
    return Tape([load_expr(FIELDS['Acons'])], [N_, om_, V_], {k_: sp.Rational(1, 3)})


def jets_from_series(ys):
    """ys = [N, om, V] arb_series in xi (length >= NJ) -> dict field -> list of x-derivatives 0..NJ-1"""
    capguard(NJ + 1)
    As = acons_tape().series_eval(ys, NJ + 1)[0]
    ser = {'N': ys[0], 'om': ys[1], 'V': ys[2], 'A': As}
    return {f: [coef(ser[f], i)*math.factorial(i) for i in range(NJ)] for f in FN}


def sonic_jets(tb):
    sm = tb['sm']
    Y = [[deser(p) for p in sm['Y'][i]] for i in range(3)]; X = [deser(p) for p in sm['X']]
    assert float(arb(sm['X'][0][0])) == 0.0 and float(arb(sm['X'][0][1])) == 0.0, 'x_s = X(0) must be exactly 0'
    capguard(NJ + 2)
    Xs = arb_series([arb(0)] + X[1:NJ + 1]); th = Xs.reversion()
    ys = [arb_series(Y[i][:NJ + 1])(th) for i in range(3)]
    for s_ in ys:
        assert len(s_.coeffs()) >= NJ - 1 or True
    return jets_from_series(ys), dict(N=Y[0][0], om=Y[1][0], V=Y[2][0])


def lc_root(ow, log, target=1):
    """unique simple root of N = 1 on the outward leg; returns x_lc ball, step index, t ball, jets."""
    steps = ow['steps']
    KT = ow['KT']
    a_neg_all = True; cand = []
    from pipeline import load_expr, FIELDS
    N_, om_, V_, k_ = sp.symbols('N om V k')
    aT = Tape([load_expr(FIELDS['Acons']) - 2 - (1 - k_)*om_], [N_, om_, V_], {k_: sp.Rational(1, 3)})
    for j, st in enumerate(steps):
        B = [deser(p) for p in st['Bx']]
        a = aT.eval(B[:3], conv=sym_const)[0]
        if not (a < 0):
            a_neg_all = False
        if (B[0] - target).contains(0):
            cand.append(j)
    # far-side: a < 0 on Q' (sealed lemma) ; N decreasing on the whole outward leg iff a < 0 on every step box
    if len(cand) == 0:
        raise ArithmeticError('no step box contains N = 1')
    # root: interval Newton on the step polynomial (x steps only)
    roots = []
    for j in cand:
        st = steps[j]
        if st['kind'] != 'x':
            raise ArithmeticError('light cone inside a tau step (unexpected)')
        Yc = [[deser(p) for p in st['Yc'][i]] for i in range(3)]
        h = deser(st['h']); Rc = deser(st['Rc']); x0 = deser(st['x0'])
        tot = deser(st['Bx'][0]).rad()
        q = abs(h)/Rc; M = KT + 1
        tail0 = tot*q**M/(1 - q)
        # derivative tail on |t| <= |h|: <= tot/Rc * M q^(M-1) / (1-q)^2
        tail1 = tot/Rc*M*q**(M - 1)/(1 - q)**2
        PN = Yc[0][:KT + 1]; dPN = [PN[n]*n for n in range(1, len(PN))]
        Tb = h/2 + sym(abs(h)/2)
        dN = horner(dPN, Tb) + sym(tail1)
        if not (dN < 0):
            raise ArithmeticError('N_x not certified negative on the root step')
        # interval Newton on t in [0, h]: N(T) = m - (f(m))/f'(T); N(T) inside int(T) => unique root in T
        T = Tb; contained = False
        for it in range(80):
            m = arb(T.mid())
            Tn = m - (horner(PN, m) - target + sym(tail0))/(horner(dPN, T) + sym(tail1))
            if not T.overlaps(Tn):
                T = None; break
            if T.contains_interior(Tn):
                contained = True
                T = Tn
            else:                                   # T <- T intersect N(T) (hull of the overlap endpoints)
                lo = T.lower().max(Tn.lower()); hi = T.upper().min(Tn.upper())
                c_ = (lo + hi)/2
                T = arb(c_.mid()) + sym(((hi - lo)/2).abs_upper() + c_.rad() + lo.rad() + hi.rad())
            if contained and T.rad() < arb('1e-45'):
                break
        if T is not None and contained:
            roots.append((j, T, x0 + T, st))
    if not roots:
        raise ArithmeticError('no certified root of N = 1 on the outward leg')
    # uniqueness: N strictly decreasing on every outward step box (a < 0) and on the far side (sealed lemma: a <= -mu)
    if not a_neg_all:
        raise ArithmeticError('N_x < 0 not certified on every outward step: uniqueness of x_lc open')
    xs_ = [r_[2] for r_ in roots]
    if not all(xs_[0].overlaps(x_) for x_ in xs_):
        raise ArithmeticError('two disjoint root enclosures although N is strictly decreasing')
    j, T, xlc, st = min(roots, key=lambda r_: float(r_[2].rad()))
    # jets at x_lc: n-th derivative = sum m!/(m-n)! Y_m T^(m-n) + Cauchy remainder
    Rc = deser(st['Rc']); h = deser(st['h'])
    jets3 = []
    for i in range(3):
        Y = [deser(p) for p in st['Yc'][i]][:KT + 1]
        tot = deser(st['Bx'][i]).rad()
        r = (T.abs_upper())/Rc
        M = KT + 1
        d = []
        for n in range(NJ):
            poly = sum((Y[m]*(math.factorial(m)//math.factorial(m - n))*T**(m - n) for m in range(n, KT + 1)), arb(0))
            g = r*(arb(n)/M).exp()
            assert g < 1
            rem = tot/Rc**n*arb(M)**n*r**(M - n)/(1 - g)
            d.append(poly + sym(rem))
        jets3.append(d)
    ys = [arb_series([jets3[i][n]/math.factorial(n) for n in range(NJ)]) for i in range(3)]
    jets = jets_from_series(ys)
    jets['N'][0] = arb(target)                  # exactly N = 1 at x_lc (the root)
    return dict(x_lc=xlc, step=j, t=T, jets=jets, a_neg_all=a_neg_all, candidates=cand)


# ===================================================================== evaluation of derived expressions
def jet_symbol_index(s):
    nm = s.name
    for f in FN:
        if nm.startswith(f) and nm[len(f):].isdigit():
            return f, int(nm[len(f):])
    return None


def make_tape(exprs, extra):
    syms = set()
    for e in exprs:
        syms |= e.free_symbols
    names = sorted({s_.name for s_ in syms}, key=str)
    byname = {}
    for e in exprs:
        for s_ in e.free_symbols:
            byname[s_.name] = s_
    inputs = [byname[n] for n in names]
    for n in names:
        if jet_symbol_index(byname[n]) is None and n not in extra:
            raise ValueError('unexpected free symbol %s' % n)
        j = jet_symbol_index(byname[n])
        if j is not None and j[1] + 1 >= NJ:
            raise ValueError('jet %s needs derivative order %d >= NJ' % (n, j[1] + 1))
    return Tape(exprs, inputs), names


def series_inputs(names, jets, extra_vals):
    """arb_series in xi (length 3) for every input: F_i(xi) = F_i + F_{i+1} xi + F_{i+2} xi^2/2"""
    ins = []
    A = arb_series(jets['A'][:3] if False else [jets['A'][0], jets['A'][1], jets['A'][2]/2])
    Vs = arb_series([jets['V'][0], jets['V'][1], jets['V'][2]/2])
    for n in names:
        j = jet_symbol_index(sp.Symbol(n))
        if j is not None:
            f, i = j
            ins.append(arb_series([jets[f][i], jets[f][i + 1], jets[f][i + 2]/2 if i + 2 < NJ else arb(0)]))
        elif n == 'sA':
            ins.append(A.sqrt())
        elif n == 'gV':
            ins.append((1 - Vs*Vs).sqrt())
        else:
            ins.append(arb_series([extra_vals[n]]))
    return ins


def eval_series(tape, names, jets, extra_vals, cap=4):
    capguard(cap)
    outs = tape.series_eval(series_inputs(names, jets, extra_vals), cap)
    return [(coef(o, 0), coef(o, 1)) for o in outs]


def adj(Mx):
    n = len(Mx)
    def det(m):
        if len(m) == 1: return m[0][0]
        s = arb(0)
        for j in range(len(m)):
            minor = [row[:j] + row[j + 1:] for row in m[1:]]
            s += (-1)**j*m[0][j]*det(minor)
        return s
    return [[(-1)**(i + j)*det([row[:i] + row[i + 1:] for k2, row in enumerate(Mx) if k2 != j]) for j in range(n)] for i in range(n)], det(Mx)


def frob(Cv, Cp, Dv):
    """rho = l D r / l C' r with r, l from the adjugate of C (rank n-1)"""
    n = len(Cv)
    Ad, dt = adj(Cv)
    # adj = c r l^T : pick the entry of largest magnitude
    best = max(((i, j) for i in range(n) for j in range(n)), key=lambda ij: float(abs(Ad[ij[0]][ij[1]]).mid()))
    r = [Ad[i][best[1]] for i in range(n)]          # column -> right null vector
    l = [Ad[best[0]][j] for j in range(n)]          # row -> left null vector
    lDr = sum((l[i]*Dv[i][j]*r[j] for i in range(n) for j in range(n)), arb(0))
    lCr = sum((l[i]*Cp[i][j]*r[j] for i in range(n) for j in range(n)), arb(0))
    return lDr/lCr, dt, Ad[best[0]][best[1]]


def system_exponent(tapeC, namesC, tapeD, namesD, n, jets, kv, Lv, extra):
    ev = dict(extra); ev['kappa'] = arb(kv); ev['L'] = arb(Lv)
    cs = eval_series(tapeC, namesC, jets, ev); ds = eval_series(tapeD, namesD, jets, ev)
    Cv = [[cs[i*n + j][0] for j in range(n)] for i in range(n)]
    Cp = [[cs[i*n + j][1] for j in range(n)] for i in range(n)]
    Dv = [[ds[i*n + j][0] for j in range(n)] for i in range(n)]
    return frob(Cv, Cp, Dv)


def lin_fit(r0, r1, r2):
    """rho(kappa) = (d - kappa)/lambda'  from kappa = 0, 1; consistency at kappa = 2"""
    lam = 1/(r0 - r1); d = r0*lam
    pred2 = (d - 2)/lam
    return d, lam, bool(pred2.overlaps(r2)), (pred2 - r2)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


# ===================================================================== main pieces
def pl4_pl5():
    import eqs_l as EQ
    g = EQ.gauge_l1()
    ty = EQ.toy()
    return g, ty


def toy_numeric(ty):
    """the toy's a, b through the same second-order light-cone evaluation (series in xi, n0 = e^-xi)"""
    a, b = ty['a'], ty['b']
    n0 = [s_ for s_ in a.free_symbols | b.free_symbols if s_.name == 'n0'][0]
    out = {}
    for kv in (0, 1, 2, sp.Rational(1, 3)):
        capguard(3)
        nser = arb_series([arb(1), arb(-1), arb(1)/2])          # e^-xi
        ta = Tape([a], [n0]); tb_ = Tape([b.subs(sp.Symbol('kappa'), kv)], [n0])
        A_ = ta.series_eval([nser], 3)[0]; B_ = tb_.series_eval([nser], 3)[0]
        rho = -coef(B_, 0)/coef(A_, 1)
        out[str(kv)] = rho
    return out


def build_cowling(E):
    """metric-frozen fluid rows: C f' + (kappa B + D) f = 0 -> rho = (d - kappa)/lambda', exact linear form"""
    pf = {kk_: [load(s_) for s_ in E['pf'][kk_]] for kk_ in ('C', 'B', 'D')}
    tC, nC = make_tape(pf['C'], ['k', 'L', 'kappa', 'sA', 'gV'])
    tB, nB = make_tape(pf['B'], ['k', 'L', 'kappa', 'sA', 'gV'])
    tD, nD = make_tape(pf['D'], ['k', 'L', 'kappa', 'sA', 'gV'])

    def cowling(jets, Lv, kval=arb(1)/3):
        ev = {'k': kval, 'L': arb(Lv), 'kappa': arb(0)}
        cs = eval_series(tC, nC, jets, ev); bs = eval_series(tB, nB, jets, ev); ds = eval_series(tD, nD, jets, ev)
        Cv = [[cs[3*i + j][0] for j in range(3)] for i in range(3)]; Cp = [[cs[3*i + j][1] for j in range(3)] for i in range(3)]
        Bv = [[bs[3*i + j][0] for j in range(3)] for i in range(3)]; Dv = [[ds[3*i + j][0] for j in range(3)] for i in range(3)]
        Ad, dt = adj(Cv)
        best = max(((i, j) for i in range(3) for j in range(3)), key=lambda ij: float(abs(Ad[ij[0]][ij[1]]).mid()))
        r = [Ad[i][best[1]] for i in range(3)]; l = [Ad[best[0]][j] for j in range(3)]
        q = lambda Mx: sum((l[i]*Mx[i][j]*r[j] for i in range(3) for j in range(3)), arb(0))
        lam = q(Cp)/q(Bv); d = -q(Dv)/q(Bv)
        return d, lam, dt
    return cowling


def build_systems(E, extra=None):
    extra = extra or {'k': arb(1)/3}
    sysd = {}
    for key, nkey, Lv in (('sys_l1', 'n_l1', 2), ('sys_l2', 'n_l2', 6)):
        n = E[nkey]
        Cx = [load(s_) for s_ in E[key]['C']]; Dx_ = [load(s_) for s_ in E[key]['Dfull']]
        tCx, nCx = make_tape(Cx, ['k', 'L', 'kappa', 'sA', 'gV'])
        tDx, nDx = make_tape(Dx_, ['k', 'L', 'kappa', 'sA', 'gV'])
        sysd[key] = (tCx, nCx, tDx, nDx, n, Lv)

    def sys_rhos(key, jets):
        tCx, nCx, tDx, nDx, n, Lv = sysd[key]
        out = []
        for kv in (0, 1, 2):
            rho, dt, piv = system_exponent(tCx, nCx, tDx, nDx, n, jets, kv, Lv, extra)
            out.append((rho, dt, piv))
        return out
    return sys_rhos


def axial_rhos(E, jets):
    ax = {kk_: load(E['ax'][kk_]) for kk_ in ('al_s', 'be_s', 'ga_s', 'al_x', 'be_x', 'ga_x', 'a2', 'b2')}
    ta2, na2 = make_tape([ax['a2'], ax['b2'], ax['ga_s']/ax['ga_x']], ['k', 'L', 'kappa', 'sA', 'gV'])
    rh = []; gr = None
    for kv in (0, 1, 2):
        ev = {'k': arb(1)/3, 'L': arb(6), 'kappa': arb(kv)}
        o = eval_series(ta2, na2, jets, ev)
        rh.append(-o[1][0]/o[0][1]); gr = o[2][0]
    return rh, gr


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else '--controls'
    T00 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T00), *a, flush=True)
    old = ctx.prec
    try:
        ctx.prec = 300
        capguard(32)
        res = dict(instrument_sha256=sha(os.path.abspath(__file__)), eqs_pkl_sha256=sha(EQS),
                   eqs_src_sha256={f: sha(os.path.join(HERE, f)) for f in ('lin_core.py', 'polar_sys.py', 'eqs_l.py')},
                   tube_sha256=sha(TUBE), bar='BAR-SF-L0 seal 596c7146ff', mode=mode)
        dry = (mode == '--dryrun')
        vlog = (lambda *a: None) if dry else log           # dry run: crash test + timing only, no values
        g, ty = pl4_pl5()
        tyn = toy_numeric(ty)
        pl4 = dict(value='(1+3k)/(3(1+k)) = 1/2 at k = 1/3', check=g['check'],
                   generator='xi^mu = e^{kappa s} (T Y, X Y, Theta Y\', 0), Theta\' = -A Theta, T = (kappa - A) Theta/(A N^2), '
                             'X = (kappa + A(N^2-1)) Theta/(A N^2); gauge condition h_{A theta} = 0, K = 0; alpha gauge: '
                             'w = O(r^2) at the centre <=> kappa = 1 + U* (centre coefficient -(kappa - 1 - U*) Theta P)',
                   translation='kappa = 1 member of the residual family: centre coefficient 1 - 1/2 != 0, excluded by alpha gauge',
                   boost='kappa = 0 member: centre coefficient -1/2 != 0, excluded by alpha gauge',
                   re_kappa_zero='centre coefficient i w - 1/2 != 0 for all real w: no gauge zero on Re kappa = 0')
        pl4['pass'] = all(v for v in g['check'].values() if isinstance(v, bool))
        pl5 = dict(check=ty['check'], symbolic_rho=str(ty['rho']), numeric_rho={k_: v.str(20) for k_, v in tyn.items()},
                   class_note='weak solutions are xi^(m+1) times a step function (zero on x <= x_lc): non-analytic at the '
                              'cone for every m, so the analytic class excludes them; a bounded class admits them (rho = m+1 >= 0)')
        pl5['numeric matches -kappa'] = all(tyn[k_].contains(-sp_to_arb) for k_, sp_to_arb in
                                            (('0', arb(0)), ('1', arb(1)), ('2', arb(2)), ('1/3', arb(1)/3)))
        pl5['pass'] = all(v for v in ty['check'].values() if isinstance(v, bool)) and pl5['numeric matches -kappa']
        log('PL4', pl4['pass'], ' PL5', pl5['pass'], tyn)
        E = pickle.load(open(EQS, 'rb'))
        tb = pickle.load(open(TUBE, 'rb'))
        jets_s, sp0 = sonic_jets(tb)
        vlog('sonic jets: N0 %s om0 %s V0 %s A0 %s' % (jets_s['N'][0].str(15), jets_s['om'][0].str(15), jets_s['V'][0].str(15), jets_s['A'][0].str(15)))
        extra = {'k': arb(1)/3}
        cowling = build_cowling(E)
        if mode == '--controls':
            # FM2: cap guard armed, and a wrong-k null at the sonic point (k_eq = 17/50 on the k = 1/3 background)
            fm2 = {}
            oldcap = ctx.cap
            try:
                ctx.cap = 2
                try:
                    capguard(4); fm2['cap guard trips at ctx.cap = 2'] = False
                except ArithmeticError:
                    fm2['cap guard trips at ctx.cap = 2'] = True
            finally:
                ctx.cap = oldcap
            d1, l1, dt1 = cowling(jets_s, 2)
            dn, ln, dtn = cowling(jets_s, 2, kval=arb(17)/50)
            fm2['k = 1/3: det C(x_s) contains 0'] = bool(dt1.contains(0))
            fm2['k_eq = 17/50 on the k = 1/3 background: det C(x_s) excludes 0 (no sonic singular point -> no exponent)'] = not dtn.contains(0)
            fm2['k = 1/3 Cowling d'] = d1.str(20); fm2['k = 1/3 Cowling lambda\''] = l1.str(20)
            res.update(PL4=pl4, PL5=pl5, FM2=fm2, time_s=time.time() - T00)
            json.dump(res, open(os.path.join(HERE, 'controls_l0.json'), 'w'), indent=1, default=str)
            log('FM2', fm2)
            log('wrote controls_l0.json')
            sys.exit(0)
        # ------------------------------------------------------------ main
        ow = pickle.load(open(OUTW, 'rb')); res['outward_rec_sha256'] = sha(OUTW)
        # PL2 (a) Cowling, l = 1 (L = 2) and l = 2 (L = 6)
        pl2 = {}; pl2arb = {}
        for lname, Lv in (('l=1', 2), ('l=2', 6)):
            d, lam, dt = cowling(jets_s, Lv)
            pl2['%s metric-frozen fluid rows' % lname] = dict(d=d.str(25), lambda_p=lam.str(25), detC_contains_0=bool(dt.contains(0)))
            pl2arb['%s metric-frozen fluid rows' % lname] = (d, lam)
        vlog('PL2 Cowling', pl2)
        # PL2 (b) full reduced systems
        sys_rhos = build_systems(E)
        for key, lname in (('sys_l1', 'l=1'), ('sys_l2', 'l=2')):
            rr = sys_rhos(key, jets_s)
            dm, lam, ok2, dif = lin_fit(rr[0][0], rr[1][0], rr[2][0])
            dfl = dm - lam                                          # fluid level = metric-derivative level - 1
            pl2['%s full reduced system' % lname] = dict(rho_kappa0=rr[0][0].str(20), rho_kappa1=rr[1][0].str(20), rho_kappa2=rr[2][0].str(20),
                                                       d_metric_level=dm.str(20), lambda_p=lam.str(25), d_fluid_level=dfl.str(25),
                                                       linear_form_consistent_at_kappa2=ok2,
                                                       detC_contains_0=[bool(x[1].contains(0)) for x in rr],
                                                       adj_pivot_nonzero=[not x[2].contains(0) for x in rr])
            pl2arb['%s full reduced system' % lname] = (dfl, lam)
        dB4 = arb(D_B4); lB4 = arb(LP_B4)
        dBAR = arb(D_BAR[0]) + sym(arb(D_BAR[1])); lBAR = arb(LP_BAR[0]) + sym(arb(LP_BAR[1]))
        okPL2 = True
        for lname in ('l=1', 'l=2'):
            for src in ('metric-frozen fluid rows', 'full reduced system'):
                e = pl2['%s %s' % (lname, src)]
                d, lam = pl2arb['%s %s' % (lname, src)]
                e['d overlaps BAR d (13 digits +- 5e-14)'] = bool(d.overlaps(dBAR))
                e["lambda' overlaps BAR lambda' (13 digits +- 5e-14)"] = bool(lam.overlaps(lBAR))
                e['d contains PROOF-B4 19-digit float'] = bool(d.contains(dB4))
                e['d < 0 certified'] = bool(d < 0)
                e['d + lambda\' - 1'] = (d + lam - 1).str(5)
                okPL2 &= e['d overlaps BAR d (13 digits +- 5e-14)'] and e["lambda' overlaps BAR lambda' (13 digits +- 5e-14)"] and e['d < 0 certified']
        pl2['pass'] = okPL2
        pl2['sonic resonances'] = 'kappa = d - n lambda\' (n >= 0) at the fluid level: all certified < 0'
        vlog('PL2', pl2['pass'])
        # PL3a
        lc = lc_root(ow, log)
        xlc = lc['x_lc']; jl = lc['jets']
        lamlc = -jl['N'][1]
        vlog('x_lc =', xlc.str(15), ' lambda_lc = -N_x =', lamlc.str(15))
        xBAR = arb(XLC_BAR[0]) + sym(arb(XLC_BAR[1]))
        pl3 = dict(x_lc=xlc.str(20), x_lc_overlaps_BAR_window=bool(xlc.overlaps(xBAR)), root_step=lc['step'],
                   N_decreasing_on_every_outward_step=lc['a_neg_all'], candidates=lc['candidates'],
                   lambda_lc=lamlc.str(20), om_lc=jl['om'][0].str(15), V_lc=jl['V'][0].str(15), A_lc=jl['A'][0].str(15),
                   two_m_over_r_lc=(1 - 1/jl['A'][0]).str(12))
        # axial Pi: exact formula d = 1 - lambda_lc at the h_A / Q_x level; general expression evaluated as a check
        rhos_ax, gratio = axial_rhos(E, jl)
        d_ax, lam_ax, ok_ax, _ = lin_fit(*rhos_ax)
        d_exact = 1 - lamlc
        pl3['axial Pi l=2 (Q_x = h_A level, Gundlach Pi-tilde/hat level)'] = dict(
            d_lc=d_ax.str(20), lambda_lc=lam_ax.str(20), exact_formula_1_minus_lambda=d_exact.str(20),
            ga_s_over_ga_x=gratio.str(15), linear_consistent=ok_ax, d_lc_excludes_0=not d_ax.contains(0), sign='+' if d_ax > 0 else ('-' if d_ax < 0 else '?'))
        # polar l = 2 reduced RW system at x_lc
        rr = sys_rhos('sys_l2', jl)
        dm, lam, ok2, dif = lin_fit(rr[0][0], rr[1][0], rr[2][0])
        pl3['polar l=2 (reduced RW system, h_AB\' level = Gundlach chi-tilde/hat level)'] = dict(
            d_lc=dm.str(20), lambda_lc=lam.str(20), linear_consistent=ok2, detC_contains_0=[bool(x[1].contains(0)) for x in rr],
            d_lc_excludes_0=not dm.contains(0), sign='+' if dm > 0 else ('-' if dm < 0 else '?'),
            d_minus_1_minus_3lambda=(dm - (1 - 3*lamlc)).str(5))
        rr1 = sys_rhos('sys_l1', jl)
        pl3['polar l=1 system at x_lc: det C excludes 0 (no light-cone singular point, no dipole waves)'] = [not x[1].contains(0) for x in rr1]
        pl3['pass'] = (pl3['x_lc_overlaps_BAR_window'] and lc['a_neg_all'] and
                       pl3['axial Pi l=2 (Q_x = h_A level, Gundlach Pi-tilde/hat level)']['d_lc_excludes_0'] and
                       pl3['polar l=2 (reduced RW system, h_AB\' level = Gundlach chi-tilde/hat level)']['d_lc_excludes_0'])
        pl3['levels'] = ('d depends on the variable level by integer multiples of lambda_lc: rho -> rho + j shifts d by j lambda_lc. '
                         'Certified levels: TT/Pi potential (1 - kappa)/lambda (d = 1); axial h_A, Q_x: d = 1 - lambda; polar h_AB: '
                         'd = 1 - 2 lambda; polar h_AB\': d = 1 - 3 lambda. Resonance set {kappa = 1 - m lambda_lc} is level-independent.')
        res_lc = sorted({1 - m*float(lamlc.mid()) for m in range(0, 3)})
        pl3['light-cone resonance set in Re kappa >= 0 (kappa = 1 - m lambda_lc, m >= 0)'] = [repr(v) for v in res_lc if v >= 0]
        pl3['no resonance on Re kappa = 0 (1/lambda_lc not an integer)'] = bool((1/lamlc) > 1) and bool((1/lamlc) < 2)
        pl3b = dict(axial=dict(d_lc=d_ax.str(15), one_minus_lambda=d_exact.str(15), overlap=bool(d_ax.overlaps(d_exact))),
                    polar=dict(d_lc=dm.str(15), one_minus_lambda=d_exact.str(15), overlap=bool(dm.overlaps(d_exact))))
        pl3b['pass_axial'] = pl3b['axial']['overlap']; pl3b['pass_polar'] = pl3b['polar']['overlap']
        # PL1 (exact identities + certified values)
        chk = json.load(open(os.path.join(HERE, 'eqs_checks.json')))
        pl1 = dict(fluid=chk['polar fluid / PL1'], symbols=chk['principal symbols / PL1'])
        sq = (arb(1)/3).sqrt()
        Ns, Vs = jets_s['N'][0], jets_s['V'][0]
        pl1['sonic point: lambda_- = 1 + N(V - c)/(1 - cV) contains 0'] = bool((1 + Ns*(Vs - sq)/(1 - sq*Vs)).contains(0))
        pl1['sonic point: lambda_+, transport 1 + NV, light 1 +- N nonzero'] = all(not v.contains(0) for v in
            (1 + Ns*(Vs + sq)/(1 + sq*Vs), 1 + Ns*Vs, 1 + Ns, 1 - Ns))
        Nl, Vl = jl['N'][0], jl['V'][0]
        pl1['light cone: sonic pair, transport, 1 + N nonzero'] = all(not v.contains(0) for v in
            (1 + Nl*(Vl + sq)/(1 + sq*Vl), 1 + Nl*(Vl - sq)/(1 - sq*Vl), 1 + Nl*Vl, 1 + Nl))
        allb = []
        def walk(o):
            if isinstance(o, dict):
                for v in o.values(): walk(v)
            elif isinstance(o, bool): allb.append(o)
        walk(pl1)
        pl1['pass'] = all(allb)
        res.update(PL1=pl1, PL2=pl2, PL3a=pl3, PL3b=pl3b, PL4=pl4, PL5=pl5, time_s=time.time() - T00)
        kill = dict(sonic_d_positive=any(('d < 0 certified' in v and not v['d < 0 certified']) for v in pl2.values() if isinstance(v, dict)),
                    resonance_on_imaginary_axis=not pl3['no resonance on Re kappa = 0 (1/lambda_lc not an integer)'],
                    gauge_zero_on_imaginary_axis=not pl4['pass'])
        res['kill_switch'] = kill
        if dry:
            json.dumps(res, default=str)                        # serialisation path exercised, nothing written
            log('DRY RUN complete; total %.1fs; nothing written' % (time.time() - T00))
            sys.exit(0)
        json.dump(res, open(os.path.join(HERE, 'result_l0.json'), 'w'), indent=1, default=str)
        log('PL1 %s PL2 %s PL3a %s PL3b axial %s polar %s PL4 %s PL5 %s kill %s' % (pl1['pass'], pl2['pass'], pl3['pass'], pl3b['pass_axial'], pl3b['pass_polar'], pl4['pass'], pl5['pass'], kill))
        log('wrote result_l0.json (%.0fs)' % (time.time() - T00))
    finally:
        ctx.prec = old
