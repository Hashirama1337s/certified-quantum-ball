"""B4 build seat 2 smokes (NON-COUNTING) on Linear2 (lin_stages_v2) + LinearTM2 (lin_tm_v2).
usage: lin_b4_smoke3.py k parts [method]
  parts: SET (setup report), RT (remainder vs r at |kappa| = 5, 20, 40), C1, C2, C3, N1, N2, N3, KS, UP (upper-path winding on a box)
  env LB_TM_R (disc radius for N1-N3), LB_RT ('kc_re,kc_im,r;...'), LB_UP ('R1')"""
import sys, os, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, acb
from core import sym
import lin_stages_v2 as LS2
import lin_tm as LT
import lin_tm_v2 as LT2
import count_b4 as CB
k = sp.Rational(sys.argv[1]); parts = sys.argv[2].split(',')
method = sys.argv[3] if len(sys.argv) > 3 else 'c1'
R_ = float(os.environ.get('LB_TM_R', '0.5'))
T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
NSUB = int(os.environ.get('LB_NSUB', '1')); TOL = float(os.environ.get('LB_TOL', '1e-30'))
MDEG = int(os.environ.get('LB_M', '24')); JDEG = int(os.environ.get('LB_J', '6')); PREC = int(os.environ.get('LB_PREC', '200'))
L = LS2.Linear2(k, log=log, nsub=NSUB, tol=TOL, prec=PREC)
TMe = LT2.LinearTM2(L, m=MDEG, rem_method=method, J=JDEG)
log('settings: nsub %d tol %g m %d J %d prec %d' % (NSUB, TOL, MDEG, JDEG, PREC))
log('tube %s: %d stored steps -> %d transport steps (%d splits), max depth %d; kappa_g = %s  Abar\'(0) = %s  method %s' % (
    os.path.basename(L.tb.get('enc_src', '?')), len(L.tb['steps']), len(L.steps), L.nsplit, max(s['depth'] for s in L.steps), L.kappa_g.str(20), L.abar_x0.str(18), method))
cnt = {'tm': 0, 'pt': 0}
def tmD(a=0, b=1, det=3):
    def f(kc, r):
        cnt['tm'] += 1; return TMe.D_tm(kc, r, a, b, det)
    return f
def box(cx, cy, hx, hy):
    return [acb(cx - hx, cy - hy), acb(cx + hx, cy - hy), acb(cx + hx, cy + hy), acb(cx - hx, cy + hy)]
def wind(name, f, corners, expect, r):
    t0 = time.time()
    try:
        w, nd, ns, wa = LT.winding_tm(TMe, f, corners, r)
        log('%s: winding %d (expect %d) %s  discs %d, segments %d, %.0fs' % (name, w, expect, 'PASS' if w == expect else 'FAIL', nd, ns, time.time() - t0))
    except ArithmeticError as e:
        log('%s: NOT CERTIFIED (%s), %.0fs' % (name, str(e)[:140], time.time() - t0))
def ptD(kap, a=0, b=1):
    cnt['pt'] += 1; return L.D(kap, a=a, b=b)
def bisect_real(a, b, lo, hi, iters):
    flo = ptD(acb(lo), a, b).real; fhi = ptD(acb(hi), a, b).real
    if not ((flo > 0 and fhi < 0) or (flo < 0 and fhi > 0)): return None, (flo, fhi)
    for it in range(iters):
        mid = arb((lo + (hi - lo)*arb('0.381966')).mid())
        fm = ptD(acb(mid), a, b).real
        if fm.contains(0):
            mid = arb((lo + (hi - lo)*arb('0.618034')).mid())
            fm = ptD(acb(mid), a, b).real
            if fm.contains(0): break
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return (lo + hi)/2 + sym((hi - lo)/2), None
if 'SET' in parts:
    big = sorted(L.steps, key=lambda s: -float((s['a0']*s['rho']).mid()))[:5]
    for s in big:
        log('  largest a0*rho: %s %s h=%s rho=%s a0=%s a1=%s depth %d' % (s['kind'], (s['a0']*s['rho']).str(3), s['h'].str(3), s['rho'].str(3), s['a0'].str(3), s['a1'].str(3), s['depth']))
    for kap in [acb(0), acb(1), acb('2.81055255'), acb(0, 20), acb(0, 40), acb(40, 40)]:
        t0 = time.time()
        try:
            d = ptD(kap); log('  point D_new(%s) = %s  (%.1fs)' % (kap.str(6), d.str(8), time.time() - t0))
        except ArithmeticError as e:
            log('  point D_new(%s) refused: %s' % (kap.str(6), str(e)[:100]))
if 'RT' in parts:
    jobs = os.environ.get('LB_RT', '0,5,0.5;0,5,1.0;0,20,0.2;0,20,0.5;0,40,0.1;0,40,0.3;40,40,0.1;40,40,0.3')
    for job in jobs.split(';'):
        re_, im_, r = [float(x) for x in job.split(',')]
        kc = acb(re_, im_); t0 = time.time()
        try:
            tm = TMe.D_tm(kc, r)
            c0 = float(abs(tm.c[0]).mid())
            log('RT %s r=%.2f: c0=%s R=%s ratio R/|c0|=%.2e  |poly|=%s  max wrap %.2e  mean K %.1f  (%.0fs)' % (
                kc.str(5), r, tm.c[0].str(6), tm.R.str(3), float(tm.R.mid())/max(c0, 1e-300), LT.TM(tm.c, arb(0), tm.r).maxnorm().str(3),
                TMe.stats.get('max_wrap', 0), TMe.stats.get('mean_K', 0), time.time() - t0))
        except ArithmeticError as e:
            log('RT %s r=%.2f: refused (%s)' % (kc.str(5), r, str(e)[:120]))
if 'C1' in parts:
    kg = float(L.kappa_g.mid())
    wind('C1a D_HKA around kappa_g (box half-width 0.05)', tmD(1, 0), box(kg, 0, 0.05, 0.05), 1, 0.08)
    wind('C1b D_new around kappa_g', tmD(0, 1), box(kg, 0, 0.05, 0.05), 0, 0.08)
    enc, _ = bisect_real(1, 0, arb(kg) - arb('0.05'), arb(kg) + arb('0.05'), 40)
    log('C1c D_HKA zero enclosure:', enc.str(15) if enc is not None else 'no sign change', ' overlaps closed-form kappa_g:', bool(enc is not None and enc.overlaps(L.kappa_g)))
if 'C2' in parts:
    for kap in [acb('0.5'), acb('1.5'), acb(2, 1), acb(4, -2)]:
        dh = ptD(kap, 1, 0); dn = ptD(kap)
        rr = dh + (kap - acb(L.kappa_g))*dn/acb(L.abar_x0)
        log('C2 at %s: residual %s contains 0: %s (|D_HKA| %s)' % (kap.str(5), rr.str(4), rr.real.contains(0) and rr.imag.contains(0), abs(dh).str(4)))
if 'C3' in parts:
    q = arb(15)/17
    log('C3 closed-form kappa_g on the k=1/2 enclosure:', L.kappa_g.str(20), ' contains 15/17:', L.kappa_g.contains(q))
    for x in ['0.78', '0.86', '0.90', '1.0']:
        log('  C3 point D_HKA(%s) = %s' % (x, ptD(acb(arb(x)), 1, 0).str(8)))
    if os.environ.get('LB_C3_WIND', '1') == '1':
        wind('C3a D_HKA around 15/17', tmD(1, 0), box(15/17, 0, 0.05, 0.05), 1, 0.08)
    enc, ends = bisect_real(1, 0, q - arb('0.05'), q + arb('0.05'), 30)
    log('C3b D_HKA zero enclosure:', enc.str(18) if enc is not None else 'no sign change %s' % [e.str(5) for e in ends], ' contains 15/17:', bool(enc is not None and enc.contains(q)))
if 'N1' in parts:
    wind('N1 shifted box [3.5,8]x[-4,4]', tmD(), box(5.75, 0, 2.25, 4), 0, R_)
if 'N2' in parts:
    kg = acb(L.kappa_g)
    def f2(kc, r):
        tm = tmD(1, 0)(kc, r)
        c0 = kc + kg
        q = arb(r)/abs(c0).lower()
        if not q < 1: raise ArithmeticError('pole too close')
        m = len(tm.c) - 1
        inv = [(-1)**j/(c0**(j + 1)) for j in range(m + 1)]
        Rinv = (1/abs(c0).lower())*q**(m + 1)/(1 - q)
        return tm*LT.TM(inv, Rinv, tm.r)
    wind('N2 D_HKA/(kappa+kappa_g) on [0.1,3.5]x[-1,1]', f2, box(1.8, 0, 1.7, 1), 2, R_)
if 'N3' in parts:
    wind('N3 e_Abar detector on [0.5,3.5]x[-1,1]', tmD(det=0), box(2.0, 0, 1.5, 1), 2, R_)
    wind('N3-ref e_V detector on [0.5,3.5]x[-1,1]', tmD(), box(2.0, 0, 1.5, 1), 1, R_)
if 'KS' in parts:
    enc, _ = bisect_real(0, 1, arb('2.7'), arb('2.9'), 40)
    log('kappa* enclosure (smoke):', enc.str(15) if enc is not None else 'no sign change', ' gamma =', (1/enc).str(15) if enc is not None else '')
if 'UPB' in parts:
    # NULL boxes via the conjugate-symmetric upper path only: box [x0,x1]x[-y,y], function F in {new, hka_over, abar}
    x0, x1, y = [float(v) for v in os.environ['LB_BOX'].split(',')]
    fn = os.environ.get('LB_F', 'new'); expect = int(os.environ.get('LB_EXPECT', '0')); rb = float(os.environ.get('LB_TM_R', '0.3'))
    kg = acb(L.kappa_g)
    def F(kc, r):
        if fn == 'new': return tmD()(kc, r)
        if fn == 'abar': return tmD(det=0)(kc, r)
        tm = tmD(1, 0)(kc, r); c0 = kc + kg; q = arb(r)/abs(c0).lower()
        m = len(tm.c) - 1
        inv = [(-1)**j/(c0**(j + 1)) for j in range(m + 1)]
        return tm*LT.TM(inv, (1/abs(c0).lower())*q**(m + 1)/(1 - q), tm.r)
    def Fpt(x):
        kap = acb(x)
        if fn == 'new': return ptD(kap)
        if fn == 'abar':
            return F(kap, 1e-6).c[0]
        return ptD(kap, 1, 0)/(kap + kg)
    SHARP = os.environ.get('LB_SHARP', '1') == '1'
    _fc = {}
    def Fnorm(z):
        """point value of F in the TM normalisation (det A(kappa) factor, lin_tm.py:196-198), cached"""
        k_ = (z.real.mid().str(30), z.imag.mid().str(30))
        if k_ not in _fc:
            dA = TMe.detA(z)
            if fn == 'new': v = dA*L.D(z)
            elif fn == 'abar': v = dA*L.D(z, detector=0)
            else: v = dA*L.D(z, 1, 0)/(z + kg)
            cnt['pt'] += 1; _fc[k_] = v
        return _fc[k_]
    path = [complex(x1, 0), complex(x1, y), complex(x0, y), complex(x0, 0)]
    t0 = time.time(); total = arb(0); raw = arb(0); nd = 0; ok = True
    edges = [int(v) for v in os.environ.get('LB_EDGES', '0,1,2').split(',')]
    for (e, a, b, r) in CB.disc_plan_path(path, lambda z: rb):
        if e not in edges: continue
        pa, pb = acb(a.real, a.imag), acb(b.real, b.imag); kc = (pa + pb)/2
        try:
            tm = F(kc, r); inc = CB.segment_arg(tm, pa, pb, kc)
        except ArithmeticError:
            inc = None
        nd += 1
        if inc is None:
            log('UPB %s: segment %s -> %s NOT CERTIFIED at r=%g' % (fn, a, b, r)); ok = False; break
        raw += inc
        if SHARP:
            inc = CB.sharpen(inc, lambda: Fnorm(pa), lambda: Fnorm(pb))
        total += inc
    if ok and edges != [0, 1, 2]:
        e0, e1 = Fpt(x1), Fpt(x0)
        log('UPB %s box [%g,%g]x[-%g,%g] edges %s: partial Delta = %s (unsharpened %s); endpoint signs %s %s  (%d discs r=%g, %.0fs)' % (
            fn, x0, x1, y, y, edges, total.str(25), raw.str(10), '+' if e0.real > 0 else '-', '+' if e1.real > 0 else '-', nd, rb, time.time() - t0))
    elif ok:
        e0, e1 = Fpt(x1), Fpt(x0)
        w = total/arb.pi(); n_ = round(float(w.mid())); near = bool((w - n_).abs_upper() < arb('0.25'))
        par = (e0.real > 0) == (e1.real > 0)
        log('UPB %s box [%g,%g]x[-%g,%g]: Delta/pi = %s (unsharpened %s) -> winding %d (expect %d) %s; within 0.25: %s; endpoint signs %s %s, parity check %s; %d discs r=%g, %.0fs' % (
            fn, x0, x1, y, y, w.str(8), (raw/arb.pi()).str(5), n_, expect, 'PASS' if (n_ == expect and near) else 'FAIL', near,
            '+' if e0.real > 0 else '-', '+' if e1.real > 0 else '-', ((-1)**n_ == 1) == par, nd, rb, time.time() - t0))
if 'UP' in parts:
    # conjugate-symmetric count on a small box [0,R1]x[-R1,R1] via the upper path only (smoke)
    R1 = float(os.environ.get('LB_UP', '4'))
    t0 = time.time(); total = arb(0); nd = 0
    D0 = ptD(acb(0)); DR = ptD(acb(R1))
    for (e, a, b, r) in CB.disc_plan(R1):
        pa, pb = acb(a.real, a.imag), acb(b.real, b.imag); kc = (pa + pb)/2
        tm = TMe.D_tm(kc, r); nd += 1
        inc = CB.segment_arg(tm, pa, pb, kc)
        if inc is None:
            log('UP: segment %s -> %s NOT CERTIFIED at r=%g' % (a, b, r)); break
        total += inc
    else:
        w = total/arb.pi()
        log('UP R1=%g: Delta/pi = %s -> N = %d (expect 1); parity sign D(0) D(R1) = %s; %d discs, %.0fs' % (
            R1, w.str(8), round(float(w.mid())), '+' if (D0.real > 0) == (DR.real > 0) else '-', nd, time.time() - t0))
log('TM evaluations %d, point evaluations %d, total %.0fs' % (cnt['tm'], cnt['pt'], time.time() - T0))
