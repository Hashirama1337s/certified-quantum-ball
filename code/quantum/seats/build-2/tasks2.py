"""build-2 QB1b items PQ6b (real off-shell null), PQ8 (seam), controls through the extended code path.
usage (memcap 8): py -3.11 tasks2.py null|seam|controls   -> pq6b.json / pq8.json / controls_tm.json
"""
import sys, os, json, math, time
from flint import arb, ctx
ctx.cap = 64
sys.path.insert(0, '<LAB>/quantum/seats/build-1')
import tm
from tm import TM
tm.install()
import bg, qb1, ctl
import charts as CH
import ajet as J
from ajet import Jet, n4
from aser import Ser, revert
import tmcell
from schart import schart, profiles_direct
HERE = os.path.dirname(os.path.abspath(__file__))
S = qb1.S
COMP = qb1.COMP


def kap():
    return arb(1)/3


# ------------------------------------------------------------------ PQ6b
def points_profiles(Bg):
    """the 28 QB1 sample points (build-1 qb1.route_points list), as (tag, x0, prof, chart)"""
    rows = []
    for xt in (-12, -9, -6, '-4.5'):
        th = arb(math.exp(float(arb(xt).mid()) - 0.212365646758945385))
        x0, prof = Bg.centre_profiles_theta(th, qb1.K)
        for ch in ('polar', 'cart'):
            rows.append(('centre chart x=%s' % xt, x0, prof, ch))
    for xt in ('-4.5', '-3', '-2', '-1.5', '-1', '-0.84', '-0.5', '-0.2489'):
        x = arb(xt)
        ybox, (cN, co, cV) = Bg.point_series('tube', x, qb1.K)
        prof, _ = bg.profiles_from_NomV(x, cN, co, cV, qb1.K)
        rows.append(('tube x=%s' % xt, x, prof, 'polar'))
        if xt in ('-3', '-1'):
            rows.append(('tube x=%s' % xt, x, prof, 'cart'))
    x0, cN, co, cV = Bg.sonic_series('bg2', arb(0), qb1.K)
    prof, _ = bg.profiles_from_NomV(x0, cN, co, cV, qb1.K)
    rows.append(('sonic x=0', x0, prof, 'polar'))
    for xt in ('0.2489', '0.5', '1', '1.357941265765732', '2', '3', '6', '9', '12'):
        x = arb(xt)
        ybox, (cN, co, cV) = Bg.point_series('outw', x, qb1.K)
        prof, _ = bg.profiles_from_NomV(x, cN, co, cV, qb1.K)
        rows.append(('outward x=%s' % xt, x, prof, 'polar'))
    return rows


def run_chart(x0, prof, ch):
    if ch == 'cart':
        return CH.cart(x0.exp(), 'x', x0, *prof, kap())
    return CH.polar(x0, *prof, kap())


def pq6b(log):
    Bg = bg.Background()
    fac = 1 + arb('1e-12')
    out = []
    for tag, x0, prof, ch in points_profiles(Bg):
        p2 = [prof[0], prof[1], [c*fac for c in prof[2]], prof[3]]     # omega -> omega (1 + 1e-12); A, N, V kept
        o = run_chart(x0, p2, ch); o0 = run_chart(x0, prof, ch)
        E = {k: o[k] for k in ('E_nn', 'E_nr', 'E_rr', 'E_T')}
        dis = {c: bool(not o['a_' + c].overlaps(o['b_' + c])) for c in COMP}
        row = dict(tag=tag, chart=ch, x=S(x0), einstein=({k: S(v) for k, v in E.items()}),
                   einstein_excludes_0=any(not v.contains(0) for v in E.values()),
                   routeA={c: S(o['a_' + c]) for c in COMP}, routeB={c: S(o['b_' + c]) for c in COMP},
                   AB_disjoint=dis, AB_disjoint_any=any(dis.values()),
                   routeB_unchanged=all(o['b_' + c].overlaps(o0['b_' + c]) for c in COMP),
                   certified_einstein_contains_0=all(o0[k].contains(0) for k in E))
        out.append(row)
        log('%-26s %-5s Eres excl0 %s  A/B disjoint %s' % (tag, ch, row['einstein_excludes_0'], [c for c, v in dis.items() if v]))
    summ = dict(n_points=len(out), einstein_excludes_0_anywhere=any(r['einstein_excludes_0'] for r in out),
                einstein_excludes_0_count=sum(r['einstein_excludes_0'] for r in out),
                AB_disjoint_anywhere=any(r['AB_disjoint_any'] for r in out),
                AB_disjoint_count=sum(r['AB_disjoint_any'] for r in out),
                AB_disjoint_components=sorted(set(c for r in out for c, v in r['AB_disjoint'].items() if v)),
                routeB_unchanged_everywhere=all(r['routeB_unchanged'] for r in out),
                certified_background_einstein_contains_0_everywhere=all(r['certified_einstein_contains_0'] for r in out),
                perturbation='omega profile (Wt = om/z^2 series, all orders) multiplied by (1 + 1e-12) after integration; F2 (A N^2 z^2), At ((A-1)/z^2), Vt unchanged')
    return dict(summary=summ, rows=out)


# ------------------------------------------------------------------ PQ8 seam
def pq8(log):
    Bg = bg.Background()
    x = arb(-9)/2                                     # exact x = -9/2
    tm.H = arb(0)
    g = tmcell.guess_theta(Bg, 'centre', x)
    profc, info = tmcell.centre_cell(Bg, x, arb(0), g)
    pc = [[tm.rng(c) for c in l] for l in profc]
    oc = schart(*pc, kap())
    rec_ybox, (cN, co, cV) = Bg.point_series('tube', x, 5)
    pt = profiles_direct(cN, co, cV, 5, Ser, kap())
    ot = schart(*pt, kap())
    # build-1 polar chart for both, too
    profc_b1 = [[tm.rng(c) for c in l] for l in profc]
    # convert direct profiles (G, A, om, V) -> build-1 (F2, At, Wt, Vt) at exact x
    def to_b1(p):
        G, A_, om, V = (Ser(l, 5) for l in p)
        ex = Ser([arb(0), arb(1)], 5).exp()*x.exp(); e2 = ex*ex
        return [(G*e2).c, ((A_ - 1)/e2).c, (om/e2).c, (V/ex).c]
    opc = CH.polar(x, *to_b1(pc), kap()); opt = CH.polar(x, *to_b1(pt), kap())
    res = dict(x='-9/2 (exact)', theta=S(tm.rng(info['L'])) if False else None,
               schart=dict(centre_b_nn=S(oc['b_nn']), tube_b_nn=S(ot['b_nn']), overlap=bool(oc['b_nn'].overlaps(ot['b_nn'])),
                           all_components_overlap=all(oc['b_' + c].overlaps(ot['b_' + c]) for c in COMP),
                           centre={c: S(oc['b_' + c]) for c in COMP}, tube={c: S(ot['b_' + c]) for c in COMP}),
               polar_build1=dict(centre_b_nn=S(opc['b_nn']), tube_b_nn=S(opt['b_nn']), overlap=bool(opc['b_nn'].overlaps(opt['b_nn']))),
               cross_chart_overlap=bool(oc['b_nn'].overlaps(opt['b_nn']) and ot['b_nn'].overlaps(opc['b_nn'])),
               theta_solve=dict(L=info['L'].str(5), residual=info['res'].str(5)))
    res.pop('theta')
    log('seam x=-9/2: centre %s  tube %s  overlap %s' % (oc['b_nn'].str(15), ot['b_nn'].str(15), res['schart']['overlap']))
    return res


# ------------------------------------------------------------------ controls (extended path)
def controls_tm(log):
    out = {}
    tm.P = 5
    # (a) build-1 controls verbatim under tm.install() (arb inputs: ArbW == flint.arb)
    import qb1 as Q
    c1 = Q.controls(log)
    out['build1_controls_under_install'] = dict(all_positive_controls_pass=c1['all_positive_controls_pass'],
                                                signflip_null_refuses=c1['signflip_null_refuses'])
    # (b) Hayward with a TM base point u = eps, |eps| <= 1/16: V_mn(u) = e^{2u} V_mn(0) (homothety, Zahn eq. 6)
    def hay(flip):
        tm.H = arb(1)/16
        J.setup(3, 4)
        u0 = TM.var(arb(0))
        d1, d2 = Jet.var(1), Jet.var(2)
        sn = (d2 + arb.pi()/2).sin()
        e = ((d1 + u0)*(-2)).exp()
        g = [[Jet() for _ in n4] for _ in n4]
        g[0][0] = e*(-2); g[1][1] = e*2; g[2][2] = e; g[3][3] = e*sn*sn
        d = J.curv(g); T1, T2 = J.zahn(d, flip=flip); V = J.Vxi(T1, T2, 0)
        sc = (u0*2).exp()
        exp_ = [arb(1)/80, arb(17)/240, arb(-7)/480]
        return [bool((V[i][i].c0() - sc*exp_[i]).contains(0)) for i in range(3)], [tm.rng(V[i][i].c0()).str(10) for i in range(3)]
    ok, vals = hay(False); okf, valsf = hay(True)
    out['C0_hayward_TM_u_in_[-1/16,1/16]'] = dict(contains=ok, values=vals, pass_=all(ok))
    out['NULL_signflip_hayward_TM'] = dict(contains=okf, values=valsf, refuses=not all(okf))
    # (c) RN M=1, Q=1/2 at r = 3 + eps, |eps| <= 1/32: closed forms (design-1) as TMs, Bach = 60 V, Route A = V
    tm.H = arb(1)/32
    J.setup(3, 4)
    R0 = TM.var(arb(3))
    d1, d2 = Jet.var(1), Jet.var(2)
    sn = (d2 + arb.pi()/2).sin()
    r = d1 + R0
    f = 1 - 2/r + (r*r).inv()*arb('0.25')
    g = [[Jet() for _ in n4] for _ in n4]
    g[0][0] = -f; g[1][1] = f.inv(); g[2][2] = r*r; g[3][3] = r*r*sn*sn
    d = J.curv(g); T1, T2 = J.zahn(d); Bc = J.bach(T1, T2)
    R = R0
    ex = [(R*R*4 - R*8 + 1)**2/(R**10*640), -1/(R**6*120), (R*R*4 - R*8 + 1)/(R**6*240)]
    rn = {}
    for nm, xi in (('xi0', 0), ('xi1/6', arb(1)/6)):
        V = J.Vxi(T1, T2, xi)
        rn[nm] = [bool((V[i][i].c0() - ex[i]).contains(0)) for i in range(3)]
    rn['Bach_eq_60V'] = [bool((Bc[i][i].c0() - ex[i]*60).contains(0)) for i in range(3)]
    VA = J.routeA(d, d['Ric'])
    rn['routeA_eq_V'] = [bool((VA[i][i].c0() - ex[i]).contains(0)) for i in range(3)]
    rn['pass_'] = all(rn['xi0'] + rn['xi1/6'] + rn['Bach_eq_60V'] + rn['routeA_eq_V'])
    out['C3_RN_TM_r_in_[95/32,97/32]'] = rn
    # (d) Minkowski (KHA N = e^{-x}) through schart with x = x_m + eps, |eps| <= 1/64: B = 0
    tm.H = arb(1)/64
    K = 5
    xs = Ser.var(TM.var(arb('0.25')), K)
    Nn = (xs*(-1)).exp()
    one = Ser.const(1, K); zero = Ser.const(0, K)
    o = schart((one*Nn*Nn).c, one.c, zero.c, zero.c, kap())
    out['C2_minkowski_schart_TM'] = dict(b={c: tm.rng(o['b_' + c]).str(5) for c in COMP},
                                         pass_=all(o['b_' + c].contains(0) for c in COMP) and o['R'].contains(0))
    # (e) Friedmann k = 1/2 through schart, z = z_m + eps (frw.py closed form), r^4 V_uu(0) = 2 z^4/243, Bach = 0
    fr = []; allp = True
    for zm in ('0.5', '1', '1.5'):
        tm.H = arb(1)/512
        zT = TM.var(arb(zm))
        q = arb(4)/9; c = q*(1 - q)
        zs = Ser.var(zT, K)
        Vz = zs*(-q); Nz = (1 + zs*zs*c)/zs; Az = (1 - zs*zs*q*q).inv(); oz = zs*zs*q*q*arb('1.5')/(1 - zs*zs*q*q)
        xz = zs.log() - (1 + zs*zs*c).log()/(2*(1 - q))
        gi = revert(xz.shift0())
        conv = lambda F_: F_.compose_into(gi)
        G = Az*Nz*Nz
        o = schart(conv(G).c, conv(Az).c, conv(oz).c, conv(Vz).c, arb(1)/2, want_routeA=False)
        pred = zT**4*2/243
        ok = (bool((o['V0_uu'] - pred).contains(0)) and o['V6_uu'].contains(0) and all(o['b_' + cc].contains(0) for cc in COMP)
              and bool((o['R'] + zT*zT*8/27).contains(0)) and all(o[e].contains(0) for e in ('E_nn', 'E_nr', 'E_rr', 'E_T'))
              and not o['xi_diff_uu'].contains(0))
        allp &= ok
        fr.append(dict(z=zm + ' +- 1/512', V0_uu=tm.rng(o['V0_uu']).str(10), pred=tm.rng(pred).str(10), V6_uu=tm.rng(o['V6_uu']).str(3), pass_=ok))
    out['C6_friedmann_k1/2_schart_TM'] = dict(rows=fr, pass_=allp)
    out['all_extended_positive_controls_pass'] = (out['C0_hayward_TM_u_in_[-1/16,1/16]']['pass_'] and rn['pass_'] and
                                                  out['C2_minkowski_schart_TM']['pass_'] and allp and
                                                  out['build1_controls_under_install']['all_positive_controls_pass'])
    out['signflip_null_refuses_TM'] = out['NULL_signflip_hayward_TM']['refuses']
    log('extended controls pass %s, sign-flip refuses %s' % (out['all_extended_positive_controls_pass'], out['signflip_null_refuses_TM']))
    return out


if __name__ == '__main__':
    mode = sys.argv[1]
    T0 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    old = ctx.prec
    try:
        ctx.prec = 256
        res = dict(null=pq6b, seam=pq8, controls=controls_tm)[mode](log)
    finally:
        ctx.prec = old
    res['time_s'] = round(time.time() - T0, 1)
    fn = dict(null='pq6b.json', seam='pq8.json', controls='controls_tm.json')[mode]
    json.dump(res, open(os.path.join(HERE, fn), 'w'), indent=1)
    log('wrote', fn)
