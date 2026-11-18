"""b1_box.py - build-5 phase 2, B1: refined sonic V0 enclosure by the E5 method (build seat, 2026-11-17).

Method = engine/cert/run_v32.py target mode (l.49-53, 77-95), imported unchanged: Problem2 (driver_v31), sealed nu_s rule
(driver_v2.choose_nu, ladder 0.004/0.002/0.001 at box radius 1e-9), MatcherV32 (R8 mean-value sonic head), residual =
MatcherV31.residual, Krawczyk test = driver2.kr_test via driver2.krawczyk_ladder.  What changes (parameters only):
  * centre: one Newton step from the E5 seed (run_v32.py:74-75) with the point-run Jacobian, Xt1 = mid(Xt0 - C R(Xt0));
  * radius ladder: LADDER below (E5: 1e-15 ... 1e-9);
  * PREC / K_t / K_s from the command line (E5: 200 / 36 / 16).  K_c 60, nu_c 0.03 as E5.  (K_s is the sonic-manifold
    order; nu_s keeps the sealed choose_nu rule.)  Probes: K_t 36 -> point residual floor 3.9e-21 (Lohner truncation);
    K_s 16 -> error bound E = 1.67e-24 (point tail eps 1.9e-32 x growth); hence K_t 50, K_s 24.
  * --nonewton: the E5 seed is used as the centre directly (run B_Ks16 measured the Newton step: [+/- 7.9e-30]).
The certified box is K(X) /\ X at the first contracting rung (run_v32.py:96 uses Kx[0] /\ Xb[0] for V0).
Null (--null): driver_v31.content_null along V0 with off = 1e-22 and rad = the B1 rung: Krawczyk must NOT contract.

usage (memcap 8):  py -3.11 b1_box.py <tag> <prec> <Kt> <Ks> [--nonewton] [--null r] [--ladder r1,r2,...]      -> b1_<tag>.json
"""
import sys, os, time, json
os.environ['LB_NGRID'] = '24'
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from flint import arb, arb_mat, ctx
import driver_v31 as V31
import driver_v32 as V32
import driver_v3 as V3
import driver_v2 as V2
import driver2 as D2
import lohner as LO
import param as PM
from param import val
from core import sym
ctx.cap = 64
assert LO.lohner_step is V3._step and V3.MODE == 'qrw', 'weighted-QR Lohner not installed'
assert PM.arb_series is V31.arb_series_exact, 'R7 not installed'
assert V32.MatcherV32.start_parts is not V31.MatcherV31.start_parts, 'R8 not installed'

HERE = os.path.dirname(os.path.abspath(__file__))
tag, PREC, KT, KS = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
NULL_R = sys.argv[sys.argv.index('--null') + 1] if '--null' in sys.argv else None
NEWTON = '--nonewton' not in sys.argv
KC, NUC = 60, '0.03'
NU_LADDER = ['0.004', '0.002', '0.001']; NU_AT = arb('1e-9')
LADDER = sys.argv[sys.argv.index('--ladder') + 1].split(',') if '--ladder' in sys.argv else ['1e-27', '1e-26', '1e-25']
E5_V0 = ('0.1124394013880983362', '3.31e-20')          # runs6/target.log:35 (K box), lin_tube_v2.py:39
k13 = sp.Rational(1, 3)
OM = '0.3256888894230684'
T0 = time.time()


def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)


def s(a):
    return a.str(radius=True) if hasattr(a, 'str') else a


def newton(M, Xt0, out):
    Rp, Ep, dp = M.residual(Xt0[0], Xt0[1], Xt0[2], OM, log=None)
    Rv = [val(r) for r in Rp]
    C = arb_mat([[arb(Rp[i].d[j].mid()) for j in range(3)] for i in range(3)]).inv()
    dX = [sum((C[i, j]*Rv[j] for j in range(3)), arb(0)) for i in range(3)]
    Xt1 = [arb((Xt0[i] - dX[i]).mid()) for i in range(3)]
    out['newton'] = dict(R0=[s(r) for r in Rv], E0=s(Ep), step=[s(d) for d in dX], Xt1=[s(x) for x in Xt1],
                         width_end_point=s(dp['width_end']), time=dp['time'])
    log('point run at E5 seed: |R| %s, E %s, Newton step %s, width_end %s (%.0fs)' % (max(abs(r) for r in Rv).str(3), Ep.str(3),
        max(abs(d) for d in dX).str(3), dp['width_end'].str(3), dp['time']))
    return Xt1


def main():
    old = ctx.prec
    try:
        ctx.prec = PREC
        V0s = arb('0.112439401388098336223763027603135130845'); W0s = arb('8.90132327775564351911721605677568608689')
        xcs = -arb('0.808668956318442829476786547218464827352').log()
        Pb = V31.Problem2(k13, prec=PREC, K_sonic=KS, K_centre=KC, K_ty=KT)
        nu, _ = V2.choose_nu(Pb, V0s + sym(NU_AT), OM, NU_LADDER, KS, log=log)
        M = V32.MatcherV32(Pb, nu, NUC, Ks=KS, x_switch='-0.05', Kt=KT)
        out = dict(tag=tag, prec=PREC, Kt=KT, Ks=KS, nu_s=nu, ladder=LADDER)
        Xt0 = [V0s, W0s, xcs]
        if NEWTON:
            Xt1 = newton(M, Xt0, out)
        else:
            Xt1 = Xt0; out['newton'] = 'skipped (--nonewton): centre = E5 seed run_v32.py:74-75'
        json.dump(out, open(os.path.join(HERE, 'b1_%s.json' % tag), 'w'), indent=1)
        if NULL_R:
            res = V31.content_null(M, Xt1, 0, OM, off='1e-22', rad=NULL_R, log=log)
            out['null'] = {k: (s(v) if hasattr(v, 'str') else v) for k, v in res.items()}
            log('NULL V0 + 1e-22 at rad %s: %s' % (NULL_R, out['null']))
        else:
            L = D2.krawczyk_ladder(M, Xt1, OM, ladder=LADDER, log=log)
            out['ladder'] = [tuple(str(x) for x in t) for t in L['tried']]
            out['point1'] = dict(R=[s(r) for r in L['R']], width_end=s(L['point']['width_end']), E=s(L['point']['Erad']))
            if L['ok']:
                V0box = L['Kx'][0].intersection(L['Xb'][0])
                e5 = arb(E5_V0[0]) + sym(arb(E5_V0[1]))
                out['B1'] = dict(r=L['r'], K=[s(k) for k in L['Kx']], X=[s(x) for x in L['Xb']], V0=s(V0box),
                                 width=s(2*V0box.rad()), inside_E5=bool(e5.contains(V0box)),
                                 W0=s(L['Kx'][1].intersection(L['Xb'][1])), xc=s(L['Kx'][2].intersection(L['Xb'][2])),
                                 Ebox=s(L['box']['Erad']), width_end_box=s(L['box']['width_end']))
                log('B1: V0 = %s width %s inside E5 %s' % (V0box.str(40, radius=True), (2*V0box.rad()).str(3), out['B1']['inside_E5']))
            else:
                out['B1'] = None
                log('B1: no contraction on the ladder')
        out['time_s'] = round(time.time() - T0, 1)
        json.dump(out, open(os.path.join(HERE, 'b1_%s.json' % tag), 'w'), indent=1)
        return 0 if (NULL_R or out.get('B1')) else 3
    finally:
        ctx.prec = old


if __name__ == '__main__':
    sys.exit(main())
