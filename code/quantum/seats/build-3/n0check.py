"""build-3 N0 cross-check: b_ab from the exact Bach Tape (qb2lib) versus the QB1 instrument's certified values
(seats/build-1/result_qb1.json: sonic point, centre c0, tube route rows).  Overlap of balls is required.
Also fixes the Bach normalisation convention (ratio printed)."""
import sys, json, pickle
sys.path.insert(0, '<LAB>/quantum/seats/build-3')
import qb2lib as Q
from flint import arb, ctx
from core import sym
import taylor as TY
from param import _c as _cc
import driver_v2 as V2
ctx.prec = 200
R1 = json.load(open('<LAB>/quantum/seats/build-1/result_qb1.json'))
tb = pickle.load(open('<LAB>/engine/cert/lin_tube2_k1_3.pkl', 'rb'))
de = lambda p: arb(p[0]) + sym(arb(p[1]))
ball = lambda p: arb(p[0]) + sym(arb(p[1]))
out = {}
# --- sonic point (tube sm, theta = 0)
sm = tb['sm']
Y = [[de(p) for p in sm['Y'][i]] for i in range(3)]; X = [de(p) for p in sm['X']]
n = 40
xp = Q.s_deriv(X)
r = Q.bach_from_series(Y[0][:n], Y[1][:n], Y[2][:n], xp, n)
ref = R1['sonic']
for i, nm in enumerate(('nn', 'nr', 'rr', 'T')):
    mine = r['b'][i][0]; qb = ball(ref['b_' + nm])
    out['sonic_b_' + nm] = dict(mine=mine.str(20), qb1=qb.str(20), overlap=bool(mine.overlaps(qb)), ratio=(mine/qb).str(8))
# --- centre c0 = lim b_nn e^{-4x}
cm = tb['cm']
Z = [[de(p) for p in cm['Y'][i]] for i in range(3)]; chi = [de(p) for p in cm['X']]
xc = de(tb['xc'])
n = 30
Pz, Wz, Uz = [z[:n + 1] for z in Z]
assert Pz[0].is_zero()
Qs = Pz[1:n + 1]                              # P/theta
fac = Q.s_add([arb(1)] + [arb(0)]*(n - 1), Q.s_up(Q.s_deriv(chi[:n + 1]))[:n])   # 1 + theta chi'
Dc = lambda s, m: Q.s_div(Q.s_up(Q.s_deriv(s))[:m], fac, m)
A = Q.series_out(Q.T_ACEN, [Pz[:n], Wz[:n], Uz[:n]], n, arb(0), Q.conv_for(arb(0)))[0]
# l' = d/dx ln N = -(theta d/dtheta ln P)/(1+theta chi') = -(1 + theta Q'/Q)/fac
tq = Q.s_up(Q.s_deriv(Qs))[:n]
l1 = Q.s_scal(Q.s_div(Q.s_add([arb(1)] + [arb(0)]*(n - 1), Q.s_div(tq, Qs, n)), fac, n), -1)
Aj = [A]; lj = [l1]
for j in range(4): Aj.append(Dc(Aj[-1], n))
for j in range(3): lj.append(Dc(lj[-1], n))
b = Q.series_out(Q.T_BACH, [Pz[:n]] + Aj[:5] + lj[:4], n, arb(0), Q.conv_for(arb(0)))
print('centre low coefficients of b_nn:', [b[0][k].str(5) for k in range(6)])
print('centre low coefficients of b_nr:', [b[1][k].str(5) for k in range(7)])
c0 = b[0][4]*(-4*xc).exp()
qb = ball(R1['centre']['c0'])
out['centre_c0'] = dict(mine=c0.str(20), qb1=qb.str(20), overlap=bool(c0.overlaps(qb)))
out['centre_low_coeffs_contain_0'] = dict(
    nn=all(b[0][k].contains(0) for k in range(4)), rr=all(b[2][k].contains(0) for k in range(4)),
    T=all(b[3][k].contains(0) for k in range(4)), nr=all(b[1][k].contains(0) for k in range(5)),
    nn_max_rad=max(b[0][k].rad() for k in range(4)).str(3), nr_max_rad=max(b[1][k].rad() for k in range(5)).str(3))
# --- tube route points: state from the record x-step (ball start c +- 1.01 setrad, jets, tail), ODE x-jets
Ftape = Q.T_F
xs = tb['X_SWITCH']; x0 = de(xs)
steps = tb['steps']
want = {}
for row in R1['routes']:
    if row['source'].startswith('tube x-step') and row['chart'] == 'polar':
        want[row['x'][0][:8]] = row
done = []
KL = 40
for st in steps:
    if st['kind'] != 'x':
        continue
    h = de(st['h'])
    for key, row in want.items():
        xr = arb(row['x'][0])
        t = xr - x0
        if (t <= 0 and t >= h) and key not in done:
            c = [de(p) for p in st['c']]; setrad = de(st['setrad'])
            y0b = [c[i] + sym(setrad*arb('1.01')) for i in range(3)]
            Rc = abs(h)/arb('0.25')
            y1, Bx, Yc = TY.step_cauchy2(Ftape, y0b, h, KL, q=arb('0.25'), conv=lambda q_: _cc(q_, arb(0)))
            tq_ = abs(t)/Rc
            y = [Q.horner(Yc[i][:KL + 1], t) + sym(Bx[i].rad()*tq_**(KL + 1)/(1 - tq_)) for i in range(3)]
            ys = Q.ode_xseries(y, 8)
            rr = Q.bach_from_series(ys[0], ys[1], ys[2], [arb(1)] + [arb(0)]*8, 8)
            res = {}
            for i, nm in enumerate(('nn', 'nr', 'rr', 'T')):
                mine = rr['b'][i][0]; qb = ball(row['B_' + nm])
                res[nm] = dict(mine=mine.str(15), qb1=qb.str(15), overlap=bool(mine.overlaps(qb)))
            out['tube_x=' + key] = res
            done.append(key)
    x0 = x0 + h
print(json.dumps(out, indent=1))
allov = all(v['overlap'] for k, v in out.items() if k.startswith('sonic')) and out['centre_c0']['overlap'] and \
    all(r_['overlap'] for k, v in out.items() if k.startswith('tube') for r_ in v.values())
out['ALL_OVERLAP'] = allov
print('ALL_OVERLAP', allov)
json.dump(out, open('<LAB>/quantum/seats/build-3/n0check.json', 'w'), indent=1)
