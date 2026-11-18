"""build-4: N6 / controls / N7 / N8 / C6 on the BG2 tube = the pinned build-3 qb2final.main (imported, patched in
process by qb2b: engine, tube path, HERE = build-4), plus C1 at x = -0.01 (forced sonic series, stages.json).
usage (memcap 8): py -3.11 final_b.py -> build-4/final.json"""
import sys, json
sys.path.insert(0, '<LAB>/quantum/seats/build-4')
import qb2b as B
from qb2b import FN, Q, arb, ctx, sym
from qb2core import NA, des
o = FN.main()
stg = json.load(open(B.B4 + 'stages.json'))
c = stg['C1_x001']
m1b = arb(o['matching']['C1']['b'][0]) + sym(arb(o['matching']['C1']['b'][1]))
old = ctx.prec
try:
    ctx.prec = 128
    yv = [des(c['C1'][i]) + m1b*des(c['hka'][i]) for i in range(4)]
    st = [des(t) for t in c['state']]
    x = arb(-1)/100
    loc = Q.ode_xseries(st, 3)
    A = Q.series_out(Q.T_ACONS, loc, 3, arb(0), Q.conv_for(arb(0)))[0]
    G2 = [A[1]/A[0], loc[0][1]/loc[0][0] + 2, loc[1][1]/loc[1][0], loc[2][1]]
    ex = [Q.phi_val(x)*g for g in G2]
    S = lambda a: [a.mid().str(20, radius=False), a.rad().str(3, radius=False)]
    o['C1']['points']['m001'] = dict(x='-0.01 (theta on the sonic manifold, X(theta) = -0.01)', solver=[S(v) for v in yv],
                                     exact=[S(v) for v in ex], overlap=bool(all(yv[i].overlaps(ex[i]) for i in range(4))))
    o['C1']['pass'] = bool(o['C1']['pass'] and o['C1']['points']['m001']['overlap'])
finally:
    ctx.prec = old
json.dump(o, open(B.B4 + 'final.json', 'w'), indent=1, default=str)
print(json.dumps({k: o.get(k) for k in ('N6', 'N-a', 'N-c', 'C8', "C8'", 'C2', 'gate_N1_C1_C8', 'N8', 'C6')}, indent=1, default=str))
print('C1', {k: v for k, v in o['C1'].items() if k != 'points'}, {p: v['overlap'] for p, v in o['C1']['points'].items()})
