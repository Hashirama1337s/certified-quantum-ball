"""build-4: the BG2 tube must overlap the E5 tube on every step.  Each record step stores the step-start set as
c (point) and setrad: the true background state at the step start lies in c +- setrad (lin_tube_v2 record()).
Check: identical step partition (kind, h overlap), and for every step and component |c_new - c_old| <= setrad_new +
setrad_old (balls overlap); also the sonic-manifold and centre-manifold coefficient balls overlap.  -> tube_compare.json"""
import pickle, json
from flint import arb, ctx
ctx.prec = 300
de = lambda p: arb(p[0]) + arb(0, arb(p[1]))
o = pickle.load(open('<LAB>/engine/cert/lin_tube2_k1_3.pkl', 'rb'))
nw = pickle.load(open('<LAB>/quantum/seats/build-4/lin_tube_bg2_k1_3.pkl', 'rb'))
res = dict(n_old=len(o['steps']), n_new=len(nw['steps']))
same_part = len(o['steps']) == len(nw['steps']) and all(a['kind'] == b['kind'] and de(a['h']).overlaps(de(b['h']) + arb(0, arb('1e-15')))
                                                         for a, b in zip(o['steps'], nw['steps']))
res['same_partition'] = same_part
bad = []; maxsr_old = arb(0); maxsr_new = arb(0); worst = 0.0
if same_part:
    for i, (a, b) in enumerate(zip(o['steps'], nw['steps'])):
        ra, rb = arb(a['setrad'][0]), arb(b['setrad'][0])
        maxsr_old = maxsr_old.max(ra); maxsr_new = maxsr_new.max(rb)
        for ca, cb in zip(a['c'], b['c']):
            d = abs(arb(ca[0]) - arb(cb[0]))
            if not d <= ra + rb:
                bad.append(i)
            worst = max(worst, float(d/(ra + rb)))
res['steps_not_overlapping'] = bad
res['max_setrad_old'] = maxsr_old.str(3); res['max_setrad_new'] = maxsr_new.str(3)
res['max_|dc|/(setrad_old+setrad_new)'] = worst
for nm in ('sm', 'cm'):
    ok = True; n = 0
    for i in range(3):
        for a, b in zip(o[nm]['Y'][i], nw[nm]['Y'][i]):
            n += 1
            if not de(a).overlaps(de(b)): ok = False
    for a, b in zip(o[nm]['X'], nw[nm]['X']):
        n += 1
        if not de(a).overlaps(de(b)): ok = False
    res[nm + '_coeffs_overlap'] = ok; res[nm + '_ncoeffs'] = n
res['sm_coeff1_rad_old_new'] = [o['sm']['Y'][0][1][1], nw['sm']['Y'][0][1][1]]
res['sm_coeff5_rad_old_new'] = [o['sm']['Y'][0][5][1], nw['sm']['Y'][0][5][1]]
res['end_overlap'] = all(de(a).overlaps(de(b)) for a, b in zip(o['end'], nw['end']))
res['end_rad_old_new'] = [[a[1] for a in o['end']], [b[1] for b in nw['end']]]
print(json.dumps(res, indent=1))
json.dump(res, open('<LAB>/quantum/seats/build-4/tube_compare.json', 'w'), indent=1)
