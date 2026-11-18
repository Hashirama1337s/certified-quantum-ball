"""build-5 probe: plain vs tight (build-4 qb2b_tight T1+T2) Linear2 stage-L cost per tube step, and the plain order-0
radius profile of the cached steps (build-4 Lbg2_nsub2.pkl)."""
import os, sys, time, copy
os.environ['QB2B_RUN1'] = '1'
sys.path.insert(0, '<LAB>/quantum/seats/build-4')
import qb2b as B
import lin_stages_v2 as LS2
L = B.engine_b(log=lambda *a: None)
prof = []
for k in []:
    st = L.steps[k]
    prof.append((k, st['kind'], '%.1e' % max(float(st['A0'][0][i, j].real.rad()) for i in range(4) for j in range(4))))
print('plain order-0 radius profile:', prof, flush=True)
tb = L.tb
for mode in ('plain', 'tight'):
    if mode == 'tight':
        import qb2b_tight
    for j0 in (690, 715):
        o = copy.copy(L); o.tb = dict(tb, steps=tb['steps'][j0:j0 + 2])
        t = time.time(); LS2.Linear2._tube_setup(o); dt = (time.time() - t)/2
        r = max(float(s['A0'][0][i, j].real.rad()) for s in o.steps for i in range(4) for j in range(4))
        print(mode, 'tube step', j0, '%.1fs/step' % dt, 'order-0 radius %.1e' % r, flush=True)
import pickle
from flint import arb
tbr = pickle.load(open(B.NEW_TUBE, 'rb'))
for j in (100, 600, 690, 715):
    print('record setrad at tube step', j, tbr['steps'][j]['setrad'][0][:8], flush=True)
