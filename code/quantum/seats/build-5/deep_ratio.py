"""build-5 PQ2c-1 evidence: order-0 and order-1 radii of the homogeneous step matrices, plain (build-4 run-1 cache) vs
tight (build-4 qb2b_tight T1+T2, rebuilt here) on every substep of tube steps J0..719 (where the e width is generated).
Also checks the T1 overlap assertion (bar kill condition).  -> deep_ratio.json"""
import os, sys, time, copy, json
os.environ['QB2B_RUN1'] = '1'
sys.path.insert(0, '<LAB>/quantum/seats/build-4')
import qb2b as B
import lin_stages_v2 as LS2
L = B.engine_b(log=lambda *a: None)
import qb2b_tight
J0 = int(sys.argv[1]) if len(sys.argv) > 1 else 480
T0 = time.time()
rows = []
for j in range(J0, 720, 8):
    o = copy.copy(L); o.tb = dict(L.tb, steps=L.tb['steps'][j:j + 8])
    LS2.Linear2._tube_setup(o)
    for s_i, st in enumerate(o.steps):
        k = 2*j + s_i
        p = L.steps[k]
        assert p['h'].overlaps(st['h'])
        r = lambda S, n: max(float(S['A0'][n][a, b].real.rad()) for a in range(4) for b in range(4))
        rows.append(dict(k=k, plain0=r(p, 0), tight0=r(st, 0), plain1=r(p, 1), tight1=r(st, 1),
                         a_plain=float(p['a0'].mid()), a_tight=float(st['a0'].mid())))
    print('tube step %d done (%.0fs)' % (j, time.time() - T0), flush=True)
    json.dump(dict(J0=J0, rows=rows, sec=round(time.time() - T0, 1)), open('deep_ratio.json', 'w'))
ratios = [x['tight0']/x['plain0'] for x in rows]
print('order-0 radius ratio tight/plain: min %.3f median %.3f max %.3f' % (min(ratios), sorted(ratios)[len(ratios)//2], max(ratios)))
