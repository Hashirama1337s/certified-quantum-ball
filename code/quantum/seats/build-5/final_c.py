"""build-5 PQ2c-3 and PQ2c-4.
PQ2c-3: replay of the matching / controls / readout (build-3 qb2final.main through build-4 qb2b, QB2B_RUN1) on the QB2b
run-1 transported record (build-4 pass_ckpt.json, stages.json; read only), in a fresh process.
PQ2c-4: replication row  Psi2 = 30 pi x_s^2 Psi_blind / A(x_V),  x_s = e^{-xc}  (centre lapse: N e^x -> e^{xc}, A -> 1, so
t_central = e^{xc} t_KHA and the sonic point x = 0 sits at r/(-t_central) = e^{-xc}); xc = E5 box of the cm record
(lin_tube_bg2_k1_3.pkl 'xc'); A(x_V) from N8; Psi_blind from seats/blind-qb2/blind_qb2_result.json.
usage (memcap 8): py -3.11 final_c.py -> final_c.json, conv_c.json"""
import os, sys, json, pickle
os.environ['QB2B_RUN1'] = '1'
sys.path.insert(0, '<LAB>/quantum/seats/build-4')
import qb2b as B
from qb2b import FN, arb, ctx, sym
B5 = '<LAB>/quantum/seats/build-5/'
o = FN.main()
json.dump(o, open(B5 + 'final_c.json', 'w'), indent=1, default=str)
old = ctx.prec
try:
    ctx.prec = 200
    tb = pickle.load(open(B.NEW_TUBE, 'rb'))
    xc = arb(tb['xc'][0]) + sym(arb(tb['xc'][1]))
    xs = (-xc).exp()
    n8 = o['N8']
    psi = arb(n8['Psi2'][0]) + sym(arb(n8['Psi2'][1])); AV = arb(n8['A_xV'][0]) + sym(arb(n8['A_xV'][1]))
    old_ = arb('3.4616358037141656915') + sym(arb('3.02e-12'))
    bl = json.load(open('<LAB>/quantum/seats/blind-qb2/blind_qb2_result.json'))
    S = lambda a: [a.mid().str(20, radius=False), a.rad().str(3, radius=False)]
    rows = {}
    for nm, val in [('headline Psi_per_unit_K_eps', bl['Psi_per_unit_K_eps'])] + list(bl['Psi_runs'].items()):
        conv = 30*arb.pi()*xs**2*arb(val)/AV
        rel = (conv - psi)/psi
        rows[nm] = dict(Psi_blind=val, converted=S(conv), rel_diff=S(rel), le_1e_9=bool(abs(rel).upper() <= arb('1e-9')))
    out = dict(x_s=S(xs), x_s_blind_float=bl['background']['x_sonic'],
               x_s_rel_diff=S((xs - arb(bl['background']['x_sonic']))/xs),
               A_xV=S(AV), A_xV_blind_float=bl['background']['A_peak'], Psi2=S(psi), rows=rows,
               PQ2c4_pass=all(r['le_1e_9'] for r in rows.values()),
               PQ2c3=dict(Psi2=S(psi), overlaps_QB2b=bool(psi.overlaps(old_)), excludes_0=bool(not psi.contains(0)),
                          rel_width=(psi.rad()/abs(psi)).str(3)))
finally:
    ctx.prec = old
json.dump(out, open(B5 + 'conv_c.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
print('N6', o['N6']); print('N-a', o['N-a']['excludes_0'], o['N-a']['e_rowzero_encloses_full_prediction']); print('N-c', o['N-c']); print('C8', o['C8']['e'], o['C8']['e_contains_0'])
