"""build-4: result_qb2b.json (PQ2b-1..PQ2b-5) and controls_qb2b.json (+ sha256).  usage: py -3.11 make_outputs_b.py"""
import json, hashlib, os, glob
H = '<LAB>/quantum/seats/build-4/'
J = lambda f: json.load(open(H + f))
F = J('final.json'); CK = J('pass_ckpt.json'); TC = J('tube_compare.json'); CE = J('centre_exact.json'); STG = J('stages.json')
TH = J('thth.json') if os.path.exists(H + 'thth.json') else None
ex = CK['extra']
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
from flint import arb
psi_old = arb('3.4616358037141') + arb(0, arb('1.15e-7'))
n8 = F.get('N8') if isinstance(F.get('N8'), dict) else {}
psi = arb(n8['Psi2'][0]) + arb(0, arb(n8['Psi2'][1])) if n8 else None
res = {
 'PQ2b-1': dict(item='sharper inward background (BG2 V0 box) + rebuilt homogeneous stage + optimised series tails',
                tube_compare=TC, N6=F.get('N6'), e_width_le_1em10=F.get('N6', {}).get('width_le_1em10'),
                Dt=dict(all_series_coeffs_contain_0=ex['dt_all_ok'], max_abs_bound_over_inward_steps=ex['max']['dt'],
                        le_1em10=float(ex['max']['dt']) <= 1e-10),
                transport=dict(substeps=CK['k'], max_remainder=ex['max']['rem'], order=ex['max']['Ks']),
                matching=F.get('matching'), det=F.get('N6_matrix_det_D_HKA2'), det_overlaps_pq7b=F.get('N6_det_overlaps_pq7b')),
 'PQ2b-2': dict(item='exact centre orders (sympy, centre_exact.py) + certified order-4 coefficient', **CE,
                N4a_numeric_low_coeffs=STG['centre_info'].get('N4a_b_low')),
 'PQ2b-3': dict(item='domain lemma (DOMAIN-LEMMA.md) + C1 at x = -6, -4.6, -2, -0.25, -0.01', C1=F.get('C1'),
                lemma_sha256=sha(H + 'DOMAIN-LEMMA.md')),
 'PQ2b-4': dict(C8=F.get('C8'), C8_thth_residual=TH, C8prime=F.get("C8'"), N_a=F.get('N-a'), N_c=F.get('N-c'), C2=F.get('C2')),
 'PQ2b-5': dict(N7=F.get('N7'), N8=F.get('N8'), C6=F.get('C6'),
                overlaps_QB2_Psi2=bool(psi is not None and psi.overlaps(psi_old)), QB2_Psi2='[3.4616358037141 +/- 1.15e-7]'),
}
mine = sorted(os.path.basename(p) for p in glob.glob(H + '*.py') + glob.glob(H + '*.sh') + glob.glob(H + '*.txt') + glob.glob(H + '*.md'))
b3 = ['qb2lib.py', 'qb2core.py', 'qb2stages.py', 'qb2pass.py', 'qb2final.py', 'bach_u.txt', 'qb2_exprs.pkl', 'n1_rows.log']
reused = ['<LAB>/engine/cert/lin_tube_v2.py', '<LAB>/engine/cert/lin_stages_v2.py',
          '<LAB>/engine/cert/lin_stages.py', '<LAB>/engine/cert/count_b4.py',
          '<LAB>/engine/cert/lin_tube2_k1_3.pkl', '<LAB>/engine/cert/lin_fields2.pkl',
          '<LAB>/engine/system.py', '<LAB>/b4speed/cacheL.py',
          '<LAB>/sideways/build-5/phase2/outward_C.pkl']
ctl = dict(C1=F.get('C1'), C2=F.get('C2'), C6=F.get('C6'), C8=F.get('C8'), C8_thth=TH, C8prime=F.get("C8'"),
           N_a=F.get('N-a'), N_c=F.get('N-c'), gate=F.get('gate_N1_C1_C8'),
           build4_sha256={f: sha(H + f) for f in mine},
           data_sha256={f: sha(H + f) for f in ('lin_tube_bg2_k1_3.pkl', 'Lbg2_nsub2.pkl', 'stages.json', 'pass_ckpt.json', 'final.json')},
           build3_imported_sha256={f: sha('<LAB>/quantum/seats/build-3/' + f) for f in b3},
           reused_sha256={f: sha(f) for f in reused})
json.dump(res, open(H + 'result_qb2b.json', 'w'), indent=1, default=str)
json.dump(ctl, open(H + 'controls_qb2b.json', 'w'), indent=1, default=str)
print('wrote result_qb2b.json, controls_qb2b.json; Psi2 overlaps QB2:', res['PQ2b-5']['overlaps_QB2_Psi2'])
