"""build-3: assemble result_qb2.json (PQ2-1..PQ2-8) and controls_qb2.json (+ instrument sha256) from final.json,
stages.json, pass_ckpt.json, n0check.json, n1_rows.log.  usage (memcap 8): py -3.11 make_outputs.py"""
import json, hashlib, os
HERE = '<LAB>/quantum/seats/build-3/'
F = json.load(open(HERE + 'final.json'))
STG = json.load(open(HERE + 'stages.json'))
CK = json.load(open(HERE + 'pass_ckpt.json'))
N0 = json.load(open(HERE + 'n0check.json'))
n1 = open(HERE + 'n1_rows.log').read()
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
mine = ['bachsym.py', 'bach_u.py', 'bach_u.txt', 'bach_exprs.txt', 'n1_rows.py', 'n1_rows.txt', 'hka_copy.py', 'qb2sym2.py',
        'qb2_exprs.pkl', 'qb2lib.py', 'qb2core.py', 'qb2stages.py', 'qb2pass.py', 'qb2final.py', 'n0check.py', 'make_outputs.py',
        'run_pass.sh', 'a3.py', 'c2check.py']
reused = ['<LAB>/engine/cert/count_b4v3.py', '<LAB>/engine/cert/count_b4.py',
          '<LAB>/engine/cert/count_b4v3_seal.py', '<LAB>/engine/cert/lin_stages_v2.py',
          '<LAB>/engine/cert/lin_stages.py', '<LAB>/engine/cert/lin_tube2_k1_3.pkl',
          '<LAB>/engine/cert/lin_fields2.pkl', '<LAB>/engine/cert/lin_mc.py',
          '<LAB>/engine/cert/taylor.py', '<LAB>/engine/cert/core.py',
          '<LAB>/engine/cert/pipeline.py', '<LAB>/engine/cert/fields.pkl',
          '<LAB>/engine/system.py', '<LAB>/engine/eqs_gen.py',
          '<LAB>/b4speed/L13_nsub2.pkl', '<LAB>/b4speed/cacheL.py',
          '<LAB>/quantum/seats/build-1/result_qb1.json']
ex = CK['extra']
res = {
 'mapping': 'PQ2-1 = N0 + N1 (source), PQ2-2 = N2, PQ2-3 = N3, PQ2-4 = N4a + N4, PQ2-5 = N5, PQ2-6 = N6, PQ2-7 = N7, PQ2-8 = N8',
 'PQ2-1': dict(N1=dict(method='undetermined coefficients (n1_rows.py), independent of r2_forced.py',
                       rows_equal_K5=[('N1 row %d == r2_forced_k13.log K5 row %d : True' % (i, i)) in n1 for i in range(4)],
                       sigma_I_equal_K5=('N1 sigma_I == K5 l.25 : True' in n1),
                       rows01_vs_srcprop='S_0 = -(1/A) srcprop l.13, S_1 = -(1/N) srcprop l.14 (derive-variable normalisation Ap = A Abar_p, Np = N Nbar_p), exact',
                       frame_forms_true=('N1 S_0 frame == A beta_nn : True' in n1 and 'N1 sigma_I frame == A(-beta_nn + N beta_nr) : True' in n1)),
               N0=dict(method='b_ab = Bach of gt (bachsym.py, exact; conformal weight -2) as certified Taylor series along every transport '
                              'substep (ball jets + mean-value form near the sonic point), sonic theta-series, centre theta-series',
                       qb1_cross_check=N0, traceless_symbolic='b_nn - b_rr - 2 b_T == 0 identically (bachsym.log)',
                       trace_note='crude complex-cell evaluation of the trace wraps (max %s); not informative, identity is exact' % ex['max']['tr'])),
 'PQ2-2': dict(N2=dict(symbolic='D_t == 0 and (div Sigma)_x == 0 identically for Bach of any gt (bachsym.log)',
                       tube_Dt_series_coeffs_contain_0=ex['dt_all_ok'], tube_max_abs_Dt_bound=ex['max']['dt'],
                       sonic=dict(contain_0=STG['sonic_info']['Dt_sonic_coeffs_contain_0'], max_abs=STG['sonic_info']['Dt_sonic_coeffs_max_abs']),
                       range='[-4.6, 0] tube cells + sonic theta-series + centre chart (x < -4.6, via b_low/N4a); outward (0, 12] NOT covered (lever 1)')),
 'PQ2-3': dict(N3=STG['sonic_info'], vectors={k: v.get('tail') for k, v in STG['sonic'].items()}),
 'PQ2-4': dict(N4a=STG['centre_info'].get('N4a_b_low'), N4a_rows=STG['centre_info'].get('N4a_rows_orders01'),
               N4={k: v for k, v in STG['centre_info'].items() if k not in ('N4a_b_low', 'N4a_rows_orders01')}),
 'PQ2-5': dict(N5=dict(substeps=CK['k'], max_remainder=ex['max']['rem'], max_order=ex['max']['Ks'], max_a=ex['max']['aaug'])),
 'PQ2-6': dict(N6=F.get('N6'), matching=F.get('matching'), det=F.get('N6_matrix_det_D_HKA2'), det_overlaps_pq7b=F.get('N6_det_overlaps_pq7b'),
               x_m=F.get('x_m'), c_V=F.get('c_V_xm')),
 'PQ2-7': dict(N7=F.get('N7')),
 'PQ2-8': dict(N8=F.get('N8')),
}
C2J = json.load(open(HERE + 'c2.json'))
ctl = dict(C1=F.get('C1'), C2=dict(final=F.get('C2'), c2check=C2J), C6=F.get('C6'), C8=F.get('C8'), C8prime=F.get("C8'"), N_a=F.get('N-a'), N_c=F.get('N-c'),
           gate=F.get('gate_N1_C1_C8'),
           instrument_sha256={f: sha(HERE + f) for f in mine if os.path.exists(HERE + f)},
           reused_sha256={f: sha(f) for f in reused if os.path.exists(f)})
json.dump(res, open(HERE + 'result_qb2.json', 'w'), indent=1, default=str)
json.dump(ctl, open(HERE + 'controls_qb2.json', 'w'), indent=1, default=str)
print('wrote result_qb2.json, controls_qb2.json')
