"""referee-3: exact conversion build-3 Psi2 <-> blind-qb2 Psi. Inputs copied from the two result JSONs."""
import json
from mpmath import mp, mpf, pi, log
mp.dps = 30
B = json.load(open('<LAB>/quantum/seats/build-3/result_qb2.json'))['PQ2-8']['N8']
L = json.load(open('<LAB>/quantum/seats/blind-qb2/blind_qb2_result.json'))
Pb, Pbr = mpf(B['Psi2'][0]), mpf(B['Psi2'][1]); AV = mpf(B['A_xV'][0]); Ab2 = mpf(B['Abar2_xV'][0]); xV = mpf(B['x_V'][0])
runs = [mpf(v) for v in L['Psi_runs'].values()]
xs = mpf(L['background']['x_sonic']); xp = mpf(L['background']['x_peak']); tm = mpf(L['background']['two_m_over_r_peak'])
print('x_V(build) - ln(x_peak/x_s) =', xV - log(xp/xs))
print('A(build) - 1/(1-2m/r)(blind) =', AV - 1/(1 - tm))
for P in runs:
    conv = 30*pi*xs**2*P/AV
    print('blind run %.17s -> build units %s ; minus build %s ; rel %s' % (P, mp.nstr(conv, 16), mp.nstr(conv - Pb, 3), mp.nstr((conv - Pb)/Pb, 3)))
print('build Psi2 radius', Pbr, ' rel', mp.nstr(Pbr/Pb, 3))
print('implied factor/pi =', mp.nstr(Pb*AV/(xs**2*runs[2])/pi, 15))
print('Abar2/(30 pi x_s^2) =', mp.nstr(Ab2/(30*pi*xs**2), 15), ' blind N60 =', runs[2])
print('null 5pi^2 instead of 30pi/A:', mp.nstr(5*pi**2*xs**2*runs[2], 12), 'vs', Pb)
print('per unit N_c eps: build delta(2m/r) coeff of s e^{2s} = Psi2/(30pi) =', mp.nstr(Pb/(30*pi), 12))
print('  same in central-proper-time u e^{2u}: Psi_blind/A =', mp.nstr(runs[2]/AV, 12), ' = Psi2/(30pi x_s^2) =', mp.nstr(Pb/(30*pi*xs**2), 12))
print('photons N_c=12, K*Psi2 = 12 Psi2/(30 pi) =', mp.nstr(12*Pb/(30*pi), 12))
