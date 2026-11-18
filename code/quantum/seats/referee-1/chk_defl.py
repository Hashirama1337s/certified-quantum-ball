import mpmath as mp
mp.mp.dps = 30
kg = mp.mpf('0.3556992037'); Ap = mp.mpf('-0.198928186801772')
for kap, Dn, Dh in [(0, '-108.52727565146470122', '194.05528271674018860'), (1, '-7.2440885428830988038', '-23.462597692195414844'), (2, '-2.4126270857981028564', '-19.942295267987575083')]:
    pred = -(kap - kg)*mp.mpf(Dn)/Ap
    print(kap, 'D_HKA pred', mp.nstr(pred, 15), 'file', Dh, 'rel diff', mp.nstr((pred - mp.mpf(Dh))/mp.mpf(Dh), 3))
# ratio D_new(x_m=-4.6)/D_new(x_m=-2) per kappa (float values from dnew_qb1.json NONCOUNTING block)
for kap, a, b in [(0, -108.52727565150171, -7.958257902778397), (1, -7.244088542878887, -0.5442148357564793), (2, -2.412627085797027, -0.18568821647419909)]:
    print('kappa', kap, 'ratio -4.6/-2 =', a/b)
mid = mp.mpf('-2.4126270857981028564'); rad = mp.mpf('1.13e-9')
print('rounded -2.41262709: |diff|+rad =', mp.nstr(abs(mid - mp.mpf('-2.41262709')) + rad, 3))
