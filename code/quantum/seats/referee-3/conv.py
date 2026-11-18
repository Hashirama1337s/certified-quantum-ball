from mpmath import mp, mpf, pi, log, exp
mp.dps=30
Pb=mpf('3.4616358037141'); Pbr=mpf('1.15e-7')
Pc=mpf('0.10751999497716666'); Pce=mpf('1e-11')
xs=mpf('0.808668956318442829477')
A=mpf('1.91434642526073')
print('raw ratio', Pb/Pc)
conv = 30*pi*xs**2*Pc/A
print('30pi xs^2 Pblind / A =', conv)
print('build                  =', Pb)
print('diff', Pb-conv, 'rel', (Pb-conv)/Pb)
print('implied factor Pb*A/(xs^2 Pc)/pi =', Pb*A/(xs**2*Pc)/pi)
print('x_V check: ln(xpeak/xs)=', log(mpf('0.62803916560721')/xs), 'build x_V=-0.252787102083')
print('A from blind 2m/r:', 1/(1-mpf('0.47762850714765')))
