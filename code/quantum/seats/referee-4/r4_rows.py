"""referee-4: frame -> coordinate source and the HKA rows 0, 1, sigma_I (PROOF-QB2-LEMMAS-v2 l.55, l.63-67).
Independent: sigma from Sigma(d_s, d_s) etc. with d_s = a r (N n - e_r), d_x = a r e_r (my own derivation, report)."""
import sympy as sp
N, A, bnn, bnr, brr = sp.symbols('N A bnn bnr brr')
# frame components of Sigma in (n, e_r); coordinate vectors d_s = sqrt(A) r (N n - e_r), d_x = sqrt(A) r e_r; r^2 factored out
Sf = sp.Matrix([[bnn, bnr], [bnr, brr]])
ds = sp.Matrix([N, -1]); dx = sp.Matrix([0, 1])
sss = A*(ds.T*Sf*ds)[0]; ssx = A*(ds.T*Sf*dx)[0]; sxx = A*(dx.T*Sf*dx)[0]
print('sigma_ss - page:', sp.expand(sss - A*(N**2*bnn - 2*N*bnr + brr)), ' sigma_sx - page:', sp.expand(ssx - A*(N*bnr - brr)),
      ' sigma_xx - page:', sp.expand(sxx - A*brr))
sg = sp.Matrix([sss, ssx, sxx])
row0 = sp.Matrix([[1/N**2, 2/N**2, 1/N**2]]); row1 = sp.Matrix([[-1/(2*N**2), -1/N**2, (N**2 - 1)/(2*N**2)]])
sI = sp.Matrix([[-1/N**2, (N**2 - 2)/N**2, (N**2 - 1)/N**2]])
print('row0 =', sp.factor((row0*sg)[0]), '| row1 =', sp.factor((row1*sg)[0]), '| sigma_I =', sp.factor((sI*sg)[0]))
# metric check: g_ss = A r^2 (1 - N^2), g_sx = -A r^2, g_xx = A r^2 from eta(n,n) = -1
eta = sp.diag(-1, 1)
print('g_ss/r^2 =', sp.expand(A*(ds.T*eta*ds)[0]), ' g_sx/r^2 =', sp.expand(A*(ds.T*eta*dx)[0]), ' g_xx/r^2 =', sp.expand(A*(dx.T*eta*dx)[0]))
