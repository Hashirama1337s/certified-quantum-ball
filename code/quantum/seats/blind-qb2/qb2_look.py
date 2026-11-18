import sys, json; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_bg import *
import numpy as np, warnings
warnings.filterwarnings('ignore')
for w0v in [17.802642822265625, 17.802644729614258]:
    x0 = 1e-4
    y0 = [1 + w0v/6*x0**2, 1 + w0v/3*x0**2, w0v, -0.5*x0]
    ev2 = lambda xx, y: 1e4 - np.abs(rhs(xx, y)).max() if (abs(y[3]) < 1 and y[1] > 0) else -1.0
    ev2.terminal = True
    r = solve_ivp(rhs, [x0, 1.2], y0, method='DOP853', rtol=1e-12, atol=1e-14, events=[ev2], dense_output=True)
    print('w0', w0v, 'end', r.t[-1], r.status)
    for xx in np.linspace(0.05, r.t[-1], 40):
        y = r.sol(xx)
        print('  x %.4f al %.6f A %.6f w %.6f v %+.6f det %+.3e' % (xx, *y, detM(xx, y)))
