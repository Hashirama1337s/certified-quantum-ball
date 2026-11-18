import sys; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_bg import *
import numpy as np
def probe(w0v, ths=(1e-2, 1e-3, 1e-4)):
    x0 = 1e-4
    y0 = [1 + w0v/6*x0**2, 1 + w0v/3*x0**2, w0v, -0.5*x0]
    evs = []
    for th_ in ths:
        f = (lambda tt: (lambda xx, y: abs(detM(xx, y)) - tt))(th_)
        f.terminal = (th_ == ths[-1]); f.direction = -1
        evs.append(f)
    r = solve_ivp(rhs, [x0, 1.2], y0, method='DOP853', rtol=1e-12, atol=1e-14, events=evs)
    out = []
    for i in range(len(ths)):
        if r.t_events[i].size:
            xe = r.t_events[i][0]; ye = r.y_events[i][0]
            out.append((xe, rhs(xe, ye)[3]))
        else:
            out.append((np.nan, np.nan))
    return out
for w0v in np.linspace(0.1, 0.8, 15):
    o = probe(float(w0v))
    print('w0 %.4f ' % w0v + '  '.join('x %.6f dv %.5g' % t for t in o), flush=True)
