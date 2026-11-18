import sys, json; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_bg import *
import numpy as np, warnings
warnings.filterwarnings('ignore')
def side(w0v):
    xe, F, r = shoot(float(w0v))
    return (-1 if (r.t[-1] < 0.9 and F < 0) else +1), r.t[-1], F
lo, hi = 17.75, 18.0
for it in range(45):
    mid = 0.5*(lo + hi)
    s, xe, F = side(mid)
    if s < 0: lo = mid
    else: hi = mid
    print(it, repr(lo), repr(hi), xe, F, flush=True)
json.dump({'w0_lo': lo, 'w0_hi': hi}, open('bisect.json', 'w'))
