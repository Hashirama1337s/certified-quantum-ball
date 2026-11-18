import sys, time; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
t0=time.time()
from qb2_bg import *
print('import', time.time()-t0, flush=True)
print('al1', sal1); print('A1', sA1)
print('M', M); print('b', b)
y=[1.0,1.0,1.0,-0.05]
t0=time.time()
for i in range(1000): rhs(0.1,y)
print('1000 rhs', time.time()-t0)
