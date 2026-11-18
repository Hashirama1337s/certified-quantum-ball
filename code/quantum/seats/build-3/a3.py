import sympy as sp
x=sp.Symbol('x'); A=sp.Function('A')(x); N=sp.Function('N')(x)
bn,br,bt=[sp.Function(n)(x) for n in ('bnn','bnr','brr')]
fss=A*(N**2*bn-2*N*br+bt); fsx=A*(N*br-bt); fxx=A*bt; kap=2
Ab=A;Nb=N
Dt = ((1/Ab - 2/(Ab*Nb**2))*sp.diff(fsx, x) + (1/Ab - 1/(Ab*Nb**2))*sp.diff(fxx, x)
      + (-kap/(Ab*Nb**2) + 2*sp.diff(Nb, x)/(Ab*Nb**3) + sp.diff(Ab, x)/(2*Ab**2*Nb**2))*fss
      + (-2*kap/(Ab*Nb**2) + 2/Ab + sp.diff(Nb, x)/(Ab*Nb) + 4*sp.diff(Nb, x)/(Ab*Nb**3) + sp.diff(Ab, x)/(Ab**2*Nb**2))*fsx
      + (-kap/(Ab*Nb**2) + 2/Ab + sp.diff(Nb, x)/(Ab*Nb) + 2*sp.diff(Nb, x)/(Ab*Nb**3) - sp.diff(Ab, x)/(2*Ab**2) + sp.diff(Ab, x)/(2*Ab**2*Nb**2))*fxx
      - sp.diff(fss, x)/(Ab*Nb**2))
Dt=sp.expand(Dt)
for f in (bn,br,bt):
    print('d/dx', f, sp.simplify(Dt.coeff(sp.diff(f,x))))
rest=sp.expand(Dt-sum(Dt.coeff(sp.diff(f,x))*sp.diff(f,x) for f in (bn,br,bt)))
for f in (bn,br,bt):
    print('alg', f, sp.factor(sp.simplify(rest.coeff(f))))
