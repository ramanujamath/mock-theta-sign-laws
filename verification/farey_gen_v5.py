import cmath, math
from math import pi, sqrt, exp, log, gcd
import farey_gen as fg
def logJ(m,q):
    s=0; n=1
    while True:
        t=q**(m*n); s+=cmath.log(1-t)
        if abs(t)<1e-18: return s
        n+=1
for (h,k,a,C) in ((1,4,pi*pi/32,sqrt(pi/2)*cmath.exp(3j*pi/8)),(1,8,pi*pi/128,sqrt(pi)/2*cmath.exp(5j*pi/16)),(3,8,pi*pi/128,sqrt(pi)/2*cmath.exp(-13j*pi/16))):
    for z in (0.02+0j,0.02-0.05j,0.01+0.03j,0.03-0.1j):
        q=cmath.exp(2j*pi*h/k)*cmath.exp(-z)
        lth=8*logJ(2,q)-3*logJ(1,q)-4*logJ(4,q)
        main=C*z**-0.5*cmath.exp(-z/8+a/z)
        diff=abs(cmath.exp(lth)-main)
        X=(1/z).real; Lam=0
        for _,m,p in fg.TH:
            g=gcd(m,k); rho=exp(-4*pi*pi*g*g*X/(m*k*k)); j=1
            while rho**j>1e-18: Lam+=abs(p)*(-math.log1p(-rho**j)); j+=1
        bnd=abs(C)*z.real**-0.5*exp(a*X)*math.expm1(Lam)
        print(f'{h}/{k} z={z}: |TH-main|={diff:.3e} <= {bnd:.3e}: {diff<=bnd}   |TH/main-1|={abs(cmath.exp(lth)/main-1):.1e}')
