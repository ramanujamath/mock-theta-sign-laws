import math, cmath, random
from math import pi, sqrt, exp, log
from fractions import Fraction as Fr
from farey_lib import *
random.seed(1)
rows=[]
for t0 in (0.02,0.01):
    K=int(sqrt(2*pi/t0))
    for _ in range(1500):
        k=random.randint(1,K); h=random.choice([x for x in range(k) if math.gcd(x,k)==1])
        phi=random.uniform(-2*pi/(k*(K+1)),2*pi/(k*(K+1)))
        z=complex(t0,-phi); q=cmath.exp(2j*pi*h/k)*cmath.exp(-z)
        hp,kp=cusp_data(h,k); rho=exp(-pi*pi*(1/z).real/(2*kp*kp))
        for a,sig,s in ((3,Fr(1,2),-1),(1,Fr(0),1)):
            A=Fr(a,8); act=abs(theta_direct(s*q**a,q**8))
            if act<1e-250: continue
            Ah=Ahat(A,sig,hp,kp)
            base=float(8*(-Fr(1,8)+A/2-A*A/2))*(-t0) - 0.5*log(4*kp*abs(z)/pi) + (float(B2(Ah))/2+1/24)*log(rho)
            lo=base+phiL(Ah,rho)+etaL(rho); hi=base+phiU(Ah,rho)+etaU(rho)
            la=log(act)
            if la<lo-1e-6 or la>hi+1e-6:
                tr=log(theta_via_transform(a,sig*1 if s==-1 else 0,h,k,z)) if True else None
                # max term size of the direct sum (roundoff floor)
                rows.append((h,k,a,round(la,2),round(lo,2),round(hi,2),round(tr,2)))
for r in rows[:12]: print('h/k=%d/%d a=%d  log|direct|=%s  bounds [%s,%s]  log|via transform|=%s'%r)
print(len(rows),'cases; transform value inside bounds in all:', all(r[4]-1e-6<=r[6]<=r[5]+1e-6 for r in rows))
