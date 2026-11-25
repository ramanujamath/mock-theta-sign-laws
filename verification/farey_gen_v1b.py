import math, cmath, random
from math import pi, sqrt, exp, log, gcd
from fractions import Fraction as Fr
import farey_gen as fg
from farey_lib import phiU, phiL, etaU, etaL
random.seed(7)
def log_theta_prod(x,Q):   # log|(x;Q)(Q/x;Q)(Q;Q)| , no cancellation
    s=math.log(abs(1-x)); n=1
    while True:
        a=Q**n; s+=math.log(abs(1-a*x))+math.log(abs(1-a/x))+math.log(abs(1-a))
        if abs(a)<1e-18*max(1,abs(x),1/abs(x)): return s
        n+=1
def rand_point(t0):
    K=int(sqrt(2*pi/t0)); k=random.randint(1,K)
    h=random.choice([x for x in range(k) if gcd(x,k)==1])
    phi=random.uniform(-2*pi/(k*(K+1)),2*pi/(k*(K+1)))
    return h,k,complex(t0,-phi)
facs=[('th',4,Fr(0),Fr(1,2)),('th',4,Fr(3,4),Fr(0)),('th',4,Fr(1,2),Fr(1,2)),('th',4,Fr(1,4),Fr(0))]
bad=0;cnt=0;worst=9
for t0 in (0.02,0.008,0.004):
    for _ in range(3000):
        h,k,z=rand_point(t0); q=cmath.exp(2j*pi*h/k)*cmath.exp(-z); X=(1/z).real; w=X/k**2; kz=k*abs(z)
        for _,M,A,s in facs:
            x=cmath.exp(2j*pi*s)*q**int(A*M); act=log_theta_prod(x,q**M)
            ff=('th',M,A,s,1); R,pr,p,qp,prod=fg.factor_data(ff,h,k)
            g=gcd(M,k); kp=k//g; hp=M*h//g; Ah=fg.ahat(A,s,hp,kp); rho=exp(-4*pi*pi*g*g*w/M)
            base=pr-0.5*log(kz)-float(qp)*t0+(pi*pi/64)*float(R)/k**2*X
            lo=base+phiL(Ah,rho)+etaL(rho); hi=base+phiU(Ah,rho)+etaU(rho)
            cnt+=1; sl=min(act-lo,hi-act); worst=min(worst,sl)
            if sl<-1e-9: bad+=1
print(f'theta factors via triple product: {cnt} checks, {bad} violations, worst slack {worst:.2e}')
