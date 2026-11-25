import pickle, math, cmath, numpy as np, sys
from math import pi, sqrt
from decimal import Decimal as Dm, getcontext
getcontext().prec=60
import farey_N0, farey_S1_N0
D=pickle.load(open('big.pkl','rb'))
nm=sys.argv[1]; n=int(sys.argv[2])
c=D[nm]
t0=pi/(16*sqrt(n)); K=int(sqrt(2*pi/t0))
mod=farey_N0 if nm=='S0' else farey_S1_N0
tot,lmi,lma,lt,*_=mod.bound_parts(n); lsup=lma-n*t0
C = sqrt(pi)/2*cmath.exp(9j*pi/16) if nm=='S0' else sqrt(pi)/2*cmath.exp(-3j*pi/16)
s = 1/16 if nm=='S0' else -7/16
# arc of 1/4 in Farey order K: neighbours
from math import gcd
fr=sorted(set((h//gcd(h,k),k//gcd(h,k)) for k in range(1,K+1) for h in range(k+1)),key=lambda x:x[0]/x[1])
i=fr.index((1,4)); lo=(fr[i-1][0]+1)/(fr[i-1][1]+4); hi=(1+fr[i+1][0])/(4+fr[i+1][1])
print('K',K,'arc',lo,hi,'phi range',2*pi*(lo-.25),2*pi*(hi-.25),' log major sup',lsup)
# high-precision direct evaluation
def fval(theta):
    # sum c_j r^j e^{i j theta} in Decimal using recurrence on complex numbers
    r=Dm(-t0).exp()
    co=Dm(repr(math.cos(theta))); si=Dm(repr(math.sin(theta)))   # theta itself given as float; consistent with z below
    # better: compute cos/sin of theta to high precision via float is insufficient; use exact float theta value -> Decimal cos
    x=Dm(theta); # Taylor
    def dcos(x):
        s=Dm(0);t=Dm(1);k=0
        while abs(t)>Dm(10)**-58: s+=t;k+=2;t=-t*x*x/(k*(k-1))
        return s
    def dsin(x):
        s=Dm(0);t=x;k=1
        while abs(t)>Dm(10)**-58: s+=t;k+=2;t=-t*x*x/(k*(k-1))
        return s
    co,si=dcos(x),dsin(x)
    wr,wi=r*co,r*si; pr,pi_=Dm(1),Dm(0); sr,si_=Dm(0),Dm(0)
    for j,cj in enumerate(c):
        if cj: sr+=cj*pr; si_+=cj*pi_
        pr,pi_=pr*wr-pi_*wi, pr*wi+pi_*wr
    return complex(float(sr),float(si_)), sr, si_
for frac in np.linspace(lo,hi,13):
    theta=2*pi*frac; phi=theta-pi/2
    z=t0-1j*phi
    main=C*z**-0.5*cmath.exp(-s*z+(pi**2/64)/z)
    fv,sr,si_=fval(theta)
    # high-precision difference: main in float has rel err 1e-16
    diff=abs(fv-main)
    print(f'phi={phi:+.5f} |f|={abs(fv):.3e} |main|={abs(main):.3e} |f-main|={diff:.3e} log ratio to sup={math.log(diff)-lsup:+.2f}  rho={math.exp(-pi**2*z.real/abs(z)**2/2):.3f}')
