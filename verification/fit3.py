import math, cmath, numpy as np
from math import gcd, pi
from al import m, Th, Pm
from radial import S, mu
def fit(f,z,ts=(0.016,0.012,0.009,0.0065,0.005)):
    rows=[];ys=[]
    for t in ts:
        v=abs(f(z*math.exp(-t)))
        rows.append([1/t,math.log(t),1]); ys.append(math.log(v))
    sol,*_=np.linalg.lstsq(np.array(rows),np.array(ys),rcond=None)
    return sol[0]*64/pi**2, sol[1]
fs={'S0':lambda q:S(q,0)[0],'S0_AL':lambda q:2*m(-q**3,q**8,-1,R=300),'S0_TH':lambda q:q*Th(1,8,q,True)*Th(2,8,q)**2/Th(3,8,q)**2,
    'S1':lambda q:S(q,2)[0],'mu2':lambda q:mu(q)[0],'mu_AL':lambda q:4*m(-q,q**4,-1,R=300)}
for name,f in fs.items():
    out=[]
    for k in (4,8,16,32,5,12):
        for h in range(1,k):
            if gcd(h,k)!=1 or h>k//2: continue
            g,b=fit(f,cmath.exp(2j*pi*h/k))
            if g>0.02: out.append(f'{h}/{k}:{g:.3f}(t^{b:+.2f})')
    print(name,' '.join(out) or 'no exponential growth found')
