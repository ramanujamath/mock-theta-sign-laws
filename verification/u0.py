import cmath, math, numpy as np
from math import pi, gcd
from al import m
RR=110
def U0(q): return 2*m(-q,q**4,-1,R=RR)
def U1(q): return -m(-q,q**4,-q**2,R=RR)
# verify forms vs series coefficients
import pickle; F=pickle.load(open('family.pkl','rb'))
q0=0.3+0.2j
RR=40
for nm,f in (('U0',U0),('U1',U1)):
    s=sum(c*q0**n for n,c in enumerate(F[nm][:400]))
    print(nm,'AL form check',f(q0),s)
RR=110
def fit(f,z,ts=(0.012,0.008,0.006,0.0045,0.0035)):
    A=[];y=[]
    for t in ts:
        v=f(z*math.exp(-t)); A.append([1/t,math.log(t),1]); y.append(math.log(abs(v)))
    s,*_=np.linalg.lstsq(np.array(A),np.array(y),rcond=None); return s
for nm,f in (('U0',U0),('U1',U1)):
    out=[]
    for k in (1,2,3,4,8,16):
        for h in range(k):
            if gcd(h,k)!=1: continue
            g,b,c=fit(f,cmath.exp(2j*pi*h/k))
            out.append((round(g*64/pi**2,3),round(b,2),f'{h}/{k}'))
    print(nm,sorted(out,reverse=True)[:6])
