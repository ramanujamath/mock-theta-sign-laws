import math, cmath
from math import gcd, pi
from al import m, Th, Pm
pieces={
 'S0_AL':lambda q:2*m(-q**3,q**8,-1,R=250),
 'S0_TH':lambda q:q*Th(1,8,q,True)*Th(2,8,q)**2/Th(3,8,q)**2,
 'S1_AL':lambda q:-(2/q)*m(-q,q**8,-1,R=250),
 'S1_TH':lambda q:(1/q)*Th(3,8,q,True)*Th(2,8,q)**2/Th(1,8,q)**2,
 'mu_AL':lambda q:4*m(-q,q**4,-1,R=250),
 'mu_TH':lambda q:-Pm(2,q,40000)**8/(Pm(1,q,40000)**3*Pm(4,q,40000)**4),
}
t1,t2=0.012,0.008
for name,f in pieces.items():
    res=[]
    for k in range(1,33):
        for h in range(k):
            if gcd(h,k)!=1: continue
            z=cmath.exp(2j*pi*h/k)
            try:
                a=abs(f(z*math.exp(-t1))); b=abs(f(z*math.exp(-t2)))
                g=(math.log(b)-math.log(a))/(1/t2-1/t1)
            except (ZeroDivisionError,ValueError,OverflowError): g=float('nan')
            res.append((g*64/pi**2,h,k))
    res.sort(reverse=True)
    print(name,' top (gamma*64/pi^2, h, k):',[(round(g,3),h,k) for g,h,k in res[:6]])
