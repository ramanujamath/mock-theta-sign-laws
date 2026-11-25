import cmath, math
from math import pi
def logJ(m,q):
    s=0;n=1
    while True:
        t=q**(m*n); s+=cmath.log(1-t)
        if abs(t)<1e-18: return s
        n+=1
def logTH(q): return 8*logJ(2,q)-3*logJ(1,q)-4*logJ(4,q)
for h,k,a in ((1,4,pi**2/32),(3,4,pi**2/32),(1,8,pi**2/128),(3,8,pi**2/128),(5,8,pi**2/128),(7,8,pi**2/128)):
    for z in (0.003,0.0015):
        lv=logTH(cmath.exp(2j*pi*h/k)*math.exp(-z))-(-0.5*math.log(z)-z/8+a/z)
        v=cmath.exp(lv)
        print(h,k,z,'|C|=%.10f arg/pi*16=%.6f'%(abs(v),16*cmath.phase(v)/pi))
print(math.sqrt(pi/2), math.sqrt(pi)/2)
