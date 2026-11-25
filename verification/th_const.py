import cmath, math
from math import pi, sqrt
def logJ(m,q):
    s=0; n=1
    while True:
        t=q**(m*n); s+=cmath.log(1-t)
        if abs(t)<1e-18: return s
        n+=1
def logTH(q): return 8*logJ(2,q)-3*logJ(1,q)-4*logJ(4,q)
for h,k,a in ((1,4,pi**2/32),(1,8,pi**2/128),(3,8,pi**2/128)):
    for z in (0.004,0.002,0.001):
        lv=logTH(cmath.exp(2j*pi*h/k)*math.exp(-z)) - (-0.5*math.log(z) - z/8 + a/z)
        v=cmath.exp(lv)
        print(f'TH at {h}/{k}, z={z}: |C|={abs(v):.8f}  arg/pi={cmath.phase(v)/pi:.8f}')
print('sqrt(pi/2)=%.8f sqrt(pi)/2=%.8f  5/16=%.5f -13/16=%.5f'%(sqrt(pi/2),sqrt(pi)/2,5/16,-13/16))
