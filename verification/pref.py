import cmath, math
from radial import S, mu
pi=math.pi
for name,f,a in [('S0',lambda q:S(q,0),pi**2/64),('S1',lambda q:S(q,2),pi**2/64),('mu2',mu,pi**2/32)]:
    print('==',name)
    prev=None
    for t in (0.04,0.02,0.01,0.005,0.0025):
        v=f(1j*math.exp(-t))[0]*math.exp(-a/t)
        s=''
        if prev: s=f' beta_est={-math.log(abs(v)/abs(prev[1]))/math.log(t/prev[0]):.4f}'
        print(f' t={t}: g={v:.6f} |g|={abs(v):.6f} arg/pi={cmath.phase(v)/pi:.6f}'+s)
        prev=(t,v)
