import cmath, math
pi=math.pi
def U0(q):
    tot=0;cur=1;n=0
    while True:
        term=cur*q**(n*n); tot+=term
        if n>5 and abs(term)<1e-18*max(1,abs(tot)): return tot
        n+=1; cur*=(1+q**(2*n-1))/(1+q**(4*n))
def U1(q):
    tot=0;cur=1/(1+q**2);n=0
    while True:
        term=cur*q**((n+1)**2); tot+=term
        if n>5 and abs(term)<1e-18*max(1,abs(tot)): return tot
        n+=1; cur*=(1+q**(2*n-1))/(1+q**(4*n+2))
for nm,f,zs,gam in (('U0',U0,(1/8,3/8),pi**2/128),('U1',U1,(1/4,),pi**2/32)):
    for h in zs:
        z=cmath.exp(2j*pi*h)
        prev=None
        for t in (0.02,0.01,0.005,0.0025,0.00125):
            g=f(z*math.exp(-t))*math.exp(-gam/t)*math.sqrt(t)
            s=f'  slope={math.log(abs(g)/abs(prev))/math.log(t/prev[1]) if False else 0:.0f}' if False else ''
            print(f'{nm} at {h}: t={t}: |g|={abs(g):.6f} arg/pi={cmath.phase(g)/pi:.6f}')
            prev=(g,t)
