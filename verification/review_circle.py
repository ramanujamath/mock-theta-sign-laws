"""Evaluate each function on |q|=e^{-t0} via exact coefficients + FFT and test the pointwise sup bounds."""
import pickle, math, cmath, numpy as np
from math import pi, sqrt, gcd, floor
import farey_N0, farey_S1_N0, farey_gen as G
D=pickle.load(open('big.pkl','rb')); F=pickle.load(open('family.pkl','rb'))
cf={'S0':D['S0'],'S1':D['S1'],'mu2':D['mu2'],'U0':F['U0'],'U1':F['U1']}
MAIN={'S0':dict(k=4,a=pi**2/64,s=1/16,C={1:sqrt(pi)/2*cmath.exp(9j*pi/16)}),
      'S1':dict(k=4,a=pi**2/64,s=-7/16,C={1:sqrt(pi)/2*cmath.exp(-3j*pi/16)}),
      'mu2':dict(k=4,a=pi**2/32,s=1/8,C={1:-sqrt(pi/2)*cmath.exp(3j*pi/8)}),
      'U1':dict(k=4,a=pi**2/32,s=1/8,C={1:0.5*sqrt(pi/2)*cmath.exp(3j*pi/8)}),
      'U0':dict(k=8,a=pi**2/128,s=1/8,C={1:sqrt(pi)/2*cmath.exp(5j*pi/16),3:sqrt(pi)/2*cmath.exp(-13j*pi/16)})}
for nm,P in MAIN.items():
    k=P['k']
    for h in list(P['C']): P['C'][k-h]=P['C'][h].conjugate()
EC={nm:G.minor_setup(nm) for nm in ('mu2','U1','U0')}
def params(nm,n):
    if nm=='S0':
        t0=pi/(16*sqrt(n)); tot,lmi,lma,lt,*_=farey_N0.bound_parts(n)
    elif nm=='S1':
        t0=pi/(16*sqrt(n)); tot,lmi,lma,lt,*_=farey_S1_N0.bound_parts(n)
    else:
        E,ch=EC[nm]; t0=(pi/8)*sqrt(float(E)/n); tot,lmi,lma,lt,_=G.bound_parts(nm,n,E,ch)
    return t0, lmi-n*t0, lma-n*t0
def farey(K):
    fs=sorted(set((h//gcd(h,k),k//gcd(h,k)) for k in range(1,K+1) for h in range(0,k+1)),key=lambda x:x[0]/x[1])
    return fs
for nm,ntest in (('S0',1500),('S0',3000),('S1',1500),('S1',3000),('mu2',300),('mu2',800),('U1',150),('U1',400),('U0',400),('U0',1200)):
    c=cf[nm]; t0,lsup_minor,lsup_major=params(nm,ntest)
    K=int(floor(sqrt(2*pi/t0)))
    nn=np.arange(len(c)); lw=-nn*t0
    a=np.array([float(x) for x in c])*np.exp(lw)
    # truncation estimate: tail magnitude bound by last terms
    L=1<<19; A=np.zeros(L); A[:len(c)]=a
    f=np.fft.ifft(A)*L          # f(theta_j) = sum a_n e^{i n theta_j}
    th=2*pi*np.arange(L)/L
    fs=farey(K); P=MAIN[nm]; k0=P['k']
    worst_minor=-1e9; worst_major=-1e9
    # mediants
    for i,(h,k) in enumerate(fs):
        if h==k and k==1: continue
        lo=(fs[i-1][0]+h)/(fs[i-1][1]+k) if i>0 else (fs[-2][0]+fs[-1][0]-1+1)/(1)  # handled below
        if i==0: lo=(fs[-2][0]-fs[-2][1]+h)/(fs[-2][1]+k)  # wrap: neighbour of 0/1 below is (K-1)/K -1
        hi=(h+fs[i+1][0])/(k+fs[i+1][1])
        j0=int(math.ceil(lo*L)); j1=int(math.floor(hi*L))
        idx=np.arange(j0,j1+1)%L
        vals=f[idx]
        if k==k0 and h in P['C']:
            phi=th[idx]-2*pi*h/k; phi=(phi+pi)%(2*pi)-pi
            z=t0-1j*phi
            main=P['C'][h]*z**-0.5*np.exp(-P['s']*z+P['a']/z)
            r=np.log(np.abs(vals-main)).max()-lsup_major; worst_major=max(worst_major,r)
        else:
            r=np.log(np.abs(vals)).max()-lsup_minor; worst_minor=max(worst_minor,r)
    tailest=max(math.log(abs(float(c[j]))+1e-300)-j*t0 for j in range(len(c)-50,len(c)))
    print(f'{nm} n={ntest} t0={t0:.5f} K={K}: max log(|f|/minor sup)={worst_minor:.2f}  max log(|f-main|/major sup)={worst_major:.2f}  (minor sup log {lsup_minor:.1f}, major {lsup_major:.1f}, trunc term log {tailest:.1f})')
