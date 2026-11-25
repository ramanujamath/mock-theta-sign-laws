import numpy as np, math
from math import gcd, pi
def nu_eff(s,a,M,h,k):
    L=2*k; j=np.arange(L)
    E=M*j*(j-1)//2 + a*j
    ph=((-s)**(j%2))*np.exp(2j*pi*((h*E)%k)/k)
    G=np.fft.ifft(ph)*L              # G_nu = sum_j ph_j e(+nu j/L), nu=0..L-1 (periodic in nu)
    b=a-M/2
    c=lambda nu: G[nu%L]*np.exp(2j*pi*nu*b/(L*M))
    tol=1e-7*L
    if abs(c(0))>tol: return 0
    for nu in range(1,3*L):
        if abs(c(nu)+c(-nu))>tol: return nu
    raise RuntimeError
def expo(s,a,M,h,k):  # coefficient of 1/t in log|Theta|, units pi^2/64
    nu=nu_eff(s,a,M,h,k); return -(pi**2*nu**2/(2*k*k*M))/(pi**2/64)
J=lambda m:(1,m,3*m)
TH={'S0':[(1,-1,1,8),(2,1,2,8),(-2,1,3,8)],'S1':[(1,-1,3,8),(2,1,2,8),(-2,1,1,8)],'mu2':[(8,)+J(2),(-3,)+J(1),(-4,)+J(4)]}
ALA={'S0':[(-1,0,8)],'S1':[(-1,0,8)],'mu2':[(-1,0,4)]}
ALB={'S0':[(-1,2,8),(-1,6,8)],'S1':[(-1,2,8),(-1,6,8)],'mu2':[(1,1,4),(1,3,4)]}
allthetas=set()
for d in (TH,):
    for v in d.values():
        for p,s,a,M in v: allthetas.add((s,a,M))
for d in (ALA,ALB):
    for v in d.values():
        for t in v: allthetas.add(t)
K=400
# 1) check nu_eff <= 6M for every theta and every cusp k<=K; record max nu_eff/M
mx={}
for (s,a,M) in allthetas:
    worst=0
    for k in range(1,K+1):
        for h in range(k):
            if gcd(h,k)==1:
                worst=max(worst,nu_eff(s,a,M,h,k))
    mx[(s,a,M)]=worst
print('max nu_eff over k<=%d (bound 6M):'%K, {f'Theta({s:+d}q^{a};q^{M})':(v,6*M) for (s,a,M),v in mx.items()})
# 2) best-form exponent at every cusp k<=K, k!=4
for name in TH:
    worst=(0,None)
    for k in range(1,K+1):
        if k==4: continue
        for h in range(k):
            if gcd(h,k)!=1: continue
            eth=sum(p*expo(s,a,M,h,k) for p,s,a,M in TH[name])
            A=max(eth,max(-expo(*z,h,k) for z in ALA[name]))
            B=max(-expo(*z,h,k) for z in ALB[name])
            e=min(A,B)
            if e>worst[0]+1e-12 or (worst[1] is None): worst=(max(e,worst[0]),f'{h}/{k}' if e>=worst[0] else worst[1])
    print(name,'max exponent over cusps k<=%d, k!=4: %.4f at %s'%(K,worst[0],worst[1]))
