# Leading growth exponent (coefficient of 1/t, in units of pi^2/64) of theta factors at cusp q = e^{2 pi i h/k} e^{-t}
# Theta(s*q^a ; q^M) = sum_n (-s)^n q^{M n(n-1)/2 + a n}.  Poisson in n = L l + j, L = 2k:
# Theta ~ (1/L) sqrt(2pi/(tM)) sum_nu G_nu exp(-2 pi^2 nu^2/(L^2 M t)),  G_nu = sum_j (-s)^j e(hE(j)/k) e(nu j/L)
import cmath, math
from math import gcd, pi
def E(j,a,M): return M*j*(j-1)//2 + a*j
def theta_exp(s,a,M,h,k,numax=None):
    L=2*k
    for nu in range(0, L):             # |nu| increasing; G_nu, G_-nu
        for sg in ((1,) if nu==0 else (1,-1)):
            G=sum((-s)**j*cmath.exp(2j*pi*(h*E(j,a,M)/k + sg*nu*j/L)) for j in range(L))
            if abs(G)>1e-8*L: return -(2*pi**2*nu**2/(L**2*M))/(pi**2/64)
    raise ValueError
# factor lists: (power, s, a, M)
def J(m):  return (1, m, 3*m)          # (q^m;q^m)_inf = Theta(q^m; q^{3m}), s=+1
TH={'S0':[(1,-1,1,8),(2,1,2,8),(-2,1,3,8)],      # q*Theta(-q;q^8)Theta(q^2;q^8)^2/Theta(q^3;q^8)^2
    'S1':[(1,-1,3,8),(2,1,2,8),(-2,1,1,8)],
    'mu2':[(8,)+J(2),(-3,)+J(1),(-4,)+J(4)]}
ALden={'S0':(1,-1,0,8),'S1':(1,-1,0,8),'mu2':(1,-1,0,4)}   # 1/Theta(-1;q^8) or 1/Theta(-1;q^4)
K=60
for name in TH:
    rows=[]
    for k in range(1,K+1):
        for h in range(k):
            if gcd(h,k)!=1: continue
            eth=sum(p*theta_exp(s,a,M,h,k) for p,s,a,M in TH[name])
            _,s,a,M=ALden[name]; eal=-theta_exp(s,a,M,h,k)
            rows.append((eth,eal,h,k))
    top=sorted(rows,reverse=True)[:6]
    print(name,'theta-quotient exponents (top):',[(round(e,4),f'{h}/{k}') for e,_,h,k in top])
    print('   crude AL bound exponents (top):',[(round(e,4),f'{h}/{k}') for _,e,h,k in sorted(rows,key=lambda r:-r[1])[:6]])
    others=max(max(e,f) for e,f,h,k in rows if k!=4)
    print('   max over all cusps k!=4 (k<=%d): %.4f'%(K,others))
