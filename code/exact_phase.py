"""Exact leading coefficient of Theta(s q^a; q^M) at q = e(h/k) e^{-t}, from the exact Poisson identity
   Theta = (1/L) sqrt(2pi/(tM)) e^{t b^2/2M} sum_nu G_nu e(nu b/(L M)) exp(-2 pi^2 nu^2/(L^2 M t)),  L=2k, b=a-M/2.
The leading coefficient c* = (1/L) sum_{|nu|=nu*} G_nu e(nu b/(LM)) is a finite sum of roots of unity."""
import cmath, math
from math import pi
from fractions import Fraction as Fr
def lead(s,a,M,h,k):
    L=2*k; b=Fr(a)-Fr(M,2)
    def G(nu):
        return sum((-s)**j*cmath.exp(2j*pi*((Fr(h*(M*j*(j-1)//2+a*j),k)+Fr(nu*j,L))%1)) for j in range(L))
    for nu in range(0,3*L):
        c=sum(G(v)*cmath.exp(2j*pi*float((Fr(v)*b/(L*M))%1)) for v in ({nu,-nu}))/L
        if abs(c)>1e-9: return nu,c
def phase_of(factors,h,k,qpref=(0,0)):
    """factors: list (p, s, a, M); qpref: (c, ) multiplier q^c -> phase of e(h c/k). returns total arg/pi, exponent in units pi^2/(64 t)"""
    ph=float(qpref[0])*2*h/k; ex=0.0; amp=1.0
    for p,s,a,M in factors:
        nu,c=lead(s,a,M,h,k)
        ph+=p*cmath.phase(c)/pi
        ex+=-p*(2*pi*pi*nu*nu/((2*k)**2*M))/(pi*pi/64)
    return ph, ex
def report(name,factors,h,k,qc=0):
    ph,ex=phase_of(factors,h,k,(qc,))
    r=Fr(ph/2).limit_denominator(64)*2   # phase/pi mod 2
    r=r-2*math.floor((r+1)/2)
    print(f'{name} at {h}/{k}: exponent {ex:.6f} (units pi^2/64t); phase/pi = {ph%2 if ph%2<=1 else ph%2-2:.12f} ~ {r}  (residual {abs(((ph-float(r))+1)%2-1):.1e})')
# S0 TQ = q Th(-q;q^8) Th(q^2;q^8)^2 / Th(q^3;q^8)^2  [qpref phase handled separately: q = e(h/k) e^{-t} -> e(h/k)]
report('S0-TQ',[(1,-1,1,8),(2,1,2,8),(-2,1,3,8)],1,4,qc=1)
report('S1-TQ',[(1,-1,3,8),(2,1,2,8),(-2,1,1,8)],1,4,qc=-1)
J=lambda p,m:(p,1,m,3*m)
TH=[J(8,2),J(-3,1),J(-4,4)]
for h,k in ((1,4),(1,8),(3,8)): report('TH',TH,h,k)
