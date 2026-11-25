"""Independent cusp exponents via Poisson/Gauss sums (no Siegel machinery).
Theta(s q^a; q^M) = sum_n (-s)^n q^{M n(n-1)/2 + a n}.  Leading growth at q=e(h/k)e^{-t}:
exp(-(pi^2/(64 t)) * 32 nu*^2/(k^2 M)) where nu* = least |nu| with nonzero combined Gauss coefficient.
Exponent returned in units pi^2/(64 t)."""
import cmath, math
from math import pi, gcd
from fractions import Fraction as Fr
from functools import lru_cache
@lru_cache(maxsize=None)
def theta_exp(s, a, M, h, k):
    L = 2*k; n0 = Fr(1,2) - Fr(a, M)
    c = []
    for j in range(L):
        E = Fr(M*j*(j-1), 2) + a*j   # integer since a integer here (or Fraction)
        c.append(((-s)**j) * cmath.exp(2j*pi*float((Fr(h)*E/k) % 1)))
    for nu in range(0, 4*L):
        tot = 0
        for v in ({nu, -nu}):
            G = sum(c[j]*cmath.exp(2j*pi*((v*j) % L)/L) for j in range(L))
            tot += G*cmath.exp(-2j*pi*float((n0*v/L) % 1))
        if abs(tot) > 1e-7*L:
            return Fr(-32*nu*nu, k*k*M)
    raise Exception('none')
def J(m, h, k):   # (q^m;q^m)_inf = Theta(q^m; q^{3m})
    return theta_exp(1, m, 3*m, h, k)
def term_exp(factors, h, k):  # factors: (p, s, a, M)
    return sum(p*theta_exp(s, a, M, h, k) for p, s, a, M in factors)
# --- S0 / S1 terms ---
TQ  = [(1,-1,1,8),(2,1,2,8),(-2,1,3,8)]
TQ1 = [(1,-1,3,8),(2,1,2,8),(-2,1,1,8)]
ALA = [(-1,-1,0,8)]; ALB1=[(-1,-1,2,8)]; ALB2=[(-1,-1,6,8)]
TH = [('J',8,2),('J',-3,1),('J',-4,4)]
def th_exp(h,k): return sum(p*J(m,h,k) for _,p,m in TH)
U0A=[(-1,-1,0,4)]; U0B=[(-1,1,3,4)]; U1A=[(-1,-1,2,4)]; U1B=[(-1,1,1,4)]
if __name__=='__main__':
    KMAX=200
    worst={}
    for k in range(1,KMAX+1):
        g=gcd(8,k); kp=k//g
        for h in range(k):
            if gcd(h,k)!=1: continue
            # S0, S1 : form A if k' odd else form B
            if kp%2:
                e0=max(term_exp(TQ,h,k), term_exp(ALA,h,k)); e1=max(term_exp(TQ1,h,k), term_exp(ALA,h,k))
            else:
                e0=e1=max(term_exp(ALB1,h,k), term_exp(ALB2,h,k))
            th=th_exp(h,k)
            u0a,u0b,u1a,u1b=(term_exp(f,h,k) for f in (U0A,U0B,U1A,U1B))
            # same class choices as farey_gen (rep index by class)
            vals={'S0':e0,'S1':e1}
            # record all available terms for gen functions
            vals['U0opts']=(u0a,u0b,max(th,u1a),max(th,u1b))
            vals['U1opts']=(max(th,u0a),max(th,u0b),u1a,u1b)
            vals['mu2opts']=(max(th,u0a),max(th,u0b))
            for nm,v in vals.items():
                if nm in('S0','S1'):
                    if k!=4: worst.setdefault(nm,[]).append((v,h,k))
            for nm,main in (('U0',8),('U1',4),('mu2',4)):
                if k!=main:
                    worst.setdefault(nm+'_bestrep',[]).append((min(vals[nm+'opts']),h,k))
            if k in (4,8):
                print('cusp',h,k,'S0',e0,'S1',e1,'TH',th,'U0 forms',u0a,u0b,'U1 forms',u1a,u1b)
    for nm,l in worst.items():
        l.sort(reverse=True); print(nm,'top:',[(str(v),h,k) for v,h,k in l[:6]])
