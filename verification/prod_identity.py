"""Proof of U0 + 2 U1 = TH as an identity of infinite products.
By HM (4.39),(4.40): U0 + 2U1 = 2[m(-q,q^4,-1) - m(-q,q^4,-q^2)].
By HM Thm 2.3 (2.7): m(x,Q,z1)-m(x,Q,z0) = z0 J^3 j(z1/z0;Q) j(x z0 z1;Q) / ( j(z0;Q) j(z1;Q) j(x z0;Q) j(x z1;Q) ),  J=(Q;Q)_inf.
We first confirm this form of (2.7) numerically, then compare exact product exponents."""
import cmath
from fractions import Fraction as Fr
from al import m
def j(x,Q,R=80): return sum((-1)**n*Q**(n*(n-1)//2)*x**n for n in range(-R,R+1))
def Jinf(Q):
    p=1
    for k in range(1,2000): p*=1-Q**k
    return p
q=0.3+0.2j; Q=q**4; x=-q; z1=-1; z0=-q**2
lhs=m(x,Q,z1,R=60)-m(x,Q,z0,R=60)
rhs=z0*Jinf(Q)**3*j(z1/z0,Q)*j(x*z0*z1,Q)/(j(z0,Q)*j(z1,Q)*j(x*z0,Q)*j(x*z1,Q))
print('(2.7) form check:', abs(lhs-rhs))
# ---- exact product algebra: element = (const, qpow, {n: e_n}) meaning const*q^qpow*prod (1-q^n)^{e_n}
L=64   # compare exponents for n<=L (all sequences have period dividing 16)
def one(): return [Fr(1),0,{}]
def mul(a,b,p=1):
    c=[a[0]*b[0]**p, a[1]+p*b[1], dict(a[2])]
    for n,e in b[2].items(): c[2][n]=c[2].get(n,0)+p*e
    return c
def lin(s,e):
    """factor (1 - s q^e), s=+-1, e integer (may be <=0)"""
    if e==0: return [Fr(1-s),0,{}]
    if e<0:  # 1 - s q^e = -s q^e (1 - s q^{-e})   (s^{-1}=s)
        f=lin(s,-e); return [-s*f[0], e+f[1], f[2]]
    if s==1: return [Fr(1),0,{e:1}]
    return [Fr(1),0,{2*e:1, e:-1}]           # 1+q^e = (1-q^{2e})/(1-q^e)
def poch(s,a,M):
    """(s q^a; q^M)_inf truncated to factors with exponent <= L (a may be <=0 only for first factor)"""
    r=one(); k=0
    while a+M*k<=L:
        r=mul(r,lin(s,a+M*k)); k+=1
    return r
def jj(s,a,M):   # j(s q^a; q^M) = (x;Q)(Q/x;Q)(Q;Q)
    return mul(mul(poch(s,a,M),poch(s,M-a,M)),poch(1,M,M))
def Jm(m_): return poch(1,m_,m_)
def trunc(e): return {n:v for n,v in e[2].items() if v!=0 and n<=L}
# 2D with x=-q, z1=-1, z0=-q^2 (Q=q^4): z0=-q^2 -> const -1, qpow 2
D=[Fr(-1)*2,2,{}]
for (s,a,p) in [(1,0,0)]: pass
D=mul(D,Jm(4),3)
D=mul(D,jj(1,-2,4))          # j(z1/z0) = j(q^{-2};q^4)
D=mul(D,jj(-1,3,4))          # j(x z0 z1) = j(-q^3)
D=mul(D,jj(-1,2,4),-1)       # j(z0)=j(-q^2)
D=mul(D,jj(-1,0,4),-1)       # j(z1)=j(-1)
D=mul(D,jj(1,3,4),-1)        # j(x z0)=j(q^3)
D=mul(D,jj(1,1,4),-1)        # j(x z1)=j(q)
TH=mul(mul(Jm(2),Jm(1),0),one())
TH=[Fr(1),0,{}]; TH=mul(TH,Jm(2),8); TH=mul(TH,Jm(1),-3); TH=mul(TH,Jm(4),-4)
print('2D : const',D[0],' qpow',D[1])
print('TH : const',TH[0],' qpow',TH[1])
same=trunc(D)==trunc(TH)
print('exponent sequences of (1-q^n), n <=',L,'agree:',same)
if not same:
    for n in range(1,L+1): print(n,trunc(D).get(n,0),trunc(TH).get(n,0))
