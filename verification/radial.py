import cmath, math
from math import gcd, pi
def S(q,shift):
    tot=0; cur=1; n=0; big=0
    while True:
        term=cur*q**(n*n+shift*n); tot+=term; big=max(big,abs(term))
        if n>5 and abs(term)<1e-18*max(1,abs(tot)): break
        n+=1; cur*=(1+q**(2*n-1))/(1+q**(2*n))
    return tot,big
def mu(q):
    tot=0; cur=1; n=0; big=0
    while True:
        term=(-1)**n*cur*q**(n*n); tot+=term; big=max(big,abs(term))
        if n>5 and abs(term)<1e-18*max(1,abs(tot)): break
        n+=1; cur*=(1-q**(2*n-1))/(1+q**(2*n))**2
    return tot,big
