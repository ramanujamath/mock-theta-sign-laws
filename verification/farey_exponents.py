# Exact cusp exponents (units pi^2 Re(1/z)/64) of the crude bounds, for ALL cusps, via Bernoulli polynomials.
# E = R/k'^2 where R depends only on (h' mod 8, k' mod 2).
from farey_lib import *
from math import gcd
def R(fl,hp,kp): return class_exponent(fl,hp,1 if kp%2 else 2,None)  # value at k'=1 or 2 representative, times... (see below)
rows=[]
for kpar in (1,2):
    for hm in range(8):
        if kpar==2 and hm%2==0: continue
        # representative with given parity: use kp=kpar (only parity enters Ahat), exponent scales 1/k'^2
        Rt=class_exponent(TQ,hm,kpar,None)*kpar**2
        Ra=class_exponent(ALA,hm,kpar,None)*kpar**2
        Rb=max(class_exponent(ALB1,hm,kpar,None),class_exponent(ALB2,hm,kpar,None))*kpar**2
        rows.append((kpar,hm,Rt,Ra,Rb))
        print(f"k' {'odd ' if kpar==1 else 'even'}  h'={hm} mod 8:  R_TQ={str(Rt):>6}  R_ALA={str(Ra):>5}  R_formB={str(Rb):>5}")
# which (h,k) give which class: k odd -> h'=0; k=2 mod 4 -> h'=4; k=4 mod 8 -> h' in {2,6}; 8|k -> h' odd
print()
best={}
for k in range(1,2001):
    for h in range(k):
        if gcd(h,k)!=1: continue
        hp,kp=cusp_data(h,k)
        Et=class_exponent(TQ,hp,kp,None); Ea=class_exponent(ALA,hp,kp,None)
        Eb=max(class_exponent(ALB1,hp,kp,None),class_exponent(ALB2,hp,kp,None))
        EA=max(Et,Ea,0); EB=max(Eb,0)
        best[(h,k)]=(min(EA,EB),EA,EB)
mx=max((v[0],hk) for hk,v in best.items() if hk[1]!=4)
print('max over k<=2000, k!=4 of min(formA,formB) exponent:',mx[0],'at',mx[1])
print('k=4:',best[(1,4)],best[(3,4)])
