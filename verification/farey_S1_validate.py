import math, cmath, random, pickle
from math import pi, sqrt, exp, log
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
from farey_lib import *
random.seed(1)
# V1': Lemma 1+2 two-sided bounds for the two NEW factors Theta(-q^3;q^8) [A=3/8,s=1/2] and Theta(q;q^8) [A=1/8,s=0]
viol=0; worst=1e9; cnt=0
for t0 in (0.02,0.01):
    K=int(sqrt(2*pi/t0))
    for _ in range(1500):
        k=random.randint(1,K); h=random.choice([x for x in range(k) if math.gcd(x,k)==1])
        phi=random.uniform(-2*pi/(k*(K+1)),2*pi/(k*(K+1)))
        z=complex(t0,-phi); q=cmath.exp(2j*pi*h/k)*cmath.exp(-z)
        hp,kp=cusp_data(h,k); rho=exp(-pi*pi*(1/z).real/(2*kp*kp))
        for a,sig,s in ((3,Fr(1,2),-1),(1,Fr(0),1)):
            A=Fr(a,8); act=abs(theta_direct(s*q**a,q**8))
            if act<1e-250: continue
            Ah=Ahat(A,sig,hp,kp)
            base=float(8*(-Fr(1,8)+A/2-A*A/2))*(-t0) - 0.5*log(4*kp*abs(z)/pi) + (float(B2(Ah))/2+1/24)*log(rho)
            lo=base+phiL(Ah,rho)+etaL(rho); hi=base+phiU(Ah,rho)+etaU(rho)
            la=log(act); cnt+=1
            if la<lo-1e-6 or la>hi+1e-6: viol+=1
            worst=min(worst,la-lo,hi-la)
print('V1 new factors: %d checks, %d violations, worst slack %.2e'%(cnt,viol,worst))
# V5': main-arc identity at 1/4 for TQ1 = q^-1 Th(-q^3)Th(q^2)^2/Th(q)^2
from farey_S1_N0 import Lambda1, C0
C=C0*cmath.exp(-3j*pi/16)
for z in (0.05+0.0j,0.05-0.1j,0.03+0.05j,0.08-0.2j):
    q=1j*cmath.exp(-z); Q=q**8
    tq=theta_direct(-q**3,Q)*theta_direct(q**2,Q)**2/theta_direct(q,Q)**2/q
    main=C*z**-0.5*cmath.exp(7*z/16+pi*pi/(64*z))
    rho=exp(-pi*pi*(1/z).real/2)
    bnd=C0*z.real**-0.5*exp(7*z.real/16)*rho**(-1/32)*math.expm1(Lambda1(rho))
    print(f'V5 z={z}: |TQ/main-1|={abs(tq/main-1):.2e}  |TQ-main|={abs(tq-main):.3e} <= bound {bnd:.3e}: {abs(tq-main)<=bnd}')
# V7': actual |c(n)-M(n)| vs bound
from farey_S1_N0 import bound_parts
getcontext().prec=120
PI=Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982')
D=pickle.load(open('big.pkl','rb'))['S1']
def dcos(r):
    a=PI*Decimal(r.numerator)/Decimal(r.denominator); s=Decimal(1); t=Decimal(1)
    for i in range(1,100): t=-t*a*a/((2*i-1)*(2*i)); s+=t
    return s
for n in (10652,20000,60000):
    N=Decimal(n)+Decimal(7)/Decimal(16); x=PI/4*N.sqrt()
    M=dcos(Fr(-3,16)-Fr(n%4,2))*(x.exp()+(-x).exp())/2/N.sqrt()
    E=abs(Decimal(D[n])-M)
    print(f'V7 n={n}: log|actual error|={float(E.ln()):.1f}  log B(n)={bound_parts(n)[0]:.1f}')
