import math, cmath, random, pickle
from math import pi, sqrt, exp, log, gcd
from fractions import Fraction as Fr
import farey_gen as fg
from farey_lib import theta_direct, B2, phiU, phiL, etaU, etaL
from u0b import U0 as U0s, U1 as U1s
from radial import mu as mus
random.seed(7)
def logJ(m,q):
    s=0; n=1
    while True:
        t=q**(m*n); s+=cmath.log(1-t)
        if abs(t)<1e-18: return s.real
        n+=1
def rand_point(t0):
    K=int(sqrt(2*pi/t0)); k=random.randint(1,K)
    h=random.choice([x for x in range(k) if gcd(x,k)==1])
    phi=random.uniform(-2*pi/(k*(K+1)),2*pi/(k*(K+1)))
    return h,k,complex(t0,-phi)
# ---- 1. two-sided factor bounds (pointwise, exact exponent) ----
facs=[('th',4,Fr(0),Fr(1,2)),('th',4,Fr(3,4),Fr(0)),('th',4,Fr(1,2),Fr(1,2)),('th',4,Fr(1,4),Fr(0)),('eta',1),('eta',2),('eta',4)]
bad=0; cnt=0; worst=9
for t0 in (0.02,0.008):
    for _ in range(1500):
        h,k,z=rand_point(t0); q=cmath.exp(2j*pi*h/k)*cmath.exp(-z)
        X=(1/z).real; w=X/k**2; kz=k*abs(z)
        for f in facs:
            if f[0]=='eta':
                m=f[1]; act=logJ(m,q); ff=('eta',m,1)
            else:
                _,M,A,s=f; x=cmath.exp(2j*pi*s)*q**int(A*M); v=theta_direct(x,q**M)
                if abs(v)<1e-280: continue
                act=log(abs(v)); ff=('th',M,A,s,1)
            R,pr,p,qp,prod=fg.factor_data(ff,h,k)
            g=gcd(ff[1],k); M=ff[1]
            base=pr-0.5*log(kz)-float(qp)*t0+(pi*pi/64)*float(R)/k**2*X
            # lower bound: replace upper products by lower
            rho=exp(-4*pi*pi*g*g*w/M)
            if ff[0]=='eta': lo=base+etaL(rho); hi=base+etaU(rho)
            else:
                _,_,A,s,_=ff; kp=k//g; hp=M*h//g; Ah=fg.ahat(A,s,hp,kp)
                lo=base+phiL(Ah,rho)+etaL(rho); hi=base+phiU(Ah,rho)+etaU(rho)
            cnt+=1
            # direct theta sum has roundoff floor ~ 1e-16 * max term; skip if act < -30 for theta (unreliable)
            reliable = ff[0]=='eta' or act>-30
            if reliable and (act<lo-1e-6 or act>hi+1e-6): bad+=1; print('VIOL',f,h,k,z,act,lo,hi)
            if reliable: worst=min(worst,act-lo,hi-act)
print(f'1) factor two-sided bounds: {cnt} checks, {bad} violations (reliable points), worst slack {worst:.2e}')
# ---- 2. pointwise bound for the whole functions (min over reps) vs direct series ----
fun={'mu2':lambda q:mus(q)[0],'U1':U1s,'U0':U0s}
for name,f in fun.items():
    mx=-1e9; viol=0; cnt=0
    for t0 in (0.03,0.015):
        for _ in range(400):
            h,k,z=rand_point(t0)
            if k==fg.MAIN[name]['k']: continue
            q=cmath.exp(2j*pi*h/k)*cmath.exp(-z); X=(1/z).real
            try: act=log(abs(f(q)))
            except Exception: continue
            best=None
            for rep in fg.reps(name):
                lv=[]
                for coeff,fs,nb in rep:
                    lb,R=fg.term_bound(fs,h,k,t0,X/k**2,coeff,nb,kz_lo=k*abs(z),kz_hi=k*abs(z))
                    lv.append(lb+(pi*pi/64)*float(R)/k**2*X)
                m=max(lv); v=m+log(sum(exp(x-m) for x in lv))
                best=v if best is None else min(best,v)
            cnt+=1; mx=max(mx,act-best)
            if act>best: viol+=1
    print(f'2) {name}: {cnt} minor-arc points, {viol} violations; max log(actual/bound) = {mx:.2f}')
# ---- 3. actual coefficient error vs bound ----
from decimal import Decimal, getcontext
getcontext().prec=120
PI=Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982')
def dcos(r):
    a=PI*Decimal(r.numerator)/Decimal(r.denominator); s=Decimal(1); t=Decimal(1)
    for i in range(1,100): t=-t*a*a/((2*i-1)*(2*i)); s+=t
    return s
D=pickle.load(open('big.pkl','rb')); F=pickle.load(open('family.pkl','rb'))
coef={'mu2':D['mu2'],'U1':F['U1'],'U0':F['U0']}
for name in ('mu2','U1','U0'):
    E,ch=fg.minor_setup(name)
    for n in (500,1000,5180,10000,19999):
        N=Decimal(n)-Decimal(1)/Decimal(8)
        if name=='U0':
            x=PI/(4*Decimal(2).sqrt())*N.sqrt(); tr=dcos(Fr(5,16)-Fr(n%8,4))+dcos(Fr(13,16)+Fr(3*(n%8),4))
        else:
            x=PI/(2*Decimal(2).sqrt())*N.sqrt()
            tr=(Decimal(2).sqrt()*dcos(Fr(-5,8)-Fr(n%4,2))) if name=='mu2' else dcos(Fr(3,8)-Fr(n%4,2))/Decimal(2).sqrt()
        M=tr*(x.exp()+(-x).exp())/2/N.sqrt()
        err=abs(Decimal(coef[name][n])-M)
        le=float(err.ln()) if err>0 else -999
        print(f'3) {name} n={n}: log|actual error|={le:.1f}  log bound={fg.bound_parts(name,n,E,ch)[0]:.1f}  log lower main={fg.log_main_lower(name,n):.1f}')
