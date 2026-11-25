# high-precision check of |f - main| near the cusp 1/4 for S0/S1 (n=3000 parameters)
import pickle, math, sys
from decimal import Decimal as Dm, getcontext
getcontext().prec=70
from review_verify import PI, dcos as dcospi
import farey_N0, farey_S1_N0
D=pickle.load(open('big.pkl','rb'))
def dcos(x):
    s=Dm(0);t=Dm(1);k=0
    while abs(t)>Dm(10)**-68: s+=t;k+=2;t=-t*x*x/(k*(k-1))
    return s
def dsin(x):
    s=Dm(0);t=x;k=1
    while abs(t)>Dm(10)**-68: s+=t;k+=2;t=-t*x*x/(k*(k-1))
    return s
def cexp(a,b): e=a.exp(); return e*dcos(b), e*dsin(b)
def cmul(x,y): return x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0]
def csqrt_inv(a,b):  # z^{-1/2}, principal, Re z>0
    r=(a*a+b*b).sqrt(); ang=Dm(math.atan2(float(b),float(a)))  # refine angle by Newton
    for _ in range(4):
        # solve tan(ang)=b/a
        c,s=dcos(ang),dsin(ang); ang-= (s*a-c*b)/(c*a+s*b)
    m=1/r.sqrt(); return cexp(-(r.sqrt()).ln()*0, -ang/2) if False else (m*dcos(-ang/2), m*dsin(-ang/2))
for nm,n in (('S0',3000),('S1',3000),('S0',1500)):
    c=D[nm]; t0=Dm(repr(math.pi))  # placeholder
    t0=PI/(16*Dm(n).sqrt())
    mod=farey_N0 if nm=='S0' else farey_S1_N0
    lsup=mod.bound_parts(n)[2]-n*float(t0)
    ph=Dm(9)/16 if nm=='S0' else Dm(-3)/16; s=Dm(1)/16 if nm=='S0' else Dm(-7)/16
    Cm=PI.sqrt()/2; C=(Cm*dcos(PI*ph), Cm*dsin(PI*ph))
    a=PI*PI/64; worst=-1e9
    for phi_f in [-0.003,-0.001,-0.0003,0,0.0002,0.0008,0.0015,0.003]:
        phi=Dm(phi_f); theta=PI/2+phi
        r=(-t0).exp(); w=(r*dcos(theta), r*dsin(theta)); p=(Dm(1),Dm(0)); S=[Dm(0),Dm(0)]
        for cj in c:
            if cj: S[0]+=cj*p[0]; S[1]+=cj*p[1]
            p=cmul(p,w)
        zr,zi=t0,-phi; den=zr*zr+zi*zi
        # exponent -s z + a/z
        er=-s*zr+a*zr/den; ei=-s*zi-a*zi/den
        m=cmul(cmul(C,csqrt_inv(zr,zi)),cexp(er,ei))
        diff=((S[0]-m[0])**2+(S[1]-m[1])**2).sqrt()
        rr=float(diff.ln())-lsup; worst=max(worst,rr)
        print(f'{nm} n={n} phi={phi_f:+.4f} |main|={float((m[0]**2+m[1]**2).sqrt()):.3e} |f-main|={float(diff):.4e} log(|f-main|/sup)={rr:+.2f}')
    print(nm,n,'worst',worst)
