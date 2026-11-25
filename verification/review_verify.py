import pickle, random, math
from decimal import Decimal as Dm, getcontext
from fractions import Fraction as Fr
getcontext().prec=130
def dpi():
    # Machin
    getcontext().prec+=10
    def arctan_inv(x):
        x=Dm(x); s=Dm(0); t=1/x; k=0; x2=x*x
        while True:
            term=t/(2*k+1)
            if term==0 or abs(term)<Dm(10)**(-getcontext().prec): break
            s+= term if k%2==0 else -term
            t/=x2; k+=1
        return s
    p=4*(4*arctan_inv(5)-arctan_inv(239)); getcontext().prec-=10; return +p
PI=dpi()
def dcos(fr):  # cos(pi*fr), fr Fraction
    fr=fr%2; x=PI*Dm(fr.numerator)/Dm(fr.denominator)
    s=Dm(0); t=Dm(1); k=0
    while abs(t)>Dm(10)**-140:
        s+=t; k+=2; t=-t*x*x/(k*(k-1))
    return s
def dcosh(x): e=x.exp(); return (e+1/e)/2
D=pickle.load(open('big.pkl','rb')); F=pickle.load(open('family.pkl','rb'))
def Mval(name,n):
    if name=='S0':
        N=Dm(n)-Dm(1)/16; x=PI*N.sqrt()/4; return dcos(Fr(9,16)-Fr(n,2))*dcosh(x)/N.sqrt(), N, x
    if name=='S1':
        N=Dm(n)+Dm(7)/16; x=PI*N.sqrt()/4; return dcos(Fr(-3,16)-Fr(n,2))*dcosh(x)/N.sqrt(), N, x
    N=Dm(n)-Dm(1)/8
    if name in('mu2','U1'):
        x=2*(PI*PI/32*N).sqrt(); c=-Dm(2).sqrt() if name=='mu2' else 1/Dm(2).sqrt()
        return c*dcos(Fr(3,8)-Fr(n,2))*dcosh(x)/N.sqrt(), N, x
    x=2*(PI*PI/128*N).sqrt()
    return (dcos(Fr(5,16)-Fr(n,4))+dcos(Fr(-13,16)-Fr(3*n,4)))*dcosh(x)/N.sqrt(), N, x
import farey_N0, farey_S1_N0, farey_gen as G
LOWER={'S0':math.cos(7*math.pi/16),'S1':math.cos(5*math.pi/16),**G.LOWER}
EC={nm:G.minor_setup(nm) for nm in ('mu2','U1','U0')}
def logB(name,n):
    if name=='S0': return farey_N0.bound_parts(n)[0]
    if name=='S1': return farey_S1_N0.bound_parts(n)[0]
    E,ch=EC[name]; return G.bound_parts(name,n,E,ch)[0]
cfs={'S0':D['S0'],'S1':D['S1'],'mu2':D['mu2'],'U0':F['U0'],'U1':F['U1']}
N0={'S0':11231,'S1':10652,'mu2':976,'U1':493,'U0':5180}
random.seed(12345)
# LOWER check over residues
for r in range(8):
    v=math.cos(5*math.pi/16-math.pi*r/4)+math.cos(13*math.pi/16+3*math.pi*r/4)
    print('U0 r=%d main coef %.4f'%(r,v))
print('U0 LOWER',G.LOWER['U0'])
for name,c in cfs.items():
    top=len(c)-1
    ns=[N0[name],N0[name]+1]+sorted(random.sample(range(N0[name],top+1),10))+[top]
    worst_ratio=-1e9; worst_marg=1e9
    for n in ns:
        M,N,x=Mval(name,n)
        err=abs(Dm(c[n])-M)
        lerr=float(err.ln()) if err>0 else -1e9
        lB=logB(name,n)
        lM=float(abs(M).ln())
        lLow=math.log(LOWER[name])+float(dcosh(x).ln())-0.5*float(N.ln())
        ok1=lerr<=lB; ok2=lB<lLow; ok3=lLow<=lM+1e-12; sgnok=(M>0)==(c[n]>0)
        worst_ratio=max(worst_ratio,lerr-lB); worst_marg=min(worst_marg,lLow-lB)
        print(f'{name} n={n}: log|c-M|={lerr:8.2f} logB={lB:8.2f} log(LOWER main)={lLow:8.2f} log|M|={lM:8.2f}  {ok1} {ok2} {ok3} sign {sgnok}')
    print(name,'max(log|c-M|-logB)=%.2f  min margin=%.5f'%(worst_ratio,worst_marg))
