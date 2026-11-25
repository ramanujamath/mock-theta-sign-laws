# Exact checks in Q(x) with minimal polynomials:
#  S0,S1: x^2 = sqrt2-1  => x^4 + 2x^2 - 1 = 0
#  mu2  : x^2 = 1/3
from fractions import Fraction as F
import math, cmath
def red(p,mp):  # p: list coeffs low->high; mp monic list low->high
    p=p[:]; d=len(mp)-1
    for i in range(len(p)-1,d-1,-1):
        c=p[i]
        if c:
            for j in range(d+1): p[i-d+j]-=c*mp[j]
    return [c for c in p[:d]]+[F(0)]*max(0,d-len(p))
def mul(a,b):
    r=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def add(*ps):
    n=max(map(len,ps)); return [sum((p[i] if i<len(p) else 0) for p in ps) for i in range(n)]
def sc(c,p): return [c*x for x in p]
mp=[F(-1),F(0),F(2),F(0),F(1)]
X=[F(0),F(1)]; one=[F(1)]; X2=mul(X,X)
# S0 bracket B = (1-x-x^2) + i x ; S1 bracket B = (1-x+x^2) - i x (times x)
for name,re_,im_,extra in [('S0',add(one,sc(-1,X),sc(-1,X2)),X,one),('S1',add(one,sc(-1,X),X2),sc(-1,X),X)]:
    # amplitude^2 claim: |extra*B|^2 = 2 (1-x)^2 (1+x)
    lhs=red(mul(mul(extra,extra),add(mul(re_,re_),mul(im_,im_))),mp)
    rhs=red(sc(2,mul(mul(add(one,sc(-1,X)),add(one,sc(-1,X))),add(one,X))),mp)
    print(name,'|amp|^2 identity holds:',lhs==rhs)
    # phase: total = pi/8 - (1/2)atan x + arg B ; claim = target. Check via W = B^2 (1 - i x): arg W = 2(target - pi/8)
    Bre,Bim=re_,im_
    W2re=add(mul(Bre,Bre),sc(-1,mul(Bim,Bim))); W2im=sc(2,mul(Bre,Bim))
    Wre=add(W2re,mul(W2im,X)); Wim=add(W2im,sc(-1,mul(W2re,X)))
    Wre,Wim=red(Wre,mp),red(Wim,mp)
    x=math.sqrt(math.sqrt(2)-1)
    ev=lambda p:sum(float(c)*x**i for i,c in enumerate(p))
    print('  arg W/pi =',cmath.phase(complex(ev(Wre),ev(Wim)))/math.pi, ' tan check Im+x^2 Re (S0) / Im - x^2 Re... ->',
          red(add(Wim,mul(X2,Wre)),mp) if name=='S0' else red(add(Wim,sc(-1,mul(add(one,sc(-1,X2)),Wre))),mp) )
# S1 exact: tan(-5pi/8)=1+sqrt2=1/x^2  <=> x^2*Im - Re = 0
re_,im_=add(one,sc(-1,X),X2),sc(-1,X)
W2re=add(mul(re_,re_),sc(-1,mul(im_,im_))); W2im=sc(2,mul(re_,im_))
Wre=red(add(W2re,mul(W2im,X)),mp); Wim=red(add(W2im,sc(-1,mul(W2re,X))),mp)
print('S1 tan identity:',red(add(mul(X2,Wim),sc(-1,Wre)),mp))
