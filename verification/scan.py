N=1200
def mul(a,b):
    r=[0]*(N+1)
    nz=[(j,y) for j,y in enumerate(b) if y]
    for i,x in enumerate(a):
        if x:
            for j,y in nz:
                if i+j>N: break
                r[i+j]+=x*y
    return r
def inv(a):
    r=[0]*(N+1); r[0]=1; nz=[(k,a[k]) for k in range(1,N+1) if a[k]]
    for n in range(1,N+1):
        s=0
        for k,v in nz:
            if k>n: break
            s+=v*r[n-k]
        r[n]=-s
    return r
def P(d):
    p=[0]*(N+1)
    for e,c in d.items():
        if e<=N: p[e]+=c
    return p
def shift(a,e):
    return [0]*e+a[:N+1-e] if e<=N else [0]*(N+1)
def gsum(expo, num, den, maxn=60, sgn=lambda n:1):
    # sum_n sgn(n) q^{expo(n)} prod_{k=1..n} num(k)/den(k) ; num/den lists of polys per step
    tot=[0]*(N+1); cur=P({0:1}); n=0
    while expo(n)<=N and n<maxn:
        if n>0:
            for p in num(n): cur=mul(cur,p)
            for p in den(n): cur=mul(cur,inv(p))
        t=shift(cur,expo(n)); s=sgn(n)
        tot=[x+s*y for x,y in zip(tot,t)]; n+=1
    return tot
F={}
# 5th order
F['f0']=gsum(lambda n:n*n,lambda k:[],lambda k:[P({0:1,k:1})])
F['f1']=gsum(lambda n:n*n+n,lambda k:[],lambda k:[P({0:1,k:1})])
# 6th order phi6 = sum (-1)^n q^{n^2} (q;q^2)_n/(-q;q)_{2n}
F['phi6']=gsum(lambda n:n*n,lambda k:[P({0:1,2*k-1:-1})],lambda k:[P({0:1,2*k-1:1}),P({0:1,2*k:1})],sgn=lambda n:(-1)**n)
# 6th psi6 = sum (-1)^n q^{(n+1)^2}(q;q^2)_n/(-q;q)_{2n+1}
t=gsum(lambda n:(n+1)**2,lambda k:[P({0:1,2*k-1:-1})],lambda k:[P({0:1,2*k:1}),P({0:1,2*k+1:1})],sgn=lambda n:(-1)**n)
F['psi6']=mul(t,inv(P({0:1,1:1})))
# 2nd order A(q)=sum q^{n+1}(-q^2;q^2)_n/(q;q^2)_{n+1}  -> positive presumably; B, mu2
F['mu2']=gsum(lambda n:n*n,lambda k:[P({0:1,2*k-1:-1})],lambda k:[P({0:1,2*k:1}),P({0:1,2*k:1})],sgn=lambda n:(-1)**n)
# 3rd order f
F['f3']=gsum(lambda n:n*n,lambda k:[],lambda k:[P({0:1,k:1}),P({0:1,k:1})])
# 8th order S0 = sum q^{n^2}(-q;q^2)_n/(-q^2;q^2)_n ; T0 = sum q^{(n+1)(n+2)}(-q^2;q^2)_n/(-q;q^2)_{n+1}
F['S0']=gsum(lambda n:n*n,lambda k:[P({0:1,2*k-1:1})],lambda k:[P({0:1,2*k:1})])
F['S1']=gsum(lambda n:n*(n+2),lambda k:[P({0:1,2*k-1:1})],lambda k:[P({0:1,2*k:1})])
import pickle; pickle.dump(F,open('F.pkl','wb'))
for name,c in F.items():
    best=None
    for M in range(1,25):
        pat=[];ok=True
        for r in range(M):
            s={(v>0)-(v<0) for n,v in enumerate(c) if n%M==r and n>N//3}
            s.discard(0)
            if len(s)!=1: ok=False;break
            pat.append('+' if 1 in s else '-')
        if ok: best=(M,''.join(pat));break
    print(name,c[:16],'period/pattern:',best)
