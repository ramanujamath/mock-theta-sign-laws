# independent recomputation of the q-series (direct definitions) up to N
import pickle
N=1200
def mul(a,b):
    r=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(N+1-i):
                if b[j]: r[i+j]+=x*b[j]
    return r
def inv(a):
    r=[0]*(N+1); r[0]=1
    for n in range(1,N+1):
        r[n]=-sum(a[k]*r[n-k] for k in range(1,n+1) if a[k])
    return r
def binom(e,s): # 1+s q^e
    p=[0]*(N+1); p[0]=1
    if e<=N: p[e]+=s
    return p
def poch(s,a,step,n):  # prod_{j=0}^{n-1}(1+s q^{a+step j})
    p=[0]*(N+1); p[0]=1
    for j in range(n): p=mul(p,binom(a+step*j,s))
    return p
def ser(expo,num,den,sign=lambda n:1):
    tot=[0]*(N+1); n=0
    while expo(n)<=N:
        t=mul(num(n),inv(den(n)))
        e=expo(n)
        for i in range(N+1-e): tot[i+e]+=sign(n)*t[i]
        n+=1
    return tot
S0=ser(lambda n:n*n,lambda n:poch(1,1,2,n),lambda n:poch(1,2,2,n))
S1=ser(lambda n:n*(n+2),lambda n:poch(1,1,2,n),lambda n:poch(1,2,2,n))
mu=ser(lambda n:n*n,lambda n:poch(-1,1,2,n),lambda n:mul(poch(1,2,2,n),poch(1,2,2,n)),lambda n:(-1)**n)
U0=ser(lambda n:n*n,lambda n:poch(1,1,2,n),lambda n:poch(1,4,4,n))
U1=ser(lambda n:(n+1)**2,lambda n:poch(1,1,2,n),lambda n:poch(1,2,4,n+1))
D=pickle.load(open('big.pkl','rb')); F=pickle.load(open('family.pkl','rb'))
for nm,a,b in (('S0',S0,D['S0']),('S1',S1,D['S1']),('mu2',mu,D['mu2']),('U0',U0,F['U0']),('U1',U1,F['U1'])):
    print(nm, a[:N+1]==b[:N+1])
