import sys,pickle
N=int(sys.argv[1])
def mulb(a,m,s):
    r=a[:]
    for i in range(N,m-1,-1): r[i]+=s*a[i-m]
    return r
def divb(a,m,s):
    r=a[:]
    for i in range(m,N+1): r[i]-=s*r[i-m]
    return r
def run(expo,step,init=lambda c:c):
    tot=[0]*(N+1); cur=[0]*(N+1); cur[0]=1; cur=init(cur); n=0
    while expo(n)<=N:
        if n>0: cur=step(cur,n)
        e=expo(n)
        for i in range(N+1-e): tot[i+e]+=cur[i]
        n+=1
    return tot
F={}
# T0: q^{(n+1)(n+2)} (-q^2;q^2)_n/(-q;q^2)_{n+1};  T1: q^{n(n+1)} same
stT=lambda c,k: divb(mulb(c,2*k,1),2*k+1,1)
F['T0']=run(lambda n:(n+1)*(n+2),stT,lambda c:divb(c,1,1))
F['T1']=run(lambda n:n*(n+1),stT,lambda c:divb(c,1,1))
F['U0']=run(lambda n:n*n,lambda c,k:divb(mulb(c,2*k-1,1),4*k,1))
F['U1']=run(lambda n:(n+1)**2,lambda c,k:divb(mulb(c,2*k-1,1),4*k+2,1),lambda c:divb(c,2,1))
V=run(lambda n:n*n,lambda c,k:divb(mulb(c,2*k-1,1),2*k-1,-1)); V=[2*x for x in V]; V[0]-=1; F['V0']=V
F['V1']=run(lambda n:(n+1)**2,lambda c,k:divb(mulb(c,2*k-1,1),2*k+1,-1),lambda c:divb(c,1,-1))
pickle.dump(F,open('family.pkl','wb'))
for k,c in F.items(): print(k,c[:20])
