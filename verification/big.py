import sys, pickle
N=int(sys.argv[1])
def mulb(a,m,s):   # a*(1+s q^m)
    r=a[:]
    for i in range(N,m-1,-1): r[i]+=s*a[i-m]
    return r
def divb(a,m,s):   # a/(1+s q^m)
    r=a[:]
    for i in range(m,N+1): r[i]-=s*r[i-m]
    return r
def run(expo,step,sgn):
    tot=[0]*(N+1); cur=[0]*(N+1); cur[0]=1; n=0
    while expo(n)<=N:
        if n>0: cur=step(cur,n)
        e=expo(n); s=sgn(n)
        for i in range(N+1-e): tot[i+e]+=s*cur[i]
        n+=1
    return tot
st8=lambda c,k: divb(mulb(c,2*k-1,1),2*k,1)
S0=run(lambda n:n*n,st8,lambda n:1)
S1=run(lambda n:n*(n+2),st8,lambda n:1)
mu=run(lambda n:n*n,lambda c,k: divb(divb(mulb(c,2*k-1,-1),2*k,1),2*k,1),lambda n:(-1)**n)
pickle.dump({'S0':S0,'S1':S1,'mu2':mu},open('big.pkl','wb'))
print(S0[:12],S1[:12],mu[:12])
