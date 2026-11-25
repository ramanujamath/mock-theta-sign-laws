N=1500
def mul(a,b):
    r=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(0,N+1-i):
                if b[j]: r[i+j]+=x*b[j]
    return r
def inv(a):  # a[0]=1
    r=[0]*(N+1); r[0]=1
    for n in range(1,N+1):
        s=0
        for k in range(1,n+1):
            if a[k]: s+=a[k]*r[n-k]
        r[n]=-s
    return r
def series(term_num, den_factor):
    # sum_n q^{n^2} * prod_{k<=n} num_k / den_k
    tot=[0]*(N+1); cur=[0]*(N+1); cur[0]=1  # running product
    n=0
    while n*n<=N:
        if n>0:
            num,den=term_num(n),den_factor(n)
            cur=mul(cur,num); cur=mul(cur,inv(den))
        for i in range(N+1-n*n): tot[i+n*n]+=cur[i]
        n+=1
    return tot
def poly(d):
    p=[0]*(N+1)
    for e,c in d.items():
        if e<=N: p[e]+=c
    return p
one=poly({0:1})
# phi(q)=sum q^{n^2}/(-q^2;q^2)_n
phi=series(lambda n:one, lambda n:poly({0:1,2*n:1}))
# chi(q)=sum q^{n^2}(-q;q)_n/(-q^3;q^3)_n ; (-q;q)_n/(-q^3;q^3)_n = prod 1/(1-q^k+q^2k)
chi=series(lambda n:one, lambda n:poly({0:1,n:-1,2*n:1}))
import pickle; pickle.dump((phi,chi),open('mt.pkl','wb'))
for name,c in [('phi',phi),('chi',chi)]:
    print(name,c[:30])
    for m in range(4):
        for r in range(m):
            pass
    # sign patterns mod 2,3,4,6
    for M in (2,3,4,6,8,12):
        pats=[]
        for r in range(M):
            vals=[c[n] for n in range(r,N+1,M) if n>200]
            s=set((v>0)-(v<0) for v in vals)
            pats.append(''.join(sorted({1:'+',-1:'-',0:'0'}[x] for x in s)))
        print(' mod',M,pats)
