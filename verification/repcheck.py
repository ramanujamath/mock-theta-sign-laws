import pickle
N=3000
def zero(): return [0]*(N+1)
def inv(a):
    r=zero(); r[0]=1; assert a[0]==1
    nz=[(k,a[k]) for k in range(1,N+1) if a[k]]
    for n in range(1,N+1):
        r[n]=-sum(v*r[n-k] for k,v in nz if k<=n)
    return r
def mul(a,b):
    r=zero(); nzb=[(j,y) for j,y in enumerate(b) if y]
    for i,x in enumerate(a):
        if x:
            for j,y in nzb:
                if i+j>N: break
                r[i+j]+=x*y
    return r
def theta(zs,zc,M):   # Theta(zs*q^zc ; q^M) = sum (-1)^n q^{M n(n-1)/2} (zs q^zc)^n
    r=zero()
    for n in range(-200,201):
        e=M*n*(n-1)//2+zc*n
        if 0<=e<=N: r[e]+=(-1)**n*zs**n
    return r
def m_series(xs,xc,M,zs,zc,shift=0):
    # m(xs q^xc, q^M, zs q^zc) * q^shift as power series (shift may be negative if result divisible)
    L=[0]*(N+1+50)
    off=50
    for r in range(-200,201):
        e0=M*r*(r-1)//2+zc*r; c0=(-1)**r*zs**r
        # 1/(1 - eps q^e), eps q^e = q^{M(r-1)} * xs q^xc * zs q^zc
        e=M*(r-1)+xc+zc; eps=xs*zs
        if e>0:
            k=0
            while e0+k*e+shift<=N:
                idx=e0+k*e+shift+off
                if idx>=0: L[idx]+=c0*eps**k
                k+=1
        elif e<0:
            # 1/(1-eps q^e) = -eps^{-1} q^{-e} /(1 - eps^{-1} q^{-e});  eps=+-1
            k=1
            while e0-k*e+shift<=N:
                idx=e0-k*e+shift+off
                if idx>=0: L[idx]+=-c0*eps**(-k) if abs(eps)==1 else None
                k+=1
        else: raise ValueError('pole')
    assert all(v==0 for v in L[:off]), 'negative powers remain'
    return mul(L[off:off+N+1],inv(theta(zs,zc,M)))
D=pickle.load(open('big.pkl','rb'))
# S1 = -q^{-1}[m(-q,q^8,-q^2)+m(-q,q^8,-q^6)]
A=m_series(-1,1,8,-1,2,shift=-1); B=m_series(-1,1,8,-1,6,shift=-1)
S1=[-(a+b) for a,b in zip(A,B)]
print('S1 rep matches through q^%d:'%(N-1), S1[:N]==D['S1'][:N])
# mu2 = 2[m(-q,q^4,q)+m(-q,q^4,q^3)]
A=m_series(-1,1,4,1,1); B=m_series(-1,1,4,1,3)
mu=[2*(a+b) for a,b in zip(A,B)]
print('mu2 rep matches through q^%d:'%(N-1), mu[:N]==D['mu2'][:N])
