import cmath, random
def theta(x,Q,R=30): return sum((-1)**n*Q**(n*(n-1)//2)*x**n for n in range(-R,R+1))
def m(x,Q,z,R=30): return sum((-1)**r*Q**(r*(r-1)//2)*z**r/(1-Q**(r-1)*x*z) for r in range(-R,R+1))/theta(z,Q)
def poch(a,q,n):
    p=1
    for j in range(n): p*=1-a*q**j
    return p
def ser(expo,term,N=200):
    return sum(q_**expo(n)*term(n) for n in range(N))
def J(mm,q,K=3000):
    p=1
    for k in range(1,K): p*=1-q**(mm*k)
    return p
random.seed(1); worst={}
for trial in range(12):
    r=random.uniform(0.2,0.75); th=random.uniform(0,6.283); q=r*cmath.exp(1j*th); q_=q
    S0=sum(q**(n*n)*poch(-q,q*q,n)/poch(-q*q,q*q,n) for n in range(60))
    S1=sum(q**(n*(n+2))*poch(-q,q*q,n)/poch(-q*q,q*q,n) for n in range(60))
    mu=sum((-1)**n*q**(n*n)*poch(q,q*q,n)/poch(-q*q,q*q,n)**2 for n in range(60))
    U0=sum(q**(n*n)*poch(-q,q*q,n)/poch(-q**4,q**4,n) for n in range(60))
    U1=sum(q**((n+1)**2)*poch(-q,q*q,n)/poch(-q*q,q**4,n+1) for n in range(60))
    Q8=q**8; Q4=q**4
    Th=lambda a,neg=False,M=8: theta((-1 if neg else 1)*q**a, q**M)
    TH=J(2,q)**8/(J(1,q)**3*J(4,q)**4)
    checks={'S0A':S0-(2*m(-q**3,Q8,-1)+q*Th(1,True)*Th(2)**2/Th(3)**2),
            'S0B':S0-(m(-q**3,Q8,-q**2)+m(-q**3,Q8,-q**6)),
            'S1A':S1-(-(2/q)*m(-q,Q8,-1)+(1/q)*Th(3,True)*Th(2)**2/Th(1)**2),
            'S1B':S1-(-(1/q)*(m(-q,Q8,-q**2)+m(-q,Q8,-q**6))),
            'U0a':U0-2*m(-q,Q4,-1),'U0b':U0-2*m(-q,Q4,q**3),
            'U1a':U1+m(-q,Q4,-q**2),'U1b':U1+m(-q,Q4,q),
            'TH':U0+2*U1-TH,'mu':mu-(2*U0-TH)}
    for k,v in checks.items(): worst[k]=max(worst.get(k,0),abs(v))
print({k:'%.1e'%v for k,v in worst.items()})
