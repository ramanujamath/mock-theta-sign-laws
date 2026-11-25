import cmath, math
from radial import S, mu
def theta(x,q,R=60):   # Θ(x;q)=Σ(-1)^n q^{n(n-1)/2} x^n
    return sum((-1)**n * q**(n*(n-1)//2) * x**n for n in range(-R,R+1))
def m(x,q,z,R=60):
    return sum((-1)**r * q**(r*(r-1)//2) * z**r/(1-q**(r-1)*x*z) for r in range(-R,R+1))/theta(z,q)
def Th(a,mm,q,neg=False): return theta((-1 if neg else 1)*q**a, q**mm)
def Pm(mm,q,K=4000):  # Θ_m = prod(1-q^{mk})
    p=1
    for k in range(1,K):
        p*=1-q**(mm*k)
        if abs(q**(mm*k))<1e-18: break
    return p
def S0_al(q): return 2*m(-q**3,q**8,-1)+q*Th(1,8,q,True)*Th(2,8,q)**2/Th(3,8,q)**2
def S1_al(q): return -(2/q)*m(-q,q**8,-1)+(1/q)*Th(3,8,q,True)*Th(2,8,q)**2/Th(1,8,q)**2
def mu_al(q): return 4*m(-q,q**4,-1)-Pm(2,q)**8/(Pm(1,q)**3*Pm(4,q)**4)
q0=0.3+0.25j
print('S0',S(q0,0)[0],S0_al(q0)); print('S1',S(q0,2)[0],S1_al(q0)); print('mu',mu(q0)[0],mu_al(q0))
