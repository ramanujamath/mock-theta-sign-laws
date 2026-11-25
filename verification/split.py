import math, cmath
from al import *
pi=math.pi
for t in (0.04,0.02,0.01):
    q=1j*math.exp(-t)
    A=2*m(-q**3,q**8,-1,R=200); T=q*Th(1,8,q,True)*Th(2,8,q)**2/Th(3,8,q)**2
    A1=-(2/q)*m(-q,q**8,-1,R=200); T1=(1/q)*Th(3,8,q,True)*Th(2,8,q)**2/Th(1,8,q)**2
    Am=4*m(-q,q**4,-1,R=200); Tm=-Pm(2,q,20000)**8/(Pm(1,q,20000)**3*Pm(4,q,20000)**4)
    g=math.exp(-pi**2/(64*t)); g2=math.exp(-pi**2/(32*t))
    print(f't={t}: S0: AL*e^-={abs(A)*g:.3e} theta*e^-={abs(T)*g:.4f} | S1: AL={abs(A1)*g:.3e} th={abs(T1)*g:.4f} | mu: AL={abs(Am)*g2:.3e} th={abs(Tm)*g2:.4f}')
