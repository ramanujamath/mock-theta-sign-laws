import pickle, cmath, math
from math import pi, sqrt
N=20000
def mulb(a,m,s):
    for i in range(N,m-1,-1): a[i]+=s*a[i-m]
def divb(a,m,s):
    for i in range(m,N+1): a[i]-=s*a[i-m]
# TH = J2^8 / (J1^3 J4^4), J_m = prod (1-q^{mk})
th=[0]*(N+1); th[0]=1
for k in range(1,N+1):
    if 2*k<=N:
        for _ in range(8): mulb(th,2*k,-1)
    for _ in range(3): divb(th,k,-1)
    if 4*k<=N:
        for _ in range(4): divb(th,4*k,-1)
pickle.dump(th,open('TH.pkl','wb'))
F=pickle.load(open('family.pkl','rb')); D=pickle.load(open('big.pkl','rb'))
U0,U1,mu=F['U0'],F['U1'],D['mu2']
print('U0 + 2U1 == TH  through q^%d:'%N, all(U0[n]+2*U1[n]==th[n] for n in range(N+1)))
print('mu2 == 2U0 - TH through q^%d:'%N, all(mu[n]==2*U0[n]-th[n] for n in range(N+1)))
# The constants of TH at the main cusps are computed (with logarithms, avoiding underflow) in th_const.py.
