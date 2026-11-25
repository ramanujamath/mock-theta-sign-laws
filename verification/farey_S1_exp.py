from fractions import Fraction as Fr
from farey_lib import Ahat, B2
# S1 factor lists (mult, A=a/8, sigma)
TQ1 =[(1,Fr(3,8),Fr(1,2)),(2,Fr(1,4),Fr(0)),(-2,Fr(1,8),Fr(0))]  # q^-1 Th(-q^3)Th(q^2)^2/Th(q)^2
ALA =[(-1,Fr(0),Fr(1,2))]
ALB1=[(-1,Fr(1,4),Fr(1,2))]; ALB2=[(-1,Fr(3,4),Fr(1,2))]
def R(fl,hm,kpar):
    s=Fr(sum(m for m,_,_ in fl),24)+sum(m*B2(Ahat(A,s_,hm,kpar))/2 for m,A,s_ in fl)
    return -32*s
CL=[(1,h) for h in range(8)]+[(2,h) for h in (1,3,5,7)]
def kmin(kpar,hm): return 3 if (kpar==1 and hm in (2,6)) else (1 if kpar==1 else 2)
print('class   R_TQ  R_ALA  R_B1  R_B2   bestE(k\'>=kmin, excluding k=4)')
Eall=0
for kpar,hm in CL:
    a=(R(TQ1,hm,kpar),R(ALA,hm,kpar)); b=(R(ALB1,hm,kpar),R(ALB2,hm,kpar))
    km=kmin(kpar,hm)
    eA=max(a)/km**2; eB=max(b)/km**2
    best=min(eA,eB); Eall=max(Eall,best)
    print((kpar,hm),[str(x) for x in a+b],'A:',eA,'B:',eB,'->',best)
print('Eall =',Eall)
print('k=4 class (odd,2):',R(TQ1,2,1),R(ALA,2,1))
