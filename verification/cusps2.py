from cusps import theta_exp, TH, ALden
from math import gcd
FB={'S0':[(-1,2,8),(-1,6,8)],'S1':[(-1,2,8),(-1,6,8)],'mu2':[(1,1,4),(1,3,4)]}   # Theta(s q^a; q^M) for z's in form B
K=120
for name in TH:
    worst=[];atfour=[]
    for k in range(1,K+1):
        for h in range(k):
            if gcd(h,k)!=1: continue
            eth=sum(p*theta_exp(s,a,M,h,k) for p,s,a,M in TH[name])
            _,s,a,M=ALden[name]; eAL=-theta_exp(s,a,M,h,k)
            A=max(eth,eAL)
            B=max(-theta_exp(s,a,M,h,k) for s,a,M in FB[name])
            if k==4: atfour.append((f'{h}/{k}',round(eth,4),round(eAL,4),round(B,4)))
            else: worst.append((min(A,B),A,B,f'{h}/{k}'))
    worst.sort(reverse=True)
    print(name,'at k=4 (theta part, AL crude, formB crude):',atfour)
    print('   worst other cusps (best, formA, formB):',[(round(b,4),round(a,4),round(c,4),s) for b,a,c,s in worst[:5]])
