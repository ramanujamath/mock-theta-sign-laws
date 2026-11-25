from al import m, Th, Pm
from radial import S, mu
import itertools
q0=0.31+0.22j
targets={'S0':S(q0,0)[0],'S1':S(q0,2)[0],'mu2':mu(q0)[0]}
print('S0 form B:', m(-q0**3,q0**8,-q0**2)+m(-q0**3,q0**8,-q0**6), targets['S0'])
# search S1 = c*q^e*[m(-q,q^8,-q^a)+m(-q,q^8,-q^b)]
for a,b in itertools.combinations_with_replacement(range(0,8),2):
    for c in (1,-1,2,-2):
        for e in (-1,0):
            v=c*q0**e*(m(-q0,q0**8,-q0**a)+m(-q0,q0**8,-q0**b))
            if abs(v-targets['S1'])<1e-10: print('S1 =',c,'q^',e,'[m(-q,q^8,-q^%d)+m(-q,q^8,-q^%d)]'%(a,b))
for a,b in itertools.combinations_with_replacement(range(0,4),2):
    for c in (1,-1,2,-2,4,-4):
        for sa in (1,-1):
          for sb in (1,-1):
            v=c*(m(-q0,q0**4,sa*q0**a)+m(-q0,q0**4,sb*q0**b))
            if abs(v-targets['mu2'])<1e-10: print('mu2 =',c,'[m(-q,q^4,%+d q^%d)+m(-q,q^4,%+d q^%d)]'%(sa,a,sb,b))
