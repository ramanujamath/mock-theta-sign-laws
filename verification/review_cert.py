import sys, math, time
import certify
from math import pi
nm=sys.argv[1]; lo=int(sys.argv[2]); hi=int(sys.argv[3])
if nm=='S0':
    import farey_N0 as A, farey_bound as FB
    mf=lambda n: FB.log_main_lower(n)-A.bound_parts(n)[0]
elif nm=='S1':
    import farey_S1_N0 as B
    mf=lambda n: B.log_main_lower(n)-B.bound_parts(n)[0]
else:
    import farey_gen as G
    E,ch=G.minor_setup(nm)
    mf=lambda n: G.log_main_lower(nm,n)-G.bound_parts(nm,n,E,ch)[0]
t=time.time(); bad=[]; mn=(1e9,None); prev=None; dec=0
for n in range(lo,hi+1):
    m=mf(n)
    if m<=1e-8: bad.append(n)
    if m<mn[0]: mn=(m,n)
    if prev is not None and m<prev-1e-9: dec+=1
    prev=m
print(nm,'certified scan',lo,hi,': n with margin<=1e-8:',bad[:20],'count',len(bad),' min margin',mn,' decreases',dec,' time %.0fs'%(time.time()-t),flush=True)
# sparse large n
for e in [x/4 for x in range(int(4*math.log10(hi)),4*13+1)]:
    n=int(10**e); print(nm,'n=%.3g margin %.3f'%(n,mf(n)),flush=True)
