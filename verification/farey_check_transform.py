# Validate the exact transformation identity for |Theta(s q^a;q^8)| at many cusps/points.
import cmath, math, random
from math import pi
from farey_lib import *
random.seed(1)
facs=[(1,Fr(1,2)),(2,Fr(0)),(3,Fr(0)),(0,Fr(1,2)),(2,Fr(1,2)),(6,Fr(1,2))]
worst=0
for (h,k) in [(0,1),(1,2),(1,3),(1,4),(3,4),(1,5),(1,8),(3,8),(1,16),(5,16),(7,12),(2,7),(5,24),(11,40),(1,10),(3,32)]:
    for trial in range(3):
        t0=random.uniform(0.02,0.06); phi=random.uniform(-1,1)*2*pi/(k*12)
        z=complex(t0,-phi)
        q=e(h/k)*cmath.exp(-z)
        for a,sg in facs:
            x=(-1 if sg else 1)*q**a
            d=abs(theta_direct(x,q**8))
            tr=theta_via_transform(a,sg,h,k,z)
            worst=max(worst,abs(math.log(d/tr)))
print('max |log(direct/transformed)| =',worst)
