import pickle, math
from decimal import Decimal, getcontext
getcontext().prec=80
F=pickle.load(open('family.pkl','rb'))
pi=math.pi
def U0pat(n): return math.cos(5*pi/16-pi*n/4)+math.cos(-13*pi/16-3*pi*n/4)
print('U0 predicted pattern:',''.join('+' if U0pat(r)>0 else '-' for r in range(8)),' min |trig| =',min(abs(U0pat(r)) for r in range(8)))
U1pat=lambda n: math.cos(3*pi/8-pi*n/2)/math.sqrt(2)
print('U1 predicted pattern:',''.join('+' if U1pat(r)>0 else '-' for r in range(4)))
PI=Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
def main(nm,n,d):
    N=Decimal(n)+Decimal(d)
    b=PI/(4*Decimal(2).sqrt()) if nm=='U0' else PI/(2*Decimal(2).sqrt())
    x=b*N.sqrt(); ch=(x.exp()+(-x).exp())/2
    tr=Decimal(U0pat(n%16)) if nm=='U0' else Decimal(U1pat(n%4))
    return tr*ch/N.sqrt()
for nm in ('U0','U1'):
    c=F[nm]
    best=None
    for num in range(-16,17):
        d=num/32
        errs=[abs(Decimal(c[n])-main(nm,n,d))/abs(main(nm,n,d)) for n in (19990,19991,19992,19993,19994,19995,19996,19997)]
        e=max(errs)
        if best is None or e<best[0]: best=(e,num)
    print(nm,'best shift delta = %d/32, max rel err at n~20000 = %.2e'%(best[1],best[0]))
