import pickle, math
from fractions import Fraction
D=pickle.load(open('big.pkl','rb'))
pi=math.pi
# main(n) = A cos(theta - pi n/2) cosh(b sqrt(N))/sqrt(N),  N = n + delta  (delta from the theta-function shift)
spec={'S0':(1,9*pi/16,pi/4,-1/16,pi/8),'S1':(1,-3*pi/16,pi/4,7/16,pi/8),'mu2':(math.sqrt(2),-5*pi/8,pi/(2*math.sqrt(2)),-1/8,pi/(4*math.sqrt(2)))}
SC=10**40
for k,(A,th,b,d,berr) in spec.items():
    c=D[k]; print('==',k,' predicted error growth exponent',round(berr,4))
    for n in (2000,5000,10000,20000,40000,59999):
        N=n+d
        # compute main in high precision via Python ints: use exp in float of log then subtract via Fraction approximations
        x=b*math.sqrt(N)
        # use decimal for accuracy
        from decimal import Decimal, getcontext
        getcontext().prec=120
        X=Decimal(x)  # float-precision argument: ~1e-16 relative -> limits error; refine below
        # better: compute x with Decimal sqrt
        Nd=Decimal(n)+Decimal(Fraction(d).numerator)/Decimal(Fraction(d).denominator)
        bd={'S0':Decimal(pi)/4,'S1':Decimal(pi)/4,'mu2':Decimal(pi)/(2*Decimal(2).sqrt())}[k]
        PI=Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651')
        bd={'S0':PI/4,'S1':PI/4,'mu2':PI/(2*Decimal(2).sqrt())}[k]
        Xd=bd*Nd.sqrt()
        ch=(Xd.exp()+(-Xd).exp())/2
        cs=Decimal(math.cos(th-pi*(n%4)/2))  # cos of a multiple of pi/16: float accurate to 1e-16 only
        # exact cos via high precision series
        ang=Decimal(th-pi*(n%4)/2)
        def dcos(a):
            # reduce using PI precisely: angles are rational multiples of pi
            return None
        r=Fraction(th/pi).limit_denominator(64)-Fraction(n%4,2)
        a=PI*Decimal(r.numerator)/Decimal(r.denominator)
        s=Decimal(1);term=Decimal(1)
        for i in range(1,80):
            term=-term*a*a/((2*i-1)*(2*i)); s+=term
        Ad=Decimal(2).sqrt() if k=='mu2' else Decimal(1)
        main=Ad*s*ch/Nd.sqrt()
        E=Decimal(c[n])-main
        rel=abs(E)/abs(main)
        print(f'  n={n}: log|E|/sqrt(n)={float(abs(E).ln())/math.sqrt(n):.4f}   |E|/|main|={float(rel):.2e}')
