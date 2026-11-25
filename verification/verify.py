import pickle, math
from fractions import Fraction
D=pickle.load(open('big.pkl','rb'))
pi=math.pi
spec={'S0':(1.0,9*pi/16,pi/4),'S1':(1.0,-3*pi/16,pi/4),'mu2':(math.sqrt(2),-5*pi/8,pi/(2*math.sqrt(2)))}
def logabs(x): # big int safe
    return math.log(abs(x)) if abs(x)<1e300 else (len(str(abs(x)))-17)*math.log(10)+math.log(int(str(abs(x))[:17]))
for k,(A,th,b) in spec.items():
    c=D[k]
    print('==',k,' predicted sign pattern:',''.join('+' if math.cos(th-pi*r/2)>0 else '-' for r in range(4)))
    for n in (100,101,102,103,1000,1001,1002,1003,10000,10001,10002,10003,59996,59997,59998,59999):
        cs=math.cos(th-pi*n/2)
        # log of main term: log(A|cos|) + b*sqrt(n) - log 2 - 0.5 log n  (cosh ~ e^x/2)
        lm=math.log(A*abs(cs))+b*math.sqrt(n)+math.log1p(math.exp(-2*b*math.sqrt(n)))-math.log(2)-0.5*math.log(n)
        ratio=math.exp(logabs(c[n])-lm)*(1 if (c[n]>0)==(cs>0) else -1)
        print(f'  n={n}: c(n)/main = {ratio:.6f}')
