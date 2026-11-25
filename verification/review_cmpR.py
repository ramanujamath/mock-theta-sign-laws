from review_exponents import theta_exp, J
import farey_gen as G, farey_lib as FL
from fractions import Fraction as Fr
from math import gcd
bad=0; cnt=0
for k in range(1,121):
    for h in range(k):
        if gcd(h,k)!=1: continue
        # farey_gen factors
        for f,mine in ((('th',4,Fr(0),Fr(1,2),-1), -theta_exp(-1,0,4,h,k)),
                       (('th',4,Fr(3,4),Fr(0),-1), -theta_exp(1,3,4,h,k)),
                       (('th',4,Fr(1,2),Fr(1,2),-1), -theta_exp(-1,2,4,h,k)),
                       (('th',4,Fr(1,4),Fr(0),-1), -theta_exp(1,1,4,h,k)),
                       (('eta',1,1),J(1,h,k)),(('eta',2,1),J(2,h,k)),(('eta',4,1),J(4,h,k))):
            R=G.factor_data(f,h,k)[0]; cnt+=1
            if Fr(R,k*k)!=mine: bad+=1; print('gen mismatch',f,h,k,R/k/k,mine)
        hp,kp=FL.cusp_data(h,k)
        for fl,mine in ((FL.TQ, theta_exp(-1,1,8,h,k)+2*theta_exp(1,2,8,h,k)-2*theta_exp(1,3,8,h,k)),
                        (FL.ALA,-theta_exp(-1,0,8,h,k)),(FL.ALB1,-theta_exp(-1,2,8,h,k)),(FL.ALB2,-theta_exp(-1,6,8,h,k)),
                        ([(1,Fr(3,8),Fr(1,2)),(2,Fr(1,4),Fr(0)),(-2,Fr(1,8),Fr(0))], theta_exp(-1,3,8,h,k)+2*theta_exp(1,2,8,h,k)-2*theta_exp(1,1,8,h,k))):
            R=FL.class_exponent(fl,hp,kp,None); cnt+=1
            # class_exponent returns -32 s / k'^2 in units pi^2 Re(1/z)/64
            if R!=mine: bad+=1; print('lib mismatch',fl,h,k,R,mine)
print('compared',cnt,'mismatches',bad)
