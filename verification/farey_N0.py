"""Uniform explicit bound (all n) and N0.  t0 = pi/(16 sqrt n), K = floor(sqrt(2 pi/t0)).
Uniform facts used (proved in farey_S0.md, Lemma 5):
   minor arcs:  Re(1/z)/k'^2 > 1/(4 pi)          => rho <= RHO_MIN := exp(-pi/8)
   k|z| <= D < sqrt(4 pi t0)
   major arc :  Re(1/z) > 1/(t0 + pi/8)          => rho <= exp(-pi^2/(2(t0+pi/8)))
   tail      :  Y >= pi/(2(sqrt(2pi/t0)+4)),  a t0/Y^2 <= (sqrt(2pi)+4 sqrt(t0))^2/16
"""
import math, sys
from math import pi, sqrt, exp, log
from fractions import Fraction as Fr
from farey_lib import *
from farey_bound import (CLASSES, NA, NB1, NB2, log_lower_single, log_upper_TQ, Lambda, C0, A_, log_main_lower)

RHO = exp(-pi/8)
# class constants (product parts only), computed once at rho = RHO (they are monotone in rho)
PROD = {}
for kpar, hm in CLASSES:
    if kpar == 2:
        PROD[(kpar,hm)] = ('B', exp(-log_lower_single(Fr(1,4),Fr(1,2),hm,kpar,RHO)),
                                exp(-log_lower_single(Fr(3,4),Fr(1,2),hm,kpar,RHO)))
    else:
        PROD[(kpar,hm)] = ('A', exp(log_upper_TQ(hm,kpar,RHO)),
                                exp(-log_lower_single(Fr(0),Fr(1,2),hm,kpar,RHO)))

def bound_parts(n):
    t0 = pi/(16*sqrt(n)); N = n - 1/16
    D = sqrt(4*pi*t0)
    preTQ = sqrt(pi/(4*t0)); preAL = sqrt(4*D/pi)
    cmin = 0; arg = None
    for key, (f, u, v) in PROD.items():
        val = (NB1(t0)*preAL*u + NB2(t0)*preAL*v) if f == 'B' else (preTQ*u + 2*NA(t0)*preAL*v)
        if val > cmin: cmin, arg = val, key
    lminor = n*t0 + pi*pi/(256*t0) + log(cmin)            # = pi sqrt(n)/8 + log cmin
    # major arc
    rmaj = exp(-pi*pi/(2*(t0+pi/8)))
    zmax = sqrt(t0*t0 + pi*t0/8)
    bmaj = C0*t0**-0.5*rmaj**(-1/32)*math.expm1(Lambda(rmaj)) \
         + 2*NA(t0)*sqrt(4*zmax/pi)*exp(-(phiL(0.5, rmaj)+etaL(rmaj)))
    lmaj = n*t0 + log(bmaj)
    # tails
    Y = pi/(2*(sqrt(2*pi/t0)+4)); at = (sqrt(2*pi)+4*sqrt(t0))**2/16
    lt = log(4*C0/(2*pi*N)) + N*t0 + at + log(2*Y**-0.5 + (2*A_/3)*Y**-1.5)
    m = max(lminor, lmaj, lt)
    tot = m + log(exp(lminor-m)+exp(lmaj-m)+exp(lt-m))
    return tot, lminor, lmaj, lt, cmin, arg

if __name__ == '__main__':
    print('rho_minor <=', RHO)
    for key,(f,u,v) in PROD.items(): print('  class',key,'form',f,' product constants %.3e %.3e'%(u,v))
    last_bad = None; first_good=None
    margins = {}
    for n in range(1, 200001):
        tot, *_ = bound_parts(n)
        mg = log_main_lower(n) - tot
        margins[n] = mg
        if mg <= 0: last_bad = n
    N0 = last_bad + 1
    print('last n <= 2e5 with bound >= 0.195*main :', last_bad, ' => N0 =', N0)
    # monotonicity of the margin beyond N0
    dec = [n for n in range(N0, 200000) if margins[n+1] < margins[n] - 1e-9]
    print('n in [N0,2e5) where margin decreases:', dec[:10], 'count', len(dec))
    for n in (N0-1, N0, 20000, 60000, 200000):
        tot, lmi, lma, lt, cmin, arg = bound_parts(n)
        print(f'n={n}: log|E|<={tot:.3f} (minor {lmi:.2f}, major {lma:.2f}, tail {lt:.2f}); log(0.195|main|)>={log_main_lower(n):.3f};'
              f' margin {log_main_lower(n)-tot:.3f}; worst class {arg}, Cmin={cmin:.3e}')
    # large-n check (up to 1e12 sparsely)
    for n in [10**j for j in range(6, 13)]:
        print('n=1e%d margin %.1f' % (round(math.log10(n)), log_main_lower(n)-bound_parts(n)[0]))
