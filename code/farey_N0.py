"""Uniform explicit bound (all n) and N0.  t0 = pi/(16 sqrt n), K = floor(sqrt(2 pi/t0)).
Uniform facts used (Lemma 5.1 of the paper):
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
