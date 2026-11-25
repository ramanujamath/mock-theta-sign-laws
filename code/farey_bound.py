"""Helper functions for the S0 bound (Lemmas 3.1-3.3, 5.2 and Proposition 4.1 of the paper)."""
import math
from math import pi, sqrt, exp, log
from fractions import Fraction as Fr
from farey_lib import *

A_ = pi*pi/64          # a = pi^2/64
C0 = sqrt(pi)/2        # |C0|

# ---- class table: (k' parity, h' mod 8) -> (form, list of terms) ----
CLASSES = [(1,h) for h in range(8)] + [(2,h) for h in (1,3,5,7)]

def NA(t0): return exp(t0)*(1+sqrt(pi/(4*t0)))/(1-exp(-3*t0))
def NB1(t0): return exp(t0/4)*(1+sqrt(pi/(4*t0)))/(1-exp(-3*t0))
def NB2(t0): return exp(t0/4)*(1+sqrt(pi/(4*t0)))/(1-exp(-t0))

def log_lower_single(A, s, hp, kp, rho):   # log lower bound of normalised |g_{(A,s)}| * |eta|
    return phiL(Ahat(A, s, hp, kp), rho) + etaL(rho)
def log_upper_TQ(hp, kp, rho):
    v = etaU(rho)
    for m, A, s in TQ:
        Ah = Ahat(A, s, hp, kp)
        v += m*phiU(Ah, rho) if m > 0 else (-m)*(-phiL(Ah, rho))
    return v

# main arc: factors of F (exponents of rho, multiplicity) for TQ at cusp 1/4: Ahat = 1/4 (g1), 1/2 (g2^2), 1/4 (g3^-2), eta
def Lambda(rho):
    exps = []
    for Ah, m in ((Fr(1,4), 1), (Fr(1,2), 2), (Fr(1,4), 2)):
        exps.append((float(Ah), m))
        for n in range(1, 3000):
            exps.append((n+float(Ah), m)); exps.append((n-float(Ah), m))
            if rho**(n-float(Ah)) < 1e-18: break
    n = 1
    while rho**n > 1e-18: exps.append((n, 1)); n += 1
    return sum(-m*math.log1p(-rho**x) for x, m in exps)

def log_main_lower(n):
    N = n - 1/16; x = pi*sqrt(N)/4
    return log(math.cos(7*pi/16)) + x + math.log1p(exp(-2*x)) - log(2) - 0.5*log(N)

