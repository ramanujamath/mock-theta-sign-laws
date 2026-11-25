"""Explicit circle-method bound and N0 for S1 (Sections 3-6 of the paper).
S1 = -(2/q) m(-q,q^8,-1) + q^{-1} Th(-q^3;q^8) Th(q^2;q^8)^2 / Th(q;q^8)^2          (form A)
   = -q^{-1} [ m(-q,q^8,-q^2) + m(-q,q^8,-q^6) ]                                     (form B)
Main term M(n) = cos(-3pi/16 - pi n/2) cosh(pi sqrt(N)/4)/sqrt(N), N = n + 7/16, |cos| >= cos(5pi/16)."""
import math
from math import pi, sqrt, exp, log
from fractions import Fraction as Fr
from farey_lib import Ahat, phiU, phiL, etaU, etaL
from farey_bound import log_lower_single, A_
C0 = sqrt(pi)/2
TQ1 = [(1, Fr(3,8), Fr(1,2)), (2, Fr(1,4), Fr(0)), (-2, Fr(1,8), Fr(0))]
CLASSES = [(1,h) for h in range(8)] + [(2,h) for h in (1,3,5,7)]
RHO = exp(-pi/8)
def log_upper_TQ1(hp, kp, rho):
    v = etaU(rho)
    for m, A, s in TQ1:
        Ah = Ahat(A, s, hp, kp)
        v += m*phiU(Ah, rho) if m > 0 else (-m)*(-phiL(Ah, rho))
    return v
PROD = {}
for kpar, hm in CLASSES:
    if kpar == 2:
        PROD[(kpar,hm)] = ('B', exp(-log_lower_single(Fr(1,4),Fr(1,2),hm,kpar,RHO)), exp(-log_lower_single(Fr(3,4),Fr(1,2),hm,kpar,RHO)))
    else:
        PROD[(kpar,hm)] = ('A', exp(log_upper_TQ1(hm,kpar,RHO)), exp(-log_lower_single(Fr(0),Fr(1,2),hm,kpar,RHO)))
# Lerch numerators for x=-q:  sum |q|^{4r(r-1)}|z|^r / min|1-q^{8r-8} x z|
G = lambda t0: 1+sqrt(pi/(4*t0))
def NA(t0):  return exp(t0)*G(t0)/(1-exp(-t0))       # z=-1 : xz=q,   denominators 1-q^{8r-7}, min |e|=1
def NB1(t0): return exp(t0/4)*G(t0)/(1-exp(-3*t0))   # z=-q^2: xz=q^3, 1-q^{8r-5}, min |e|=3
def NB2(t0): return exp(t0/4)*G(t0)/(1-exp(-t0))     # z=-q^6: xz=q^7, 1-q^{8r-1}, min |e|=1
QINV = lambda t0: exp(t0)          # |q^{-1}|
def Lambda1(rho):                  # F at cusp 1/4: Ahat = 3/4 (m=1), 1/2 (m=2), 3/4 (m=-2), eta (m=1)
    tot = 0.0
    for Ah, m in ((0.75,1),(0.5,2),(0.75,2)):
        tot += -m*math.log1p(-rho**Ah)
        n = 1
        while rho**(n-Ah) > 1e-18:
            tot += -m*(math.log1p(-rho**(n+Ah)) + math.log1p(-rho**(n-Ah))); n += 1
    n = 1
    while rho**n > 1e-18: tot += -math.log1p(-rho**n); n += 1
    return tot
def bound_parts(n):
    t0 = pi/(16*sqrt(n)); N = n + 7/16
    D = sqrt(4*pi*t0); preTQ = sqrt(pi/(4*t0)); preAL = sqrt(4*D/pi)
    cmin = 0; arg = None
    for key, (f, u, v) in PROD.items():
        if f == 'B': val = QINV(t0)*preAL*(NB1(t0)*u + NB2(t0)*v)
        else:        val = exp(7*t0/16)*preTQ*u + 2*QINV(t0)*NA(t0)*preAL*v
        if val > cmin: cmin, arg = val, key
    lminor = n*t0 + pi*pi/(256*t0) + log(cmin)
    rmaj = exp(-pi*pi/(2*(t0+pi/8))); zmax = sqrt(t0*t0 + pi*t0/8)
    bmaj = C0*t0**-0.5*exp(7*t0/16)*rmaj**(-1/32)*math.expm1(Lambda1(rmaj)) \
         + 2*QINV(t0)*NA(t0)*sqrt(4*zmax/pi)*exp(-(phiL(0.5, rmaj)+etaL(rmaj)))
    lmaj = n*t0 + log(bmaj)
    Y = pi/(2*(sqrt(2*pi/t0)+4)); at = (sqrt(2*pi)+4*sqrt(t0))**2/16
    lt = log(4*C0/(2*pi*N)) + N*t0 + at + log(2*Y**-0.5 + (2*A_/3)*Y**-1.5)
    m = max(lminor, lmaj, lt)
    return m + log(exp(lminor-m)+exp(lmaj-m)+exp(lt-m)), lminor, lmaj, lt, cmin, arg
def log_main_lower(n):
    N = n + 7/16; x = pi*sqrt(N)/4
    return log(math.cos(5*pi/16)) + x + math.log1p(exp(-2*x)) - log(2) - 0.5*log(N)
