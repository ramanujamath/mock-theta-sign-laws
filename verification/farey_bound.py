"""Assemble the explicit bound |E(n)| <= Bound(n) and find N0.
All inequalities are those of Lemmas 1-7 in farey_S0.md."""
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

def minor_constant(t0, K):
    """returns (log Cmin, details): sup over minor arcs of |S0| <= Cmin * exp(pi^2/(256 t0))."""
    D = sqrt(K*K*t0*t0 + (2*pi/(K+1))**2)          # max of k|z| on arcs of order K
    w0 = t0/D**2                                    # min of Re(1/z)/k'^2
    rho = exp(-pi*pi*w0/2)
    preTQ = sqrt(pi/(4*t0)); preAL = sqrt(4*D/pi)
    Eall = Fr(1, 4)
    worst = -1e300; det = {}
    for kpar, hm in CLASSES:
        # exponents: actual per-term exponents are R/k'^2 <= Eall (checked in farey_exponents.py); use Eall
        if kpar == 2:      # form B
            t1 = log(NB1(t0)*preAL) - log_lower_single(Fr(1,4), Fr(1,2), hm, kpar, rho)
            t2 = log(NB2(t0)*preAL) - log_lower_single(Fr(3,4), Fr(1,2), hm, kpar, rho)
            val = log(exp(t1) + exp(t2))
        else:              # form A
            t1 = log(preTQ) + log_upper_TQ(hm, kpar, rho)
            t2 = log(2*NA(t0)*preAL) - log_lower_single(Fr(0), Fr(1,2), hm, kpar, rho)
            val = log(exp(t1) + exp(t2))
        det[(kpar, hm)] = val
        worst = max(worst, val)
    return worst, rho, det

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

def major_constant(t0, K):
    phimax = 2*pi/(4*(K+1))
    X = t0/(t0*t0 + phimax**2)
    rho = exp(-pi*pi*X/2)
    zmax = sqrt(t0*t0 + phimax**2)
    L = Lambda(rho)
    b_tq = C0*t0**-0.5 * rho**(-1/32) * math.expm1(L)
    # AL_A at cusp 1/4: Ahat of (0,1/2) is 1/2, exponent 0
    b_al = 2*NA(t0)*sqrt(4*zmax/pi)*exp(-(phiL(0.5, rho)+etaL(rho)))
    return b_tq + b_al, rho

def tail(n, t0, K):
    N = n - 1/16
    Y = pi/(2*(K+4))
    return 4*C0*exp(N*t0)/(2*pi*N)*exp(A_*t0/Y**2)*(2*Y**-0.5 + (2*A_/3)*Y**-1.5)

def log_bound(n, lam, Kfac=1.0):
    t0 = lam*pi/(16*sqrt(n))
    K = max(4, int(Kfac*sqrt(2*pi/t0)))
    lc, rho, _ = minor_constant(t0, K)
    bmaj, rmaj = major_constant(t0, K)
    lminor = n*t0 + lc + pi*pi/(256*t0)
    lmaj = n*t0 + log(bmaj)
    lt = log(tail(n, t0, K))
    m = max(lminor, lmaj, lt)
    return m + log(exp(lminor-m) + exp(lmaj-m) + exp(lt-m)), dict(t0=t0, K=K, rho=rho, rmaj=rmaj, lc=lc, lminor=lminor, lmaj=lmaj, lt=lt)

def log_main_lower(n):
    N = n - 1/16; x = pi*sqrt(N)/4
    return log(math.cos(7*pi/16)) + x + math.log1p(exp(-2*x)) - log(2) - 0.5*log(N)

def best(n):
    bestv = None
    for lam in [0.6+0.05*i for i in range(30)]:
        for Kf in (0.7, 0.85, 1.0, 1.2, 1.4):
            v, d = log_bound(n, lam, Kf)
            if bestv is None or v < bestv[0]: bestv = (v, lam, Kf, d)
    return bestv

if __name__ == '__main__':
    import sys
    for n in (500, 1000, 2000, 3000, 5000, 10000, 20000, 60000):
        v, lam, Kf, d = best(n)
        print(f"n={n:6d} log|E|<= {v:9.3f}  log(0.195 main)= {log_main_lower(n):9.3f}  margin={log_main_lower(n)-v:8.3f}"
              f"  lam={lam:.2f} Kf={Kf} K={d['K']} rho_minor={d['rho']:.3f} logCmin={d['lc']:.2f} parts(minor,major,tail)=({d['lminor']:.1f},{d['lmaj']:.1f},{d['lt']:.1f})")
