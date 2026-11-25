"""Core objects for the explicit circle method (Sections 3 and 5 of the paper).

Notation.  q = e(h/k) e^{-z},  Re z = t0 > 0,  tau = h/k + i z/(2 pi),  T = 8 tau.
Every theta factor is Theta(s q^a; q^8) with s = e(sigma), sigma in {0,1/2}; A = a/8.

Identity (Jacobi triple product, exact):
   |Theta(s q^a; q^8)| = |q|^{8(-1/8 + A/2 - A^2/2)} |g_{(A,sigma)}(T)| |eta(T)|
with the Siegel function g_{(A,B)}(T) = -Q^{B2(A)/2} e(B(A-1)/2) (1-qz) prod_{n>=1}(1-Q^n qz)(1-Q^n/qz),
Q = e(T), qz = e(A T + B), B2(x) = x^2 - x + 1/6.
Transformation: cusp 8h/k = h'/k' (reduced), gamma = (alpha beta; k' -h') in SL2(Z), T' = gamma T:
   |eta(T)| = |k'T-h'|^{-1/2} |eta(T')|,  |k'T-h'| = 4k'|z|/pi,
   |g_{(A,B)}(T)| = |g_{(A,B) gamma^{-1}}(T')|,   (A,B)gamma^{-1} = (-A h' - B k', B alpha - A beta),
   |g| depends only on the characteristic mod Z^2,
   rho := |e(T')| = exp(-pi^2 Re(1/z) / (2 k'^2)).
"""
import cmath, math
from math import pi, gcd
from fractions import Fraction as Fr

def B2(x): return x*x - x + Fr(1, 6)

def e(x): return cmath.exp(2j*pi*x)

# factor lists: (multiplicity, A, sigma)  [A = a/8]
TQ  = [(1, Fr(1,8), Fr(1,2)), (2, Fr(1,4), Fr(0)), (-2, Fr(3,8), Fr(0))]   # q*Theta(-q)Theta(q^2)^2/Theta(q^3)^2
ALA = [(-1, Fr(0), Fr(1,2))]          # 1/Theta(-1;q^8)
ALB1 = [(-1, Fr(1,4), Fr(1,2))]       # 1/Theta(-q^2;q^8)
ALB2 = [(-1, Fr(3,4), Fr(1,2))]       # 1/Theta(-q^6;q^8)

def qpow(fl):   # exponent p with |prod| = |q|^p * ... (TQ includes the extra q)
    return sum(m*8*(-Fr(1,8)+A/2-A*A/2) for m,A,s in fl)

def cusp_data(h, k):
    g = gcd(8, k); kp = k//g; hp = 8*h//g
    assert gcd(hp, kp) == 1
    return hp, kp

def Ahat(A, B, hp, kp):
    x = -A*hp - B*kp
    return x - math.floor(x)

def class_exponent(fl, hp, kp, extra_eta):
    """exponent E (units pi^2 Re(1/z)/64) of the rho-power in the bound for prod of factors in fl
       (eta power = sum of multiplicities, unless extra_eta given)."""
    etap = sum(m for m,_,_ in fl) if extra_eta is None else extra_eta
    s = Fr(etap, 24) + sum(m*B2(Ahat(A, s_, hp, kp))/2 for m, A, s_ in fl)
    # rho^s = exp(-pi^2 X s/(2k'^2)) = exp(pi^2 X/64 * (-32 s/k'^2))
    return -32*s/kp**2

# ---- product factors (normalised, rho-power removed) ----
def phiU(Ah, rho, nmax=4000):
    v = math.log1p(rho**Ah) if Ah > 0 else math.log(2.0)
    for n in range(1, nmax):
        a = rho**(n+Ah); b = rho**(n-Ah)
        v += math.log1p(a) + math.log1p(b)
        if b < 1e-18: break
    return v            # log of upper bound
def phiL(Ah, rho, nmax=4000):
    v = math.log1p(-rho**Ah) if Ah > 0 else math.log(2*math.sin(pi/8))
    for n in range(1, nmax):
        a = rho**(n+Ah); b = rho**(n-Ah)
        v += math.log1p(-a) + math.log1p(-b)
        if b < 1e-18: break
    return v            # log of lower bound
def etaU(rho):
    v = 0.0; n = 1
    while rho**n > 1e-18: v += math.log1p(rho**n); n += 1
    return v
def etaL(rho):
    v = 0.0; n = 1
    while rho**n > 1e-18: v += math.log1p(-rho**n); n += 1
    return v

def log_prod_bound(fl, hp, kp, rho):
    """log of the upper bound of the normalised product part of prod |g|^m |eta|^{sum m} (magnitude)."""
    etap = sum(m for m,_,_ in fl)
    v = etap*(etaU(rho) if etap > 0 else -etaL(rho)) if etap != 0 else 0.0
    if etap < 0: v = -etap*(-etaL(rho))
    for m, A, s in fl:
        Ah = Ahat(A, s, hp, kp)
        v += m*phiU(Ah, rho) if m > 0 else (-m)*(-phiL(Ah, rho))
    return v
