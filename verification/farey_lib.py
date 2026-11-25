"""Core objects for the explicit circle method for S0.

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

# ---- exact numerics (for validation) ----
def theta_direct(x, Q):
    lQ = cmath.log(Q); lx = cmath.log(x); tot = 0j
    for n in range(-6000, 6001):
        ex = n*(n-1)/2
        lm = ex*lQ.real + n*lx.real
        if lm < -700: continue
        tot += (-1)**(n & 1) * cmath.exp(ex*lQ + n*lx)
    return tot

def egcd(a, b):
    if b == 0: return (a, 1, 0)
    g, x, y = egcd(b, a % b); return (g, y, x - (a//b)*y)

def gamma_for(hp, kp):
    # need alpha*d - beta*c = 1 with c=kp, d=-hp :  -alpha*hp - beta*kp = 1
    g, x, y = egcd(hp, kp)      # x*hp + y*kp = g = 1
    alpha, beta = -x, -y
    assert alpha*(-hp) - beta*kp == 1
    return alpha, beta

def siegel_abs(A, B, Tp, nmax=3000):
    A = A - math.floor(A); B = B - math.floor(B)
    Q = e(Tp); qz = e(A*Tp + B)
    v = abs(Q)**(float(B2(Fr(A).limit_denominator(1000)))/2) * abs(1-qz)
    for n in range(1, nmax):
        a = Q**n
        v *= abs(1-a*qz)*abs(1-a/qz)
        if abs(a) < 1e-18 * max(1, abs(qz), 1/abs(qz)): break
    return v
def eta_abs(Tp):
    Q = e(Tp); v = abs(Q)**(1/24); n = 1
    while abs(Q)**n > 1e-18: v *= abs(1-Q**n); n += 1
    return v

def theta_via_transform(a, sigma, h, k, z):
    """|Theta(s q^a; q^8)| computed on the transformed side."""
    hp, kp = cusp_data(h, k); al, be = gamma_for(hp, kp)
    T = 8*(h/k + 1j*z/(2*pi))
    Tp = (al*T + be)/(kp*T - hp)
    A = Fr(a, 8); B = Fr(sigma)
    Ap = -A*hp - B*kp; Bp = B*al - A*be
    q_abs = math.exp(-z.real)
    p = 8*(-Fr(1,8) + A/2 - A*A/2)
    return q_abs**float(p) * abs(kp*T - hp)**-0.5 * eta_abs(Tp) * siegel_abs(float(Ap), float(Bp), Tp)
