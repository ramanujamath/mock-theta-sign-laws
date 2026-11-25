"""Numerical validation of every inequality used in farey_S0.md."""
import cmath, math, random, pickle
from math import pi, sqrt, exp, log, gcd
from fractions import Fraction as Fr
from farey_lib import *
from farey_bound import NA, NB1, NB2, Lambda, C0, A_, log_main_lower
from farey_N0 import bound_parts, PROD
from radial import S as S_series
random.seed(7)

def theta_bounds(a, sg, h, k, z):
    hp, kp = cusp_data(h, k); A = Fr(a, 8); s = Fr(sg, 2)  # sg=1 <=> s=-1 <=> sigma=1/2
    X = (1/z).real; rho = exp(-pi*pi*X/(2*kp*kp))
    Ah = Ahat(A, s, hp, kp)
    p = float(8*(-Fr(1,8)+A/2-A*A/2))
    base = exp(-z.real*p)*sqrt(pi/(4*kp*abs(z)))*rho**(float(B2(Ah))/2+1/24)
    return base*exp(phiL(Ah, rho)+etaL(rho)), base*exp(phiU(Ah, rho)+etaU(rho))

def jtp_abs(x, Q):
    # |Theta(x;Q)| via the Jacobi triple product (x;Q)(Q/x;Q)(Q;Q), summed in logs
    L = 0.0; n = 0
    while True:
        Qn = Q**n
        L += math.log(abs(1-Qn*x)) + math.log(abs(1-Qn*Q/x)) + math.log(abs(1-Qn*Q))
        if abs(Qn) < 1e-19: break
        n += 1
    return math.exp(L)

def farey_arcs(K):
    F = sorted({Fr(h, k) for k in range(1, K+1) for h in range(0, k+1)})
    arcs = []
    for i, x in enumerate(F[:-1]):
        y = F[i+1]; med = Fr(x.numerator+y.numerator, x.denominator+y.denominator)
        arcs.append((x, med, y))
    return F, arcs

def random_point_on_arc(K, t0):
    """random theta in [0,1), returns (h,k,z) for its Farey arc of order K."""
    F = sorted({Fr(h, k) for k in range(1, K+1) for h in range(0, k+1)})
    th = random.random()
    # nearest Farey fraction in the mediant sense
    import bisect
    i = bisect.bisect_left(F, Fr(th).limit_denominator(10**9))
    cand = [F[max(i-1,0)], F[min(i,len(F)-1)]]
    x, y = (cand[0], cand[1]) if cand[0] <= th else (cand[0], cand[0])
    if x == y: hk = x
    else:
        med = (x.numerator+y.numerator)/(x.denominator+y.denominator)
        hk = x if th <= med else y
    h, k = hk.numerator % hk.denominator, hk.denominator
    if k == 1: h = 0
    phi = 2*pi*(th - float(hk))
    return h, k, complex(t0, -phi)

print('== V1: |Theta| lower/upper bounds at random points of Farey arcs (direct summation) ==')
viol = 0; cnt = 0; rows = []
for t0, K in ((0.02, 17), (0.01, 25)):
    for trial in range(250):
        h, k, z = random_point_on_arc(K, t0)
        q = e(h/k)*cmath.exp(-z)
        for a, sg in ((1,1),(2,0),(3,0),(0,1),(2,1),(6,1)):
            x = (-1 if sg else 1)*q**a
            d = abs(theta_direct(x, q**8))
            if d < 1e-10: d = jtp_abs(x, q**8)     # direct sum is roundoff-limited there
            lo, up = theta_bounds(a, sg, h, k, z)
            cnt += 1
            if not (lo*(1-1e-6) <= d <= up*(1+1e-6)): viol += 1
            rows.append((log(d/lo), log(up/d)))
print(f'  {cnt} checks, violations: {viol};  min log(actual/lower) = {min(r[0] for r in rows):.2e}, min log(upper/actual) = {min(r[1] for r in rows):.2e}')
print(f'  median slack: lower {sorted(r[0] for r in rows)[len(rows)//2]:.2f}, upper {sorted(r[1] for r in rows)[len(rows)//2]:.2f} (natural log)')

print('== V2: geometry: min Re(1/z)/k\'^2 over Farey arcs vs 1/(4pi) ==')
for n in (100, 1000, 10000, 60000):
    t0 = pi/(16*sqrt(n)); K = int(sqrt(2*pi/t0))
    F = sorted({Fr(h, k) for k in range(1, K+1) for h in range(0, k+1)})
    w = 1e9
    for i, x in enumerate(F):
        lm = Fr(F[i-1].numerator+x.numerator, F[i-1].denominator+x.denominator) if i > 0 else x
        rm = Fr(F[i+1].numerator+x.numerator, F[i+1].denominator+x.denominator) if i < len(F)-1 else x
        ph = 2*pi*float(max(x-lm, rm-x)); kk = x.denominator
        w = min(w, t0/(kk*kk*(t0*t0+ph*ph)))
    print(f'  n={n}: t0={t0:.5f} K={K}  min Re(1/z)/k^2 = {w:.5f}  >= 1/(4pi) = {1/(4*pi):.5f}: {w >= 1/(4*pi)}')

import sys
if "--v1" in sys.argv: sys.exit()
print("== V3: S0 on the circle vs assembled sup bound (minor arcs) ==")
for n in (400, 1000, 2000):
    t0 = pi/(16*sqrt(n)); K = int(sqrt(2*pi/t0))
    tot, lmi, lma, lt, cmin, arg = bound_parts(n)
    supbound = exp(lmi - n*t0)          # sup bound for |S0| on minor arcs
    worst = 0; worstpt = None; worstmaj = 0
    for trial in range(400):
        h, k, z = random_point_on_arc(K, t0)
        q = e(h/k)*cmath.exp(-z)
        val = abs(S_series(q, 0)[0])
        if k == 4: worstmaj = max(worstmaj, val); continue
        if val > worst: worst, worstpt = val, (h, k)
    # also near the k=8,16 cusps where the bound is tight in exponent
    for (h, k) in ((5,8),(3,8),(1,16),(3,16),(1,12)):
        q = e(h/k)*exp(-t0); val = abs(S_series(q, 0)[0])
        if val > worst: worst, worstpt = val, (h, k)
    print(f'  n={n}: t0={t0:.4f}: max|S0| sampled on minor arcs = {worst:.3e} at {worstpt};  bound = {supbound:.3e}  (log ratio {log(supbound/worst):.1f})')

print('== V4: Bessel identity (1/2pi i) int z^{-1/2} e^{Nz+a/z} dz = cosh(pi sqrt N/4)/sqrt(pi N) ==')
import numpy as np
for N in (50.0, 200.0):
    t0 = pi/(8*sqrt(N)); Ymax = 4000.0
    y = np.linspace(-Ymax, Ymax, 4_000_001); zz = t0+1j*y
    f = zz**-0.5*np.exp(N*zz + A_/zz)
    # damp the oscillatory tail with a smooth cutoff and integration by parts estimate negligible
    w = np.exp(-(y/(0.8*Ymax))**8)
    val = np.trapezoid(f*w, y)/(2*pi)
    print(f'  N={N}: numeric {val.real:.8e} (imag {val.imag:.1e}),  formula {math.cosh(pi*sqrt(N)/4)/sqrt(pi*N):.8e}')

print('== V5: major arc: TQ(i e^{-z}) vs C0 z^{-1/2} e^{-z/16} e^{pi^2/(64z)}, and |TQ-Main| bound ==')
C0c = C0*cmath.exp(9j*pi/16)
for t0, phi in ((0.05, 0.0), (0.03, 0.02), (0.02, -0.05), (0.01, 0.1)):
    z = complex(t0, -phi); q = 1j*cmath.exp(-z)
    tq = q*theta_direct(-q, q**8)*theta_direct(q**2, q**8)**2/theta_direct(q**3, q**8)**2
    main = C0c*z**-0.5*cmath.exp(-z/16 + pi*pi/(64*z))
    X = (1/z).real; rho = exp(-pi*pi*X/2)
    bnd = C0*abs(z)**-0.5*exp(-t0/16)*rho**(-1/32)*math.expm1(Lambda(rho))
    print(f'  z={z}: |TQ/main - 1| = {abs(tq/main-1):.3e};  |TQ-main| = {abs(tq-main):.3e} <= bound {bnd:.3e}')

print('== V6: Lerch numerator bounds ==')
for t0 in (0.05, 0.01, 0.003):
    mx = [0, 0, 0]
    for trial in range(60):
        q = cmath.exp(-t0 + 2j*pi*random.random())
        Q = q**8
        def num(x, zz):
            return sum((-1)**(r & 1)*Q**(r*(r-1)//2)*zz**r/(1-Q**(r-1)*x*zz) for r in range(-int(8/sqrt(t0)), int(8/sqrt(t0))))
        vals = [abs(num(-q**3, -1)), abs(num(-q**3, -q**2)), abs(num(-q**3, -q**6))]
        mx = [max(u, v) for u, v in zip(mx, vals)]
    print(f'  t0={t0}: sampled max |N_A|,|N_B1|,|N_B2| = {mx[0]:.3e},{mx[1]:.3e},{mx[2]:.3e};  bounds {NA(t0):.3e},{NB1(t0):.3e},{NB2(t0):.3e}')

print('== V7: actual |c(n)-M(n)| vs rigorous bound (log scale) ==')
from decimal import Decimal, getcontext
getcontext().prec = 150
PI = Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955')
def dcos(r):  # cos(pi*r), r Fraction
    a = PI*Decimal(r.numerator)/Decimal(r.denominator); s = Decimal(1); t = Decimal(1)
    for i in range(1, 120): t = -t*a*a/((2*i-1)*(2*i)); s += t
    return s
c = pickle.load(open('big.pkl', 'rb'))['S0']
for n in (100, 1000, 5000, 11231, 20000, 40000, 60000):
    Nd = Decimal(n) - Decimal(1)/16; x = PI/4*Nd.sqrt()
    M = dcos(Fr(9,16) - Fr(n % 4, 2))*(x.exp()+(-x).exp())/2/Nd.sqrt()
    E = abs(Decimal(c[n]) - M)
    tot = bound_parts(n)[0]
    print(f'  n={n}: log|E(n)| actual = {float(E.ln()):9.3f}   rigorous bound {tot:9.3f}   log(0.195|main|) {log_main_lower(n):9.3f}')
