# u_cusps.py -- cusp growth exponents for U0, U1 (Gordon-McIntosh 8th order)
#   python3 u_cusps.py search|exact|cusps|radial
# Units: exponents = coefficient of 1/t in log|f(zeta e^{-t})|, in units of pi^2/64.
import sys, math, cmath, pickle, itertools
from math import gcd, pi
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
F = pickle.load(open('family.pkl', 'rb'))
X = (-1, 1)          # x = -q, base Q = q^4 throughout: m(-q, q^4, z)

# ---------------------------------------------------------------- numeric helpers
def theta_n(x, q, R=70): return sum((-1)**n * q**(n*(n-1)//2) * x**n for n in range(-R, R+1))
def m_n(x, q, z, R=70):
    return sum((-1)**r * q**(r*(r-1)//2) * z**r / (1 - q**(r-1)*x*z) for r in range(-R, R+1)) / theta_n(z, q)

def search():
    for q0 in (0.31+0.22j, -0.2+0.35j, 0.1-0.4j):
        T = {nm: sum(c*q0**n for n, c in enumerate(F[nm][:600])) for nm in ('U0', 'U1')}
        Z = [(e, a) for e in (1, -1, 1j, -1j) for a in range(4) if not (e == -1 and a == 3)]
        M = {z: m_n(-q0, q0**4, z[0]*q0**z[1]) for z in Z}
        for nm in ('U0', 'U1'):
            hits = []
            for z1, z2 in itertools.combinations_with_replacement(Z, 2):
                for s2 in (1, -1):
                    if z1 == z2 and s2 == -1: continue
                    v = M[z1] + s2*M[z2]
                    for c in (1, -1, 2, -2, .5, -.5, 1j, -1j, (1+1j)/2, (1-1j)/2, (-1+1j)/2, (-1-1j)/2):
                        if abs(c*v - T[nm]) < 1e-9: hits.append((c, z1, s2, z2))
            print(q0, nm, 'c*[m(z1)+s*m(z2)], z=(eps,a) means eps*q^a:', hits)

# ---------------------------------------------------------------- exact power-series check
N = 3000
def lerch_num(xs, xc, M, zs, zc):
    """sum_r (-1)^r Q^{r(r-1)/2} z^r / (1 - Q^{r-1} x z),  Q=q^M, x=xs q^xc, z=zs q^zc (zs,xs = +-1)."""
    off = 4000; L = [0]*(N+1+off)
    for r in range(-200, 201):
        e0 = M*r*(r-1)//2 + zc*r; c0 = (-1)**abs(r) * zs**abs(r)
        e = M*(r-1) + xc + zc; eps = xs*zs
        if e > 0:
            k = 0
            while e0 + k*e <= N:
                if e0 + k*e + off >= 0: L[e0 + k*e + off] += c0 * eps**k
                k += 1
        elif e < 0:
            k = 1
            while e0 - k*e <= N:
                if e0 - k*e + off >= 0: L[e0 - k*e + off] += -c0 * eps**k
                k += 1
        else:
            assert eps != 1, 'exact pole: denominator identically 0'
            if 0 <= e0 + off and e0 <= N: L[e0 + off] += Fraction(c0, 1 - eps)   # constant denominator 1-eps (=2)
    assert all(v == 0 for v in L[:off])
    return L[off:off+N+1]
def theta_s(zs, zc, M):
    r = [0]*(N+1)
    for n in range(-200, 201):
        e = M*n*(n-1)//2 + zc*n
        if 0 <= e <= N: r[e] += (-1)**abs(n) * zs**abs(n)
    return r
def mul_sparse(a, b):
    r = [0]*(N+1); nzb = [(j, y) for j, y in enumerate(b) if y]
    for i, x in enumerate(a):
        if x:
            for j, y in nzb:
                if i+j > N: break
                r[i+j] += x*y
    return r
def exact():
    # f = c*m(-q,q^4,z)  <=>  f * Theta(z;q^4) == c * LerchNumerator  (Theta has nonzero constant term)
    for nm, c, (zs, zc) in (('U0', 2, (-1, 0)), ('U0', 2, (1, 3)), ('U1', -1, (-1, 2)), ('U1', -1, (1, 1))):
        lhs = mul_sparse(F[nm][:N+1], theta_s(zs, zc, 4))
        rhs = [c*v for v in lerch_num(-1, 1, 4, zs, zc)]
        ok = all(lhs[n] == rhs[n] for n in range(N))
        print(f'{nm} = {c}*m(-q,q^4,{zs:+d}q^{zc}) exact through q^{N-1}: {ok}')

# ---------------------------------------------------------------- cusp exponents
def nu_eff(s, a, M, h, k):          # same as lemma.py
    L = 2*k; j = np.arange(L)
    E = M*j*(j-1)//2 + a*j
    ph = ((-s)**(j % 2)) * np.exp(2j*pi*((h*E) % k)/k)
    G = np.fft.ifft(ph)*L
    b = a - M/2
    c = lambda nu: G[nu % L]*np.exp(2j*pi*nu*b/(L*M))
    tol = 1e-7*L
    if abs(c(0)) > tol: return 0
    for nu in range(1, 3*L):
        if abs(c(nu)+c(-nu)) > tol: return nu
    raise RuntimeError
def expo(s, a, M, h, k): return -32*nu_eff(s, a, M, h, k)**2/(k*k*M)

def cls(s, b):
    """Theta(s q^b; q^4): reduce to class (s, a), a in {0,1,2}, |Theta| equal up to poly/monomial factors."""
    a = b % 4
    if a == 3: a = 1          # Theta(Q/x)=Theta(x), s=+-1
    return (s, a)
CL = [(1, 1), (1, 2), (-1, 0), (-1, 1), (-1, 2)]   # (1,0) is identically zero
ZOK = [(-1, 0), (1, 1), (-1, 1), (1, 2), (-1, 2), (1, 3)]   # z=s q^a: Theta(z;Q)!=0, xz not in Q^Z

def forms(z1):
    """U = c*m(x,Q,z1).  For each z0: U = c*m(x,Q,z0) + c*[m(z1)-m(z0)],
       m(z1)-m(z0) = z0 J(Q)^3 Th(z1/z0)Th(x z0 z1)/(Th(z0)Th(z1)Th(x z0)Th(x z1)).
       Returns list (label, crude-m class, numerator classes, denominator classes, J power)."""
    out = []
    for z0 in ZOK:
        s1, b1 = z1; s0, b0 = z0
        lab = f'm(z={s0:+d}q^{b0})'
        num = [cls(s1*s0, b1-b0), cls(-s0*s1, 1+b0+b1)]
        den = [cls(s0, b0), cls(s1, b1), cls(-s0, 1+b0), cls(-s1, 1+b1)]
        if z0 == z1 or (1, 0) in num:      # quotient identically zero -> theta-free single-m form
            out.append((lab, cls(s0, b0), None, None))
        else:
            out.append((lab + '+thq', cls(s0, b0), num, den))
    return out
FORMS = {'U0': forms((-1, 0)), 'U1': forms((-1, 2))}

def cusp_table(K=400):
    res = {}; maxnu = {c: 0 for c in CL}; maxnuJ = 0
    for k in range(1, K+1):
        for h in range(k):
            if gcd(h, k) != 1: continue
            ex = {}
            for c in CL:
                nu = nu_eff(c[0], c[1], 4, h, k); maxnu[c] = max(maxnu[c], nu)
                ex[c] = -32*nu*nu/(k*k*4)
            nuJ = nu_eff(1, 4, 12, h, k); maxnuJ = max(maxnuJ, nuJ); eJ = -32*nuJ**2/(k*k*12)
            for nm, fl in FORMS.items():
                best = None
                for lab, cm, num, den in sorted(fl, key=lambda f: f[2] is not None):   # theta-free forms first (win ties)
                    e = -ex[cm]                      # crude m bound: 1/Theta(z;Q) (Lerch numerator polynomial)
                    if num is not None:
                        eq = 3*eJ + sum(ex[c] for c in num) - sum(ex[c] for c in den)
                        e = max(e, eq)
                    if best is None or e < best[0] - 1e-12: best = (e, lab)
                res[(nm, h, k)] = best
    return res, maxnu, maxnuJ

def cusps(K=400):
    res, maxnu, maxnuJ = cusp_table(K)
    print('max nu_eff over k<=%d for Theta(s q^a;q^4) classes (lemma bound 24):' % K, maxnu, ' J4=Theta(q^4;q^12) (bound 72):', maxnuJ)
    dom = {'U0': (8, 0.5), 'U1': (4, 2.0)}
    lines = []
    for nm in ('U0', 'U1'):
        kd, ed = dom[nm]
        dvals = sorted({(f'{h}/{k}', round(v[0], 6), v[1]) for (n, h, k), v in res.items() if n == nm and k == kd})
        others = sorted(((v[0], f'{h}/{k}', v[1]) for (n, h, k), v in res.items() if n == nm and k != kd), reverse=True)
        top = others[:8]
        # which form achieves the best at each non-dominant cusp
        used = {}
        for e, c, lab in others: used[lab] = used.get(lab, 0) + 1
        s = (f'{nm}: dominant cusps (k={kd}): {dvals}\n'
             f'   worst non-dominant cusps k<={K}: {[(round(e,6), c, lab) for e, c, lab in top]}\n'
             f'   strictly below {ed}: {others[0][0] < ed - 1e-9};  forms used (best): {used}\n'
             f'   k=1,2,4 values: {[(f"{h}/{k}", round(res[(nm,h,k)][0],6), res[(nm,h,k)][1]) for k in (1,2,4,8) for h in range(k) if gcd(h,k)==1]}')
        # theta-free single form alone:
        print(s); lines.append(s)
    # exponents using ONLY the single theta-free forms U0=2m(-q,q^4,q^3), U1=-m(-q,q^4,q)
    for nm, c, kd in (('U0', (1, 1), 8), ('U1', (1, 1), 4)):
        w = (-1, None)
        for k in range(1, K+1):
            if k == kd: continue
            for h in range(k):
                if gcd(h, k) == 1:
                    e = -expo(1, 3 if nm == 'U0' else 1, 4, h, k)
                    if e > w[0] + 1e-12: w = (e, f'{h}/{k}')
        s = f'{nm} single theta-free form z=q^{3 if nm=="U0" else 1} only: worst non-dominant exponent {w[0]:.6f} at {w[1]}'
        print(s); lines.append(s)
    for nm, ed in (('U0', .5), ('U1', 2.)):
        s = f'{nm} lemma bound k>{K}: any single-m form, crude bound = 1/Theta(z;q^4), 32*nu^2/(4k^2) <= 32*24^2/(4k^2) = 4608/k^2 < {4608/K**2:.4f} (dominant {ed})'
        print(s); lines.append(s)
    open('results/u_cusps_out.txt', 'w').write('\n'.join(lines) + '\n')

# ---------------------------------------------------------------- radial cross-check (high precision, defining series)
from decimal import Decimal as D, getcontext
getcontext().prec = 90
def dpi():
    getcontext().prec += 2
    three = D(3); lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s; n, na = n+na, na+8; d, da = d+da, da+32; t = (t*n)/d; s += t
    getcontext().prec -= 2
    return +s
PI = dpi()
def dcos(x):
    getcontext().prec += 2; i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s; i += 2; fact *= i*(i-1); num *= x*x; sign *= -1; s += num/fact*sign
    getcontext().prec -= 2; return +s
def dsin(x): return dcos(PI/2 - x)
class C:
    __slots__ = ('r', 'i')
    def __init__(s, r, i=D(0)): s.r = D(r); s.i = D(i)
    def __add__(s, o): o = o if isinstance(o, C) else C(o); return C(s.r+o.r, s.i+o.i)
    __radd__ = __add__
    def __mul__(s, o): o = o if isinstance(o, C) else C(o); return C(s.r*o.r-s.i*o.i, s.r*o.i+s.i*o.r)
    __rmul__ = __mul__
    def __truediv__(s, o):
        d = o.r*o.r+o.i*o.i; return C((s.r*o.r+s.i*o.i)/d, (s.i*o.r-s.r*o.i)/d)
    def ab(s): return (s.r*s.r+s.i*s.i).sqrt()
def qpt(h, k, t):
    ang = 2*PI*h/k; rad = (-D(t)).exp(); return C(rad*dcos(ang), rad*dsin(ang))
def powc(q, n):
    r = C(1); b = q
    while n:
        if n & 1: r = r*b
        b = b*b; n >>= 1
    return r
def U_series(nm, q, t):
    """U0 = sum q^{n^2}(-q;q^2)_n/(-q^4;q^4)_n ;  U1 = sum q^{(n+1)^2}(-q;q^2)_n/(-q^2;q^4)_{n+1}"""
    nmax = int(math.sqrt(230/t)) + 10
    q2 = q*q; q4 = q2*q2
    if nm == 'U0': cur = C(1); sh = 0
    else: cur = C(1)/(C(1)+q2); sh = 1
    tot = C(0); big = D(0)
    qodd = q            # q^{2n-1}
    qden = q4 if nm == 'U0' else q2*q4   # q^{4n} or q^{4n+2}
    for n in range(nmax):
        if n > 0:
            cur = cur*(C(1)+qodd)/(C(1)+qden); qodd = qodd*q2; qden = qden*q4
        term = cur*powc(q, (n+sh)**2); tot = tot+term; big = max(big, term.ab())
    return tot, big
def radial():
    ts = (0.01, 0.0075, 0.005, 0.0035, 0.0025)
    pred = {}
    try:
        for ln in open('results/u_cusps_out.txt'): pass
    except Exception: pass
    out = []
    for nm, cl in (('U0', ((1, 8), (3, 8), (1, 4), (1, 16), (3, 16), (1, 12), (1, 2), (0, 1), (1, 3), (1, 5))),
                   ('U1', ((1, 4), (1, 12), (5, 12), (1, 8), (3, 8), (1, 20), (1, 2), (0, 1), (1, 3), (1, 16)))):
        for h, k in cl:
            A = []; y = []; cancel = 0
            for t in ts:
                v, big = U_series(nm, qpt(h, k, t), t)
                A.append([1/t, math.log(t), 1]); y.append(float(v.ab().ln()))
                cancel = max(cancel, float((big/v.ab()).log10()) if v.ab() > 0 else 99)
            sol, *_ = np.linalg.lstsq(np.array(A), np.array(y), rcond=None)
            # also fit using only the two smallest t with beta fixed from full fit
            s = f'{nm} at {h}/{k}: gamma = {sol[0]*64/pi**2:+.4f} (pi^2/64 units), beta = {sol[1]:+.3f}, max digits cancelled {cancel:.1f}'
            print(s, flush=True); out.append(s)
    open('results/u_cusps_radial_out.txt', 'w').write('\n'.join(out) + '\n')

if __name__ == '__main__':
    for a in sys.argv[1:]:
        {'search': search, 'exact': exact, 'cusps': cusps, 'radial': radial}[a]()
