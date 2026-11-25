"""Explicit circle-method bounds for mu2, U1, U0 (base q^4 Appell-Lerch sums + eta quotient TH).
Identities (verified exactly through q^20000, th_check.py):
    mu2 = 2 U0 - TH,   U1 = (TH - U0)/2,   U0 = TH - 2 U1,   TH = J2^8/(J1^3 J4^4).
Single-m forms (verified exactly through q^2999 by the U0/U1 cusp agent, u_cusps.py):
    U0 = 2 m(-q,q^4,-1) = 2 m(-q,q^4,q^3),   U1 = -m(-q,q^4,-q^2) = -m(-q,q^4,q).
Lemmas are those of farey_S0.md, generalised from q^8 to q^4 / eta(m tau):
  * |Theta(s q^a;q^M)| = |q|^{M(-1/8+A/2-A^2/2)} |g_{(A,sigma)}(M tau)| |eta(M tau)|, A=a/M.
  * At cusp h/k: g=gcd(M,k), k'=k/g, h'=Mh/g; |eta(M tau)| = (2 pi g/(M k|z|))^{1/2} |eta(T')|,
    rho = |e(T')| = exp(-4 pi^2 g^2 Re(1/z)/(M k^2)); Ahat = {-A h' - sigma k'}.
  * Lemma 2 product bounds (phiU/phiL/etaU/etaL from farey_lib).
  * Farey arcs of order K = floor(sqrt(2 pi/t0)) on |q| = e^{-t0}: Re(1/z)/k^2 > 1/(4 pi), t0 <= k|z| < sqrt(4 pi t0),
    and on the arc of h/k: Re(1/z) >= 1/(t0 + 2 pi/k^2), arc half-width >= 2 pi/(k(sqrt(2pi/t0)+k)).
"""
import math
from math import pi, sqrt, exp, log, gcd
from fractions import Fraction as Fr
from farey_lib import B2, phiU, phiL, etaU, etaL

def ahat(A, s, hp, kp):
    x = -A*hp - s*kp
    return x - math.floor(x)

# ---------------- objects ----------------
# factor: ('eta', m, p)  for J_m^p ;  ('th', M, A, sigma, p) for Theta(e(sigma) q^{A M}; q^M)^p
TH = (('eta', 2, 8), ('eta', 1, -3), ('eta', 4, -4))
def G(t0): return 1 + sqrt(pi/(2*t0))          # sum_r exp(-2 t0 (r-delta)^2) <= G
# Appell-Lerch single-m forms: (coefficient, theta-denominator factor, numerator bound function)
U0_FORMS = {
  'U0[z=-1]':  (2.0, ('th', 4, Fr(0), Fr(1,2), -1),   lambda t0: exp(t0/2)*G(t0)/(1-exp(-t0))),
  'U0[z=q^3]': (2.0, ('th', 4, Fr(3,4), Fr(0), -1),   lambda t0: exp(t0/8)*G(t0)/(1-exp(-4*t0))),
}
U1_FORMS = {
  'U1[z=-q^2]': (1.0, ('th', 4, Fr(1,2), Fr(1,2), -1), lambda t0: G(t0)/(1-exp(-t0))),
  'U1[z=q]':    (1.0, ('th', 4, Fr(1,4), Fr(0), -1),   lambda t0: exp(t0/8)*G(t0)/(1-exp(-2*t0))),
}

def factor_data(f, h, k):
    """returns (E*k^2 contribution [units pi^2 Re(1/z)/64], log prefactor const, sum p, qpow P, product-bound function of w)"""
    if f[0] == 'eta':
        _, m, p = f; g = gcd(m, k)
        c = Fr(1, 24); M = m; Ah = None
        qp = Fr(-m*p, 24)
    else:
        _, M, A, s, p = f; g = gcd(M, k); kp = k//g; hp = M*h//g
        Ah = ahat(A, s, hp, kp); c = B2(Ah)/2 + Fr(1, 24)
        qp = p*M*(-Fr(1,8) + A/2 - A*A/2)
    R = -Fr(256)*p*g*g*c/M                       # E = R/k^2
    pref = (p/2)*log(2*pi*g/M)                   # times (k|z|)^{-p/2}
    def prod(w):                                 # w = lower bound of Re(1/z)/k^2
        rho = exp(-4*pi*pi*g*g*w/M)
        if f[0] == 'eta':
            return p*etaU(rho) if p > 0 else (-p)*(-etaL(rho))
        return p*(phiU(Ah, rho)+etaU(rho)) if p > 0 else (-p)*(-(phiL(Ah, rho)+etaL(rho)))
    return R, pref, p, qp, prod

from functools import lru_cache
@lru_cache(maxsize=None)
def _static(factors, h, k, w):
    R = Fr(0); lp = 0.0; sp = 0; qp = Fr(0); lprod = 0.0
    for f in factors:
        r, pr, p, q_, prod = factor_data(f, h, k)
        R += r; lp += pr; sp += p; qp += q_; lprod += prod(w)
    return R, lp, sp, qp, lprod
def term_bound(factors, h, k, t0, w, coeff=1.0, numer=None, Eglob=None, kz_lo=None, kz_hi=None):
    """log upper bound of |coeff * numer * prod factors| on an arc; returns (logbound_without_exp, R)
       the exponential factor exp(pi^2/64 * E * Re(1/z)) is returned separately via R (E = R/k^2)."""
    R, lp, sp, qp, lprod = _static(tuple(factors), h, k, w)
    # (k|z|)^{-sp/2}
    kz = kz_lo if sp > 0 else kz_hi
    lp += (-sp/2)*log(kz)
    lq = max(0.0, -float(qp))*t0                 # |q|^qp <= e^{max(0,-qp) t0}
    ln = log(numer(t0)) if numer else 0.0
    return log(abs(coeff)) + ln + lp + lq + lprod, R

# ---------------- functions as choices of representation ----------------
def reps(name):
    """list of alternative representations; each a list of (coeff, factors, numer)"""
    th = (1.0, TH, None)
    u0 = [[(c*2/2, [f], nb)] for c, f, nb in [(v[0], v[1], v[2]) for v in U0_FORMS.values()]]
    u1 = [[(c, [f], nb)] for c, f, nb in [(v[0], v[1], v[2]) for v in U1_FORMS.values()]]
    if name == 'mu2':      # mu2 = -TH + 2 U0
        return [[(1.0, TH, None)] + [(2*c, fs, nb) for c, fs, nb in r] for r in u0]
    if name == 'U1':       # U1 = TH/2 - U0/2   or  U1 directly
        return [[(0.5, TH, None)] + [(0.5*c, fs, nb) for c, fs, nb in r] for r in u0] + u1
    if name == 'U0':       # U0 directly  or  TH - 2 U1
        return u0 + [[(1.0, TH, None)] + [(2*c, fs, nb) for c, fs, nb in r] for r in u1]

MAIN = {  # main cusps, constant |C|, a (exponent), shift, k
  'mu2': dict(k=4, hs=(1, 3), C=sqrt(pi/2), coef=1.0, a=pi*pi/32),
  'U1':  dict(k=4, hs=(1, 3), C=sqrt(pi/2), coef=0.5, a=pi*pi/32),
  'U0':  dict(k=8, hs=(1, 3, 5, 7), C=sqrt(pi)/2, coef=1.0, a=pi*pi/128),
}
LOWER = {  # |M(n)| >= LOWER * cosh(2 sqrt(aN))/sqrt(N)
  'mu2': sqrt(2)*math.cos(3*pi/8),
  'U1':  math.cos(3*pi/8)/sqrt(2),
  'U0':  min(abs(math.cos(5*pi/16-pi*r/4)+math.cos(13*pi/16+3*pi*r/4)) for r in range(8)),
}

def classes():
    """class key (k mod 8, h mod 4) -> list of representative (h,k), k<=48"""
    cl = {}
    for k in range(1, 49):
        for h in range(k):
            if gcd(h, k) == 1:
                cl.setdefault((k % 8, h % 4), []).append((h, k))
    return cl

def class_table(name):
    """for each class and representation: per-term R (must be constant across reps) ; returns dict"""
    cl = classes(); out = {}
    for key, lst in cl.items():
        per = []
        for ri, rep in enumerate(reps(name)):
            Rs = []
            for coeff, fs, nb in rep:
                vals = set()
                for (h, k) in lst:
                    r = term_bound(fs, h, k, 0.01, 0.1, kz_lo=1, kz_hi=1)[1]
                    vals.add(r)
                assert len(vals) == 1, (name, key, ri, vals)
                Rs.append(vals.pop())
            per.append(Rs)
        kmin_all = min(k for h, k in lst)
        out[key] = (per, kmin_all, lst)
    return out

def kmin_excl(name, lst):
    km = MAIN[name]['k']
    ks = sorted(set(k for h, k in lst if k != km))
    return ks[0] if ks else None

def minor_setup(name):
    """choose rep per class minimising max exponent; return Eall and choices"""
    tab = class_table(name); choice = {}; Eall = Fr(0)
    for key, (per, _, lst) in tab.items():
        km = kmin_excl(name, lst)
        best = None
        for ri, Rs in enumerate(per):
            E = max(max(R, 0) for R in Rs)/km**2
            if best is None or E < best[0]: best = (E, ri)
        choice[key] = (best[1], km, lst[0])
        Eall = max(Eall, best[0])
    return Eall, choice

def bound_parts(name, n, Eall, choice):
    t0 = (pi/8)*sqrt(float(Eall)/n)
    D = sqrt(4*pi*t0); w = 1/(4*pi)
    # ---- minor arcs: sup over classes of sum of term bounds (each term exp <= (pi^2/64) Eall / t0)
    cmin = -1e300; worst = None
    for key, (ri, km, (h, k)) in choice.items():
        tot = []
        for coeff, fs, nb in reps(name)[ri]:
            lb, R = term_bound(fs, h, k, t0, w, coeff, nb, kz_lo=t0, kz_hi=D)
            tot.append(lb)
        m = max(tot); v = m + log(sum(exp(x-m) for x in tot))
        if v > cmin: cmin, worst = v, key
    lminor = n*t0 + (pi*pi/64)*float(Eall)/t0 + cmin
    # ---- major arcs
    P = MAIN[name]; k = P['k']; a = P['a']; C = P['C']; N = n - 1/8
    X = 1/(t0 + 2*pi/k**2)                         # Re(1/z) lower bound on main arcs
    zmax = sqrt(t0*t0 + 2*pi*t0/k**2)
    # TH - main: |C| t0^{-1/2} e^{a X}(e^Lambda - 1), Lambda over eta factors (|p| sum -log(1-rho_m^n))
    Lam = 0.0
    for _, m, p in TH:
        g = gcd(m, k); rho = exp(-4*pi*pi*g*g*X/(m*k*k)); j = 1
        while rho**j > 1e-18: Lam += abs(p)*(-math.log1p(-rho**j)); j += 1
    b_th = P['coef']*C*t0**-0.5*exp(a*X)*math.expm1(Lam)
    # the other terms of the representation used on the main arc: pick rep with TH first, minimal bound
    best_other = None
    for rep in reps(name):
        if rep[0][1] is not TH: continue
        others = rep[1:]; lv = []
        for coeff, fs, nb in others:
            for h in P['hs'][:1]:
                lb, R = term_bound(fs, h, k, t0, X/k**2, coeff, nb, kz_lo=t0, kz_hi=k*zmax)
                Ek = float(R)/k**2
                lv.append(lb + (pi*pi/64)*max(Ek, 0)/t0)
        mm = max(lv); val = mm + log(sum(exp(x-mm) for x in lv))
        if best_other is None or val < best_other: best_other = val
    lb_th = log(b_th); mm = max(lb_th, best_other)
    lmaj = n*t0 + mm + log(exp(lb_th-mm) + exp(best_other-mm))
    # tails: 2 per main arc
    Y = 2*pi/(k*(sqrt(2*pi/t0)+k)); at = a*t0/Y**2
    ntails = 2*len(P['hs'])
    lt = log(ntails*P['coef']*C/(2*pi*N)) + N*t0 + at + log(2*Y**-0.5 + (2*a/3)*Y**-1.5)
    m = max(lminor, lmaj, lt)
    return m + log(exp(lminor-m)+exp(lmaj-m)+exp(lt-m)), lminor, lmaj, lt, worst

def log_main_lower(name, n):
    a = MAIN[name]['a']; N = n - 1/8; x = 2*sqrt(a*N)
    return log(LOWER[name]) + x + math.log1p(exp(-2*x)) - log(2) - 0.5*log(N)

if __name__ == '__main__':
    import sys
    names = sys.argv[1:] or ['mu2', 'U1', 'U0']
    for name in names:
        Eall, choice = minor_setup(name)
        print(f'== {name}: minor-arc exponent Eall = {Eall} (units pi^2/64), main exponent = {MAIN[name]["a"]/(pi*pi/64):.3f}')
        tab = class_table(name)
        for key in sorted(choice):
            ri, km, rep = choice[key]
            print(f'   class {key}: rep {ri}, kmin {km}, term R values {[str(r) for r in tab[key][0][ri]]}')
        last_bad = None; prev = None; dec = 0
        for n in range(1, 200001):
            mg = log_main_lower(name, n) - bound_parts(name, n, Eall, choice)[0]
            if mg <= 0: last_bad = n
            elif last_bad is not None and prev is not None and n > last_bad + 1 and mg < prev - 1e-9: dec += 1
            prev = mg
        N0 = (last_bad or 0) + 1
        print(f'   N0 = {N0}; margin decreases beyond N0 at {dec} values of n <= 2e5')
        for n in (N0, N0+200, 20000, 60000, 200000, 10**6, 10**9, 10**12):
            tot, a, b, c, wk = bound_parts(name, n, Eall, choice)
            print(f'   n={n}: log B={tot:.2f} (minor {a:.1f}, major {b:.1f}, tail {c:.1f}); log lower main={log_main_lower(name, n):.2f}; margin {log_main_lower(name, n)-tot:.3f}; worst class {wk}')
