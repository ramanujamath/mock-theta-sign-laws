"""Certified (outward-rounded) product constants, patched into farey_lib before the bound modules are imported."""
import math
from decimal import Decimal as Dm, getcontext
getcontext().prec = 50
import farey_lib
BUMP = 1 + 1e-13          # rho is used as an UPPER bound; U increasing, L decreasing in rho -> bump up is conservative
def up(x):   return math.nextafter(math.nextafter(float(x), math.inf), math.inf)
def down(x): return math.nextafter(math.nextafter(float(x), -math.inf), -math.inf)
def _rho(rho): return Dm(repr(rho*BUMP))
def _pw(r, e): return (r.ln()*Dm(repr(float(e)))).exp() if e > 0 else Dm(1)
def phiU(Ah, rho, nmax=None):
    r = _rho(rho); Ah = float(Ah)
    v = (1 + _pw(r, Ah)).ln() if Ah > 0 else Dm(2).ln()
    n = 1
    while True:
        a = _pw(r, n+Ah); b = _pw(r, n-Ah)
        v += (1+a).ln() + (1+b).ln()
        if b < Dm('1e-30'):
            v += (a + b)*r/(1-r)/(1-r) + Dm('1e-40')   # tail sum_{m>n} (r^{m+Ah}+r^{m-Ah}) <= (a+b) r/(1-r), extra /(1-r) slack
            return up(v)
        n += 1
def phiL(Ah, rho, nmax=None):
    r = _rho(rho); Ah = float(Ah)
    v = (1 - _pw(r, Ah)).ln() if Ah > 0 else (2 - Dm(2).sqrt()).sqrt().ln()   # 2 sin(pi/8) = sqrt(2 - sqrt 2)
    n = 1
    while True:
        a = _pw(r, n+Ah); b = _pw(r, n-Ah)
        v += (1-a).ln() + (1-b).ln()
        if b < Dm('1e-30'):
            v -= 2*(a + b)*r/(1-r) + Dm('1e-40')        # -log(1-x) <= 2x for x<=1/2
            return down(v)
        n += 1
def etaU(rho):
    r = _rho(rho); v = Dm(0); n = 1
    while True:
        a = _pw(r, n); v += (1+a).ln()
        if a < Dm('1e-30'): return up(v + a*r/(1-r) + Dm('1e-40'))
        n += 1
def etaL(rho):
    r = _rho(rho); v = Dm(0); n = 1
    while True:
        a = _pw(r, n); v += (1-a).ln()
        if a < Dm('1e-30'): return down(v - 2*a*r/(1-r) - Dm('1e-40'))
        n += 1
for f in (phiU, phiL, etaU, etaL): setattr(farey_lib, f.__name__, f)
# sanity: compare with float versions
import importlib
# ---- caching on an upward-rounded rho grid (conservative: U increasing, L decreasing in rho) ----
from functools import lru_cache
_GRID = 10**7
def _g(rho): return math.ceil(rho*_GRID)/_GRID
_cU = lru_cache(maxsize=None)(lambda Ah, r: phiU.__wrapped__(Ah, r) if hasattr(phiU,'__wrapped__') else None)
_raw = dict(phiU=phiU, phiL=phiL, etaU=etaU, etaL=etaL)
@lru_cache(maxsize=None)
def _phiU(Ah, r): return _raw['phiU'](Ah, r)
@lru_cache(maxsize=None)
def _phiL(Ah, r): return _raw['phiL'](Ah, r)
@lru_cache(maxsize=None)
def _etaU(r): return _raw['etaU'](r)
@lru_cache(maxsize=None)
def _etaL(r): return _raw['etaL'](r)
def phiU_c(Ah, rho, nmax=None): return _phiU(float(Ah), _g(rho))
def phiL_c(Ah, rho, nmax=None): return _phiL(float(Ah), _g(rho))
def etaU_c(rho): return _etaU(_g(rho))
def etaL_c(rho): return _etaL(_g(rho))
farey_lib.phiU, farey_lib.phiL, farey_lib.etaU, farey_lib.etaL = phiU_c, phiL_c, etaU_c, etaL_c
