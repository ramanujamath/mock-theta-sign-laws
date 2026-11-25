import sys, math, time
import certify                      # patch first
import farey_lib
from math import pi, exp
print('check certified vs float at rho=e^{-pi/8}: phiU(1/4)', certify.phiU(0.25, exp(-pi/8)), ' phiL(0)', certify.phiL(0, exp(-pi/8)))
SLACK = 1e-8
def scan(name, margin_fn, nmax=200000):
    t = time.time(); last_bad = None; prev = None; dec = 0
    for n in range(1, nmax+1):
        mg = margin_fn(n) - SLACK
        if mg <= 0: last_bad = n
        elif last_bad is not None and prev is not None and n > last_bad+1 and mg < prev - 1e-9: dec += 1
        prev = mg
    print(f'{name}: certified N0 = {last_bad+1}; margin non-monotone points beyond N0: {dec}; margins at N0: {margin_fn(last_bad+1):.6f}, at 2e5: {margin_fn(nmax):.2f}  ({time.time()-t:.0f}s)')
which = sys.argv[1:] or ['S0','S1','mu2','U1','U0']
if 'S0' in which:
    import farey_N0 as A, farey_bound as FB
    scan('S0', lambda n: FB.log_main_lower(n) - A.bound_parts(n)[0])
if 'S1' in which:
    import farey_S1_N0 as B
    scan('S1', lambda n: B.log_main_lower(n) - B.bound_parts(n)[0])
import farey_gen as G
for nm in ('mu2','U1','U0'):
    if nm in which:
        E, ch = G.minor_setup(nm)
        scan(nm, lambda n, nm=nm, E=E, ch=ch: G.log_main_lower(nm, n) - G.bound_parts(nm, n, E, ch)[0])
