"""Main computation of the paper (Proposition 6.1): certified scan of the explicit bound B(n) of Theorem 5.3
against the main-term lower bound, for N0 <= n <= 2*10^5, and sparse values beyond (Lemma 6.2).
Usage: python3 certify_run.py [S0 S1 mu2 U1 U0]"""
import sys, math, time
import certify                       # must be imported first: installs certified product bounds
SLACK = 1e-8
def scan(name, margin_fn, nmax=200000):
    t = time.time(); last_bad = None; prev = None; dec = 0
    for n in range(1, nmax+1):
        mg = margin_fn(n) - SLACK
        if mg <= 0: last_bad = n
        elif last_bad is not None and prev is not None and n > last_bad+1 and mg < prev - 1e-9: dec += 1
        prev = mg
    N0 = (last_bad or 0) + 1
    print(f'{name}: N0 = {N0}; margin at N0 = {margin_fn(N0):.3e}; decreases on [N0, {nmax}]: {dec}; margin at {nmax}: {margin_fn(nmax):.2f}  ({time.time()-t:.0f} s)')
    print('    margins at 10^6..10^12:', [round(margin_fn(10**j), 1) for j in range(6, 13)])
which = sys.argv[1:] or ['S0', 'S1', 'mu2', 'U1', 'U0']
if 'S0' in which:
    import farey_N0 as A, farey_bound as FB
    scan('S0', lambda n: FB.log_main_lower(n) - A.bound_parts(n)[0])
if 'S1' in which:
    import farey_S1_N0 as B
    scan('S1', lambda n: B.log_main_lower(n) - B.bound_parts(n)[0])
import farey_gen as G
for nm in ('mu2', 'U1', 'U0'):
    if nm in which:
        E, ch = G.minor_setup(nm)        # also checks that every cusp class has a constant exponent R/k^2
        print(f'{nm}: largest exponent away from the dominant cusps = {E} (units pi^2/64)')
        scan(nm, lambda n, nm=nm, E=E, ch=ch: G.log_main_lower(nm, n) - G.bound_parts(nm, n, E, ch)[0])
