"""Exact integer coefficients of S0, S1, mu2, U0, U1 from their defining q-series, and the sign-law check
(Main Theorem, finite range).  Usage: python3 coeffs.py [N]   (default N = 12000 > max N0 = 11231)."""
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
def mulb(a, m, s):            # a *= (1 + s q^m)
    r = a[:]
    for i in range(N, m-1, -1): r[i] += s*a[i-m]
    return r
def divb(a, m, s):            # a /= (1 + s q^m)
    r = a[:]
    for i in range(m, N+1): r[i] -= s*r[i-m]
    return r
def series(expo, step, sign=lambda n: 1, init=lambda c: c):
    """sum_n sign(n) q^{expo(n)} P_n(q), with P_0 = init(1) and P_n = step(P_{n-1}, n)"""
    tot = [0]*(N+1); cur = [0]*(N+1); cur[0] = 1; cur = init(cur); n = 0
    while expo(n) <= N:
        if n > 0: cur = step(cur, n)
        e = expo(n); s = sign(n)
        for i in range(N+1-e): tot[i+e] += s*cur[i]
        n += 1
    return tot
step8 = lambda c, k: divb(mulb(c, 2*k-1, 1), 2*k, 1)                      # (1+q^{2k-1})/(1+q^{2k})
F = {
 'S0':  series(lambda n: n*n, step8),
 'S1':  series(lambda n: n*(n+2), step8),
 'mu2': series(lambda n: n*n, lambda c, k: divb(divb(mulb(c, 2*k-1, -1), 2*k, 1), 2*k, 1), lambda n: (-1)**n),
 'U0':  series(lambda n: n*n, lambda c, k: divb(mulb(c, 2*k-1, 1), 4*k, 1)),
 'U1':  series(lambda n: (n+1)**2, lambda c, k: divb(mulb(c, 2*k-1, 1), 4*k+2, 1), init=lambda c: divb(c, 2, 1)),
}
PATTERN = {'S0': '-++-', 'S1': '+--+', 'mu2': '--++', 'U1': '++--', 'U0': '-++-+--+'}
EXPECTED = {'S0': {0: 1, 4: 0, 6: 0, 8: 0, 12: 0, 14: 0, 20: 0, 30: 0}, 'S1': {1: 0, 2: 0, 9: 0}, 'mu2': {0: 1},
            'U1': {0: 0, 2: 0}, 'U0': {0: 1, 3: 0, 5: 0, 8: 0, 12: 0, 16: 0, 19: 0, 24: 0, 28: 0, 40: 0, 60: 0}}
ok_all = True
for name, c in F.items():
    p = PATTERN[name]; L = len(p)
    exc = {n: c[n] for n in range(N+1) if not ((c[n] > 0) if p[n % L] == '+' else (c[n] < 0))}
    ok = exc == EXPECTED[name]; ok_all &= ok
    print(f'{name}: pattern {p} for 0 <= n <= {N}; exceptions {exc}; matches Main Theorem: {ok}')
print('ALL OK' if ok_all else 'MISMATCH')
