"""Cusp-exponent tables for S0 and S1 (Section 3 of the paper), exact in rational arithmetic.
For a product of factors Theta(s q^a; q^8)^p the growth exponent on the arc of h/k is E = R/k'^2
(units pi^2 Re(1/z)/64), k' = k/gcd(8,k), where R depends only on (h' mod 8, k' mod 2).
Form A is used for odd k', form B for even k'.  Output: the class table and the largest exponent
away from the dominant cusps 1/4, 3/4 (k' = 1 in class (odd, h' = 2 or 6) occurs only for k = 4)."""
from fractions import Fraction as Fr
from farey_lib import Ahat, B2
def R(fl, hm, kpar):
    s = Fr(sum(m for m, _, _ in fl), 24) + sum(m*B2(Ahat(A, s_, hm, kpar))/2 for m, A, s_ in fl)
    return -32*s                   # E = R/k'^2
FACT = {
 'S0': dict(TQ=[(1, Fr(1,8), Fr(1,2)), (2, Fr(1,4), Fr(0)), (-2, Fr(3,8), Fr(0))]),
 'S1': dict(TQ=[(1, Fr(3,8), Fr(1,2)), (2, Fr(1,4), Fr(0)), (-2, Fr(1,8), Fr(0))]),
}
ALA = [(-1, Fr(0), Fr(1,2))]; ALB1 = [(-1, Fr(1,4), Fr(1,2))]; ALB2 = [(-1, Fr(3,4), Fr(1,2))]
CL = [(1, h) for h in range(8)] + [(2, h) for h in (1, 3, 5, 7)]
def kmin(kpar, hm):                 # smallest admissible k' in the class, excluding the dominant cusps (k = 4)
    return 3 if (kpar == 1 and hm in (2, 6)) else kpar
for name, d in FACT.items():
    worst = Fr(0)
    print(f'== {name}')
    for kpar, hm in CL:
        if kpar == 1: terms = [R(d['TQ'], hm, 1), R(ALA, hm, 1)]
        else:         terms = [R(ALB1, hm, 2), R(ALB2, hm, 2)]
        E = max(max(terms), 0)/kmin(kpar, hm)**2
        worst = max(worst, E)
        print(f"  k' {'odd ' if kpar == 1 else 'even'} h'={hm}: R terms {[str(t) for t in terms]} -> max exponent {E}")
    print(f'  largest exponent away from 1/4, 3/4: {worst};  at 1/4: TQ {R(d["TQ"], 2, 1)}, AL {R(ALA, 2, 1)}')
