Ancillary files for "Sign laws for the coefficients of second- and eighth-order mock theta functions" (Ramanuja).

Python 3, standard library only. Run from this directory.

  python3 coeffs.py          Exact integer coefficients of S0, S1, mu, U0, U1 for n <= 12000 and the check of
                             Table 1 in that range (Main Theorem for n < N0).
  python3 exponents_S.py     Exponent classes for S0 and S1 (Proposition 3.6).
  python3 exact_phase.py     Exact phases of the constants C_zeta (Proposition 4.2).
  python3 certify_run.py     Certified verification of Proposition 6.1 for all five functions (N0 and margins),
                             including the exponent classes for mu, U1, U0, and the sparse values for Lemma 6.2.

Library modules (imported by certify_run.py):
  certify.py         outward-rounded 50-digit evaluation of the products of Lemma 3.4 (must be imported first)
  farey_lib.py       Siegel-function data: B2, Ahat, product bounds (Lemmas 3.2-3.4)
  farey_bound.py     helper functions for S0 (Lemma 5.2, Proposition 4.2, main-term lower bound)
  farey_N0.py        the bound B(n) of Theorem 5.4 for S0
  farey_S1_N0.py     the bound B(n) of Theorem 5.4 for S1
  farey_gen.py       the bound B(n) of Theorem 5.4 for mu, U1, U0 (base q^4, eta quotient J2^8/(J1^3 J4^4))

Expected output of certify_run.py (a few minutes):
  S0:  N0 = 11231    S1: N0 = 10652    mu: N0 = 976    U1: N0 = 493    U0: N0 = 5180
