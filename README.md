# mock-theta-sign-laws

Computer-assisted proof of periodic sign laws for the coefficients of the eighth-order mock theta functions S₀, S₁, U₀, U₁ (Gordon–McIntosh) and the second-order function μ (McIntosh).

**Paper:** *Sign laws for the coefficients of second- and eighth-order mock theta functions*, by Ramanuja. See [`paper/signlaws.pdf`](paper/signlaws.pdf) (source: [`paper/signlaws.tex`](paper/signlaws.tex)).

## Main Theorem

For every n ≥ 0 the sign of the n-th coefficient c(n) is as follows, except at the listed n.

| Function | Sign of c(n) | Exceptions |
| --- | --- | --- |
| S₀ | −, +, +, − for n ≡ 0, 1, 2, 3 (mod 4) | c(0) = 1; c(n) = 0 for n = 4, 6, 8, 12, 14, 20, 30 |
| S₁ | +, −, −, + for n ≡ 0, 1, 2, 3 (mod 4) | c(n) = 0 for n = 1, 2, 9 |
| μ | −, −, +, + for n ≡ 0, 1, 2, 3 (mod 4) | c(0) = 1 |
| U₁ | +, +, −, − for n ≡ 0, 1, 2, 3 (mod 4) | c(n) = 0 for n = 0, 2 |
| U₀ | −, +, +, −, +, −, −, + for n ≡ 0, …, 7 (mod 8) | c(0) = 1; c(n) = 0 for n = 3, 5, 8, 12, 16, 19, 24, 28, 40, 60 |

The proof combines asymptotics with exact constants, for example

c_{S₀}(n) ~ cos(9π/16 − πn/2) · cosh(π√(n − 1/16)/4) / √(n − 1/16),

with an explicit circle-method error bound that forces the sign law for n ≥ N₀. The cutoffs are N₀ = 11231 (S₀), 10652 (S₁), 976 (μ), 493 (U₁) and 5180 (U₀). An exact computation covers all smaller n.

## Code

The directory [`code/`](code) contains the programs needed for the proof (Python 3, standard library only). Run them from `code/`:

| Script | Purpose |
| --- | --- |
| [`coeffs.py`](code/coeffs.py) | exact integer coefficients and the sign check for n ≤ 12000 |
| [`exponents_S.py`](code/exponents_S.py) | cusp-exponent classes for S₀ and S₁ (Prop. 3.6) |
| [`exact_phase.py`](code/exact_phase.py) | exact phases of the main-term constants (Prop. 4.2) |
| [`certify_run.py`](code/certify_run.py) | certified verification of Prop. 6.1: N₀ and margins for all five functions |

Library modules: [`certify.py`](code/certify.py), [`farey_lib.py`](code/farey_lib.py), [`farey_bound.py`](code/farey_bound.py), [`farey_N0.py`](code/farey_N0.py), [`farey_S1_N0.py`](code/farey_S1_N0.py), [`farey_gen.py`](code/farey_gen.py).

```
cd code
python3 coeffs.py        # ALL OK
python3 certify_run.py   # S0: N0 = 11231, S1: 10652, mu: 976, U1: 493, U0: 5180  (a few minutes)
```

## Verification

The folder [`verification/`](verification) contains independent checks and exploratory computations that are not part of the proof. They include direct and high-precision validation of every bound, an independent recomputation of the cusp exponents, the referee report, and the searches that led to the results. Its [README](verification/README.md) gives usage and the results of each script, and [`verification/results/`](verification/results) holds the outputs of a full run.

## Contact

Ramanuja, Kachhubi, P.O. Tamuria, Lakhnaur Block, Madhubani, Bihar 847410, India. ramanujamath@yahoo.com
