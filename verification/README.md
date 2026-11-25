# Verification and exploration scripts

This folder holds the independent checks, cross-validations and exploratory computations behind the paper. **None of it is needed for the proof**; the proof uses only the scripts in [`../code`](../code). These scripts check the proof's ingredients in different ways: by direct summation, high precision, independent re-derivation, and an independent referee pass.

The folder is self-contained. It has its own copies of the library modules (`farey_lib.py`, `farey_bound.py`, `farey_N0.py`, `farey_S1_N0.py`, `farey_gen.py`, `certify.py`, `certify_run.py`). These are the untrimmed versions, which include helpers used only for validation, such as direct theta summation. Python 3; some scripts also need `numpy`.

## Usage

```
cd verification
sh run_all.sh                    # regenerate data, run everything, write results/*.txt  (~20 min on 2 cores)
```

Or run the steps by hand. First generate the data files (about 8 MB of pickles, not committed):

```
python3 big.py 60000             # exact coefficients of S0, S1, mu2 for n <= 60000      -> big.pkl     (20 s)
python3 family.py 20000          # exact coefficients of U0, U1, T0, T1, V0, V1, n <= 20000 -> family.pkl (5 s)
python3 mock.py                  # third-order phi, chi, n <= 1500                        -> mt.pkl
python3 scan.py                  # other orders, sign-period scan                          -> F.pkl
python3 th_check.py              # TH = J2^8/(J1^3 J4^4) coefficients                      -> TH.pkl     (3 min)
```

Then run any script with `python3 <script>`. The arguments are in `run_all.sh`. The outputs of a full run are saved in [`results/`](results), one file per script; `results/run_all.log` lists exit codes and timings. In that run `th_check.py` exited with an error: after its identity checks it also evaluated the TH constants and hit a floating-point underflow. That part now lives in `th_const.py`, and `results/th_check.txt` is from the corrected script.

Some scripts import `al.py` or `radial.py`. These print a three-line self-check (lines beginning `S0 (`, `S1 (`, `mu (`) at the top of the output; those lines are not part of the script's own result.

## Results

### Finite checks and identities

| Script | What it checks | Result |
| --- | --- | --- |
| `review_finite.py` | sign laws from `big.pkl`, `family.pkl` | exceptions exactly as in Table 1, for n ≤ 60000 (S₀, S₁, μ) and n ≤ 20000 (U₀, U₁) |
| `review_series.py` | the five series recomputed independently from their definitions | match the pickles through q¹²⁰⁰ |
| `th_check.py` | U₀ + 2U₁ = TH and μ = 2U₀ − TH as power series | both true through q²⁰⁰⁰⁰ |
| `prod_identity.py` | U₀ + 2U₁ = TH via the change-of-z formula, as an infinite-product identity | formula (2.7) matches to 1.3·10⁻¹⁵; constant, q-power and exponents of (1 − qⁿ) for n ≤ 64 agree |
| `al.py`, `altforms.py` | Appell–Lerch forms A and B of S₀, S₁, μ at a sample point; search for alternative forms | agree to ~10⁻¹⁵; finds S₁ = −q⁻¹[m(−q,q⁸,−q²) + m(−q,q⁸,−q⁶)] and the μ forms |
| `repcheck2.py` | the S₁ form B and μ = 2[m(−q,q⁴,q) + m(−q,q⁴,q³)] as exact power series | true through q²⁹⁹⁹ |
| `review_forms.py` | all ten representations (S₀/S₁ forms A, B; U₀, U₁ forms; TH; μ) at 12 random q | max deviation 2.3·10⁻¹⁴ |

### Main-term constants

| Script | What it checks | Result |
| --- | --- | --- |
| `exact.py` | saddle-point amplitude and phase identities for S₀, S₁ in exact algebra over ℚ(x), x⁴ + 2x² − 1 = 0 | amplitude identities true; phases 7π/8 and −5π/8 confirmed by the "tan identity" lines (the intermediate "tan check" printed for S₁ uses a different expression; the following "S1 tan identity" line is the correct check and is 0) |
| `exact_phase.py` | phases of C_ζ from exact finite Gauss sums | 9π/16 (S₀), −3π/16 (S₁), 3π/8, 5π/16, −13π/16 (TH), residual ≤ 4.4·10⁻¹⁶ |
| `th_const.py`, `review_phase.py` | TH(ζe⁻ᶻ) / (z^(−1/2) e^(−z/8+a/z)) at small z | \|C\| = √(π/2) at i and √π/2 at the 8th roots; phases as above, constant in z to 8–10 digits |
| `radial.py`, `pref.py` | radial growth of S₀, S₁, μ at q = ζe⁻ᵗ | exponent π²/64t (S₀, S₁), π²/32t (μ) at ±i only; amplitudes √π/2, √π/2, √(π/2); phases 0.5625π, −0.1875π, −0.625π |
| `u0b.py` | the same for U₀ (8th roots) and U₁ (±i) | U₀: √π/2 with phases 5π/16, −13π/16; U₁: √(π/8) with phase 3π/8 |
| `verify.py`, `resid.py` | c(n)/M(n) and c(n) − M(n) for S₀, S₁, μ | ratio → 1 (0.99972 for μ at n = 60000); relative error 10⁻⁴² (S₀), 9·10⁻⁴³ (S₁), 5·10⁻⁶⁰ (μ) at n = 59999 |
| `ufit.py` | asymptotic formulas for U₀, U₁ | shift n − 1/8; relative error 3·10⁻¹⁵ (U₀), 2·10⁻¹⁶ (U₁) at n ≈ 20000; U₀ trig factor ≥ 0.2759 |

### Cusp exponents

| Script | What it checks | Result |
| --- | --- | --- |
| `cusps.py`, `cusps2.py`, `lemma.py` | exponents via Poisson / Gauss sums for k ≤ 400; the ν_eff ≤ 6M lemma | max away from main cusps: 1/4 (S₀, S₁), 1/2 (μ); observed ν_eff ≤ M |
| `farey_exponents.py`, `farey_S1_exp.py` | exact Bernoulli-polynomial class tables (Siegel route) | max 1/4 away from k = 4, all k |
| `review_exponents.py`, `review_cmpR.py` | independent recomputation of all factor exponents | 52 632 exponents (k ≤ 120), 0 mismatches; maxima 1/4 (S₀, S₁), 1/2 (μ), 2/9 (U₁), 1/8 (U₀) |
| `u_cusps.py` | U₀, U₁ representations and exponents, k ≤ 400, with high-precision radial fits | U₀: 1/2 at the 8th roots, ≤ 1/8 elsewhere; U₁: 2 at ±i, ≤ 2/9 elsewhere; fits agree (e.g. 0.2214 vs 2/9 at 5/12) |
| `fit3.py`, `scan2.py`, `split.py`, `u0.py` | early least-squares fits of radial growth | exploratory; the fits on the m-forms in `u0.py` are numerically unreliable and are superseded by `u0b.py` and `u_cusps.py` |

### Circle-method bounds

| Script | What it checks | Result |
| --- | --- | --- |
| `farey_check_transform.py` | Siegel/η transformation vs direct summation | max \|log(direct/transformed)\| = 1.1·10⁻⁷ |
| `farey_validate.py` | S₀: V1 theta bounds (3000 pts), V2 Farey geometry, V3 sup bound, V4 Bessel identity, V5 main arc, V6 Lerch numerators, V7 actual error vs bound | 0 violations; min Re(1/z)/k² = 0.0805 ≥ 1/(4π); bound e^31.6 above sampled \|S₀\|; Bessel to 3·10⁻⁷; main-arc bound holds; Lerch bounds 4–70× above samples; actual error e^36.0 vs bound e^76.2 at N₀ |
| `farey_S1_validate.py`, `farey_S1_v1dbg.py` | the same for S₁ | 24 apparent factor violations from direct summation, all roundoff (transformed values lie inside the bounds in every case); main arc holds; error e^35.1 vs bound e^75.1 at N₀ |
| `farey_gen_validate.py`, `farey_gen_v1b.py` | μ, U₁, U₀: base-q⁴ factor bounds, whole-function bounds, actual errors | 36 000 factor checks via the triple product, 0 violations; ~2250 random minor-arc points, all below the bound by ≥ e^1.8; e.g. μ error e^72.7 vs bound e^93.0 at n = 19999 |
| `farey_gen_v5.py` | main-arc identity for TH | holds at all complex test points; the "False" rows at real z are double-precision roundoff (\|TH/main − 1\| ≈ 10⁻¹³, true value ~10⁻⁵⁰) |
| `review_circle.py`, `review_circle2.py` | pointwise sup bounds on the whole circle via FFT of exact coefficients | minor arcs e^−9.6 to e^−31.6 below the bound; remaining major-arc excesses at n = 3000 (S₀, S₁) are FFT roundoff, resolved by `review_center.py` |
| `review_center.py`, `review_major.py` | 70-digit evaluation of \|f − main\| near q = i | e^−6.15 (S₀) and e^−7.24 (S₁) below the bound at n = 3000; the single large ratio in `review_major` at the cusp centre is double-precision roundoff (relative 10⁻¹³) |
| `review_verify.py` | \|c − M\| ≤ B < L·cosh/√N ≤ \|M\| at 13 n per function, sign of M | all true; worst log(\|c − M\|/B) = −15.6 |
| `review_cert.py` | certified margins near N₀ and up to n = 10¹³ | e.g. S₀ margin 1.24·10⁶ and μ margin 1.76·10⁶ at n = 10¹³ |
| `certify_run.py` | certified N₀ scan (copy of the proof run) | N₀ = 11231 (S₀), 10652 (S₁), 976 (μ), 493 (U₁), 5180 (U₀) |

### Exploration

| Script | Result |
| --- | --- |
| `mock.py` | third order: φ has signs + + − − (mod 4), zeros n = 2, 8, 10, 26; χ has signs + + + − − − (mod 6), 21 zeros up to n = 80 |
| `scan.py` | sign periods: f₀, f₁, φ₆, ψ₆ alternate; μ, S₀, S₁ have period 4 |
| `big.py`, `family.py` | T₀, T₁ alternate; V₀, V₁ positive; U₀ period 8, U₁ period 4 |

## Notes

- `repcheck.py` is the first version of `repcheck2.py`. It mishandled signs for negative r, so it is superseded and not run.
- `review_report.md` is the referee report from the independent review pass. It found no fatal error; its "serious" items have been resolved (see the paper, Sections 2 and 6).
- `farey_S0.md` is the working write-up of the S₀ circle-method bound; `u_cusps_report.txt` is the report of the U₀/U₁ cusp analysis.
