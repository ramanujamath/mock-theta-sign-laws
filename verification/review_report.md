# Referee report: sign laws for S0, S1, mu2, U0, U1 via an explicit circle method

Verdict: I found **no FATAL error**. Each bound I re-derived holds, and each numerical claim I re-ran reproduces.
The remaining gaps are about provenance and documentation: some inputs are cited or numerical, some certification runs have no saved output, and the argument for n > 2·10^5 is only a sketch.
None of the existing files was modified. My scripts are `review_*.py`, and their outputs are `review_cert_*.out`.

## Issues, ranked

### SERIOUS (repairable)
1. **Certification runs are missing for 4 of the 5 functions.**
   - `cert_mu2.txt`, `cert_U1.txt` and `cert_U0.txt` hold only the sanity line. They contain no N0 result, so those runs died or were killed.
   - There is no `cert_S0.txt`. Only `cert_S1.txt` has a certified N0 (10652).
   - I re-ran the outward-rounded (certify-patched) scans myself (`review_cert.py`). They confirm N0 = 11231, 976, 493 and 5180 exactly.
   - The certified log-margins at N0 are: S0 1.03e-4, S1 4.27e-4, mu2 7.35e-4, U1 1.03e-2, U0 1.66e-3. At N0−1 every margin is negative.
   - The margins increase monotonically on [N0, 60000].
   - Fix: re-run `certify_run.py` and keep the logs.
2. **The q-series ↔ Appell–Lerch inputs are not proved.**
   - Forms A and B of S0 and S1 (farey_S0.md §0, I1) were checked only numerically at one point (`al.py`, `altforms.py`).
   - The U0/U1 single-m forms were checked as power series only to q^2999.
   - mu2 = 2U0 − TH was checked exactly to q^20000. U0 + 2U1 = TH is reduced to Hickerson–Mortenson (2.7) in `prod_identity.py`.
   - I re-checked all 10 identities at 12 random q with 0.2 < |q| < 0.75 (`review_forms.py`). The maximum discrepancy is 2e-14.
   - Fix: cite Gordon–McIntosh and Hickerson–Mortenson (4.39)/(4.40) for each form, or prove them, for example with a Sturm-type bound for the finite q-expansion check.

### MINOR
3. **The argument for n > 2·10^5 (S0 §8, item [S]) is only a sketch.**
   - The dense scans reach 2·10^5 (float) and 60000 (certified). Beyond that there are sparse checks only.
   - B(n) is an explicit smooth function of n, with no floor(K) dependence in `farey_N0`, `farey_S1_N0` and `farey_gen`, so the analytic monotonicity argument is easy. It should be written out.
   - In my certified sparse checks the margin grows: about 65–212 at 1e5 and about 9e5–2.3e6 at 1e13.
4. **Main-term phases are fixed numerically [N].**
   - The phases are e^{9πi/16} (S0), e^{−3πi/16} (S1), and for TH 3π/8 at 1/4, 5π/16 at 1/8, −13π/16 at 3/8.
   - I confirmed them independently in two ways:
     - At the cusp, |f − main| ≈ 1.40 while |main| ≈ 7·10^19, a relative agreement of 2e-20 (`review_center.py`).
     - |c(n) − M(n)| is e^{16}–e^{40} below B(n) at 63 values of n.
   - A one-line multiplier-system argument (the constant is a root of unity of bounded order) should be added.
5. **Some float pieces are not certified.**
   - `Lambda`, `Lambda1` and the TH Λ in `farey_gen` use plain floats, are truncated at 1e-18, and have no tail bound. All exp/log/sqrt in the assembly are also plain floats.
   - This is harmless here: the major-arc term is about e^{-100} below the minor-arc term, and float error (about 1e-13) is far below the margins (≥ 1e-4). It is still not interval arithmetic.
6. **`farey_gen.bound_parts` bounds the major-arc "other" terms only at h = hs[0].**
   - This is valid, because complex conjugation maps h → k−h, and for k = 8 the cusps 5/8 and 1/8 are in the same class (h mod 4). I also confirmed |C| is the same at all four 1/8-type cusps.
   - The justification should be documented.
7. **Wording.** The claim lists n = 0 as a "zero" for S0 and U0, but c(0) = 1 there. It is an exception, not a zero.

## What I verified as correct (independently)
- **Finite checks.** I re-ran the exception scan from `big.pkl` and `family.pkl`, and every claimed exception list matches exactly. I recomputed all 5 series independently from their definitions to q^1200, and they agree with the pickles.
- **Lemma 1.** I re-derived the q-power M(−1/8 + A/2 − A²/2). The q-powers are:
  - TQ: q^{1/16}
  - TQ1: q^{−7/16}
  - 1/Θ(−1;q^8): |q|^1
  - 1/Θ(−q^{2,6};q^8): |q|^{1/4}
  - TH: q^{1/8}

  The ≤1 powers are correctly dropped and the ≥1 powers (S1's q^{−1}, e^{7t0/16}, and the e^{t0}, e^{t0/4}, e^{t0/8}, e^{t0/2} numerator factors) are included.
- **Transformation.** I re-derived |k'T − h'|, ρ = exp(−4π²g²Re(1/z)/(Mk²)), the η prefactor (2πg/(Mk|z|))^{1/2}, and Â. The prefactor bounds go in the correct direction: t0 for positive η-power, √(4πt0) for negative.
- **Cusp exponents.** I computed them by an independent Gauss-sum/Poisson method (`review_exponents.py`, no Siegel machinery) at all cusps with k ≤ 120. All 52,632 factor exponents agree with `farey_lib` and `farey_gen`.
  - Over all cusps with k ≤ 200, the maximum off-main exponent is 1/4 (S0, S1), 1/2 (mu2), 2/9 (U1) and 1/8 (U0). These equal the Eall values used.
  - The class (k mod 8, h mod 4) determines R for all k, not only k ≤ 48. k = 1 and k = 2 are included.
- **Lemma 2.** It is correct, including the Â = 0 case. For base q^4, B̂ ∈ (1/4)Z∖Z, so |1−e(B̂)| ≥ √2 > 2 sin(π/8), and the bound is conservative.
- **Lerch numerators.** I re-derived every one by hand:
  - The minimum denominator exponents are 3,3,1 (S0), 1,3,1 (S1), 1 and 4 (U0, including the constant 2 at r=0) and 1,2 (U1).
  - The shifts and Gaussian sums give e^{t0}, e^{t0/4}, e^{t0/2}, e^{t0/8}, and 1 + √(π/c) with c = 4t0 (base q^8) or c = 2t0 (base q^4).
- **Farey geometry.** All four inequalities hold.
- **Major arcs.**
  - |F−1| ≤ e^Λ−1 holds.
  - Monotonicity holds: the minimum ρ-exponent exceeds the e^{a/z} exponent, for S0/S1 (1/4 > 1/32) and for TH (π²/4 > π²/32 at k=4, π²/16 > π²/128 at k=8).
- **Bessel and tail lemma.** I re-derived the integration by parts. The bound e^{Nt0}e^{at0/Y²}(2Y^{−1/2} + (2a/3)Y^{−3/2})/(2πN) is correct, and the tail counts are correct (4 for S0/S1/mu2/U1, 8 for U0). N = n−1/16, n+7/16, n−1/8 is used consistently.
- **certify.py.** The tail bounds are valid. Rounding ρ up to the 1e-7 grid, with BUMP, is conservative everywhere, including for tiny ρ (which rounds up to 1e-7). The patch is applied before the bound modules are imported.
- **Task 3.** At 10 random n in [N0, top] plus N0, N0+1 and top for every function, |c−M| ≤ B and B < LOWER·cosh/√N ≤ |M|. The sign of M matches the law. The worst log(|c−M|/B) is −15.6.
- **Pointwise sup bounds.** I checked these on the circle (FFT of exact coefficients, Farey-classified, 10 (function, n) cases):
  - Minor arcs are e^{-9.6} to e^{-31.6} below their bounds.
  - Major arcs are below their bounds wherever float precision can resolve them. The apparent FFT "violations" at n = 3000 are float roundoff: a 70-digit direct evaluation gives |f−main| ≈ 1.4, about e^{-6} below the bound.
