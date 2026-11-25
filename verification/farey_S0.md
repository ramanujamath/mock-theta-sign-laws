# Explicit circle-method bound for the sign law of S0

**Result.** Let c(n) be the coefficients of S0, N = n − 1/16, and

  M(n) = cos(9π/16 − πn/2) · cosh(π√N/4)/√N.

For every n ≥ 1 the bound |c(n) − M(n)| ≤ B(n) holds, where B(n) is the explicit function in
`farey_N0.py::bound_parts`. Roughly, B(n) ≈ 1.1·10^15 · (n/11231)^{3/4} · e^{(π/8)√n}, so θ = 1/2.
We also have B(n) < cos(7π/16)·cosh(π√N/4)/√N ≤ |M(n)| for every **n ≥ N0 = 11231**. Since the sign law was checked exactly
for n ≤ 60000, it then holds for all n, **subject to the status labels below**.

Status labels used in this file:
- [R] rigorous, with the proof given here.
- [R-std] rigorous, given a standard cited theorem.
- [N] the constant or fact comes from a numerical computation, as stated.
- [S] sketch only.

Scripts (all new; no existing file was changed):
- `farey_lib.py`: Siegel/eta machinery, product bounds, and exact Bernoulli exponents.
- `farey_exponents.py`: cusp exponents for all cusps.
- `farey_bound.py`: per-n optimised bound.
- `farey_N0.py`: uniform bound, the N0 scan, and checks up to 10^12.
- `farey_check_transform.py`: checks the transformation identity.
- `farey_validate.py`: checks V1 to V7.

---

## 0. Inputs taken from the project (not re-proved here)

* (I1) Form A: S0 = 2m(−q³,q⁸,−1) + TQ, with TQ := qΘ(−q;q⁸)Θ(q²;q⁸)²/Θ(q³;q⁸)².
  Form B: S0 = m(−q³,q⁸,−q²) + m(−q³,q⁸,−q⁶).
  Both identities were **checked numerically** in `al.py`/`altforms.py`, but I did not prove them. They are presumably
  the Gordon–McIntosh Appell–Lerch forms. **Status: [N] input.**
* (I2) The exact finite check: the sign law holds for 0 ≤ n ≤ 60000 (`big.pkl`).

## 1. Theta functions as Siegel functions [R]

Write q = e(τ), with e(x) = e^{2πix}. Put T = 8τ and Q = e(T). For s = e(σ), σ ∈ {0, 1/2}, put A = a/8. Let B₂(x) = x² − x + 1/6. The Siegel function is
g_{(A,B)}(T) = −Q^{B₂(A)/2} e(B(A−1)/2) (1−q_z) ∏_{n≥1}(1−Qⁿq_z)(1−Qⁿ/q_z), with q_z = e(AT+B).

The Jacobi triple product Θ(x;Q) = (x;Q)∞(Q/x;Q)∞(Q;Q)∞ gives, for 0 ≤ A < 1:

**Lemma 1.** |Θ(s q^a; q⁸)| = |q|^{8(−1/8 + A/2 − A²/2)} · |g_{(A,σ)}(8τ)| · |η(8τ)|.

*Proof.* Take the triple product and compare the factors. The Q-power is −B₂(A)/2 − 1/24 = −A²/2 + A/2 − 1/8. ∎

The q-powers that occur are these:
- TQ: q^{1/16}.
- Θ(−1;q⁸): |q|^{−1}.
- Θ(−q²;q⁸) and Θ(−q⁶;q⁸): |q|^{−1/4}.

In every bound below they contribute factors ≤ 1.

## 2. Transformation to the cusp [R-std]

Let q = e(h/k)e^{−z} with Re z = t0. Then τ = h/k + iz/2π. Let g = gcd(8,k), k' = k/g and h' = 8h/g; then gcd(h',k') = 1.
Take γ = (α β; k' −h') ∈ SL₂(Z) and T' = γT. Then:

* |k'T − h'| = 4k'|z|/π, and **ρ := |e(T')| = exp(−π² Re(1/z)/(2k'²))**.
* |η(T)| = |k'T−h'|^{−1/2}|η(T')|. This is the Dedekind η transformation.
* |g_a(T)| = |g_{aγ⁻¹}(T')|, and |g_a| depends only on a mod Z². This is the Siegel-function transformation from Kubert–Lang, *Modular Units*, Ch. 2, K1–K2.
  Here (A,B)γ⁻¹ = (−Ah' − Bk', Bα − Aβ), so **Â := {−Ah' − Bk'}** depends only on (h' mod 8, k' mod 2).

Numerical check (`farey_check_transform.py`): the identity |Θ| = (transformed expression) agrees with direct summation to a log error of
5·10⁻¹⁵, over 16 cusps and 6 theta factors at random z. I did not re-derive the K1/K2 laws; I cite them.

**Lemma 2 (two-sided product bounds) [R].** Let (Â,B̂) ∈ ((1/8)Z/Z)² \ {0}, which is the case here because (A,B) ∉ Z² and γ ∈ SL₂(Z).

 ρ^{B₂(Â)/2} · L(Â,ρ) ≤ |g_{(Â,B̂)}(T')| ≤ ρ^{B₂(Â)/2} · U(Â,ρ),

where:
- U = (1+ρ^Â) ∏_{n≥1}(1+ρ^{n+Â})(1+ρ^{n−Â}).
- L = ℓ ∏_{n≥1}(1−ρ^{n+Â})(1−ρ^{n−Â}), with ℓ = 1−ρ^Â if Â > 0, and ℓ = 2 sin(π/8) if Â = 0 (then B̂ ∈ (1/8)Z \ Z).
- If Â = 0, use 2 in place of 1+ρ^Â in U.

Likewise ρ^{1/24}∏(1−ρⁿ) ≤ |η(T')| ≤ ρ^{1/24}∏(1+ρⁿ). U is increasing in ρ and L is decreasing in ρ. Since Â ∈ (1/8)Z, every exponent n − Â is at least 1/8.

*Proof.* Use |q_z| = ρ^Â and the triangle inequalities |1−w| ∈ [1−|w|, 1+|w|]. ∎

Check V1: there were 3000 random checks of six factors on Farey arcs, with t0 = 0.02 and 0.01. **There were 0 violations.**
The worst slack was −3·10⁻⁷ in log, which is roundoff in the direct sum. These equality cases come from the bound becoming exact as ρ → 0.
The median slack was 0.8 (lower) and 0.6 (upper).

## 3. Exact cusp exponents for all cusps [R]

Multiply the ρ-powers of Lemma 2 by the powers ρ^{±1/24} from η. For a product of factors this gives exp((π²/64)·Re(1/z)·E), with
**E = R(h' mod 8, k' mod 2)/k'²**, where R is an exact rational number (`farey_exponents.py`):

| k' | h' mod 8 | R_TQ | R_{AL,A} = 1/Θ(−1;q⁸) | R_{form B} |
|---|---|---|---|---|
| odd | 0 | 0 | 0 | 0 |
| odd | 1, 7 | −7/4 | 0 | 1 |
| odd | 2, 6 | **1** | 0 | 4 |
| odd | 3, 5 | 1/4 | 0 | 1 |
| odd | 4 | −12 | 0 | 0 |
| even | 1, 7 | −15/4 | 4 | 1 |
| even | 3, 5 | 9/4 | 4 | 1 |

The cases map to k as follows:
- k odd ⇒ h' ≡ 0.
- k ≡ 2 (mod 4) ⇒ h' ≡ 4.
- k ≡ 4 (mod 8) ⇒ h' ∈ {2,6}, with k' = k/4 odd.
- 8 | k ⇒ h' = h is odd.

So class (odd; 2,6) with k' = 1 is exactly k = 4.

**Lemma 3 [R].** Choose the form by k': use form A when k' is odd and form B when k' is even. Then every term of the bound has E ≤ 1/4 at every cusp h/k with k ≠ 4.
The value 1/4 is attained at k = 8 (h = 3, 5) by TQ, and at k = 16 (k' = 2) by form B. At k = 4, TQ has E = 1 and AL_A has E = 0.

This settles the cusp-exponent question exactly for all k, with no k ≤ 400 cutoff. It agrees with `lemma.py`, which found a maximum of 0.25.

## 4. Lerch numerators [R]

Let N(x,z) := Σ_r (−1)^r q^{4r(r−1)} z^r/(1 − q^{8(r−1)}xz). The bound |1−q^e| ≥ 1−e^{−|e|t0} holds for both signs of e, and Σ_r e^{−c(r−δ)²} ≤ 1 + √(π/c). These give:
- |N_A| = |N(−q³,−1)| ≤ e^{t0}(1+√(π/4t0))/(1−e^{−3t0}).
- |N(−q³,−q²)| ≤ e^{t0/4}(1+√(π/4t0))/(1−e^{−3t0}).
- |N(−q³,−q⁶)| ≤ e^{t0/4}(1+√(π/4t0))/(1−e^{−t0}).

Check V6: the sampled maxima are 4 to 70 times below these bounds.

## 5. Circle method geometry [R]

Fix n. Put t0 = π/(16√n) and K = ⌊√(2π/t0)⌋, and use Farey arcs of order K on |q| = e^{−t0}.
On the arc of h/k we have z = t0 − iφ with |φ| ≤ 2π/(k(K+1)). Since K²t0² ≤ 2πt0 and (2π/(K+1))² < 2πt0:

**Lemma 4.**
- On every arc: k|z| < √(4πt0), and Re(1/z)/k'² ≥ Re(1/z)/k² > 1/(4π). Hence ρ ≤ ρ* := e^{−π/8} = 0.6752 uniformly in n.
- On the arcs at 1/4 and 3/4: Re(1/z) > 1/(t0+π/8). Each arc extends at least Y ≥ π/(2(√(2π/t0)+4)) on each side.

Check V2: the minimum of Re(1/z)/k² over all arcs is 0.0842, 0.0821, 0.0810 and 0.0805 for n = 10², 10³, 10⁴ and 6·10⁴. The bound 1/(4π) is 0.0796.

The choice t0 = π/(16√n) balances n t0 against (1/4)π²/(64 t0). This gives e^{(π/8)√n}, the same rate as the error observed numerically.

## 6. Minor arcs (k ≠ 4) [R]

Combine Lemmas 1–4. The constants (π/(4k'|z|))^{1/2} ≤ (π/4t0)^{1/2} and (4k'|z|/π)^{1/2} ≤ (4√(4πt0)/π)^{1/2} also enter. On every minor arc:

 |S0(q)| ≤ C_min(t0) · exp(π²/(256 t0)).

C_min is the maximum over the 12 classes of:
- form A: (π/4t0)^{1/2}·U_TQ(ρ*) + 2|N_A|bound·(4D/π)^{1/2}·L₀(ρ*)⁻¹;
- form B: the analogue with N_B1, N_B2 and L_{1/4}, L_{3/4}.

The product constants at ρ* are as follows (`farey_N0.py` output):
- Form A classes: U_TQ ranges from 2·10¹¹ to 5.3·10¹³, and L₀⁻¹ = 3.6·10⁴.
- Form B classes: L_{1/4}⁻¹ = L_{3/4}⁻¹ = 5.0·10⁴.

The worst class is k' odd with h' ≡ 3, 5 (the cusps 3/8 and 5/8), where C_min ≈ 1.1·10¹⁵ at n = N0.

**These constants are computed by summing explicit convergent products in floating point. [N, routine]**
They are values of explicit monotone functions, not fitted numbers. All the tails are geometric with ratio ≤ 0.675, and they are truncated at 10⁻¹⁸.

Check V3: at n = 400, 1000 and 2000 I sampled 400 random points plus the cusps 3/8, 5/8, 1/16, 3/16 and 1/12.
The maximum of |S0| sits at 5/8, as predicted, and lies **e^{31.6} below** the sup bound. That gap is the price of the crude products at ρ ≈ 0.675.

## 7. Major arcs (h/k = 1/4, 3/4) [R, apart from one numerically fixed phase]

At 1/4 we have k' = 1, h' = 2 and T' = iπ/(4z), so Q' = e^{−π²/(2z)}. The Â values of the three TQ factors are 1/4, 1/2 and 1/4. Lemmas 1–2 and holomorphy give the exact identity

 TQ(i e^{−z}) = C₀ z^{−1/2} e^{−z/16} e^{π²/(64z)} F(z),
 with F = ∏(1−u_j)^{m_j}, |u_j| = ρ^{e_j}, e_j ∈ {1/4, 1/2, 3/4, …}, F → 1, and |C₀| = √π/2.

**The phase of C₀, namely e^{9πi/16}, is fixed numerically.** Check V5 gives |TQ/main − 1| = 4·10⁻¹¹ at z = 0.05, which matches the predicted size ρ^{1/4}. By the Siegel/η multiplier system the phase is a root of unity of bounded order, so a numerical match determines it. I did not derive it symbolically. **[N]**

**Lemma 5 [R].** |TQ − C₀z^{−1/2}e^{−z/16}e^{π²/64z}| ≤ (√π/2) t0^{−1/2} · ρ^{−1/32}(e^{Λ(ρ)} − 1), where Λ(ρ) = Σ_j |m_j|(−log(1−ρ^{e_j})).
The right side is increasing in ρ, because it is a series in ρ^{1/8} with nonnegative coefficients and minimal exponent 1/4 − 1/32 > 0. So it is bounded by its value at ρ = exp(−π²/(2(t0+π/8))).
For AL_A at 1/4 we have Â = 1/2 and the exponent is exactly 0. It is bounded by 2|N_A|(4|z|/π)^{1/2}/(L_{1/2}·η_L).

Check V5: the actual |TQ − main| is below the bound at four test points. The bound is tight to within a factor of 1.3–8.

**Lemma 6 (Bessel and tails) [R].** We have (1/2πi)∫_{t0−i∞}^{t0+i∞} z^{−1/2}e^{Nz+a/z}dz = cosh(2√(aN))/√(πN), with a = π²/64.
This is the I_{−1/2} formula; check V4 confirms it to 3·10⁻⁷. So the two arcs give exactly M(n), minus four tails.
Integrating by parts in y, each tail is at most e^{Nt0}/(2πN) · e^{at0/Y²}(2Y^{−1/2} + (2a/3)Y^{−3/2}). Lemma 4 gives at0/Y² ≤ (√(2π)+4√t0)²/16.

## 8. Assembly [R, given I1 and the numerically fixed phase]

The Farey arcs partition the circle. So, with N = n − 1/16:

 |c(n) − M(n)| ≤ e^{nt0}·[ C_min e^{π²/(256t0)} + B_maj ] + 4(√π/2)·Tail = B(n),

and n t0 + π²/(256 t0) = (π/8)√n. Each piece of C_min(t0)·t0^{3/2} is an increasing function of t0, which gives C_min(n) ≤ C_min(n*)(n/n*)^{3/4}. So

 **|E(n)| ≤ 2.4·10¹¹ · n^{3/4} · e^{(π/8)√n}·(1 + e^{−100})  for n ≥ 2·10⁵,**

and the exact B(n) is used for smaller n. This is θ = 1/2 in the requested normalisation. For n ≥ 2·10⁵ the major-arc and tail pieces are below the minor-arc piece by factors e^{−114} and e^{−131}, and these ratios decrease further.

**N0.** I scanned every n from 1 to 2·10⁵ (`farey_N0.py`). The last n with B(n) ≥ cos(7π/16)cosh(π√N/4)/√N is n = 11230. The log margins are:

| n | 11229 | 11230 | 11231 | 11232 | 11300 | 20000 | 60000 | 2·10⁵ |
|---|---|---|---|---|---|---|---|---|
| log(0.195·main) − log B(n) | −0.0035 | −0.0017 | +0.0001 | +0.0019 | +0.12 | 13.5 | 53.3 | 131.8 |

The margin increases at every n in [11231, 2·10⁵), as checked by the scan. For n > 2·10⁵, the bound in the display above gives
margin ≥ π√N/4 − π√n/8 − (3/4)log(n/n*) − (1/2)log N − const. Its derivative is positive. The value is 347.7 at 10⁶ and 3.9·10⁵ at 10¹², as checked by the sparse evaluation. **[S]**: this is an elementary calculus argument that I have not written out in full.

**So N0 = 11231 (11300 with a margin of 0.1).** With the per-n optimised t0 and K of `farey_bound.py`, the crossing is at about n ≈ 1.05·10⁴. Both are well below 60000.

## 9. Validation summary

| check | what | result |
|---|---|---|
| T | transformation identity, 16 cusps × 6 factors | log error 5e−15 |
| V1 | Lemma 2 + 1 bounds vs actual |Θ|, 3000 points | 0 violations |
| V2 | Farey geometry Re(1/z)/k² > 1/(4π) | holds, min 0.0805 |
| V3 | sup bound vs sampled |S0| on minor arcs | bound/actual = e^{31.6} |
| V4 | Bessel identity | 3e−7 relative |
| V5 | major arc identity and Lemma 5 | ratio error 4e−11; bound holds |
| V6 | Lerch numerator bounds | bound/actual 4 to 70 |
| V7 | actual log|E(n)| vs log B(n) | 36.0 vs 76.2 at N0; 87.7 vs 131.2 at 60000 |

V7 shows that the true error is about e^{40} smaller than B(n). The whole loss is the product constants of Lemma 2 at ρ ≈ 0.675; the exponent itself is sharp.

## 10. Caveats

1. The Form A/B identities (I1) are numerical inputs. A citation or a proof is needed, for example Gordon–McIntosh's Appell–Lerch representations of the eighth-order functions.
2. Two transformation laws are cited, not re-proved: the Siegel-function law |g_a∘γ⁻¹| = |g_{aγ⁻¹}| with mod-Z² invariance, and the η law. Both are standard, and both are confirmed numerically to 5e−15.
3. The phase e^{9πi/16} of C₀ is determined numerically. Only the modulus √π/2 is derived.
4. Constants such as U, L, Λ and C_min are evaluated in double precision from explicit products; there is no interval arithmetic. The margin at n = 11231 is only 1e−4 in log, so for a certified statement use N0 = 11300, which has a margin of 0.12. Any N0 ≤ 60000 suffices.
5. The argument for n > 2·10⁵ uses an elementary monotonicity argument (Section 8), which I give in sketch form; the dense scan only goes up to 2·10⁵.
