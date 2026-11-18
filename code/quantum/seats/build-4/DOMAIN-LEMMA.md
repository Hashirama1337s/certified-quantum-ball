# DOMAIN-LEMMA (build-4, QB2b item 3): Ψ2 uses only data on x ≤ 0

**Status: PROVED (to referee).** The argument uses only PROOF-QB2-LEMMAS-v2 Lemmas (b)-(d), the four-row system
(n1_rows.log, which equals K5), and the Frobenius exponent bound ρ(2) ≤ -1.9085 (CERTIFIED, PROOF-B4-v4 l.81 and
REFEREE-QB2 1(b)). Control C1 is re-run at x = -6, -4.6, -2, -0.25 and -0.01 (controls_qb2b.json).

## Statement
Let 𝔇₋ = (background N, ω, V on (-∞, 0]; the source b_ab on (-∞, 0]). Then
Ψ2 = Ā_2(x_V)/A(x_V) is a function of 𝔇₋ alone. Changing the background or the source on x > 0 leaves Ψ2 unchanged,
provided 𝔇₋ is unchanged and the data stay analytic at x = 0.

## Proof
1. **The system on x < 0.** h2 solves Q Y' = (P0 + 2P1)Y + S on (-∞, 0) with the identity I(Y) = σ_I. Here
   Q = diag(1, 1, M_x), and P0, P1, Q, S, σ_I are explicit in (N, ω, V, ω', V') and b_ab at the same x. These are the
   four HKA rows of n1_rows.log, which equal K5 exactly. On x < 0, every coefficient at x is a function of 𝔇₋ at x.
2. **The sonic condition is one-sided.** Near x = 0⁻ the system is regular-singular (Lemma (b), PROOF-B4 Fact 1.3).
   - Every solution on (-δ, 0) has the form y = y_an(x) + c·|x|^{ρ(2)}·g(x), with y_an and g analytic at 0. The real
     exponent satisfies ρ(2) ≤ -1.9085.
   - Hence the singular branch is unbounded as x → 0⁻ when c ≠ 0.
   - So "analytic at 0" is equivalent to "bounded as x → 0⁻", which is a property of the solution restricted to x < 0.
   - The forced family S_2^F (Lemma (b)) is therefore the set of solutions on (-δ, 0) bounded at 0⁻ with J ≡ 0. Its
     Taylor data at 0 are the one-sided limits of the coefficients, fixed by 𝔇₋ by continuity.
   - The gauge N̄_p(0) = 0 is a one-sided limit as well.
3. **The centre condition** (Lemma (c)) concerns x → -∞, so it uses 𝔇₋ only.
4. **Matching** (Lemma (d)) happens at x_m = -4.6 < 0. The 4×4 system [s^HKA, ρ1, ρ2, e_V] (b, c, d, e) = y_c^p - y_s^p
   uses values at x_m of solutions constructed in steps 2 and 3. Its determinant D_HKA(2) ≠ 0 (pq7b and N6 here), so b
   is a function of 𝔇₋.
5. **Readout.** x_V is the unique zero of V in U = [-0.30, -0.20] ⊂ (-4.6, 0) (N7). Also
   A(x_V) = Acons(N, ω, V)(x_V), and Ā_2(x_V) is the first component of y_s^p(x_V) + b s^HKA(x_V). Both are values on
   the sonic leg (x_m, 0). ∎

## Consequences
- Item 1 of the lead's levers (domain cut) is valid for Ψ2, N6, N7 and N8. No outward (x > 0) background, source or
  transport enters them.
- The analytic continuation of h2 to x > 0 is determined by Ψ2's data, not the reverse. Control C1 at x = 0.5 and 3
  would test an outward instrument that Ψ2 does not use. It is replaced by the near-sonic point x = -0.01, evaluated
  on the forced sonic series, plus the earlier points.
- The far-side items (T1, the h2 profile on (0, 12], and N2 on (0, 12]) remain open as stated in DESIGN-QB2-v2
  sec. 5. They are not needed for Ψ2.
