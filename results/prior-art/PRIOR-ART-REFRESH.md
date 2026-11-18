# PRIOR-ART-REFRESH (gate-2): one-loop Bach-source back-reaction on Evans-Coleman (k = 1/3), certified

Prior-art scout report, 2026-11-18. Headline gate, strongest tier. This refreshes <LAB>/quantum/PRIOR-ART-QB.md
(gate-1, 2026-11-17). It does not redo gate-1's sources; it re-uses gate-1's downloaded copies in <LAB>/quantum/src/
for new quotes only.

- Query log: <LAB>/quantum/seats/gate-2/QUERIES.log. It holds 90 queries (34 arXiv, 45 INSPIRE, 11 Zenodo), with positive controls on every database. `lab.py searchcheck` printed: "search floor OK: every database with a zero-hit query has a passing positive control".
- New sources downloaded: <LAB>/quantum/seats/gate-2/src/. They are gr-qc_9712084.pdf/.txt and two arXiv abstract files, arxiv-abstracts-1.xml and arxiv-abstracts-2.xml.
- **Newest material covered:**
  - arXiv: every gr-qc submission dated 2026-11-05 to 2026-11-18 was title-screened (373 entries; newest 2611.23200, 2026-11-17). A cross-list window (hep-th, gr-qc, math-ph, astro-ph.CO with collapse / self-similar / critical solution / trace anomaly) returned 37 entries.
  - INSPIRE: citers up to 2611.21356 (2026-11-16).
  - Zenodo: records up to 2026-11-18 (record 24338153).

**ID correction for the brief:** gr-qc/9805008 is not Brady-Ottewill. INSPIRE gives it as recid 470078, "The Late time singularity inside nonspherical black holes". Brady-Ottewill is **gr-qc/9804058**, INSPIRE recid **469529**, 17 citers, and all 17 titles were read.

## In plain words

- **Nothing new since gate-1 takes the claim.**
  - In the last two weeks nobody has posted a quantum correction to the radiation (Evans-Coleman) critical solution, on arXiv, INSPIRE or Zenodo.
  - Since 2026-11-05 the closest group (Tomasevic, Wu and collaborators at UW) has posted only classical large-D work (2611.07007).
- **(a) Which way the quantum push goes for radiation: still OPEN.**
  - Nobody gives a sign for radiation.
  - Wu's newest paper (large-D scalar toy, not radiation) says the quantum push "can favor or oppose trapping".
  - Tomasevic-Wu still list radiation as a job for future simulations.
- **(b) A certified (rigorous-numerics) quantum back-reaction on a true critical solution: OPEN.**
  - All one-loop back-reaction results we found are on toy self-similar spacetimes (Roberts, Garfinkle, large-D families). Their own authors say these are not the true critical solutions.
  - Rigorous work on true critical solutions is all classical.
- **(c) The centre form of the Bach source: PARTIAL.**
  - Two pieces are already published as general principles:
    - quantum stress terms grow like 1/t^4 at fixed x, which follows from self-similarity (Zahn, Brady-Ottewill);
    - at a smooth centre a one-loop stress tensor must be isotropic and carry no radial flux (Wu 2610.20705 writes these conditions out).
  - Not found anywhere: the Bach tensor's actual centre limit on a self-similar fluid, its coefficient (our c0 = -173.10...), or the statement that it is finite there. Earlier semiclassical sources were the opposite: their tensors diverged or were singular at the centre.
- **One old lead closed:** Chiba-Siino (1997, paywalled) is described by two open citers as a **2-dimensional** model. Its "divergence at the origin" is therefore not about the 4D Bach term.

## Aliases searched (leak-proofing)
- **Object:**
  - Evans-Coleman solution; radiation-fluid critical solution; k = 1/3 or p = rho/3 CSS perfect-fluid critical solution;
  - Koike-Hara-Adachi solution; Type II CSS fluid threshold.
- **Operation:**
  - semiclassical / one-loop / hbar back-reaction (backreaction); vacuum polarization; RSET;
  - trace / conformal / Weyl anomaly; Weyl-squared (c-) anomaly term;
  - anomaly-induced effective action; Riegert action; Page conformal transformation law;
  - Zahn's state-independent V tensor; "quantum growing mode"; quadratic-gravity / Weyl^2 EFT correction (same Bach tensor in the field equations).
- **Output:**
  - threshold shift; mass gap; trapping / apparent horizon; compactness 2m/r; compaction function (PBH community); quantum Lyapunov exponent kappa_q = 2.
- **Centre item:**
  - Bach tensor at a regular centre / regular center / origin / r = 0; centre smoothness / isotropy of a stress tensor; homothetic Killing vector.

## Step 0: our own files
- **Runs** (all with `--exclude <LAB>/quantum`, so the current project's notes are dropped):
  - `lab.py selfcite "Evans-Coleman|radiation fluid critical|critical solution" "semiclassical|one-loop|trace anomaly|Bach" "back-reaction"`: 5 hits in 3 files.
  - `"Bach tensor|Bach source|conformal anomaly" "self-similar|homothetic|regular centre|regular center" "centre expansion"`: no hits.
  - `"Brady-Ottewill|Brady and Ottewill|Tomasevic|Zahn" "critical collapse" "quantum correction"`: 3 hits in 2 files.
- **What the hits are:**
  - Our lab copy of the review 2507.07636 (<LAB>/seats/src/2507.07636.txt l.1520-1563). Already quoted in gate-1.
  - <RECON>/seats/BIGGER-IDEA.md l.58-59 and l.138-139, and OUR-EDGE.md l.171. These are our own earlier recon notes, with no outside claim.
  - None of our files claims or cites the result under another name.
- **Reference sweep:**
  - In the review 2507.07636, section "D. Quantum effects", the semiclassical threshold list is: Chiba-Siino; Ayal-Piran; Strominger-Thorlacius; Kiem hep-th/9407100; Zhou-Mueller-Kirsten-Yang PRD 51 R314; Peleg-Bose-Parker; Bose-Parker-Peleg; Ashtekar-Pretorius-Ramazanoglu 1012.0077.
  - The review places all of them in "1+1 dimensional models which serve as toy models for spherical symmetry" (l.1530-1538), except Brady-Ottewill ("3+1-dimensional but perturbative", l.1541).
  - The 2D works are in the class gate-1 already rated no overlap.
- **Suggested alias lines for <TOOLS>/lab/tools/selfcite-aliases.txt** (the lead adds them):
  - `bach tensor | bach source | weyl-squared anomaly | weyl anomaly | c-anomaly | conformal anomaly | trace anomaly`
  - `regular centre | regular center | centre smoothness | center smoothness | spherical center | origin r = 0`
  - `compactness | compaction function | 2m/r | trapping | apparent horizon`

## Main table

Labels: OPENED = full text read (copy in src/ of gate-1 or gate-2). ABSTRACT = abstract only. SNIPPET = search snippet only.

| item | verdict | best source | label | quoted sentence | locator |
|---|---|---|---|---|---|
| (a) sign of the one-loop push on a true critical solution for radiation | **OPEN** | Wu, arXiv:2610.20705v1 (large-D DSS/CSS scalar, not radiation) | OPENED (gate-1 copy) | "we use the radial constraint to relate the effective energy enclosed by each sphere to the change in the outgoing expansion, showing how the backreaction can favor or oppose trapping." | <LAB>/quantum/src/2610.20705.plain.txt l.1760 |
| (a) cont.: radiation named, sign left to simulation | OPEN (support) | Tomasevic-Wu, arXiv:2509.03587v3 | OPENED (gate-1 copy) | "Quantum vacuum polarization, which is similarly captured for a radiation fluid by the conformal anomaly, adds a second, quantum growing mode ... a quantitative determination of both [Delta p] and M_gap for radiation-fluid collapse requires dedicated semiclassical numerical simulations." | src/tw.txt l.406, l.417 |
| (a) cont.: newest papers by the same group are classical or off-topic | no overlap | 2611.07007 (classical large-D DSS); 2609.26875 (charged shell, third law); 2605.03015 (Nariai nucleation) | OPENED (2611.07007, gate-1 copy) / ABSTRACT (others) | 2609.26875: "A near-extremal black hole with an AdS2 x X throat is first prepared and then driven to extremality in finite time by a charged shell." | gate-2 src/arxiv-abstracts-1.xml |
| (b) certified (rigorous-numerics) back-reaction on any true critical solution | **OPEN** | Tomasevic-Wu, arXiv:2509.03587v3 (closest: closed form on toy CSS) | OPENED (gate-1 copy) | "we now apply it to two analytically tractable CSS Einstein-scalar spacetimes: the Garfinkle solution in D = 2 + 1 [53] and the Roberts solution in D = 3+1 [54-56]. Neither is the exact numerical critical solution of its corresponding collapse problem" | src/tw.txt l.80 |
| (b) cont.: rigour on true critical solutions is classical only | context | Tomasevic-Wu, arXiv:2509.03584v3 | OPENED (gate-1 copy) | "an explicit analytic form of the true critical solution remains unknown, except attempts for analytic approximations [281] and existence proofs [22]." (the existence proof is classical; [22] not opened here, UNVERIFIED as to content) | src/tw2.txt l.2698 |
| (b) cont.: new numerical RSET-on-compact-object method (not rigorous, not critical) | neighbour | Saffin & Xie, arXiv:2601.05129 (PRD 113, 125025), cites Brady-Ottewill | ABSTRACT (INSPIRE record 3098875) | "we compute the quantum scalar fields and stress tensor in boson star spacetimes within the framework of semiclassical gravity. Divergences are regularized using Pauli-Villars fields, and accurate numerical results are obtained through spectral methods." | INSPIRE literature/3098875 abstract |
| (b) cont.: 4D anomaly-induced collapse with back-reaction (null shell, not critical) | neighbour / named obstruction | Lowe & Thorlacius, arXiv:2605.00780v1 | ABSTRACT | "This behavior arises from a nonlinear instability of the higher-derivative semiclassical equations and is generic in models with anomaly-induced quantum corrections." | gate-2 src/arxiv-abstracts-1.xml |
| (c) exact centre form of the Bach source on a regular self-similar centre | **PARTIAL** | Wu, arXiv:2610.20705v1 (general centre-smoothness conditions for one-loop stress tensors) | OPENED (gate-1 copy) | "Across a regular center, the outward radial direction reverses, so T^0^1(D) changes sign and must vanish at least linearly with r. By contrast, T^1^1(D) - T^i^i(D) is unchanged by this reversal and vanishes at the center, where the radial and angular directions become equivalent." | src/2610.20705.plain.txt l.2226 |
| (c) cont.: 1/t^4 scaling of the state-independent term is structural | PARTIAL (covered part) | Zahn, arXiv:2509.10360v3 eq. 5-6; Brady-Ottewill gr-qc/9804058 eq. 3.8 | OPENED (gate-1 copies) | Zahn: "under the assumption (4) on the state, the expectation value of the stress tensor behaves under "self-similarity time" shift as T(u + lambda, y) = e^{2 lambda} T(u, y) + lambda [factor lost in extraction; 1/(4 pi^2) per draft l.181] V(u, y)"; "For conformal coupling (xi = 1/6), V reduces to (a multiple of) the Bach tensor" | src/zahn.txt l.135-139, l.347; src/bo.txt l.144-166 |
| (c) cont.: what is NOT covered | OPEN part | (none found) | - | No source evaluates B on any regular-centre CSS fluid, or gives its centre limit, isotropy, coefficient or finiteness. BO and Zahn never evaluate B at y = 0 (grep of src/bo.txt and src/zahn.txt for origin/centre/center: only coordinate remarks, BO l.56). | QUERIES.log: inspire `fulltext:"Bach tensor" and fulltext:"regular centre"` = 1 (static Einstein-Weyl, 2512.14148); `... "regular center" and "spherically symmetric"` = 3 (static/regular-BH reviews); arXiv `abs:"Bach tensor" AND (origin OR center OR centre OR regular)` = 3 (Riemannian geometry, CCC matching) |
| (c) cont.: prior semiclassical centre results are obstructions, not forms | context | Wu 2610.20705 | OPENED | "The simplest effective action induced by the trace anomaly realizes J = 1 ... However, its stress tensor is singular at the spherical center, where the classical DSS geometry is smooth at every finite T." | src/2610.20705.plain.txt l.199 |
| (c) cont.: Chiba-Siino centre divergence is a 2D-model statement (lead closed) | context | Gundlach review gr-qc/9712084v3 (early citer's own words) | OPENED (gate-2 src) | "Chiba and Siino [22] write down their own 2-dimensional toy model, and add a quantum stress tensor that is determined by the trace anomaly and stress-energy conservation. They note that the quantum stress tensor diverges at r = 0." | gate-2 src/gr-qc_9712084.txt l.1980-1983 |
| (c) cont.: same, second open citer | context | Tomasevic-Wu arXiv:2509.03584v3 | OPENED (gate-1 copy) | "The study by Chiba and Siino [68] is the closest one to ours since the same two-dimensional model and one-loop effective action are employed ... It is argued that regularity at the center of critical spacetime demands the disappearance of echoes" | src/tw2.txt l.869 |

## Other sources checked this round (no overlap)
- **INSPIRE citers of Brady-Ottewill (recid 469529, 17 hits).** The new titles relative to gate-1's discussion:
  - Saffin-Xie 2601.05129 (boson stars, table above).
  - Anderson-Mull gr-qc/9806069 (ABSTRACT): "a gap exists in the allowed sizes of all zero temperature static spherically symmetric black holes in semiclassical gravity when only conformally invariant fields are present". Static, not collapse.
  - Tomasevic PhD thesis 2021, INSPIRE 1971627 (record only; it predates the 2025 work).
  - Reviews, and 0808.3024 (topology-change toy).
- **INSPIRE citers of TW 2509.03587 (recid 2966624, 9) and TW 2509.03584 (recid 2966695, 11).** New since gate-1's discussion: 2609.26875 and 2605.03015 (both ABSTRACT, off-topic), plus 2602.10185 and 2601.14358 (classical large-D DSS).
- **Other citer trails:**
  - Zahn 2509.10360 (recid 2969255, 5) and Wu 2606.17135 (2): nothing new.
  - Wu 2610.20705 (0) and 2611.07007 (0): no citers.
  - Jiang-Jiang 2607.04386 (recid 3177545, 2 citers: 2611.21356, Higgs vacuum inside Schwarzschild; 2608.03538, black-hole/white-hole transitions): not collapse.
- **Chiba-Siino citers (recid 444860, 50 hits).** Mostly 2D dilaton / dimensional-reduction anomaly papers 1998-2004 and the 2025-26 TW/Wu set. None evaluates a 4D fluid critical solution.
- **Quadratic-gravity collapse (the Weyl^2 term gives the same Bach tensor in the field equations):**
  - Held-style order-reduced collapse 2307.15163 (ABSTRACT): "the gravitational collapse of a massless scalar field within Quadratic Gravity treated as a dimension-four operator Effective Field Theory". Not critical, not fluid.
  - Starobinsky R^2 critical collapse 2212.14805 (ABSTRACT): R^2 only, scalar field.
  - The arXiv query `abs:"critical collapse" AND (EFT OR higher curvature OR higher derivative)` returned 0, and `abs:"quadratic gravity" AND self-similar AND collapse` returned 0.
- **Semiclassical stars:** Arrechea-Barcelo-Carballo-Rubio-Garay 2110.15680 (ABSTRACT). It uses a Polyakov-approximation RSET on static constant-density stars and finds the "critical configuration already displays a pressure that diverges exactly at the center". Static, and 2D-derived.
- **Interval / computer-assisted with GR self-similar:**
  - arXiv: 2 hits, both PDE blow-up outside GR.
  - INSPIRE fulltext: the review and 2112.15241 (shock waves).
  - `fulltext:"computer-assisted" and "Evans-Coleman"`: the review and gr-qc/9607010 (classical RG).
- **Zenodo:**
  - Controls: "critical collapse" (66) and "Evans-Coleman" (1, our own 24327088).
  - The semiclassical, Bach, trace-anomaly and "quantum critical collapse" queries returned 0 or only off-topic records. No Brady-Ottewill or Tomasevic-Wu related deposits.
- **Web search (extended)** "quantum back-reaction radiation fluid critical collapse self-similar trace anomaly Bach tensor 2026" (SNIPPET only). It returned only the known set: 2610.04960, 2509.03584 (OSTI copy of JHEP04(2026)151), 2509.03587 (PRD DOI 10.1103/7rkc-hmxd), 2610.20705, 2606.17135, 2509.10360, gr-qc/9712084.

## Normalisations
- **Sign convention for the centre coefficient.**
  - Our draft gives B_(nn)(r=0,t) = c0/t^4 with c0 = [-173.1004443922995663 +/- 2.99e-13] (paper.md l.123).
  - Wu's conditions (l.2226) are written for a D-dimensional s-wave tensor in mixed components, T^1_1 - T^i_i -> 0, i.e. radial = angular pressure at the centre.
  - Our B_(rr) = B_(theta theta) = c0/(3 t^4) is the same isotropy statement in a frame basis.
  - Wu prints no number for any Bach-type tensor, so there are no digits to compare. The overlap is structural only.
- **Rates.** kappa_q = 2 (BO, Zahn), and "D - 2" (Wu), which equals 2 at D = 4. Our kappa* = 2.81 for k = 1/3 is unchanged from gate-1.
- **Trapping observable.**
  - Wu 2610.20705 tracks "the change in the outgoing expansion" via the radial constraint, which in spherical symmetry is equivalent to tracking 2m/r.
  - Our Psi2 is the first-order change of max_U 2m/r.
  - Same observable family; different background (large-D scalar vs k = 1/3 fluid) and different source (2D-reduced vs 4D local Bach).

## Gap that is genuinely open (one sentence)
At this search floor (arXiv to 2611.23200, INSPIRE to 2611.21356, Zenodo to 2026-11-18), no source evaluates the 4D state-independent one-loop (Bach / Weyl^2-anomaly) source on the Evans-Coleman k = 1/3 background, gives its finite centre limit, or derives a sign for the first-order change in peak compactness 2m/r of any true critical solution, and no semiclassical quantity has been certified with rigorous numerics on any true critical solution.

## Wording notes for the lead (from this refresh; the lead decides)
- **Chiba-Siino.** The draft (paper.md l.424-427) says "We cannot rule on their centre-regularity argument."
  - Two open citers describe it as a 2-dimensional model: Gundlach gr-qc/9712084 l.1980-1983, and TW 2509.03584 l.869 ("the same two-dimensional model and one-loop effective action").
  - The draft can cite those sentences and keep the paper itself labelled "abstract only".
- **Centre item (Theorem 1, item 3).** Do not present isotropy, zero flux or the 1/t^4 law as new. Wu 2610.20705 l.2226 states the smoothness and isotropy conditions in general, and the scaling is Zahn eq. 5-6 / BO eq. 3.8. The new content is the evaluation on Evans-Coleman: finiteness of the local Bach tensor and the certified c0.
- **"Every field" scope.** The brief says spin <= 1 conformal fields. Zahn's theorem covers a scalar of any xi; BO's law covers conformal fields. This is unchanged from gate-1. The state-dependent e^{2u}T(y) term (Zahn eq. 6) is what Psi1 carries. Psi1 is not claimed.

## Most dangerous neighbour
- **Chih-Hung Wu (UW Seattle), with Tomasevic and Turiaci.**
  - Wu 2610.20705 already has the three ingredients:
    - the centre-smoothness conditions;
    - a covariant one-loop action smooth at a regular centre;
    - the radial-constraint test of "favor or oppose trapping".
  - TW name radiation-fluid collapse as "the important next step".
  - Moving from large-D scalar to a 4D k = 1/3 fluid with the local anomaly term is a short step for them. They plan to take it by simulation.
- **Saffin (Nottingham).** Saffin & Xie 2601.05129 now compute full RSETs (Pauli-Villars, spectral) on strongly curved compact spacetimes, citing Brady-Ottewill. Berczi-Saffin-Zhou 2111.11400 already ran semiclassical Choptuik collapse.
- **Lowe-Thorlacius 2605.00780.** 4D anomaly-induced collapse with back-reaction, null shell.
- No race framing is intended: the overlap is noted so that we do not duplicate their work.

## Methods worth borrowing
1. **Wu's centre-smoothness criteria** (2610.20705 l.2225-2237): flux vanishes like r, anisotropy like r^2. They give an independent positive control on our certified centre values. Check that B_(nr) = O(r) and B_(rr) - B_(theta theta) = O(r^2) from the centre series.
2. **Wu's radial-constraint reading of trapping** (l.1760): enclosed effective energy gives the change in outgoing expansion. It is a second route to sgn Psi2 that avoids the maximiser x_V.
3. **Saffin-Xie Pauli-Villars + spectral RSET** (2601.05129, ABSTRACT only). A possible later route to the state-dependent Psi1 term. Not needed for Psi2.

## Gate lines

The brief's short form is below. Note that `--file seats/gate-2/PRIOR-ART-REFRESH.md` is relative to the project dir <LAB>/quantum, not <LAB>. The full-path form follows it. Not run by this seat.

```
lab.py gate <LAB> "(a) sign of the one-loop quantum push on a true critical solution for radiation (k=1/3)" OPEN --file seats/gate-2/PRIOR-ART-REFRESH.md --source "arXiv:2610.20705"
lab.py gate <LAB> "(b) certified (rigorous-numerics) semiclassical back-reaction on any true critical solution" OPEN --file seats/gate-2/PRIOR-ART-REFRESH.md --source "arXiv:2509.03587"
lab.py gate <LAB> "(c) exact centre form of the Bach source on a regular self-similar centre (smoothness/isotropy and 1/t^4 scaling known; Bach centre limit and coefficient not found)" PARTIAL --file seats/gate-2/PRIOR-ART-REFRESH.md --source "arXiv:2610.20705"
```

```
python "<MEMCAP>" 1 python "<TOOLS>/lab/tools/lab.py" gate <LAB>/quantum "(a) sign of the one-loop quantum push on a true critical solution for radiation (k=1/3)" OPEN --file <LAB>/quantum/seats/gate-2/PRIOR-ART-REFRESH.md --source "arXiv:2610.20705" --queries <LAB>/quantum/seats/gate-2/QUERIES.log
python "<MEMCAP>" 1 python "<TOOLS>/lab/tools/lab.py" gate <LAB>/quantum "(b) certified (rigorous-numerics) semiclassical back-reaction on any true critical solution" OPEN --file <LAB>/quantum/seats/gate-2/PRIOR-ART-REFRESH.md --source "arXiv:2509.03587" --queries <LAB>/quantum/seats/gate-2/QUERIES.log
python "<MEMCAP>" 1 python "<TOOLS>/lab/tools/lab.py" gate <LAB>/quantum "(c) exact centre form of the Bach source on a regular self-similar centre (smoothness/isotropy and 1/t^4 scaling known; Bach centre limit and coefficient not found)" PARTIAL --file <LAB>/quantum/seats/gate-2/PRIOR-ART-REFRESH.md --source "arXiv:2610.20705" --queries <LAB>/quantum/seats/gate-2/QUERIES.log
```

**Ledger note on (b).** Gate-1 recorded PARTIAL for the wider object "certified or rigorous semiclassical back-reaction on a critical collapse solution". That wording counted TW's closed-form toy results. The object (b) here is narrower: "on any true critical solution". For that object nothing was found, so it is OPEN. The two rows are consistent.
