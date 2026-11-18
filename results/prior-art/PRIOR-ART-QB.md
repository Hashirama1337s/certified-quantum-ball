# PRIOR-ART-QB: the quantum (one-hbar) correction to the Evans-Coleman critical solution

Prior-art scout report, 2026-11-17. Project: <LAB>/quantum (SCOPE.md). Query log: <LAB>/quantum/QUERIES-QB.log
(`lab.py searchcheck` result: "search floor OK: every database with a zero-hit query has a passing positive control").
Downloaded sources (PDF + extracted text) are in <LAB>/quantum/src/. Line locators below refer to those text files.
Newest material covered: arXiv to 2611.10320 (submitted 2026-11-06); INSPIRE citers to 2611.07007 (2026-11-04); Zenodo to 2026-11-18.

## In plain words

- **Shape of the quantum correction: known since 1998.**
  - Brady and Ottewill worked out what the quantum push on any continuously self-similar collapse looks like, perfect fluids included. It grows like e^{2u}.
  - So the quantum growth rate is 2. The classical rate for radiation is 2.81.
  - They never put numbers on the actual radiation (Evans-Coleman) solution. They assumed the quantum push opposes black-hole formation, and they wrote that this assumption "can not be addressed without a complete calculation".
  - For radiation, where kappa = 2.81 > 2, they wrote "we can say less".
- **Toy models: done in 2025-26.**
  - Tomasevic and Wu, then Wu alone, solved the one-loop problem in closed form on toy self-similar scalar-field spacetimes (Garfinkle, Roberts).
  - In those models the quantum push lowers the threshold. It favours black holes, the opposite of Brady-Ottewill's assumption.
  - They name radiation-fluid collapse as "the important next step", to be done by "dedicated semiclassical numerical simulations".
- **Our target: not found.** No source we found computes the quantum correction on the true Evans-Coleman background. None gives its sign for radiation, and none certifies any semiclassical number on any true critical solution.
- **The sign has no agreed answer.**
  - Brady-Ottewill assume the push opposes black-hole formation.
  - Ayal-Piran find evaporation acts like a smaller initial amplitude, which also opposes it.
  - Tomasevic-Wu find the threshold lowered, which favours it.
  - Wu's 2610.20705 says the direct quantum contribution "can favor or oppose trapping".

## Aliases searched (leak-proofing)
- **The object:** Evans-Coleman solution; radiation-fluid critical solution; p = rho/3 (k = 1/3) CSS perfect-fluid critical solution; Koike-Hara-Adachi solution.
- **The operation:**
  - semiclassical back-reaction; vacuum polarization;
  - renormalized stress(-energy) tensor (RSET);
  - trace / conformal anomaly; anomaly-induced effective action; Polyakov (s-wave, 2D dilaton) approximation; Page approximation (conformal transformation of the RSET);
  - one-loop correction; hbar correction.
- **The output:**
  - quantum growing mode; quantum Lyapunov exponent (kappa_q = 2 = D - 2);
  - threshold shift; quantum mass gap; Type II to Type I transition.
- **The geometry:** homothetic Killing vector; continuous self-similarity (CSS); naked singularity / Cauchy horizon of self-similar collapse.

## Step 0: our own files
- **Runs:**
  - `lab.py selfcite "Evans-Coleman semiclassical" "hbar correction" "critical collapse"` returned no hits (2101 files scanned).
  - The alias runs (semiclassical|quantum back-reaction|Polyakov|trace anomaly x critical collapse|Choptuik|self-similar|Evans-Coleman) gave 5 hits in 3 files.
  - The author-name run gave 38 hits in 19 files. Its relevant hits are <RECON>/seats/BIGGER-IDEA.md l.57-83 and l.123-138, and <TOOLS>/seats/strategy-1/STRATEGY.md l.92-146.
- **What our files hold:** both files are our own earlier gates of this lead, both rated OPEN for radiation. Neither claims the result. They name the same sources re-opened below: Brady-Ottewill, Tomasevic-Wu, Zahn, Wu 2606.17135 and 2610.20705, and Lai-Wen 2610.04960 (title only there).
- **The review:** our lab copy of the critical-collapse review, 2507.07636 (<LAB>/seats/src/2507.07636.txt), has a section "D. Quantum effects" (l.1520-1549). It is quoted in the table.
- **Suggested alias lines for selfcite-aliases.txt** (for the lead to add):
  - `semiclassical | quantum back-reaction | quantum backreaction | vacuum polarization | renormalized stress tensor | renormalised stress-energy tensor | rset | trace anomaly | conformal anomaly | polyakov | anomaly-induced | one-loop | hbar correction`
  - `quantum growing mode | quantum lyapunov exponent | quantum mass gap | threshold shift`
  - `critical collapse | critical phenomena in gravitational collapse | threshold of black hole formation | choptuik scaling | type ii critical`

## Main table

Labels: OPENED = full text read (copy in src/). ABSTRACT = abstract only. SNIPPET = title or search snippet only.

| item | verdict | best source | label | quoted sentence | locator |
|---|---|---|---|---|---|
| 1. Any O(hbar) correction to Evans-Coleman / radiation-fluid critical solution | **PARTIAL** | Brady & Ottewill, gr-qc/9804058 (PRD 58, 024006) | OPENED | "We focus attention on models of gravitational collapse in which black-hole formation turns on at infinitesimal mass, and the critical solution exhibits continuous self-similarity: perfect fluid collapse and a class of Brans-Dicke models belong in this category." and "kappa_q = 2 is easily determined by computing the Einstein tensor to linear order" (Greek letters lost in extraction) | src/bo.txt l.32, l.230 |
| 1 (cont.): what is NOT covered | PARTIAL | same | OPENED | "(ii) When kappa_c >= kappa_q we can say less about the critical point." and "This conclusion relies heavily on the assumption that quantum matter tends to oppose black-hole formation. The validity of this assumption can not be addressed without a complete calculation of the renormalized stress-energy tensor in the dynamical spacetimes of near critical collapse." | src/bo.txt l.267, l.278 |
| 1 (cont.): radiation named, not computed | PARTIAL | Tomasevic & Wu, arXiv:2509.03587v3 | OPENED | "Quantum vacuum polarization, which is similarly captured for a radiation fluid by the conformal anomaly, adds a second, quantum growing mode" ... "a quantitative determination of both [Delta p] and M_gap for radiation-fluid collapse requires dedicated semiclassical numerical simulations ... We view this as the important next step" | src/tw.txt l.406, l.417 |
| 1 (cont.): general structure, any coupling | PARTIAL | Zahn, arXiv:2509.10360v3 | OPENED | "<T>(u,y) = e^{2u} T(y) + (u/4pi) V(y)" (eq. 6); "Essentially the same result on the stress tensor was obtained by Brady and Ottewill [22] under the additional assumption of conformal coupling"; "It is striking that for large u the leading coefficient of e^{2u} is state-independent and of geometric nature." No perfect-fluid background is evaluated: his three examples are a self-similar Minkowski patch, Hayward and critical Roberts (src/zahn.txt l.358, l.382, l.522). | src/zahn.txt l.136-139, l.153 |
| 1 (cont.): review's summary of BO | context | Gundlach-Hilditch-Martin-Garcia review 2507.07636 | OPENED (lab copy) | "In Brady and Ottewill [230] a 3+1-dimensional but perturbative approach is taken. The quantum effects then give rise to an additional unstable mode with lambda = 2. If this is larger than the classical positive Lyapunov exponent lambda_0, it will become the dominant perturbation" | <LAB>/seats/src/2507.07636.txt l.1541-1549 |
| 2. Certified or rigorous semiclassical back-reaction on any critical collapse solution | **PARTIAL** | Tomasevic & Wu, arXiv:2509.03584v3 (JHEP) + 2509.03587v3 | OPENED | Analytic one-loop back-reaction on exact toy CSS spacetimes: "we now apply it to two analytically tractable CSS Einstein-scalar spacetimes: the Garfinkle solution in D = 2 + 1 [53] and the Roberts solution in D = 3+1 [54-56]. Neither is the exact numerical critical solution of its corresponding collapse problem". Their horizon tracing is numerical ("numerical horizon tracing", Wu 2606.17135 abstract). | src/tw.txt l.80 |
| 2 (cont.): nothing certified on a true critical solution | PARTIAL | TW 2509.03584 App. D | OPENED | "to our knowledge, a complete evaluation of T_{mu nu} under physically specified regularity/flux conditions across the relevant domains has not been carried out yet." | src/tw2.txt l.6310 |
| 2 (cont.): rigidity theorem, different background | neighbour | Jiang & Jiang, arXiv:2607.04386v2 | ABSTRACT | "Uniform estimates of the complete mode sums then establish a rigidity theorem ... at conformal coupling, all leading projections are positive." (Schwarzschild interior, not collapse) | arXiv abstract |
| 3. Sign of the quantum push for radiation | **OPEN** (sign uncomputed; literature disagrees on scalar models) | Wu, arXiv:2610.20705v1 | OPENED | "Our preliminary analysis shows that the direct quantum contribution can favor or oppose trapping." (DSS large D) | src/2610.20705.plain.txt l.200 |
| 3 (cont.): sign rule exists, applied only to scalar toys | method | TW 2509.03584 Sec. 6 | OPENED | "There are two equivalent ways to determine the sign of the threshold shift" ... "sgn(Delta p) = sgn(f_q/f_c)" ... "The relative ordering of kappa_c and kappa_q does not affect the sign of Delta p." | src/tw2.txt l.4196-4218, l.4270 |
| 3 (cont.): scalar toy result = favours BH | context | TW 2509.03587v3 | OPENED | "quantum corrections lower the critical point" | src/tw.txt l.275-278 |
| 3 (cont.): opposite reading, 1997 | context | Ayal & Piran, gr-qc/9704027 | OPENED | "In the intermediate regime, evaporation is equivalent to a decrease of the initial amplitude." | src/gr-qc_9704027.plain.txt abstract (l.6) |
| 3 (cont.): Wu restates radiation as open | context | Wu, arXiv:2606.17135v2 | OPENED | "A realistic treatment requires simulations of radiation-fluid collapse with quantum corrections" | src/2606.17135.plain.txt l.3103 |
| Obstruction: exact CSS broken by hbar | context | Chiba & Siino, Mod. Phys. Lett. A 12 (1997) 709, DOI 10.1142/S0217732397000741 (no arXiv; paywalled, not fetched) | ABSTRACT (INSPIRE 444860) | "It is found that the demand of regularity at the origin implies the disappearance of the echo. It is also found that semiclassical equations of motion do not admit continuously self-similar solutions." (scalar field, trace-anomaly method) | INSPIRE abstract |

## Other sources checked (no overlap with the claim)
- **2D dilaton gravity, not 4D fluid:**
  - Strominger & Thorlacius hep-th/9312017 (OPENED abstract page): "In certain two-dimensional models, collapsing matter forms a black hole if and only if the incoming energy flux exceeds the Hawking radiation rate".
  - Bose, Parker & Peleg hep-th/9606152 (OPENED abstract page).
  - Peleg, Bose & Parker gr-qc/9608040 (OPENED abstract page): "quantum effects produce a mass gap not present classically".
- **Scalar field, numerical quantum/semiclassical:**
  - Berczi-Saffin-Zhou 2111.11400 (OPENED abstract page): "Choptuik scaling in critical collapse survives in the semiclassical simulation".
  - Benitez et al. 2002.04044 (loop-quantum, SNIPPET title).
  - Xin-Zhang 2601.18191 (loop-quantum massive scalar, SNIPPET title).
- **2D CSS models:** Moitra 2211.01394 (OPENED abstract page) treats 2D CSS dilaton gravity "including a certain model of semi-classical gravity". No 4D fluid.
- **Large-D DSS, Polyakov:** Lai & Wen 2610.04960 (OPENED abstract), Einstein-Klein-Gordon: "the minimally coupled Polyakov model does not yield a generic smooth, exactly self-similar completion of the classical center."
- **Naked-singularity emission in self-similar collapse (Cauchy horizon, not the growing mode):**
  - Barve-Singh-Vaz-Witten gr-qc/9805095 (OPENED abstract): "the outgoing radiation flux diverges on the Cauchy horizon" (2D, self-similar dust).
  - Miyamoto-Harada gr-qc/0312080 (OPENED abstract + Sec. VII): power "proportional to the inverse square of the remaining time to the Cauchy horizon". Applied to Larson-Penston (GRLP) perfect fluid, not to Evans-Coleman.
- **Frolov-Thorlacius-Giddings astro-ph/9412046** ("Quantum aspects of gravity", GR14 workshop report, OPENED). Its full-text hit on Evans-Coleman plus trace anomaly is a classical review paragraph ("Evans & Coleman examined the spherical collapse of a relativistic perfect fluid (satisfying p = rho/3)"). No semiclassical fluid computation.
- **Zenodo:** no semiclassical critical-collapse deposit. The relevant hits are classical:
  - our own record 24327088;
  - Choptuik computer-assisted spectral count 23270506 and data 23222363;
  - "The exterior of the Choptuik spacetime and a characteristic completion" 24226610.

## Which stress-tensor model the literature uses on perfect-fluid / self-similar backgrounds, and named obstructions
1. **4D conformal fields + Page conformal transformation law (Brady-Ottewill 1998).**
   - BO, eq. 3.9: <T> = e^{4u}<T~>(y) + 8 beta u e^{4u} B(y) + e^{4u} S(y). The scaling is written in BO's conformal frame; in physical units it gives growth e^{2u}, kappa_q = 2. B is the Bach tensor.
   - Named obstructions:
     - BO's own words: the sign is an assumption.
     - TW 2509.03584 App. D (src/tw2.txt l.6313-6320) criticise the Page approximation: "the Page approximation cannot consistently describe this perturbative process either, as it assumes a specific thermal equilibrium structure incompatible with critical collapse".
     - TW also say the state is unspecified: "The authors do not specify the quantum state beyond assuming it is time-independent".
2. **4D Hadamard point-splitting, any coupling xi (Zahn 2025).**
   - Leading term e^{2u}(T(y) + (u/4pi) V(y)), with V state-independent and local in curvature, and "independent of the renormalization scheme".
   - Named limits:
     - The state must "respect self-similarity", which Zahn argues for "(though not rigorously prove)" (src/zahn.txt l.131).
     - Back-reaction ends self-similarity.
     - TW App. D dispute whether self-similar states are physical and whether V dominates at finite u.
3. **2D s-wave anomaly / Polyakov-type dilaton-coupled effective action (Tomasevic-Wu 2025-26; Wu 2026; Lai-Wen 2026).**
   - The state is selected by regularity at the centre and at the past light cone.
   - Named obstructions:
     - "dimensional reduction and quantization do not commute in general" (src/tw.txt l.409).
     - The tensor must be regular at both the centre and the past light cone (TW 2509.03584 Fig. 6 caption, src/tw2.txt l.1214).
     - Lai-Wen: "Spherical reduction introduces a quantum coupling that grows toward the center, so the periodic first-order expansion is nonuniform".
4. **Ad hoc 2D tensor divided by 4 pi r^2 (Ayal-Piran 1997).** Obstruction, verbatim: "Unfortunately this tensor diverges at the origin." Regularising it by a factor Q(r) violates the Bianchi identity "unless Q = Const" (src/gr-qc_9704027.plain.txt l.145-151). BO also note it "is not exactly conserved".
5. **Cauchy-horizon divergence** (Hiscock-Williams-Eardley, as reported by Barve et al.; Barve et al.; Miyamoto-Harada). The quantum flux diverges on the Cauchy horizon of self-similar naked-singularity collapse. This sits beyond the region where the growing mode and the apparent horizon are computed. It does not block a linear-order interior calculation, but it does block any claim about the global exterior.
6. **Breakdown of exact CSS** (Chiba-Siino, ABSTRACT): "semiclassical equations of motion do not admit continuously self-similar solutions". This is consistent with treating hbar as a perturbation that grows like e^{2u}, but it rules out an "exact quantum CSS solution" claim.

## Normalisations (both sides in the same convention)
- **Convention:** BO, Zahn and TW all write perturbations as e^{kappa u}, with u = -ln(t_s - t) up to a constant (BO: "(tau - tau_s) proportional to e^{-u}", src/bo.txt section II). This is the same convention as our kappa*.
- **Quantum rate:** kappa_q = 2 in 4D (BO; Zahn; TW "kappa_q = D - 2").
- **BO's classical rates:** "kappa_c decreases monotonically from about 9.46 when k = 0 to 1.22 when k = 0.899 passing through kappa_c = 2 at k ~ 0.53" (src/bo.txt l.260), and beta = 1/kappa_c.
- **Ours:** kappa* = 2.810552549 (+/- 3.83e-10), beta = 0.3558019224 (+/- 6.97e-11), at k = 1/3 < 0.53. So kappa_c > kappa_q: BO case (ii), "we can say less".
- **TW's view:** "The relative ordering of kappa_c and kappa_q does not affect the sign of Delta p." In their framework the sign question stays well posed for radiation. This contrast is the doctrine #73 tension already logged in STRATEGY.md.
- **Exact digits:** INSPIRE full text has no hit for "0.3558" together with "semiclassical" outside the classical reviews (0711.4620, gr-qc/0210101) and a 1996 newsletter. No source prints a quantum coefficient for k = 1/3.

## Gap that is genuinely open (one sentence)
At the search floor, no source evaluates the one-hbar quantum source on the actual Evans-Coleman (k = 1/3) background, in any named stress-tensor model (BO/Page-Bach, Zahn's state-independent V, or a TW-style s-wave anomaly). No source therefore gives the forced kappa = 2 response, the sign of f_q at the trapping extremum, or the resulting threshold-shift sign, and nothing semiclassical has been certified on any true critical solution.

## Most dangerous neighbour
- **Tomasevic & Wu (UW / Perimeter).** Three papers in 2025-26 (2509.03584, 2509.03587, 2606.17135) plus Wu 2610.20705. They name radiation-fluid collapse as "the important next step", by simulation.
- **Zahn (Leipzig).** His Hadamard structure theorem makes V on any CSS background a local computation.
- **Jiang & Jiang (2607.04386).** They already prove rigidity of the leading RSET and compute the sign of its null projections on the Schwarzschild interior. The same move on a CSS fluid background is one step away for them.
- No race framing is intended. Ours would be a certified linear-order coefficient on a proved background, a different kind of result from a simulation. The overlap is real, and the owner already weighed it (SCOPE.md).

## Methods worth borrowing
1. **Zahn's local V tensor** (src/zahn.txt eq. 5-6 and eq. 11-12 basis): V_{mu nu} = (1/720 - (6 xi - 1)^2/288) dT(1) - (1/240) dT(2) (Zahn eq. 19, src/zahn.txt l.334-343, where the PDF text keeps only the denominators 720, 288, 240; the full expression as quoted in our RECON note BIGGER-IDEA.md l.60 from the HTML version).
   - It is state-free and local in curvature, so it can be evaluated with intervals on the certified E5 series.
   - At xi = 1/6 it reduces to the anomaly / Bach-type piece, which matches BO's 8 beta u B term.
2. **TW's sign criterion** sgn(Delta p) = sgn(f_q/f_c) at the extremum x* of the trapping function (src/tw2.txt l.4196-4218, l.4323 footnote 45 on where f_c > 0). This turns "which way the quantum push points" into a single certified sign.
3. **The review's dimensionful-parameter expansion** (2507.07636 Sect. II F, lab copy l.458-505). The Planck length enters as a relevant perturbation (L_P^2 e^{2u}), the opposite of a mass term.
4. **Positive controls:** BO's own structure (kappa_q = 2 and the u e^{2u} secular term), and TW's closed-form Roberts results (2509.03584 Sec. 5).

## Gate lines (OPEN is not refused: searchcheck passed on QUERIES-QB.log)
```
python "<MEMCAP>" 1 python "<TOOLS>/lab/tools/lab.py" gate <LAB>/quantum "O(hbar) correction to the Evans-Coleman / radiation-fluid critical solution (structure known, Evans-Coleman profile and amplitude not computed)" PARTIAL --file <LAB>/quantum/PRIOR-ART-QB.md --source "arXiv:gr-qc/9804058" --queries <LAB>/quantum/QUERIES-QB.log
python "<MEMCAP>" 1 python "<TOOLS>/lab/tools/lab.py" gate <LAB>/quantum "certified or rigorous semiclassical back-reaction on a critical collapse solution (analytic only on toy CSS Garfinkle/Roberts; nothing certified on a true critical solution)" PARTIAL --file <LAB>/quantum/PRIOR-ART-QB.md --source "arXiv:2509.03584" --queries <LAB>/quantum/QUERIES-QB.log
python "<MEMCAP>" 1 python "<TOOLS>/lab/tools/lab.py" gate <LAB>/quantum "sign of the one-hbar quantum push (threshold shift / trapping) for the radiation-fluid critical solution" OPEN --file <LAB>/quantum/PRIOR-ART-QB.md --source "arXiv:2610.20705" --queries <LAB>/quantum/QUERIES-QB.log
```
The same lines in the brief's short form:
```
lab.py gate <LAB>/quantum "O(hbar) correction to the Evans-Coleman / radiation-fluid critical solution" PARTIAL --file PRIOR-ART-QB.md --source "arXiv:gr-qc/9804058"
lab.py gate <LAB>/quantum "certified or rigorous semiclassical back-reaction on any critical collapse solution" PARTIAL --file PRIOR-ART-QB.md --source "arXiv:2509.03584"
lab.py gate <LAB>/quantum "sign of the quantum push for radiation" OPEN --file PRIOR-ART-QB.md --source "arXiv:2610.20705"
```
