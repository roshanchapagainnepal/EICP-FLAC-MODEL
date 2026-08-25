# EICP model in FLAC 8.0 — implementation notes

Implementation of Gai & Sánchez (2019), *An elastoplastic mechanical constitutive
model for microbially mediated cemented soils*, Acta Geotechnica 14:709–726,
as a FISH user-defined model (UDM) for FLAC 8.0.

## Programme status: COMPLETE

All four validation and sensitivity studies are finished. **16 FLAC runs.**

| Study | Cases | Status |
|---|---|---|
| Regression | t01 elastic, t02 Cam-Clay reduction | PASS |
| Fig. 10 validation | m_c = 0, 1.2, 2.4, 5.3 % | COMPLETE |
| p_c calibration | 900, 800 kPa (675 skipped) | COMPLETE — p_c = 800 kPa |
| μ sensitivity | μ = 1.5, 4.5, 8.5, 12.5 | COMPLETE |
| `a` sensitivity | a = 200, 300, 400 kPa | COMPLETE |
| Table 6 confinement | σ₃ = 100, 200, 400 kPa | COMPLETE |

**Headline result.** The model reproduces every qualitative behaviour the paper
reports — strength and stiffness increasing with calcite content, earlier and
sharper peaks, bond degradation driving post-peak softening, residuals
converging toward the untreated host sand, and the dilative-to-compressive
transition with confinement. Quantitatively it underpredicts peak strength by
up to 16 % using published parameters, with p_c calibrated once against the
untreated curve and never re-tuned.

**Parameter provenance.** κ, λ, M, D_s, ν, η, a and μ are the published values
from Tables 1 and 2, unchanged. Only p_c was calibrated, because the paper does
not report it. The `a` and μ sensitivity studies are reported sensitivities,
not recalibrations — both parameters are restored to their published values in
`01_parameters.fis`.

## How to run

```
call scripts\01_parameters.fis
call scripts\02_m_eicp.fis
model m_eicp
prop dens 1700.0 m_e 0.723 m_pc 450000.0 m_chi 1.0 m_mc 2.4 m_kb 1.0e10
ini sxx -200000.0 syy -200000.0 szz -200000.0
ini_eicp_r
print $eicp_err          ; must be 0.0
```

`01_parameters.fis` must be called before `02_m_eicp.fis`.
`ini_eicp_r` must be called after `PROPERTY` and `INITIAL`, before the first `STEP`.

## Corrections made to the original brief

Five points in the original specification were wrong or inconsistent and would
have produced silently incorrect results. All are verified against the manuals
in this folder.

### 1. `ln()` and `log()` were inverted

The FISH intrinsic table (FISH Reference p. 2-53/54) defines `ln(a)` as the
**natural** logarithm and `log(a)` as **base-ten**. The manual's own
`CAMCLAY.DAT` uses `lnp = ln(sp)`. The sub-loading law Eq. 6 requires the
natural log; using `log()` would introduce a factor of 2.3026.

### 2. Elastic moduli cannot be precomputed in case 1

FISH Reference §2.8.5 states that zone stresses are **undefined in modes 1, 3
and 4**, and that mode 3 is called *before* mode 1. Because K = (1+e)p′/κ is
stress-dependent, `m_e1`/`m_e2`/`m_g2` are updated at the end of each `zsub`
cycle in case 2 and stored in `f_prop` so case 3 can report `cm_max`.
`ini_eicp` seeds them from the initial stress state — this is the analogue of
`camclay_ini_p` in the built-in Cam-Clay model.

### 3. `m_mc` is in PERCENT, not a fraction

Eq. 3 is p_b = a·χ·m_c with a = 200 kPa. A calcite content of 2.4 % entered as
a fraction (0.024) gives p_b = 4.8 kPa, which is negligible next to the
confining pressure — the cementation mechanism would effectively vanish.
Entered as percent it gives p_b = 480 kPa, which is the order of magnitude
plotted in paper Fig. 12b.

### 4. `p_c` and `R` are not independent

The sub-loading surface passes through the current stress point by
construction. Setting F_SUB = 0 at the initial state and solving Eq. 5 gives

    R₀ = (q₀²/M² + p₀²) / (p₀ (p_c + p_b))

which for an isotropic start reduces to R₀ = p′₀/(p_c + p_b). The brief's
`ini_pc = 150000` with `ini_R = 0.01` describes a sample at p′₀ = 1.5 kPa, not
the 100–200 kPa of the tests being modelled. Starting a run that way projects
the stress onto a surface two orders of magnitude too small on the first step.

`ini_eicp_r` computes R₀ consistently and sets `$eicp_err = 1.0` if R₀ > 1
(which means p_c is too small for the applied confining stress). **Always check
`$eicp_err` after initialising.**

Use plain `ini_eicp` instead when R is deliberately held fixed — for example
R = 1 to collapse the model to Cam-Clay in test 2.

### 5. Source of the UDM structure

`WRITING NEW CONSTITUTIVE MODELS.pdf` documents **C++ DLL** models only; it
ends at Example 2.4 and never mentions FISH. The FISH UDM specification is
FISH Reference §2.8, Examples 2.26–2.32. The structural notes in the brief were
correct, but §2.8.7 adds a constraint: `friend` functions must be defined
*after* the UDM that declares them, so a master file using `friend` must be
called *before* its helpers.

Minor: the 12-character limit on `f_prop` names is a display limit on plots and
print headings, not a compile error.

## Design decisions

- **Sign convention.** FLAC is tension-positive; p′ and q inside the model are
  compression-positive. Conversion happens at the top and bottom of case 2.
- **2D.** FLAC 2D provides `zs12` only, so J₂ = ½(s₁₁²+s₂₂²+s₃₃²) + s₁₂².
- **Hardening law, Eq. 10 — CONFIRMED.** Version B implemented and verified
  against paper p. 713. D_s is dimensionless. p_c does **NOT** multiply the
  D_s term:

      dp_c = p_c (1+e)/(λ−κ) dε_v^p  +  D_s (1+e)/(λ−κ) dε_q^p

  p_c scales the volumetric term only; the D_s term is a separate additive
  increment. Implemented at `scripts/02_m_eicp.fis:236`.

  The D_s term is numerically small by design — this is correct behaviour,
  not a bug. With D_s = 0.04 and (1+e)/(λ−κ) ≈ 3.5 it contributes a fraction
  of a pascal per step against a p_c of hundreds of kPa. D_s is described in
  Table 1 as a dimensionless soil dilation parameter associated with dilation
  at failure, and carries no units.
- **Damage integration.** dχ = −μχ dε_q^p is stiff, so it is integrated
  exponentially as χ ← χ·exp(−μ dε_q^p), which keeps χ strictly positive.
- **Return mapping.** Follows `modelcamclay.cpp` lines 250–300 — an elastic
  predictor with a quadratic solve for the plastic multiplier, hardening lagged
  to the end of the zone cycle. This is FLAC's standard explicit scheme.
- **Bulk modulus bound `m_kb`.** K grows with p′; without a cap the timestep
  becomes unstable. The built-in Cam-Clay carries the same guard
  (`bulk_bound`) and aborts if it is exceeded; here K is clamped instead.

## Test suite

| Test | Isolates | Acceptance criterion |
|---|---|---|
| `t01_elastic.dat` | Elastic predictor, K = (1+e)p′/κ | e vs ln p′ straight, slope −κ = −0.003; `m_ind` stays 0 |
| `t02_isotropic.dat` | Yield surface, return map, p_c hardening | Slope breaks from −κ to −λ = −0.25 at p′ = 150 kPa; `m_ind` → 1 |
| `t03_triaxial.dat` | Full model incl. bonding and damage | Paper Figs. 10–12: peak q rises with m_c, dilation increases, χ falls, R → 1 near 7 % strain |

`t02` calls `eicp_mcc`, which zeroes D_s and η to collapse the model exactly to
modified Cam-Clay. Run the tests in order — t01 and t02 must pass before any
discrepancy in t03 can be attributed to the MICP terms.

## p_c Sensitivity Study — COMPLETE

> **Outcome: p_c = 800 000 Pa is the CALIBRATED value.** The untreated baseline
> is validated against Fig. 10a of Gai & Sánchez (2019) at σ₃ = 200 kPa: peak
> q = 420 kPa at ~5 % axial strain, against a target of 410–420 kPa at
> 4.5–5.0 %. This single value carries forward to every m_c case unchanged,
> per the paper's "reference soil" methodology — it must not be re-tuned per
> calcite content.


p_c is an initial state variable, not a measured parameter, and Gai & Sánchez do
not report it. It is therefore calibrated against Fig. 10a. It is set from a
single line — `$ini_pc0` in `scripts/01_parameters.fis` — which `set_case` in
`t03_triaxial.dat` pushes into `m_pc`.

### Calibration targets (read from Fig. 10a, σ₃ = 200 kPa)

| Quantity | Target |
|---|---|
| Untreated peak q | 410–420 kPa |
| Strain at untreated peak | 4.5–5.0 % |
| Untreated large-strain residual | 340–350 kPa |

The residual is an independent confirmation of **M = 1.09**. Critical state
requires q = M·p′ on the drained path p′ = σ₃ + q/3:

    q = 1.09 (200 + q/3)  →  q = 218 / 0.6367 = 342 kPa

This falls inside the observed 340–350 kPa band, so M needs no adjustment. The
calibration is therefore about p_c alone — it sets the *peak*, while M already
sets the *residual* correctly.

### Trial values

| Trial | p_c |
|---|---|
| A | 675 kPa |
| B | 900 kPa |
| C | 1125 kPa |

### Independent check on the required p_c

At the peak, with p_b = 0 (untreated) and R at its ceiling of 1, the stress must
lie on the yield surface while following the drained path, giving

    p0 = (q²/M² + p'²) / p'    with   p' = 200 + q/3   [kPa]

| Target peak q | p′ at peak | Required p₀ |
|---|---|---|
| 400 kPa | 333 kPa | **737 kPa** |
| 415 kPa | 338 kPa | **767 kPa** |
| 420 kPa | 340 kPa | **777 kPa** |

**Note a discrepancy:** this calculation gives 737 kPa for q = 400 kPa, not the
810 kPa quoted when the study was set up. 810 kPa would follow if R were about
0.91 rather than 1.0 at the peak. Worth resolving, since it shifts where the
answer sits between trials A and B.

Two effects push the *initial* p_c **above** these numbers:

1. R < 1 at the peak in practice, and p₀ = R(p_c + p_b), so p_c = p₀/R.
2. At peak the sample is on the dry side and dilating, so dε_v^p < 0 and p_c has
   been *decreasing* from its initial value.

**Expectation for Trial A:** 675 kPa is below even the most favourable estimate
(737 kPa), so Trial A should **undershoot** the 410–420 kPa target. That is
still a useful bracket point — it is the lower bound of the study.

### Results

Trial A was skipped. All runs: untreated (m_c = 0, e₀ = 0.723), σ₃ = 200 kPa,
velocity −1.0e-7, 2 200 000 steps = 22 % axial strain.

| Trial | p_c | Peak q | Strain at peak | Residual at 22 % | Verdict |
|---|---|---|---|---|---|
| Target (Fig. 10a) | — | 410–420 kPa | 4.5–5.0 % | 340–350 kPa | — |
| A | 675 kPa | skipped | — | — | predicted undershoot |
| B | 900 kPa | 460 kPa | 4–5 % | 390 kPa | peak ~10 % high |
| **D** | **800 kPa** | **420 kPa** | **~5 %** | **375–380 kPa** | **PASS — calibrated** |

**Trial B:** shape correct, peak ~10 % high. Strain at peak already landed in
the target window, so the *form* of the response was right and only the
magnitude needed trimming.

**Trial D:** peak q = 420 kPa at ~5 % strain, hitting the target band, with
correct post-peak softening to 375–380 kPa at 22 %. R₀ = 0.25 saturating to
0.95–0.97 by 7–8 % strain. p_b = 0 and `unbal` = 0 throughout.
**p_c = 800 000 Pa adopted as the calibrated value.**

### η = 60 confirmed

R rose from R₀ = 0.22 to ~0.96 by 7–8 % strain, against the paper's report of
R → 1 at ~7 % (Fig. 14). **η = 60 from Table 1 is correct and needs no
adjustment.** This matters: η governs the strain scale of the response while
p_c governs its magnitude, so confirming η leaves p_c as the only free
quantity. An earlier concern that η was wrong came from runs that had only
reached 1.2 % strain — that was a drive-length artefact, not a model defect.

### Calibrated p_c is between 675 and 900 kPa

Independent check of the interpolation, using p₀ = (q²/M² + p′²)/p′ with
p′ = 200 + q/3:

| | q | p′ | required p₀ |
|---|---|---|---|
| Trial B achieved | 460 kPa | 353 kPa | 857 kPa |
| Target | 415 kPa | 338 kPa | 767 kPa |

Ratio 767/857 = 0.894, so p_c ≈ 900 × 0.894 = **805 kPa**. Next trial set to
**800 kPa**, which this estimate supports independently of the linear
interpolation between trials.

### Note on the residual

Trial B's 390 kPa at 22 % is above the 340–350 kPa target, but this is **not a
separate calibration problem**. At q = 390 kPa, p′ = 330 kPa and q/p′ = 1.18,
still above M = 1.09 — the sample has not yet reached critical state and is
still softening toward it (p_c was still declining at 22 %, 920 → 710 kPa).
The asymptote is M-controlled and independent of p_c: q_crit = 342 kPa. Do not
attempt to correct the residual with p_c; lowering p_c to 800 kPa will reduce
the 22 % residual as a side effect simply because there is less dilation to
unwind.

## m_c = 2.4 % Validation — treated sand

σ₃ = 200 kPa, e₀ = 0.715, p_c = 800 kPa (calibrated on untreated data, **not**
re-tuned), a = 200 kPa and μ = 6.5 taken unchanged from Table 2. Velocity
−1.0e-7, 2 200 000 steps = 22 % axial strain.

**This run is a one-calibration-parameter prediction.** p_c was calibrated from
the untreated data; every parameter specific to the cementation (a, μ) is the
published Table 2 value, unchanged. Nothing was fitted to the treated response.

| Quantity | Model | Paper (Fig. 10) | |
|---|---|---|---|
| Peak q | 570 kPa | ~680 kPa | **−16 %** |
| Strain at peak | 3.4–3.5 % | ~2–3 % | close |
| Residual at 22 % | 375 kPa | — | — |
| p_b decay | 480 → 95 kPa (−80 %) | — | mechanism active |
| p_c decay | 800 → 580 kPa | — | dilation correct |

### Qualitative validation: PASS

All five signature behaviours of MICP-treated sand reproduce correctly:

1. treated curve above untreated throughout
2. higher peak than untreated (570 vs 420 kPa)
3. earlier peak than untreated (3.4 % vs 5 %)
4. more post-peak softening than untreated
5. p_b degrades by 80 %, driving that extra softening

`unbal` = 0 throughout; R evolution matches the untreated case.

**Qualitative consistency observation — not a validation.** Both curves
approach a similar residual at 22 % strain: 375 kPa treated against
375–380 kPa untreated. This is *consistent with* the paper's remark that a
treated specimen's residual strength tends toward the untreated one at large
deformation, but it cannot be claimed as independent validation, for two
reasons. First, p_b is still ~95 000 Pa at 22 % — the bonding has degraded by
80 % but is not fully removed, so the two materials are not yet equivalent.
Second, Fig. 10 presents experimental data only to 12 % strain, so a
comparison at 22 % lies beyond the direct experimental range.

### Quantitative finding: 16 % underprediction of peak strength

**Recorded as a finding, not corrected by adjusting parameters.** `a` = 200 kPa
is the published Table 2 value; re-tuning it to close the gap would convert an
independent prediction into a curve fit and destroy the validation's value.

**The corresponding gap in `a` cannot be inferred directly from the q
shortfall.** R, p_c, χ, the dilation response and the stress path all evolve in
a coupled way during the test, so there is no valid closed-form
back-calculation from a peak-strength difference to a parameter value. Any
earlier estimate along those lines should be disregarded. Determining what
value of `a` would close the gap is a job for the sensitivity study listed
below, run case by case — not for algebra.

### Context for the discrepancy

- Gai & Sánchez state their aim was to capture the main tendencies rather than
  precisely reproduce the experiments, and they deliberately used one
  "reference soil" parameter set across all tests rather than fitting each.
- p_c is not reported in the paper and had to be calibrated here, so one
  degree of freedom was introduced that the original authors set differently.
- The paper integrates the model with Sloan et al. (1987) sub-stepping at the
  "point integration level"; this implementation uses FLAC's explicit
  single-step return mapping with lagged hardening. This is a **plausible
  contributing factor, not a confirmed explanation.** Establishing it as
  causal would require a numerical convergence check — re-running with
  progressively smaller strain increments per step and demonstrating that the
  peak converges toward the paper's value. Until that is done, the integration
  scheme should be cited as a candidate factor only.

## Figure 10 Validation Set — COMPLETE

All four Table 4 cases at σ₃ = 200 kPa, p_c = 800 kPa (calibrated once on the
untreated case and **never re-tuned**), a = 200 kPa and μ = 6.5 unchanged from
Table 2. Velocity −1.0e-7, 2 200 000 steps = 22 % axial strain.

| m_c | e₀ | Model peak q | Paper peak q | Ratio | Strain at peak |
|---|---|---|---|---|---|
| 0 % | 0.723 | 420 kPa | 420 kPa | **1.00** | 5.0 % |
| 1.2 % | 0.718 | 490 kPa | ~510 kPa | **0.96** | 4.0 % |
| 2.4 % | 0.715 | 570 kPa | ~680 kPa | **0.84** | 3.4 % |
| 5.3 % | 0.709 | 745 kPa | ~825 kPa | **0.90** | 3.0 % |

m_c = 0 % is the calibration point, so its ratio of 1.00 is by construction and
is not an independent result. The other three are predictions.

Supporting state variables, m_c = 5.3 %: p_b 1060 → 190 kPa (82 % degradation),
χ 0.99 → 0.18–0.19, R 0.108 → 0.95–0.97 by 7–8 % strain, residual 371 kPa at
22 %, `unbal` = 0 throughout.

### Qualitative trends: all correct

Across the full set, peak strength rises with m_c, the peak occurs
progressively earlier (5.0 → 4.0 → 3.4 → 3.0 %), post-peak softening becomes
more pronounced, and all four residuals converge toward a similar value
(370–380 kPa). Every trend the paper reports is reproduced.

### Finding: the model/paper ratio is non-monotone

The ratio falls from 0.96 to 0.84 and then partially recovers to 0.90. **This
is inconsistent with a simple "a is too small" explanation**, which would
predict a monotone decline as the p_b contribution grows with m_c. The pattern
instead suggests coupled interaction between p_b, p_c decay and χ degradation
at high cementation.

Three caveats belong with this finding and should be stated alongside it:

1. **Reading precision.** The paper's peaks were read from a plotted figure.
   An uncertainty of a few per cent in those readings is comparable to the
   0.84 → 0.90 recovery itself, so the non-monotonicity is suggestive rather
   than established.
2. **The m_c = 5.3 % case is compromised.** Gai & Sánchez report localised
   shear-band behaviour in that specimen and place its modelling out of scope.
   Their own model deviates there. The experimental value is affected by
   localisation the model does not represent, so the 0.90 ratio is the least
   reliable of the four.
3. **No mechanism is inferred here.** As established for the m_c = 2.4 % case,
   the coupled evolution of R, p_c, χ, dilation and the stress path means no
   parameter value can be back-calculated from peak-strength ratios. The `a`
   and μ sensitivity studies are the valid route to attribution.

## μ Sensitivity Study — COMPLETE

Reproduces Gai & Sánchez Figs. 13 and 14. m_c = 2.4 %, e₀ = 0.715,
**σ₃ = 100 kPa** (not the 200 kPa of the Fig. 10 set), p_c = 800 kPa,
a = 200 kPa. Velocity −1.0e-7, 2 200 000 steps = 22 % axial strain. Only
`$mu_p` varied. `unbal` = 0 throughout all four trials.

| μ | Peak q | Strain at peak | χ at 22 % | p_b at 22 % | p_c at 22 % |
|---|---|---|---|---|---|
| 1.5 | 527 kPa | 3.0–3.2 % | 0.65 | 310 kPa | 163 kPa |
| 4.5 | 501 kPa | 2.8–3.0 % | 0.28 | 135 kPa | 320 kPa |
| 8.5 | 480 kPa | 2.5 % | 0.09 | 50 kPa | 362 kPa |
| 12.5 | 470 kPa | 2.6 % | 0.03 | 13 kPa | 388 kPa |

### Discrepancy in the paper: μ = 16.5 vs 12.5

**Table 5 lists the fourth case as μ = 16.5, but Figs. 13 and 14 are plotted
for μ = 12.5.** These trials used **12.5** so the results are comparable with
the plotted figures. This inconsistency is in the published paper, not in this
implementation, and should be stated explicitly in the write-up rather than
silently resolved.

### Findings

1. **Peak strength is insensitive to μ.** An 8.3-fold increase in μ moves the
   peak only 527 → 470 kPa, an 11 % range. This is consistent with the paper's
   own statement that μ has practically no effect on initial stiffness and
   only a modest effect on peak strength.

2. **χ and p_b are highly sensitive, and nonlinearly so.** The same 8.3-fold
   change in μ produces a 22-fold change in χ (0.65 → 0.03) and drives p_b
   from 310 kPa to 13 kPa. The damage law amplifies rather than tracks its
   input, which follows from the exponential integration χ ← χ·exp(−μ dε_q^p).

3. **p_c rises monotonically with μ — the compensation mechanism.** p_c at
   22 % goes 163 → 320 → 362 → 388 kPa, i.e. *opposite* in direction to p_b.
   Faster bond degradation shrinks p_b, which shrinks the total surface
   p₀ = R(p_c + p_b), which moves the stress point relatively less far onto
   the dry side, so the sample dilates less and p_c is eroded less. **The two
   hardening parameters partially compensate**, and this is what makes peak
   strength so much less sensitive to μ than χ is (finding 1). μ does not act
   on p_b alone; it reaches p_c indirectly through the dilation response.

   This is the same coupling that makes parameter values impossible to
   back-calculate from peak-strength comparisons — here it is visible directly
   in the data rather than argued in the abstract.

4. **Post-peak softening plateaus above μ ≈ 4.5.** The post-peak drop is
   263 / 250 / 255 kPa for μ = 4.5 / 8.5 / 12.5 — essentially flat. Softening
   is bounded below by the untreated critical-state residual: once the bonding
   is substantially gone, the material cannot soften past the host sand, so
   further increases in μ have nothing left to remove.

5. **μ = 12.5 approaches complete debonding.** p_b = 13 kPa at 22 % strain is
   97 % degradation, so the treated response has effectively collapsed onto
   the untreated one at large strain.

### Housekeeping after the study

`$mu_p` restored to **6.5** and σ₃ restored to **200 kPa** in both places
(`ini sxx/syy/szz` and the servo target `$s_target`). Every Fig. 10 result
depends on those values.

## `a` Sensitivity Study — COMPLETE

m_c = 2.4 %, e₀ = 0.715, σ₃ = 200 kPa, p_c = 800 kPa, μ = 6.5. Only `$a_pa`
varied. Velocity −1.0e-7, 2 200 000 steps = 22 % axial strain. `unbal` = 0
throughout all three trials.

**This is a reported sensitivity, not a recalibration.** a = 200 kPa is the
published Table 2 value and remains the validation case for every Fig. 10
result. It has been restored in `01_parameters.fis`.

| a | R₀ | Peak q | vs paper (~680 kPa) | p_c at 22 % | p_b at 22 % | χ at 22 % |
|---|---|---|---|---|---|---|
| 200 kPa | 0.156 | 570 kPa | **−16 %** | 580 kPa | 95 kPa | 0.18 |
| 300 kPa | 0.132 | 640 kPa | **−6 %** | 540 kPa | 130 kPa | 0.18 |
| 400 kPa | 0.11 | 710 kPa | **+4 %** | 496 kPa | 170 kPa | 0.17 |

### Finding: effective a ≈ 320–340 kPa

The paper's peak is bracketed between a = 300 and a = 400 kPa, giving an
effective value of roughly **320–340 kPa — some 60–70 % above the published
200 kPa**. This is the quantity that could not be obtained analytically; the
sensitivity run is the only valid route to it.

The response is strongly **sublinear**: doubling a from 200 to 400 kPa raises
the peak only 570 → 710 kPa (25 %), because p_b is one component of
p₀ = R(p_c + p_b) and the compensation below absorbs part of every increase.

### The p_c compensation mechanism, confirmed from both directions

This is the most robust mechanistic result of the study. The same coupling was
reached independently through two different parameters, moving in opposite
directions:

| Study | Parameter change | p_b at 22 % | p_c at 22 % |
|---|---|---|---|
| μ | 1.5 → 12.5 (bonds lost faster) | 310 → 13 kPa ↓ | 163 → 388 kPa **↑** |
| a | 200 → 400 kPa (more bonding) | 95 → 170 kPa ↑ | 580 → 496 kPa **↓** |

Mechanism: p_b and p_c both enter the surface size p₀ = R(p_c + p_b). More
bonding enlarges p₀, pushing the stress point further onto the dry side, so
the sample dilates harder and p_c is eroded faster — and vice versa. The two
hardening parameters partially cancel.

**Control:** χ stayed at 0.17–0.18 across all three a-trials, as it must —
Eq. 4 makes χ a function of plastic deviatoric strain and μ only, with no
dependence on a. This confirms the p_c shift came from the dilation coupling
and not from any change in the damage law.

**Implication.** This is the direct evidence for the reviewer's point that no
parameter value can be back-calculated from peak-strength comparisons. The
coupling is not a theoretical caveat — it is visible in the data, in two
independent parameter directions, with a clean control.

## Table 6 Confinement Series — COMPLETE

Reproduces Fig. 15. p_c = 800 kPa, a = 200 kPa, μ = 6.5, all published or
previously calibrated values, unchanged. Velocity −1.0e-7, 2 200 000 steps
= 22 % axial strain. `unbal` = 0 throughout all three trials.

| σ₃ | e₀ | m_c | R₀ | Peak / plateau q | q_crit | p_c at 22 % | Side |
|---|---|---|---|---|---|---|---|
| 100 kPa | 0.723 | 0.9 % | 0.10 | 405 kPa at 3.0–3.5 % | 171 kPa | 390 kPa ↓ | **dry** |
| 200 kPa | 0.718 | 1.2 % | 0.192 | 490 kPa at 4.0 % | 342 kPa | 632 kPa ↓ | **dry** |
| 400 kPa | 0.715 | 1.4 % | 0.370 | ~590 kPa, still rising | 684 kPa | 1040 kPa ↑ | **wet** |

χ fell to 0.16 / 0.18 / 0.22 and p_b lost 83 / 81 / 78 % across the three runs.

**Consistency check.** The σ₃ = 200 kPa trial reproduced the Fig. 10 m_c = 1.2 %
run exactly — 490 kPa at 4.0 % strain — confirming no parameter drift across
the μ and `a` sensitivity studies that ran between them.

### Finding: the model reproduces the wet/dry side transition

At 100 and 200 kPa the specimen sits on the **dry** side of critical state: it
dilates, reaches a peak, and softens. At 400 kPa it sits on the **wet** side:
it compresses and hardens continuously with no peak. This is exactly the
transition Fig. 15 reports, and it emerges from the model rather than being
imposed.

The mechanism is R₀ = σ₃/(p_c + p_b). Raising σ₃ from 100 to 400 kPa while
p_b barely changes (m_c only moves 0.9 → 1.4 %) drives R₀ from 0.10 to 0.370,
so the stress point starts far closer to the surface apex and never travels
onto the dry side.

### The clearest signature: the direction of p_c

**p_c falls at 100 and 200 kPa but rises at 400 kPa.**

Sign of dε_v^p is set by ∂F/∂p = M²(2p′ − p₀). On the dry side p′ < p₀/2, so
the term is negative, the specimen dilates and p_c erodes. On the wet side
p′ > p₀/2, the term is positive, the specimen compacts and p_c hardens. The
p_c history is therefore a direct readout of which side the specimen is on —
a cleaner diagnostic than the volumetric strain curve, which mixes elastic and
plastic contributions.

### Caveats

1. **22 % strain is not sufficient at 400 kPa.** q_crit = 684 kPa and the
   model is at ~590 kPa still rising, so this case has not reached critical
   state. The 100 kPa case is also incompletely converged (230 kPa residual
   against q_crit = 171 kPa). **Residuals are not comparable across the series
   at a fixed 22 % strain**, and the write-up should say so rather than
   tabulating them as though they were.
2. **This is not a controlled confinement study.** Table 6 varies σ₃, e₀ *and*
   m_c simultaneously, and the paper notes the higher-confinement samples were
   incidentally prepared with more calcite. The series tests whether the model
   tracks the combined experimental conditions — which is the right validation
   question — but it cannot isolate a confining-pressure effect. Isolating it
   would need a synthetic series holding e₀ and m_c fixed, which has no
   experimental counterpart.

## Still to verify

Nothing remains for the validation programme. Optional extensions:

- **Table 7 loading paths** (Figs. 17, 18): radial extension and constant-p′
  tests at σ₃ = 100 kPa, m_c = 1.25 %. These need a modified servo — constant
  p′ requires driving both stress components together — so they are a
  meaningful extension rather than a parameter change.
- **Numerical convergence check.** Re-running with smaller strain increments
  per step would establish whether the explicit single-step return mapping
  contributes to the 16 % peak shortfall. Until done, the integration scheme
  remains a candidate factor only, not a confirmed explanation.
- **Longer runs at 100 and 400 kPa** if comparable residuals are wanted.
