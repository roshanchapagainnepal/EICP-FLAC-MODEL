# EICP model in FLAC 8.0 — implementation notes

Implementation of Gai & Sánchez (2019), *An elastoplastic mechanical constitutive
model for microbially mediated cemented soils*, Acta Geotechnica 14:709–726,
as a FISH user-defined model (UDM) for FLAC 8.0.

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

## p_c Sensitivity Study

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

## Still to verify

- Velocity magnitudes and step counts in the test files are estimates. Watch
  the `unbal` history and the axial strain, and tune so that loading stays
  quasi-static and reaches ~12 % axial strain in t03.
- `t03` currently uses p_c = 450 kPa for all cases. The paper does not state
  p_c directly; this value was chosen so that R(p_c+p_b) reaches roughly
  930 kPa as R → 1, consistent with the reported peak strengths.
