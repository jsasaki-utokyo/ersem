# Modification history of this ERSEM fork (entry point)

This file is the entry point for users of this fork (other projects, other sessions): what was changed, why, how
it is switched, and what a configuration must set to use it. The newest entry is first. Per-module notes stay in
`doc/*.md`; this file links them.

---

## 2026-11-28: TB-FVCOM box-ersem aggregation sinking of diatoms (new default-off module, one default-off producer option)

Made by the box-ersem work (`~/Github/TB-FVCOM/box-ersem`) from its registered design `docs/aggregation_sinking.md` (gpt-6.2-sol review, PASS in round 6; section 5.1). Nothing changes unless a configuration declares the new instances or sets the new key; with the module compiled in and no instance declared the model is bit identical to the previous build (gate I-G1, box-ersem `analysis/results/aggsink/gates/`).

| Item | Switch (default) | What | Record (`~/Github/TB-FVCOM/box-ersem/docs/`) |
|---|---|---|---|
| new module `phyto_aggregation` (`src/phyto_aggregation.F90`, model name `ersem/phyto_aggregation`) | opt-in: only a configuration that declares an instance uses it | shear coagulation and break-up as a mass transfer between two pools of one functional type (single cells P1 and aggregates P1A) on c, n, p, s and chlorophyll: `k_agg = k0 + (8/pi) alpha_c gamma_c phi1` (single-cell volume fraction `phi1 = v_C C`), `k_dis = min(k_d <eps>/nu / gamma_b^2, k_dis_max)`; donor and receiver sources are equal and opposite, so every conserved quantity is unchanged. Turbulence input: the dissipation moments `box_shear_sqrt_dissipation`, `box_shear_dissipation` (`shear_source` 1) or one constant `box_shear_rate` (`shear_source` 2), linked by the host by standard variable name; `c_star` > 0 replaces the single-cell carbon of the collision term by a fixed value. Needs a host that supplies those variables (the TB-FVCOM box driver does; FVCOM-FABM and fabm0d do not) | `aggregation_sinking.md` sections 4.3, 5.1; gates I-G3, I-G4 |
| `primary_producer` | `isw_agg_pool` (0 = off: nothing read, nothing registered, former arithmetic) | 1 = zero-stock-safe producer for the aggregate pool: no background carbon (`c0` forced to 0), the empty-stock ratios (p/c, n/c, Chl/c, Fe/c, also in `get_sinking_rate`) guarded so that an empty pool has exactly zero biological sources and keeps the prescribed sinking velocity `rm`; the arithmetic of a stock above zero is unchanged | `aggregation_sinking.md` 5.1; gates I-G1, I-G2 |

Not changed: grazers do not eat the aggregate pool (the configuration generator adds no prey coupling); the configuration generator is box-ersem `analysis/aggsink/make_agg_yaml.py`.

---

## 2026-11-28: MUSE rev53 ERSEM text in the MUSE source copies (not in this repository)

Made by the nippon-steel / muse rev53 work (`~/Github/muse/docs/carb_s10m/rev53/`; tank closures, adopted state s53, binary `fix53-v5a`); records `design.md`, `design_b.md`, `stage3.md`. The ERSEM modules are edited ONLY in the MUSE source copies
`~/Github/muse/docs/carb_s10m/rev53/source/ersem/src/`; nothing was added to this repository by rev53 and the box-ersem build does not contain them. Every switch defaults to the former behaviour (all off: byte identical to the pinned rev52 binary).

| Module (MUSE copy) | Switch / key (default) | What | Relation to this repository |
|---|---|---|---|
| `primary_producer`, `microzooplankton`, `mesozooplankton`, `benthos_shared_laws` | `isw_o2lim` (0) and its `o2lim_*` keys | O2 limitation of rest respiration and feeding, unmet part to lysis / mortality | the port of ersem 41cfcee (box-ersem), re-marked `jsasaki 2026-11-26: rev53 E1/E2` |
| `primary_producer`, `pom_decay` | none | light-response guard (E4, ersem bf61176 + da38386) and shared `CMass` (E7, ersem 82d67a2) | ports of the commits named |
| `bacteria_docdyn` | `o2law_mask` (0; bits 1 Q10 law, 2 uptake limiter, 4 saturation-dependent yield, 8 rest respiration with withheld maintenance), `o2law_K` (3.0), `o2law_n` (1.0); diagnostics `o2law_bits`, `o2law_F`, `o2law_rest_unl` and others | pelagic bacteria B1 O2 law (rev53 part B) | NEW in the MUSE copy only; not ported to this fork; port only with its gates (muse `docs/carb_s10m/rev53/stage2_merged.md`) |

Adoption evidence: `~/Github/muse/docs/carb_s10m/rev53/stage3.md` (audited PASS), `~/Github/muse/docs/MODIFICATION_HISTORY.md` (2026-11-28 entry).

---

## 2026-11-22 to 2026-11-23: TB-FVCOM box-ersem changes (FABM box driver, Tokyo Bay)

Made by the box-ersem work (`~/Github/TB-FVCOM/box-ersem`); records in its `docs/`. One line each: switch and default.

| Commit | Module | Switch (default) | What | Record (`~/Github/TB-FVCOM/box-ersem/docs/`) |
|---|---|---|---|---|
| 41cfcee | `primary_producer`, `microzooplankton`, `mesozooplankton`, `benthos_shared_laws` | `isw_o2lim` (0 = off: nothing read, nothing registered) | O2 limitation of plankton aerobic metabolism (P1-P6, Z4-Z6): thresholded Hill factor `o2_limit` shared with the MUSE benthos; the withheld rest respiration goes to lysis / mortality | `plankton_o2_limitation.md` (review rounds 1-3 PASS) |
| 13a6a05 | new `phyto_migration` (models `phyto_migration_target`, `phyto_migration_move`) | opt-in: only a configuration that declares the instances uses it | phytoplankton diel vertical migration by swimming; production yaml unchanged with the module compiled in (771 / 771 variables identical, job 7653670) | `phytoplankton_types.md` 11.3 |
| bf61176, da38386 | `primary_producer` | `iswTemp` 3 = Rosso CTMI (default 1 unchanged; 2 unchanged apart from a logged warning); the light-response guard `sum /= 0` is a defect fix without switch | iswTemp 3 (valid for Topt > (Tmin + Tmax)/2); the guard removes NaN / division by zero at a zero gross rate (sum = -0 at a root of the iswTemp 2 cubic) and keeps every finite result, bit identical in the Q10 gate | `ersem_issues.md` section 1; `doc/phytoplankton_temperature_ctmi.md` |
| ea83637 | `sulfur_cycle`, new `sulfur_laws` | `isw_law` (0 = legacy, untouched) | 1 = the shared water / bed sulfur laws of MUSE (mass action, per-reaction Q10), 32 Strang sub-passes of exact extents, sulfate anomaly SO4a, conserved total_reactive_sulfur; needs a host that supplies the step length AND integrates with Euler (`ode_method = 1`), as the TB-FVCOM box driver does; FVCOM-FABM (3-D) does not supply the step length (box-ersem `sulfur_unification.md`, "Host dependence"; corrected 2026-11-25 at box-ersem's request; fabm0d not verified). **Never enable it together with the rev47 `common47` sulfur of the MUSE ERSEM copies**: they are two implementations of the same law (muse `docs/carb_s10m/rev52/box_ersem_intake.md` E6) | `sulfur_unification.md` |
| 82d67a2 | `pom_decay` | none (defect fix; changes results of configurations that use `pom_decay`) | DIC of direct oxidation uses `ersem_shared` CMass 12.011 instead of a local 12.0 (12.0 created carbon) | `m1_deposition.md` section 6.1 |

---

## 2026-11-25: MUSE rev52 model-consistency fixes (tank closures, Nippon Steel mesocosm paper)

**Why it matters to every ERSEM configuration of this fork.** An equation inventory of the coupled MUSE (bed) +
ERSEM (water) model found inconsistencies, several of which sit in ERSEM defaults and therefore affect any
configuration that keeps them (e.g. TB-FVCOM box-ersem). The full record is in the muse repository:
`~/Github/muse/docs/carb_s10m/rev52/triage.md` (25 candidates, 12 defects), `design.md` (fixes, sources, gates;
Stage 1 review PASS), `stage2.md` (gates and start-state audit PASS). Every code change is behind a switch whose
default reproduces the former behaviour bit for bit (verified: gate c0, muse rev52).

### 1. Code changes and their switches

| Module | Switch (default) | `1` means | Status in this repository | Record |
|---|---|---|---|---|
| `src/dom_decay.F90` | `isw_law` (0 = former) | pure Q10 temperature law `th_q10(T, q10, Tref)` (no heat cut; the law the MUSE bed uses for DOC decay) and DIC produced with `ersem_shared` CMass 12.011 instead of the local 12.0 (which lost 0.09 % of the oxidised carbon) | **in this repository** (e5282f5, ad03fc7) | rev52 triage T01, T02 |
| `src/primary_producer.F90` | `isw_no3red` (0 = former) | nitrate assimilated by P1-P4 returns 2 O2 per N (nitrate to amine is 8 e-): closes the electron balance of the producers under the CH2O convention (section 2); carbon, N and TA untouched; diagnostic `ledger_no3red_O2o` registered only when on | **in this repository** (a9bc3af; ported from the muse rev52 source copy that the s52 binary executed, merged into this fork's current file) | rev52 design 2.4, M2/M3 |
| `src/carbonate_engine.F90` | none (new public function) | `carbonate_engine_kw(T, S, Pbar, opt_pH_scale, opt_k_carbonic, opt_total_borate)`: K_W on the pH scale of the solve (used by the MUSE bed; fixes K_W on the SWS scale used with a total-scale pH) | **in this repository** (a7ae58a; identical to the muse rev52 copy apart from the edit marks) | rev52 triage T06 |

**Port gate (2026-11-25, after the adoption of s52).** `tests/isw_no3red_gate/` (`gate.pjsub`, `gate_check.py`;
decision of the passing run in `tests/isw_no3red_gate/gate_result_g15.txt`): the FABM 0-D driver built twice from
one frozen FABM (f1c4ac1b) / GOTM (4a4754bd) copy with this repository at fc61123 (before the port) and a7ae58a (after),
ifx Release, test case `testcases/fabm-ersem-15.06-L4-noben-docdyn-iop.yaml` on `ersem-setups/0d-aquarium` (one year,
RK4, every variable in double precision; the benthic 15.06 test cases no longer initialise with this fork). Results
(job 7685855): floor (base twice) identical; switch absent and `isw_no3red: 0` both **bit identical** to the base
(365 / 365 variables, 2929 records); `isw_no3red: 1`: `ledger_no3red_O2o` = 2 x `fN3PIn` exactly at every record of
P1-P4, the term active (year-mean 0.06-0.17 mmol O2 m-3 d-1 per producer); one Euler step from the same state (the
year-run state of 2008-09-27 15:00, where all of P1-P4 take up nitrate) under constant forcing: only the state O2_o
differs, and **O2 gained = dt x 2 x nitrate uptake** to round-off (0.0478846870129 mmol O2 m-3, difference 4.6e-15,
round-off bound 2.3e-13).
gpt-6.2-sol review of the port: `doc/port52-review-20261125.md`, 7 rounds, **PASS** in round 7 (job 7686024); the
port itself matched the executed rev52 code in every round; the findings concerned the gate (P3 not exercised, wind,
commit guard, documentation) and were fixed, gate rerun PASS (job 7685855, `tests/isw_no3red_gate/gate_result_g15.txt`).

Not ported, muse-only (exist only in the muse source copies since rev47-rev51, not in this repository): the photon
cap of gross production and its accepted-extent counters (`receipts51`, `accepted_gross51`, `accepted_resp51`, in
rev52 moved to a child model `counters51` so that grazing no longer scales them). Whether they belong in this fork
is an open decision.

### 2. Configuration conventions every user should check (no code change needed)

The defaults below are ERSEM 15.06 values. They are inconsistent with the CH2O convention (one O2 per C for every
aerobic oxidation of organic carbon and for photosynthesis), which the MUSE bed uses; mixing them gives a water column
that consumes 20 % more O2 per C than the bed for the same reaction, and producers whose O2 balance does not close.

| Module / key (units mmol O2 per mg C) | ERSEM default | mol O2 : mol C | CH2O value used in MUSE rev52 |
|---|---|---|---|
| `dom_decay` `ur_O2` | 0.1 | 1.201 | 1/12.011 = 0.0832570144 |
| `bacteria*` `ur_O2` (B1) | 0.1 | 1.201 | 0.0832570144 |
| `primary_producer` `uB1c_O2`, `urB1_O2` (P1) | 0.1 | 1.201 | 0.0832570144 |
| same, P2 | 0.08 | 0.961 | 0.0832570144 |
| same, P3, P4 | 0.11 | 1.321 | 0.0832570144 |
| `microzooplankton` / `mesozooplankton` `urB1_O2` (Z4-Z6) | 0.1 | 1.201 | 0.0832570144 |

Also changed in the MUSE tank configuration (configuration-specific, check your own basis before copying):
`dom_decay` `isw_law: 1` and `K_O2: 3` mmol m-3 (one O2 half-saturation for DOC oxidation in water and bed; OMEXDIA,
Soetaert et al. 1996); `primary_producer` `isw_no3red: 1`; P1, P2 `q10: 2.0` (3.2 had no documented basis; 2.0 is the
ERSEM reference configuration value, Eppley 1972, Bissinger et al. 2008).

### 3. Effect seen in the MUSE tank closures (for orientation only)

With the fixes and the former fitted conditions (no refit), the objective of the eight closures changed by +0.1 % to
+12.6 % (dark closures most, through more O2 left in the water and less re-oxidation); after the refit seven of eight
closures reproduce the observations better than before (refit audit pending). Numbers: muse `rev52/stage2.md`,
`rev52/refit.md`.
