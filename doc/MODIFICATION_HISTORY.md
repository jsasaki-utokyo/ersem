# Modification history of this ERSEM fork (entry point)

This file is the entry point for users of this fork (other projects, other sessions): what was changed, why, how
it is switched, and what a configuration must set to use it. The newest entry is first. Per-module notes stay in
`doc/*.md`; this file links them.

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
| `src/primary_producer.F90` | `isw_no3red` (0 = former) | nitrate assimilated by P1-P4 returns 2 O2 per N (nitrate to amine is 8 e-): closes the electron balance of the producers under the CH2O convention (section 2); carbon, N and TA untouched; diagnostic `ledger_no3red_O2o` registered only when on | **pending port** (executed in the muse rev52 source copy `~/Github/muse/docs/carb_s10m/rev52/source/ersem/src/primary_producer.F90`) | rev52 design 2.4, M2/M3 |
| `src/carbonate_engine.F90` | none (new public function) | `carbonate_engine_kw(T, S, Pbar, opt_pH_scale, opt_k_carbonic, opt_total_borate)`: K_W on the pH scale of the solve (used by the MUSE bed; fixes K_W on the SWS scale used with a total-scale pH) | **pending port** (same source copy) | rev52 triage T06 |

The port of the two pending items into this repository, with the same default-off behaviour and a bit-identity gate
against the current head, is scheduled after the MUSE rev52 refit is adopted; this entry will then name the commits.

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
