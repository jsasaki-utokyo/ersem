# ERSEM Sulfur-Process Review (Literature-Based)

Date: 2026-02-15
Target repository: `ersem`
Scope: Review of current sulfur implementation against the stated objective: approximate sediment redox succession (oxic degradation -> denitrification zone -> sulfate reduction zone) under a non-vertically-resolved benthic architecture.

## Findings (ordered by severity)

1. High: Example configuration in `testcases/sulfur_cycle_example.yaml` is internally inconsistent with `pelagic_base` capabilities.
- Evidence in code: `testcases/sulfur_cycle_example.yaml:29` uses `composition: e` with `model: ersem/pelagic_base`, but `src/pelagic_base.F90:145` only adds constituent `h`; `src/pelagic_base.F90:214` falls through unknown constituents to fatal error.
- Impact: the provided sulfur example is likely not directly runnable as written, which blocks reproducible validation of the sulfur implementation.

2. High: Sulfate-reduction alkalinity source is missing in Layer 3 (explicit TODO), biasing carbonate chemistry under sulfidic conditions.
- Evidence in code: sulfate reduction source is applied at `src/benthic_sulfur_cycle.F90:315`; alkalinity addition is explicitly missing at `src/benthic_sulfur_cycle.F90:461`.
- Why it matters: sulfate reduction is an alkalinity-generating process; omitting it can underpredict benthic TA and pH buffering.
- Literature basis: alkalinity generation during sulfate reduction is established in sediment biogeochemistry (e.g., Johnston et al., 2014 GCA; Norbisrath et al., 2023 BG).
- Inference note: the +2 TA per mol sulfate-reduced relation is inferred from standard stoichiometric balancing.

3. High: Nitrate-coupled sulfide oxidation is hard-wired to denitrification (N2 production), with no DNRA branch.
- Evidence in code: `src/benthic_sulfur_cycle.F90:441` consumes `H2S_2`; `src/benthic_sulfur_cycle.F90:442` consumes NO3; `src/benthic_sulfur_cycle.F90:444` routes product N exclusively to `G4n` (N2-equivalent), and no NH4 pathway is represented.
- Why it matters: modern sediment studies report that sulfur-nitrogen coupling can shift toward DNRA under oscillatory redox or sulfidic conditions, retaining fixed N as ammonium.
- Literature basis: Bourceau et al. (2023) showed co-occurring sulfate and nitrate respiration with strong sulfate-reduction/DNRA linkage; Jones et al. (2017) and related estuarine studies report sulfide-stimulated DNRA.
- Model consequence: likely overestimation of N loss (as N2) and underestimation of NH4 retention in environments with intermittent nitrate supply.

4. Medium: FeS precipitation is represented as an effectively unlimited first-order sink, without reactive-iron state limitation.
- Evidence in code: fixed-rate removal in layers at `src/benthic_sulfur_cycle.F90:416`, `src/benthic_sulfur_cycle.F90:420`, `src/benthic_sulfur_cycle.F90:424`, plus pelagic scavenging at `src/benthic_sulfur_cycle.F90:431`.
- Why it matters: FeS/pyrite burial is a real sulfur sink, but capacity is controlled by iron mineral availability and redox history.
- Literature basis: marine sulfur-cycle reviews describe only a fraction of sulfide being trapped as FeS/FeS2, with strong coupling to iron cycling (Jorgensen et al., 2019 Frontiers).

5. Medium: Active sulfur pools are not mass-conservative by design (no explicit sulfate pool; direct sinks).
- Evidence in code: pelagic S0 sink removes mass at `src/sulfur_cycle.F90:147` and `src/sulfur_cycle.F90:161`; benthic S0 burial removes from active pools at `src/benthic_sulfur_cycle.F90:336` and `src/benthic_sulfur_cycle.F90:454`.
- Why it matters: acceptable as a simplification for phenomenology, but it prevents sulfur budget closure and can hide compensation errors during tuning.

6. Medium: Strict vertical separation is implemented, while field systems can show overlap and rapid switching.
- Evidence in code: Layer 3 sulfate reduction is always active by construction (`src/benthic_sulfur_cycle.F90:310`), and the module assumes fixed zonation logic (`src/benthic_sulfur_cycle.F90:17`).
- Literature basis: canonical zonation is valid as a first-order pattern (Jorgensen, 1982), but permeable/oscillatory sediments can exhibit simultaneous sulfate and nitrate respiration and pathway switching (Bourceau et al., 2023).
- Model consequence: good first-order seasonal behavior, but limited fidelity during transient redox oscillations.

7. Low: Layer-2 elemental sulfur (`S0_2`) has production but no explicit in-layer consumption term.
- Evidence in code: production at `src/benthic_sulfur_cycle.F90:443`; no explicit `S0_2` sink in this module.
- Interpretation: transport can still redistribute/remove `S0_2` through the benthic dissolved-matter framework, but local oxidation/disproportionation chemistry is not explicitly represented here.
- Literature basis: elemental sulfur often forms a distinct shallow sediment intermediate pool and is actively transformed (e.g., stratification observations and sulfur-disproportionation studies).

## What is scientifically well aligned with your objective

- The implemented reaction placement is directionally correct for your target process chain:
  - Oxic sulfide oxidation near the top (`src/benthic_sulfur_cycle.F90:327`).
  - Nitrate-coupled sulfide oxidation in the suboxic zone (`src/benthic_sulfur_cycle.F90:339`).
  - Sulfate-reduction source in deeper anoxic layer (`src/benthic_sulfur_cycle.F90:310`).
- The architecture leverages ERSEM's existing 3-layer benthic transport machinery (`src/benthic_column_dissolved_matter.F90:282`), which is a pragmatic choice given the stated non-multilayer constraint.

## Overall judgment

The current sulfur implementation is a reasonable first-order approximation of the intended vertical redox succession in a 3-layer ERSEM framework. It is suitable for qualitative/phenomenological studies of sulfide emergence and oxidation barriers. However, for quantitative sediment biogeochemistry (especially nitrogen partitioning and carbonate response), the present formulation has material structural biases: missing sulfate-reduction alkalinity, no DNRA branch, and unbounded FeS sink parameterization.

## References used

- Jorgensen, B. B. (1982). Mineralization of organic matter in the sea bed — the role of sulfate reduction. Nature 296, 643-645. https://www.nature.com/articles/296643a0
- Jorgensen, B. B., Findlay, A. J., & Pellerin, A. (2019). The Biogeochemical Sulfur Cycle of Marine Sediments. Frontiers in Microbiology, 10:849. https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.00849/full
- Bourceau, O. M., et al. (2023). Simultaneous sulfate and nitrate reduction in coastal sediments. ISME Communications. https://pmc.ncbi.nlm.nih.gov/articles/PMC9992702/
- Jones, Z. L., et al. (2017). Sulfide-Induced Dissimilatory Nitrate Reduction to Ammonium Supports Anammox. Appl. Environ. Microbiol. https://pmc.ncbi.nlm.nih.gov/articles/PMC5514666/
- Finster, K., Liesack, W., & Thamdrup, B. (1998). Elemental sulfur and thiosulfate disproportionation by Desulfocapsa sulfoexigens. Appl. Environ. Microbiol. https://pmc.ncbi.nlm.nih.gov/articles/PMC124681/
- Norbisrath, M., et al. (2023). Alkalinity and nitrate dynamics reveal dominance of anammox in a hyper-turbid estuary. Biogeosciences. https://bg.copernicus.org/articles/20/4307/2023/
- Johnston, S. G., et al. (2014). Alkalinity capture during microbial sulfate reduction and implications for acidification. Geochimica et Cosmochimica Acta. https://doi.org/10.1016/j.gca.2014.01.009
- Sulpis, O., et al. (2022). RADIv1: a non-steady-state early diagenetic model. Geosci. Model Dev. https://gmd.copernicus.org/articles/15/2105/2022/
