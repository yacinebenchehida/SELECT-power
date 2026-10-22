# Related work for SPOT (Selection Power Over Time)

Literature scan, October 2026. Citations and DOIs were checked against PubMed unless marked otherwise. This is a working list for the paper's introduction and discussion, not a final reference list.

## Bottom line

No published tool provides an interactive, installation-free power calculator for detecting selection at a candidate locus from allele frequencies sampled at arbitrary time points, with the test threshold calibrated by a drift-aware neutral null. Existing design work is either:

- **(a)** lab evolve-and-resequence (E&R) studies with pool-seq, done with command-line forward simulators (MimicrEE2, forqs, SLiM frameworks), where users must build their own test and null; or
- **(b)** one-off simulation studies tied to specific scenarios (human aDNA, epidemics), which assume perfect genotypes and are diploid only.

Analytical theory exists (Lynch & Ho 2020) but has no interface and ignores genotype-calling error. The many inference methods analyse data that already exist; none offers sample-size planning.

SPOT's niche is a front-end planning tool for museum, ancient DNA and monitoring designs. It covers haploid, diploid and haplodiploid organisms, models low-coverage genotype calls, uses an empirical drift-aware null, and shows directly how much a naive test is inflated.

## 1. Power and study design for temporal selection

| Reference | What it does | Difference from SPOT |
|---|---|---|
| Kofler & Schlötterer 2014, *Mol Biol Evol* 31:474. [doi:10.1093/molbev/mst221](https://doi.org/10.1093/molbev/mst221) | Genome-wide forward simulations comparing E&R designs (replicates, population size, duration, LD) for power and false positives | Lab E&R, pool-seq, command line; diploid; no genotype calling |
| Baldwin-Brown, Long & Thornton 2014, *Mol Biol Evol* 31:1040. [doi:10.1093/molbev/msu048](https://doi.org/10.1093/molbev/msu048) | Simulation power to map QTL in E&R as a function of population size, replicates, generations and coverage | QTL in lab populations; a paper, not a tool |
| Kessner & Novembre 2015, *Genetics* 199:991. [doi:10.1534/genetics.115.175075](https://doi.org/10.1534/genetics.115.175075) | forqs simulator, applied to power analysis of artificial selection / E&R | Command-line simulator; user builds the pipeline |
| Vlachos & Kofler 2018, *PLoS Comput Biol* 14:e1006413. [doi:10.1371/journal.pcbi.1006413](https://doi.org/10.1371/journal.pcbi.1006413) | MimicrEE2: genome-wide forward simulator of E&R, intended to support design | No built-in test, null or power; lab focus |
| Vlachos & Kofler 2019, *Mol Biol Evol* 36:2890. [doi:10.1093/molbev/msz183](https://doi.org/10.1093/molbev/msz183) | Optimising power in E&R: selection regimes | Concerns the selection regime, not sampling of natural populations |
| Vlachos et al. 2019, *Genome Biol* 20:169. [doi:10.1186/s13059-019-1770-8](https://doi.org/10.1186/s13059-019-1770-8) | Benchmark of E&R selection tests (CMH, χ², LLS, CLEAR, WFABC, …) | Compares methods; not a design tool |
| Lou, Therkildsen & Messer 2020, *G3* 10:3213. [doi:10.1534/g3.120.401287](https://doi.org/10.1534/g3.120.401287) | SLiM framework for QTL detection power in artificial selection | Polygenic E&R; scripted |
| Topa et al. 2015, *Bioinformatics* 31:1762. [doi:10.1093/bioinformatics/btv014](https://doi.org/10.1093/bioinformatics/btv014) | Gaussian-process test (BBGP), including simulations of design choices (time points, depth) | Mainly inference; E&R pool-seq |
| Shim et al. 2016, *G3* 6:893. [doi:10.1534/g3.115.023200](https://doi.org/10.1534/g3.115.023200) | CP-WFABC; simulations suggest at least four time points are needed to detect a change in selection | Design advice is a by-product of inference |
| Lynch & Ho 2020, *Genome Biol Evol* 12:443. [doi:10.1093/gbe/evaa056](https://doi.org/10.1093/gbe/evaa056) | Analytical sampling variances of selection and Ne estimators for temporal designs | Theory only; no simulator or UI; perfect genotypes. **The closest theoretical treatment of how many time points to take and when.** |
| Anchieri et al. 2026, *Genome Biol Evol* 18:evaf234. [doi:10.1093/gbe/evaf234](https://doi.org/10.1093/gbe/evaf234) (verified on PMC) | SLiM study of realistic aDNA sampling strategies for four inference methods | **The closest aDNA design study**; scripts only, human-specific, diploid, no planner |
| Santander & Moltke 2025, *Nat Commun* 16:6814. [doi:10.1038/s41467-025-61574-8](https://doi.org/10.1038/s41467-025-61574-8) (verified on publisher page) | SimOutbreakSelection: simulation pipeline to optimise sampling for epidemic-driven selection | **The closest "design tool"**; command line; epidemic scenario; ranking-based power; perfect genotypes |
| Clark et al. 2023, *Mol Ecol Resour* 25:e13789. [doi:10.1111/1755-0998.13789](https://doi.org/10.1111/1755-0998.13789) | Review of temporal genomics, with qualitative design guidance | Qualitative only |
| Jensen & Leigh 2022, *Ecol Evol* 12:e9340. [doi:10.1002/ece3.9340](https://doi.org/10.1002/ece3.9340) | Review of temporal genomics for climate-change responses | Qualitative |
| Dehasque et al. 2020, *Evol Lett* 4:94. [doi:10.1002/evl3.165](https://doi.org/10.1002/evl3.165) | Perspective on inferring selection from aDNA; limitations | No power quantification |
| Malaspinas 2016, *Mol Ecol* 25:24. [doi:10.1111/mec.13492](https://doi.org/10.1111/mec.13492) | Review of time-series sweep methods from an aDNA perspective | Review |
| Waples et al. 2025, *Mol Ecol Resour* 25:e14057. [doi:10.1111/1755-0998.14057](https://doi.org/10.1111/1755-0998.14057) | MaxTemp: maximising precision of the temporal Ne estimate | Ne precision, not selection power |

No web or interactive tool for power to detect selection was found; the only interactive drift–selection simulators found are teaching tools.

## 2. Inference methods for selection from time series

All of these are inference software or methods, not design tools.

- Bollback, York & Nielsen 2008, *Genetics* 179:497. [doi:10.1534/genetics.107.085019](https://doi.org/10.1534/genetics.107.085019) — diffusion HMM.
- Malaspinas et al. 2012, *Genetics* 192:599. [doi:10.1534/genetics.112.140939](https://doi.org/10.1534/genetics.112.140939) — HMM; allele age and s.
- Mathieson & McVean 2013, *Genetics* 193:973. [doi:10.1534/genetics.112.147611](https://doi.org/10.1534/genetics.112.147611) — HMM/EM; spatial structure.
- Steinrücken, Bhaskar & Song 2014, *Ann Appl Stat* 8:2203. [doi:10.1214/14-AOAS764](https://doi.org/10.1214/14-AOAS764) — spectral HMM; general diploid selection.
- Feder, Kryazhimskiy & Plotkin 2014, *Genetics* 196:509. [doi:10.1534/genetics.113.158220](https://doi.org/10.1534/genetics.113.158220) — frequency increment test (FIT); discusses how the sampling interval affects power.
- Foll et al. 2014, *PLoS Genet* 10:e1004185. [doi:10.1371/journal.pgen.1004185](https://doi.org/10.1371/journal.pgen.1004185) — time-sampled influenza.
- Foll, Shim & Jensen 2015, *Mol Ecol Resour* 15:87. [doi:10.1111/1755-0998.12280](https://doi.org/10.1111/1755-0998.12280) — WFABC.
- Terhorst, Schlötterer & Song 2015, *PLoS Genet* 11:e1005069. [doi:10.1371/journal.pgen.1005069](https://doi.org/10.1371/journal.pgen.1005069) — multi-locus E&R.
- Ferrer-Admetlla et al. 2016, *Genetics* 203:831. [doi:10.1534/genetics.115.184598](https://doi.org/10.1534/genetics.115.184598) — ApproxWF.
- Schraiber, Evans & Slatkin 2016, *Genetics* 203:493. [doi:10.1534/genetics.116.187278](https://doi.org/10.1534/genetics.116.187278) — Bayesian path augmentation.
- Gompert 2016, *Mol Ecol* 25:121. [doi:10.1111/mec.13323](https://doi.org/10.1111/mec.13323) — fluctuating selection.
- Taus, Futschik & Schlötterer 2017, *Mol Biol Evol* 34:3023. [doi:10.1093/molbev/msx225](https://doi.org/10.1093/molbev/msx225) — poolSeq R package.
- Iranmehr et al. 2017, *Genetics* 206:1011. [doi:10.1534/genetics.116.197566](https://doi.org/10.1534/genetics.116.197566) — CLEAR.
- Paris, Servin & Boitard 2019, *G3* 9:4073. [doi:10.1534/g3.119.400778](https://doi.org/10.1534/g3.119.400778) — comparison of Wright–Fisher approximations.
- He et al. 2020, *Genetics* 216:521. [doi:10.1534/genetics.120.303463](https://doi.org/10.1534/genetics.120.303463) — two linked loci.
- He et al. 2023, *Mol Biol Evol* 40:msad008. [doi:10.1093/molbev/msad008](https://doi.org/10.1093/molbev/msad008) — aDNA genotype likelihoods; time-varying selection.
- He et al. 2023, *Mol Ecol Resour* 23:1226. [doi:10.1111/1755-0998.13790](https://doi.org/10.1111/1755-0998.13790) — linkage and epistasis.
- Mathieson 2020, bioRxiv. [doi:10.1101/2020.11.17.387761](https://doi.org/10.1101/2020.11.17.387761) — bmws; time-varying s.
- Mathieson & Terhorst 2022, *Genome Res* 32:2057. [doi:10.1101/gr.276862.122](https://doi.org/10.1101/gr.276862.122) — selection in Bronze Age Britain.
- Cheng & Steinrücken 2025, *G3* 15:jkaf123. [doi:10.1093/g3journal/jkaf123](https://doi.org/10.1093/g3journal/jkaf123) — diplo-locus. **The closest software building blocks** (inference plus simulator), but command line, diploid, no power module.
- Fine & Steinrücken 2025, *PLoS Genet* 21:e1011769. [doi:10.1371/journal.pgen.1011769](https://doi.org/10.1371/journal.pgen.1011769) — EM; general diploid selection.
- Whitehouse & Schrider 2023, *Genetics* 224:iyad084. [doi:10.1093/genetics/iyad084](https://doi.org/10.1093/genetics/iyad084) — Timesweeper (CNN).
- Buffalo & Coop 2019, *Genetics* 213:1007. [doi:10.1534/genetics.119.302581](https://doi.org/10.1534/genetics.119.302581) — temporal covariance, genome-wide.
- Buffalo & Coop 2020, *PNAS* 117:20672. [doi:10.1073/pnas.1919039117](https://doi.org/10.1073/pnas.1919039117) — temporal covariance, genome-wide.
- Akbari et al. 2024, bioRxiv. [doi:10.1101/2024.09.14.613021](https://doi.org/10.1101/2024.09.14.613021) — large human aDNA scan.

## 3. Drift-inflated false positives and low-coverage error

These support SPOT's drift-aware null and its "naive false-positive rate" output.

- Waples 1989, *Genetics* 121:379. [doi:10.1093/genetics/121.2.379](https://doi.org/10.1093/genetics/121.2.379) — temporal method; separates sampling variance from drift variance.
- Jorde & Ryman 2007, *Genetics* 177:927. [doi:10.1534/genetics.107.075481](https://doi.org/10.1534/genetics.107.075481) — unbiased temporal F.
- Wang 2001, *Genet Res* 78:243. [doi:10.1017/S0016672301005286](https://doi.org/10.1017/S0016672301005286) — pseudo-likelihood temporal Ne.
- Jónás et al. 2016, *Genetics* 204:723. [doi:10.1534/genetics.116.191197](https://doi.org/10.1534/genetics.116.191197) — temporal Ne with pool-seq two-step sampling.
- Spitzer, Pelizzola & Futschik 2020, *Ann Appl Stat* 14:202. [doi:10.1214/19-AOAS1301](https://doi.org/10.1214/19-AOAS1301) (verified via institutional record and arXiv) — χ²/CMH anti-conservative under drift; adjusted tests.
- Wiberg et al. 2017, *Methods Ecol Evol* 8:1899. [doi:10.1111/2041-210X.12810](https://doi.org/10.1111/2041-210X.12810) — CMH false positives; quasibinomial GLMs.
- Han, Sinsheimer & Novembre 2014, *Mol Biol Evol* 31:723. [doi:10.1093/molbev/mst229](https://doi.org/10.1093/molbev/mst229) — bias from low-coverage genotype calls.
- Barton et al. 2025, *Nature* (Matters Arising). [doi:10.1038/s41586-024-08496-5](https://doi.org/10.1038/s41586-024-08496-5) — the Black Death selection signals (Klunk et al. 2022, [doi:10.1038/s41586-022-05349-x](https://doi.org/10.1038/s41586-022-05349-x)) are consistent with drift plus sampling and coverage differences; reply by Vilgalys et al. 2025, [doi:10.1038/s41586-024-08497-4](https://doi.org/10.1038/s41586-024-08497-4). The author list and pages of Barton et al. are not fully verified. **A high-profile example of the problem SPOT makes visible.**

## Not verified

- An R package named "ACER" implementing the Spitzer et al. tests.
- A dedicated study of genotype-error effects on false positives in museum or aDNA temporal selection scans. None was found; He et al. 2023, Anchieri et al. 2026 and Barton et al. 2025 come closest.
