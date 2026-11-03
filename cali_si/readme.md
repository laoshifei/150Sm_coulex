# Si array and S3 source energy calibration

## Saved-peak-center recalibration

The latest all-channel center-based analysis is stored in
`parameters/source_0457_peak_centers_20261102/` and
`parameters/source_0460_peak_centers_20261102/`. It uses the already fitted ADC
response centers to recompute `E_base_keV = c0 + c1 * ADC` for both source runs.
The original local-deviance cutoff is not used to veto an otherwise numerically
usable center; optimizer, covariance, center and boundary checks remain.

```bash
python3 cali_si/recalibrate_from_centers.py --run 0457 --primary PRIMARY.json --window120 WINDOW_120.json --window160 WINDOW_160.json --previous-status PREVIOUS_CHANNELS.csv --output NEW_CENTER_CALIBRATION
```

The program requires an unused output directory and never installs production
parameters. `all_channels.csv` retains all channels and flags;
`candidate_coefficients.csv` contains candidates with at least two numerically
usable centers; `window_checked_coefficients.csv` adds window and multistart
checks. Failed-center all-peak diagnostic lines are explicitly separate.
The 0457 package also contains the complete low-beam wiring remapping.

Of 304 source-channel instances, 279 have primary candidates: 255 pass the
working window checks and 24 need review. The remaining 25 retain diagnostic
curves. All 112 S3 instances pass the new center/window checks; their two-point
zero residual does not validate linearity. Existing response-fit discrepancies
remain unchanged. Notebook 0242 shows every spectrum, calibration curve,
coefficient and failed/review state. Formal `pars` and build files are unchanged.

## Source-response fitting entry

Run `make -C cali_si` from the project root. The C++/ROOT executable is
`cali_si/calibrate`; the original `auto_cali.cpp`, `get_par.cpp`, and
`compute_average.cpp` remain the workflow entry files. Shared numerical code is
in `include/sm150_source_calibration.h`. ROOT Spectrum/Minuit2, OpenSSL and
nlohmann/json are required.

```bash
cali_si/calibrate source --run 460 --input HIST_0460.root --seed pars/cali_si_r460.txt --output NEW_0460
cali_si/calibrate source --run 457 --input HIST_0457.root --seed pars/cali_si_r457_high_current.txt --output NEW_0457
cali_si/calibrate remap-low --input NEW_0457/calibration.json --output NEW_0457_LOW
cali_si/calibrate average --input NEW_0460/calibration.json --keys 10400,10401,10403,10407,10413,10607,10707 --output NEW_ESTIMATES
```

Inputs must be explicit and outputs must not already exist. The source input
is a modern `raw2hist` file (`ch10400/adc_clean`, etc.). Old
`cid1_sid04/h_cid1_sid04_ch00` ADC histograms require `--legacy`.
The thin ROOT macro prints usage with no arguments and does no deletion:

```cpp
.x cali_si/macro_get_par.cpp(460,"HIST_0460.root","pars/cali_si_r460.txt","NEW_0460")
```

## Source model and interpretation

The fitted base curve is `E_keV = c0 + c1*ADC`, with c2=0. Hardware gains,
in-beam attenuation and S3 front/back normalization are separate operations.
Source-fit coefficients must not incorporate those multipliers a second time.

The original source assignment uses three nominal alpha peaks on slots 2–8
and two on slots 9–12. The default configuration now treats these features as
Pu-239, Am-241 and Cm-244 multiplets: evaluated **alpha emission** energies and
fixed relative probabilities, not total transition energies including recoil.
Source identity is inferred from the old nominal lines; a source certificate
and detector/source energy-loss assessment are still absent. References and
printed pages are recorded for each line in `source_config.json`:

- [Pu-239, LNHB, section 4, p. 8, 2008 evaluation](https://www.lnhb.fr/nuclides/Pu-239_tables.pdf)
- [Am-241, LNHB, section 4, p. 7, 2010 evaluation](https://www.lnhb.fr/nuclides/Am-241_tables.pdf)
- [Cm-244, LNHB, section 4, p. 3, 2010 evaluation](https://www.lnhb.fr/nuclides/Cm-244_tables.pdf)

The default `joint_alpha_fit` path is implemented in `joint_source.h`. It fits
all source features in a channel simultaneously, so the tail of an upper
feature also contributes below the next feature. Each multiplet shares a
Gaussian or Gaussian plus left-tail response. Relative
weights and component separations are fixed, using the seed-plus-affine
matching scale to convert separations to ADC. Very weak emissions are omitted.
The full fit extends 140 keV below the first matched feature and above the last;
local residual summaries use separate feature windows. The tail fraction is
allowed in [0,1], with explicit nested endpoint reductions and three initial
fractions. A positive linear background is shared by the entire fit. This model
does not assert that all observed shoulders are nuclear
lines or that source self-absorption/dead-layer effects have been corrected.

The fitter integrates actual bin boundaries, minimizes Poisson deviance,
keeps all fit attempts, and rejects unreliable or shape-bound centers.
Equal-weight linear calibration and the correlated joint-center covariance are saved;
two points determine a line exactly and provide **no independent linearity
test**. The source range around 5–6 MeV cannot certify the in-beam 15–100 MeV
extrapolation. See the Ge README for shared fit/output details.

`joint_response: double_left_tail` selects a two-exponential left response
convolved with a Gaussian for controlled diagnostics. Its extra shape freedom
did not consistently resolve the experimental discrepancies; it is not the
default. `joint_alpha_fit: false` retains the local-fitting interface, including
the configurable response limits and multistart records. Historical results
must use their frozen configuration, not the current defaults.

The `conditional_center_chi2` field retains the legacy diagonal diagnostic;
it is not a goodness-of-fit statistic that includes all center correlations.
The coefficient covariance does propagate the joint center correlations.

## Failure handling, remapping and estimates

`calibration.json` and `parameters.csv` include every requested channel.
`candidate_parameters.txt` contains eligible channels only; it is explicitly
partial when `candidate_text_is_complete` is false. There is no zero-filled
failure row or automatic production installation.

`remap-low` swaps complete 0457 coefficient rows between slots 2/9 and 3/10,
records both electronics keys and leaves hardware gains separate. It is a
wiring remapping, not a second independent source fit.

`average` requires explicit targets and at least three eligible donors on the
same board, excluding all targets. It accepts valid zero intercepts and uses
address keys rather than input row order. Output is a separately labelled
estimate with donor keys/scatter and a null fit chi-square. It never copies a
mean chi-square or substitutes its estimates into fitted rows.

## 2026-11-02 joint-fit candidates

Both full source histograms were refitted for all 304 channel instances.
The primary joint fit gives 84/96 Si and 41/56 S3 candidates for 0457, and
76/96 Si and 44/56 S3 candidates for 0460. Repeating the fits with 120/160 keV
outer windows and checking centers leaves 69 Si + 34 S3 and 65 Si + 44 S3
recommended source candidates, respectively. The 5 keV center-stability
threshold is a working review criterion, not a claimed calibration accuracy.

Local historical packages are `parameters/source_0457_joint_20261102/` and
`parameters/source_0460_joint_20261102/`. `all_channels.csv` retains every
channel; `conditional_parameters.txt` contains primary-fit candidates;
`candidate_parameters.txt` contains candidates passing the additional checks.
These tables remain partial. The 0457 package includes explicitly remapped
low-beam rows; no donor means are substituted for missing measurements.

The even S3 rings (source keys 11200–11211) still have structured Am-241
residuals in both source runs. Si10311 also remains unresolved. All such
spectra and attempts are shown in the appended section of notebook 0238.
Matching-response synthetic recovery and independent bin/count/calibration
checks are separate from failed model-mismatch stress tests, which are kept.
Source-response systematics remain outside the conditional covariance.
No formal parameters or existing build ROOT files were replaced.

### Si10311 three-center calibration check

The local historical `parameters/si10311_three_peak_20261102/` separately stores two partial,
linear candidates calculated from all three fitted Si10311 centers, including
the two peaks rejected by the previous local-deviance rule. Both primary
fits have numerically usable centers. For 0457 the source-center RMS changes
from 30.343 to 5.811 keV; for 0460 it changes from 29.038 to 5.897 keV.
The middle-source residual remains about +8.2 to +8.3 keV.

This changes the ADC-to-energy line, not the fitted ADC counts or response
residuals. It is a source-range candidate with response-model review, rather
than a replacement for the full calibration table. The 0460 alternate-window
fits have shape-bound parameters and are flagged sensitivity diagnostics.
Notebook 0238 section 11 contains the coefficients, before/after spectra,
window comparisons and source provenance. Production parameters are unchanged.

## Historical 2026-10-31 result

Source runs 0457 and 0460 were read in full from E; all 152 mapped particle
channels were examined. The revised fit provides 83 and 78 eligible
**Si-array** candidates respectively, including 5 and 7 two-point candidates.
All 56 S3 channels in each source still lack two reliable anchors: the
approximately 5486 keV feature is inadequately described even though counts
are abundant. The 5805 keV multiplet is generally better described. No complete
S3 replacement table is justified by this pass. Several historical same-board
replacement channels also remain unresolved.

Compact packages live in `parameters/source_*_20261031/`, with explicit
missing/review states. Notebook 0238 contains all channel spectra, parameters
and fit diagnostics; full local evidence lives under
`codex/analysis/source_calibration_modernization_20261031/`. Formal `pars`,
S3 attenuation/front-back parameters and existing event ROOT files are unchanged.
