# Ge source efficiency calibration

## Selected scope: all-Ge only

For subsequent combined-Ge analysis, use the retained original all-Ge
coefficients recorded in [parameters/all_ge_reference.json](parameters/all_ge_reference.json):
`A=14.8, B=4.618, C=0, D=10.8, E=-0.7352, F=-0.1101, G=15`.
Energy is in keV and logarithms are natural. The yield-program normalization is
`f(E)/max(f, 50--2048 keV)`, not the 334-keV normalization used in comparison plots.
This records the user's group choice without adopting a refitted candidate or
changing the existing yield macros. Separate ring curves remain historical
comparisons and are not required for subsequent all-Ge-only use.

## Radware refit, 2026-11-02

The original Radware function has now been refitted with `C=0` and `G=15`,
matching the legacy macro. `refit_radware.py` reads saved source peak areas;
it does not decode events, repeat addback, or change production parameters.
The full frozen results remain local in `parameters/radware_0722_20261102/`.
The repository includes its `efficiency_points.csv` as a reproducible input;
publishing the fitting code does not adopt the fitted candidate coefficients.

The package includes all/ring3/ring4/ring5, free relative Ba/Eu normalization
and the original 356/344 matching as a conditional comparison, plus fits with
different low-energy thresholds. All starts and limits are retained in the
local analysis record. The main display uses 15 points, including an 81-keV
refit on a wider window. That peak region remains model dependent and ring4
has substantial residuals: this package is a **candidate, not an installed
replacement**. The source addback uses coarse timestamps and is not identical
to the production build's timing/grouping.

Run from the project root, with a new output directory:

```bash
python3 ge_eff/refit_radware.py --points ge_eff/parameters/radware_0722_20261102/efficiency_points.csv --output NEW_DIR --mode clover2_ab100 --starts 8
```

Notebook 0237 section 11 contains the new curves, coefficients, source-point
residuals and old/new comparisons. Sections 1-9 retain the earlier polynomial
candidate; section 10 shows the original saved Radware curves. Neither package
is automatically selected for production.

## Current source-data workflow

The explicit workflow reads the decoded 0722 Ba-133 / Eu-152 source run. It
produces candidate relative efficiencies for all 20 selected electronics
channels and ring 3 / 4 / 5, with the original counts and fit diagnostics retained.
The candidate package is `parameters/source_0722_20261030/`. It is **not installed
as the production efficiency**. In particular, the ring-4 fit and the strongly
energy-dependent quality selection of Ge304 need further review.

* `source_lines.json`: evaluated photon probabilities, uncertainties, source
  references, peak ranges and neighboring features.
* `recalibrate.py`: bin-integrated Poisson peak models and source energy
  calibration. `--energy-model left`, `gaussian`, or `both` selects the energy
  anchor model explicitly. It does not select the efficiency-area model.
* `source_spectra.cpp`: full C++ source-event reading, source energy spectra and
  explicit Clover addback diagnostics. Build with `make -C ge_eff`.
* `fit_efficiency.py`: compare Gaussian / left-tail / two-sided-tail peak fits
  on identical bins; retain all attempts; fit two free source normalizations
  and the relative-efficiency shape. NumPy, SciPy, pandas and uproot are needed.

Run from the project root, replacing `NEW_DIR` with a new output directory:

```bash
make -C ge_eff
raw2hist/histogram --run 722 --input /mnt/e/Sm150/rootfile/data_C1_0722.root --output NEW_DIR/source_adc.root --mode raw --detector ge
python3 ge_eff/recalibrate.py energy --input NEW_DIR/source_adc.root --seed pars/cali_ge.txt --output NEW_DIR/energy_fit --energy-model left
ge_eff/source_spectra /mnt/e/Sm150/rootfile/data_C1_0722.root NEW_DIR/energy_fit/energy_calibration.json NEW_DIR/source_energy.root
python3 ge_eff/fit_efficiency.py --adc NEW_DIR/source_adc.root --energy NEW_DIR/source_energy.root --seed pars/cali_ge.txt --output NEW_DIR/efficiency_fit
```

To reproduce the published candidate's exact energy scale, pass the frozen
`ge_eff/parameters/source_0722_20261030/source_energy.json` to `source_spectra`
instead of recalibrating the energy scale. The first-stage code and line-data
snapshots are retained in the local analysis record. All tools refuse existing
scientific output destinations.

The source clean selection excludes pileup, out-of-range flags and ADC=65535.
Normal and failed CFD are combined; there is no particle coincidence, random
subtraction, in-beam gain correction or Doppler correction. Source energy
spectra have 1-keV bins. Single-crystal peak fits use the original 2-ADC bins
and map their physical edges, without resampling their counts.

Clover 212--215 diagnostics use native coarse time `8*ts`, E>20 keV per crystal,
and forward cluster spans of 50 / 100 / 200 ns. Only single crystals or two
distinct adjacent crystals are retained; opposite, repeated and higher
multiplicity combinations are rejected. These source clusters do not reproduce
the production build's all-Ge grouping and corrected-time selection exactly.
The ROOT `dt_*` histograms have **ns** on the horizontal axis (8-ns bins); the
first output's generic ROOT axis title says energy, so plots explicitly supply
the correct time units.

The primary relative-efficiency model is
`exp(c1*x + c2*x*x)`, with `x=ln(E_keV/334)` and `epsilon(334 keV)=1`.
Ba and Eu normalizations are fitted separately. The supported source-line range
is 244.6974--1408.013 keV; 81 / 122-keV diagnostics are not primary curve inputs.
The 963.390 / 964.079-keV unresolved pair is integrated jointly. The source
activity, live time and geometry have not been certified, and true coincidence
summing is not corrected. Reported covariance is conditional and does not
include model selection, all probability correlations or a complete systematic
budget. Curve differences or successful numerical fits alone do not establish
that the candidate should replace the previous efficiency.

The local reading notebook is
`jupyter/codex_0237_150sm_ge_efficiency_recalibration.ipynb`; detailed counts,
fits, validation and source hashes are under
`codex/analysis/ge_eff_recalibration_20261030/`. These local research artifacts
are separate from the domain code and parameter package.

## Historical-method recovery and quantization display (2026-11-02)

The saved original efficiency remains the production choice. The independent
local package `parameters/legacy_recovery_0722_20261102/` contains recovered peak
areas, all seven Radware coefficients, source configuration, reconstructed fit
windows and complete fit diagnostics. It is **not an exact reproduction** of
the missing historical source ROOT, source calibration, fit-window file or
random seed. The original individual peak-area table is unavailable.

`recover_legacy_spectra.cpp` reads the full decoded source file and reconstructs
the known legacy all-Ge grouping and Clover rules. The source-only acceptance
is explicit because the preserved in-beam builder requires particle hits. The
primary recovery retains flagged events and the historical +/-0.5-keV energy
randomization. `reproduce_legacy.py` runs instrumented copies of the original
ROOT macros, preserving the peak model, source intensities, 356/344 source
normalization, sqrt(area) errors and fixed Radware C=0, G=15. Optional
`--stable-tail --retry-status` uses an algebraically equivalent stable tail and
two fixed extra starts only for peaks with nonzero optimizer status. All
attempts and failed statuses are retained.

From the project root, with new output paths and ROOT available:

```bash
make -C ge_eff
ge_eff/recover_legacy_spectra /mnt/e/Sm150/rootfile/data_C1_0722.root ge_eff/parameters/legacy_recovery_0722_20261102/source_config.json NEW_DIR/recovered.root
python3 ge_eff/reproduce_legacy.py --spectra NEW_DIR/recovered.root --ranges ge_eff/parameters/legacy_recovery_0722_20261102/fit_range_reconstructed.txt --output NEW_DIR/fits --stable-tail --retry-status
```

The recovered 57 peak areas include two nonzero fit statuses and substantial
peak-model residuals. Historical sqrt(area) errors do not describe the full
background, line-shape and common-normalization uncertainty. Similar efficiency
curves do not establish exact peak-area reproduction or better calibration.

For source-spectrum display, `source_spectra` now supports an optional fixed-seed
ADC dither. It adds Uniform(-0.5,0.5) **ADC codes before calibration**, separately
for each crystal and before addback. This differs from the historical
energy-unit jitter used in the recovery above. Threshold membership uses the
undithered energy; counts, quality selection, timing and cluster membership are
preserved. The default invocation without this option retains its old behavior.

```bash
ge_eff/source_spectra /mnt/e/Sm150/rootfile/data_C1_0722.root ge_eff/parameters/source_0722_20261030/source_energy.json NEW_DIR/source_adc_dither.root --adc-dither-seed 722
```

This removes integer-code rebinning dips under a uniform within-code assumption;
it does not repair bad ADC codes or recover measured sub-code information.
Notebook 0237 figures 1 and 3--7 use the new display histograms, while existing
fit inputs and results are preserved. Notebook 0241 provides the before/after
spectra, recovered peak fits, complete area/parameter tables and limitations.

## Legacy macros and original notes

The two original ROOT macros below are preserved as historical implementations.
Their built source ROOT, `fit_range.txt` and source-specific energy calibration
are missing at the legacy paths in this workspace. Do not load them as a way
to inspect the new candidate; they open input files globally. The historical
intensities below are provenance, not the evaluated data used by the new code.

Efficiency calibration using 133Ba and 152Eu $\gamma$ source.

## energy and intensity

### 133Ba

| energy [keV] | intensity |
| --- | --- |
| 81.00 |34.06 | 
| 302.85|18.33 |
|356.02|62.05|
|383.85|8.94|

### 152Eu

| energy [keV] | intensity |
| --- | --- |
|121.78 |28.58| 
|244.70|7.58|
|344.28|26.5|
|411.12|2.234|
|443.96|3.148|
|778.90|12.94|
|867.37|4.245|
|964.08|14.60|
|1112.1|13.64|
|1408.0|21.00|

## info

+ input file: data0722_build_200ns_jump_300ns.root, fit_range.txt
+ output file: terminal out of some certain peak effciency and error, some figure with fit results.

### How to run

#### eff_cali_ring_ge.cpp

```
> root
> .L eff_cali_ring_ge.cpp
> eff_cali_ring_ge(3,2) // (3,2) is one example, need to change
```

#### eff_cali_all_ge.cpp

```
> root
> .L eff_cali_all_ge.cpp
> eff_cali_all_ge()
```

### details

+ `eff_cali_ring_ge.cpp` and `eff_cali_all_ge.cpp` code.
+ In each `cpp` file, including the folowing functions.
  + Double_t gpeakexregion(Double_t *v, Double_t *par), Gaussian + Left-Skew + Step + Background.
  + Double_t bgstep(Double_t *v, Double_t *par), Step + Pol1.
  + Double_t gpeak(Double_t *v, Double_t *par), Gaussian + Left-Skew.
  + double eff_fun(double *x, double *p), radware Ge detector efficiency function.
  + double *eff_fit(map<double, double> &m,  map<double, double> &m_err), fitting.
  + void get_ring_ge_h(int r), get histgoram of certain ge detector.
  + double eff_cali_ring_ge_single_peak(TH1D *h, double energy, double l_bg, double l_peak, double u_peak, double u_bg, double ll_ex, double lr_ex, double hl_ex, double hr_ex), get peak yield used for fitting.
  + void eff_cali_ring_ge(int r), fit certain ge detector.



 
