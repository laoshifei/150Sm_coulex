# Ge source efficiency calibration

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



 
