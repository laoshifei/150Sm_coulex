# Ge source energy calibration and in-beam gain checks

Build with `make -C cali_ge` from the project root (ROOT with Spectrum and
Minuit2, OpenSSL, and nlohmann/json are required). The maintained executable is
`cali_ge/calibrate`. `auto_cali.cpp` and `get_par.cpp` retain the original source
calibration roles; `compute_k.cpp` handles the separate 511 keV gain step.
Shared numerical and input/output code is in
`include/sm150_source_calibration.h`.

The adopted Ge base polynomial and per-run gain tables reside in
[`include/parameter`](../include/parameter/README.md).
`pars/cali_ge.txt` and `pars/run_k/` remain compatible links. The 2026-11-02
relocation preserved every value; the source and gain candidates below remain
unadopted. Runtime selection is maintained in `pars/active.json` under
`ge_source_calibration` and `ge_gain`.

```bash
cali_ge/calibrate source --run 377 --input NEW_HIST.root --seed pars/cali_ge.txt --output NEW_SOURCE_DIR
cali_ge/calibrate gain --run 472 --input INBEAM_HIST.root --seed pars/cali_ge.txt --output NEW_GAIN_DIR
```

Both output directories must be new. `--keys 204,205` selects explicit electronics
addresses. No command installs production parameters, appends to an old table,
deletes old results, or automatically submits all experimental runs.

Source mode reads `ch00204/adc_clean` etc. from `raw2hist` output and checks the
run metadata. For the original `cid0_sid02/h_cid0_sid02_ch04` ADC histograms,
use `--legacy` explicitly. There is no guessed fallback to another input.
The thin ROOT macros also require explicit paths and a compiled executable:

```cpp
.x cali_ge/macro_get_par.cpp(377,"NEW_HIST.root","pars/cali_ge.txt","NEW_SOURCE_DIR")
.x cali_ge/macro_compute_k.cpp(472,"INBEAM_HIST.root","pars/cali_ge.txt","NEW_GAIN_DIR")
```

## Method and outputs

`source_config.json` defines the Ba-133/Eu-152 lines, anchors, search range,
matching limits and quality gates. Peak matching is deterministic and allows
an affine displacement of the seed energy. The seed is prior matching
information, not independent validation.

Local fits integrate a Gaussian or Gaussian plus exponential left tail over
each ADC bin, together with a positive linear background. They minimize raw
Poisson deviance, preserve every attempt and compare AICc only on identical
bins. Integer ADC bin edges are shifted by -0.5. Centers with unreliable
covariance, insufficient significance, or constrained mean/core width/upper
tail shape are excluded as anchors. Quality limits are working diagnostics,
not confidence levels.

The energy curve is `E_keV = c0 + c1*ADC + c2*ADC^2`. Anchors have equal weight;
the saved covariance propagates conditional peak-center errors through the
linear polynomial solve and omits model/source systematic uncertainties.
The tabulated center chi-square is a diagnostic, not the minimized calibration
objective. Each curve records its actual anchor interval; extrapolation is not
validated by monotonicity alone.

- `calibration.json`: all requested channels, matching, fits, counts, residuals,
  covariance, status and provenance, including failures.
- `parameters.csv`: full channel status and available coefficients.
- `candidate_parameters.txt`: **only eligible rows**, in the historical seven
  column format. It may be incomplete and must not replace a complete `pars`
  table. Check `completion.json` and `candidate_text_is_complete`.
- `gain.json`, `gain.csv`, `k_candidate.txt`: corresponding gain diagnostics.

Gain mode requires modern decoded-input metadata and an exact base-table hash.
It fits **energy_base**, which has no previous run-dependent k, and computes
`k = 511 / fitted_center`. It never fits `energy_final` to infer a total k.
A new source polynomial and a k derived for the old polynomial are not a
compatible calibration chain without reprocessing.

## 2026-10-31 validation and candidates

Full source 0377 was read anew; the previously verified full 0722 ADC histogram
was reused. All 37 Ge addresses in the build mapping were checked, including
addresses outside the 20-channel efficiency selection. This is not a claim
that all 37 have acceptable source response. Final candidates retain 29/37 and
30/37 respectively; other channels remain visible for review. Three existing
10-million-event segmented histograms (0472, 0617, 0680) give only representative
511 gain checks, not a new 173-run gain calibration.

Compact candidate packages are in `parameters/source_*_20261031/` and
`parameters/gain_checks_20261031/`. They are not installed. Notebook 0238 and
the full fits/validation are maintained locally under `jupyter/` and
`codex/analysis/source_calibration_modernization_20261031/`.
The source energy snapshot used by efficiency notebook 0237 remains frozen in
`ge_eff/parameters/source_0722_20261030/`; this work does not relabel or replace it.
