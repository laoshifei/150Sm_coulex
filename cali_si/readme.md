# Si array and S3 source energy calibration

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

Each multiplet shares a Gaussian or Gaussian plus left-tail response. Relative
weights and component separations are fixed, using the seed-plus-affine
matching scale to convert separations to ADC. Very weak emissions are omitted.
The local half window is at most 140 keV and is clipped to avoid adjacent source
features. This model does not assert that all observed shoulders are nuclear
lines or that source self-absorption/dead-layer effects have been corrected.

The fitter integrates actual bin boundaries, minimizes Poisson deviance,
keeps all fit attempts, and rejects unreliable or shape-bound centers.
Equal-weight linear calibration and conditional center covariance are saved;
two points determine a line exactly and provide **no independent linearity
test**. The source range around 5–6 MeV cannot certify the in-beam 15–100 MeV
extrapolation. See the Ge README for shared fit/output details.

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

## 2026-10-31 result

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
