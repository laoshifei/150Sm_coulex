# Current build interface (2026-10-28)

Use `parameters/angle_table_baseline.json` for the retained baseline: original
Ge/particle angles, Si exit-energy factor 1.0, and S3 factor 0.5. The beta and
per-Ge angle trials did not establish a consistent overall improvement that
warrants replacing this baseline. The executable uses the file explicitly
supplied through `--parameters`; no trial is activated automatically.

The local-only `parameters/summary_20261029/` directory contains the parameter
reference and historical five-/seven-Ge trial exports. Those research exports
are not required by the executable and are not part of this source commit.

Compilation requires ROOT, a C++17 compiler, the nlohmann JSON headers, and
OpenSSL development libraries. The local `pars/angle_info/` and `pars/beta_info/`
tables and the input ROOT's matching `dataNNNN_build_cpp_complete.json` receipt
must be supplied with the experimental data; they are not bundled here.

Build with `make -C doppler` from the project root. The same executable supports
the current narrow/wide build schema through explicit options:

```bash
./doppler/doppler --project "$PWD" \
  --parameters doppler/parameters/angle_table_baseline.json \
  --input rootfile/data0617_build_wide.root \
  --output /tmp/data0617_doppler_wide.root
```

The output must not exist. Each row of `pairs` represents one accepted Ge–Si or
Ge–S3 pair. Both detector families and all retained Si hits are preserved. No
timing correction, energy calibration, gate, or deduplication is applied again.
Input SHA256 and pair counts are checked against the build receipt.

`e_raw` is calibrated laboratory gamma energy; `e_dc_r` is the recoil-corrected
energy, both in keV. `particle_MeV`, detector indices, source identities, CFD flags,
`dt_ns`, `alpha`, and `random_side` are retained. `alpha` is the individual
prompt/random window-width ratio. The adopted spectrum policy is loaded from
`doppler/parameters/random_subtraction.json`: high-beam Si has
an additional factor 1.09; the other groups have factor 1. The signed `weight` is
1 for prompt and `-alpha * extra_factor` for random. Weighted variances use the
sum of squared weights. Never subtract the raw random spectrum with a universal
coefficient or apply the extra factor twice.

This tracked policy is an identical runtime copy of the project's local
`codex/resources/spectrum_processing/random_subtraction.json` policy. Keep the
two consistent if the adopted coefficients change. The executable only needs
the tracked copy; the entire policy and its SHA256 are embedded in each output.

The baseline preserves historical recoil exit-energy factors (Si: 1, S3: 0.5)
and the historical mass convention. It changes pair retention, not the baseline
physics. `exclusive_max_energy` records an optional historical-style selection
diagnostic. It is not applied to the main spectra.

`angle_info.cpp` owns the common correction model. Angle offsets also update the
interpolated recoil angle and exit energy; endpoint extrapolation is flagged.
`particle_energy_fraction > 0` enables a table-anchored momentum perturbation,
which is an experimental sensitivity candidate, not a validated new particle
calibration. Parameters and all input table hashes are embedded in the output.

For parameter scans of already extracted pair files, see
[`doppler_optimize/readme.md`](../doppler_optimize/readme.md). Pilot results and
their adoption limits are recorded locally in
`codex/analysis/doppler_modernization_20261028/`.

The positional `RUN WIN JUMP` command and documentation below describe the
historical file format and selection. They do not process the current build files.


# doppler

## how to run ?

```
./doppler run_number coin_win random_win / ./batch.sh
```

## info

+ form build rootfile get doppler correction results;
+ input file: ../rootfile/data%04d_build_%dns_jump_%dns.root
+ output file: ../rootfile/data%04d_doppler_%dns_jump_%dns.root
+ * parmeter file: run.txt, angle information and beta information.

## evaluation/d2hist

```
./batch_d2hist.sh
```

### info

+ from doppler rootfile to get all histograms, very slow.

## evaluation/d2hist_simple

```
./batch_d2hist_simple.sh
```

### info

+ from doppler rootfile to get part of histograms, such as `ge all different spider rings`, `ge ring3,4,5 different spider rings`,`ge all different s3 rings` and `ge ring3,4,5 different s3 rings`.

## evaluation/d2th2

```
> root
> .L d2th2.cpp
> d2th2("data0472_0605_doppler_100ns_jump_300ns")
```

### info

+ from doppler rootfile to get $\gamma$-$\gamma$ coincidenve th2 histogram.

## evaluation/get_no_bg_hist

```
> root
> .L get_no_bg_hist.cpp
> get_no_bg_hist("doppler_100ns_p1")
```

### info

+ from `doppler_100ns_p1_hist.root/doppler_100ns_p2_hist.root` to get `doppler_100ns_p1_no_bg_hist.root/doppler_100ns_p2_no_bg_hist.root`
+ include event* and bg* histograms

## evaluation/get_rebin_hist

```
> root
> .L get_rebin_hist.cpp
> get_rebin_hist("doppler_100ns_p1_no_bg")
```

### info

+ add 4 s3 rings together.

## evaluation/get_hist_aa

```
> root
> .L get_hist_aa.cpp
> get_hist_aa("doppler_100ns_p1_no_bg")
```

### info

+ add all spider and s3 together to seewhere I can get.

