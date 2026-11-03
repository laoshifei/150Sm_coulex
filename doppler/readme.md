# Current build interface (2026-11-03)

The executable defaults to the adopted configuration selected by `pars/active.json`.
Its payload is `include/parameter/doppler.json`, reached through
`pars/doppler/current.json`. It stores final beta values (`v/c`) for each
period and particle ring: high1/low/high2, Si rings1-8 and S3 rings1-24.
The adopted scales have already been included in these values. Current event
processing reads beta directly and applies no further energy or beta scale.
`include/parameter/doppler_beta.csv` lists all96 entries.

All original Ge/particle angles are retained. In particular, the Ge214 phi
trial was not adopted. Source scales and their conversion are recorded in
`include/parameter/doppler_beta_provenance.json`, which is not loaded for event
correction. Corrected energies remain identical to the adopted source model.
An explicit `--parameters` still selects a separate file; no trial is activated automatically.
`parameters/angle_table_baseline.json` remains the historical all-period
Si=1.0/S3=0.5 configuration for explicit comparisons.

Schema3 uses `beta_by_period`, with all three periods and both detector
families explicitly present. Each ring stores beta, particle/recoil polar
angles and expected projectile energy. The build manifest determines the
period. The correction cache includes the period, including in mixed-run
scans. Output metadata records the period and `precomputed_ring_beta` mode.
No scale can be supplied together with this table. Particle polar-angle scans
and event-energy blends require a separately generated compatible model;
they cannot reuse a fixed ring-beta table. Input geometry/kinematics hashes
are checked to prevent mixing beta values with changed source tables.

Schemas1/2 remain available for explicit historical comparisons. Older
executables reject schema3 instead of silently applying a second scale.

The local-only `parameters/summary_20261029/` directory contains the parameter
reference and historical five-/seven-Ge trial exports. Those research exports
are not required by the executable and are not part of this source commit.

Compilation requires ROOT, a C++17 compiler, the nlohmann JSON headers, and
OpenSSL development libraries. The local `pars/doppler/geometry/` and `pars/doppler/kinematics/`
tables and the input ROOT's matching `dataNNNN_build_cpp_complete.json` receipt
must be supplied with the experimental data; they are not bundled here.
See [the runtime parameter index](../pars/README.md) for the central manifest,
checksums and compatibility paths. The default random policy is selected from
`pars/coincidence/`; the numerical policy is unchanged.

Build with `make -C doppler` from the project root. The same executable supports
the current narrow/wide build schema through explicit options:

```bash
./doppler/doppler --project "$PWD" \
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
prompt/random window-width ratio. The adopted spectrum policy is selected by
the `random_subtraction` entry in `pars/active.json`. As of 2026-11-03, its
payload is `include/parameter/random_subtraction.json`, with additional factor
1 for all four intensity/particle groups. The signed `weight` is
1 for prompt and `-alpha * extra_factor` for random. Weighted variances use the
sum of squared weights. Never subtract the raw random spectrum with a universal
coefficient or apply the extra factor twice.

The older tracked `doppler/parameters/random_subtraction.json` remains an
explicit legacy fallback, including the former high-beam Si factor 1.09. It is
not the current default while the central selection above is present. Do not
mix results from these two policies. The actual loaded policy and its SHA256
are embedded in each output.

The historical mass convention is retained. Changing only the Doppler factor
does not require rebuilding the accepted Ge-particle pairs; existing Doppler
outputs retain their own saved parameters until explicitly reprocessed.
`exclusive_max_energy` records an optional historical-style selection
diagnostic. It is not applied to the main spectra.

`angle_info.cpp` owns the common correction model. Angle offsets also update the
interpolated recoil angle and exit energy; endpoint extrapolation is flagged.
`particle_energy_fraction > 0` enables a table-anchored momentum perturbation,
which is an experimental sensitivity candidate, not a validated new particle
calibration. Parameters and all input table hashes are embedded in the output.

`make -C doppler test` compares stored beta with its frozen source model,
covering all particle rings/sectors and interleaved periods. It also checks
unchanged geometry and rejection of duplicate scales or incomplete tables.
It uses local parameter files and does not read ROOT event data.

`make -C doppler export_beta` builds the parameter exporter. Its interface is
`doppler/export_beta PROJECT SOURCE_JSON NEW_JSON`; it refuses to overwrite an
existing output. It uses the same correction model and saves beta with
round-trip double precision. It does not read events or activate the export.

For parameter scans of already extracted pair files, see
[`doppler_optimize/readme.md`](../doppler_optimize/readme.md). Pilot results and
their adoption limits are recorded locally in
`codex/analysis/doppler_modernization_20261028/`.

The positional `RUN WIN JUMP` command and documentation below describe the
historical file format and selection. They do not process the current build files.
That interface reads its historical energy factor from
`pars/doppler/angle_table_baseline.json`, rather than a number embedded in C++.


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

