# 150Sm_coulex

Data processing and analysis for the 32S + 150Sm Coulomb-excitation experiment:
Ge gamma-ray detectors, a Si array, and the double-sided S3 detector.

This overview describes the maintained processing entries and the adopted
configuration checked on **2026-11-03**. Runtime choices come from
`pars/active.json`; each generated ROOT file retains the parameters used to
produce it. Changing a configuration does not update existing ROOT files.

## Workflow and program entries

Production processing starts from decoded C2 ROOT hits. Calibration and timing
studies establish the parameters; `build_all/build` applies those parameters
directly to decoded hits, produces narrow/wide coincidence files, and supplies
the input for timing checks and Doppler correction.

| Directory | Maintained entry or purpose |
| --- | --- |
| [raw2hist](raw2hist/README.md) | `histogram`: ADC and calibrated single-channel spectra from decoded C2 or explicitly supplied C1 data |
| [cali_ge](cali_ge/readme.md) | `calibrate`: Ge source calibration and run-dependent gain diagnostics |
| [cali_si](cali_si/readme.md) | `calibrate`: Si array and S3 alpha-source calibration |
| [shift_resolution](shift_resolution/) | Source/in-beam peak shifts and Ge energy-resolution diagnostics |
| [get_ts_offset](get_ts_offset/main.cpp) | `analysis`: particle-channel timing measurements and offset studies |
| [ge_time_walk](ge_time_walk/README.md) | `build` and `fit`: Ge timing diagnostics, energy-dependent fits, and complete timing-package maintenance |
| [build_s3](build_s3/readme.md) | `build`: particle-only S3 calibration events after particle timing correction |
| [s3_fb_cor](s3_fb_cor/readme.md) | Front/back A/B fitting, review, and complete parameter export |
| [s3_energy_attenuation](s3_energy_attenuation/readme.md) | Run-dependent ring/pixel q estimation and confirmation |
| [s3_coin_ana](s3_coin_ana/README.md) | `ana`: S3 two-particle diagnostics using the complete S3-only event tree |
| [build_all](build_all/main.cpp) | `build`: production Ge-particle coincidences with separate narrow and wide selections |
| [check_ts_align](check_ts_align/check_build.cpp) | `check_build`: time-difference spectra and consistency checks from production build files |
| [beta_info](beta_info/readme.md) | Kinematics and target energy-loss inputs used to derive the recoil correction |
| [doppler](doppler/readme.md) | `doppler`: calibrated laboratory and recoil-corrected gamma energies for accepted pairs |
| [doppler_optimize](doppler_optimize/readme.md) | Explicit parameter scans and comparisons on saved pairs |
| [ge_eff](ge_eff/readme.md) | Source peak areas, efficiency studies, and the retained all-Ge efficiency reference |

`real_cali_sort_ge_si/` retains the separate calibrated-hit sorting workflow.
It is not a prerequisite for the current direct-C2 production build. Older
macros under `evaluation/` and positional `RUN WIN JUMP` commands may require
historical ROOT schemas; use the maintained entries above for current files.

## Runtime parameters

The local [parameter index](pars/README.md) documents the selection in
[pars/active.json](pars/active.json). [pars/manifest.json](pars/manifest.json)
records file identities and source paths. The adopted numerical payloads are
stored directly in [include/parameter](include/parameter/README.md):

| Parameter | Current payload |
| --- | --- |
| Ge base energy calibration | [cali_ge.txt](include/parameter/cali_ge.txt) |
| Ge run-dependent gains | [run_k/](include/parameter/run_k/) |
| Si/S3 source calibration | [cali_si_r460.txt](include/parameter/cali_si_r460.txt), [cali_si_r457_low_current.txt](include/parameter/cali_si_r457_low_current.txt), [cali_si_r457_high_current.txt](include/parameter/cali_si_r457_high_current.txt) |
| Si/S3 adopted and retained channel values | [particle_source_parameters.csv](include/parameter/particle_source_parameters.csv) |
| Complete S3 A/B/q correction | [s3_coefficients.json](include/parameter/s3_coefficients.json), [A/B table](include/parameter/s3_ab_parameters.csv), [q table](include/parameter/s3_q_parameters.csv) |
| Complete Si/S3/Ge time correction | [time_alignment.json](include/parameter/time_alignment.json) |
| Narrow/wide prompt and random windows | [coincidence_windows.json](include/parameter/coincidence_windows.json), [readable table](include/parameter/coincidence_windows.csv) |
| Additional random-subtraction factors | [random_subtraction.json](include/parameter/random_subtraction.json) |
| Final Doppler beta values and recoil directions | [doppler.json](include/parameter/doppler.json), [beta table](include/parameter/doppler_beta.csv), [provenance](include/parameter/doppler_beta_provenance.json) |
| Numerical payload checksums | [manifest.json](include/parameter/manifest.json) |

The older `pars/` numerical paths and `include/parameters/` paths remain
compatibility links. They do not represent independent copies to edit.
Detector masks, run lists, geometry, and kinematics retain their `pars/` layout.
Domain directories maintain the generating programs and frozen source packages;
installing a candidate does not select it for production.

The current selector is `doppler_tabulated_beta_20261103`, with the adopted
`energy_parameters_20261102` energy settings, `joint_zero_20261028` timing,
and `narrow_review_20261103` coincidence windows. These names identify versions,
not extra corrections to apply in sequence.

The strict run list is [pars/runs/run.txt](pars/runs/run.txt):

| Period | Run-number range | Listed runs | Production S3 rings |
| --- | --- | ---: | --- |
| High beam I (`high1`) | 0472-0605 | 101 | 1-8 |
| Low beam (`low`) | 0616-0674 | 40 | 1-24 |
| High beam II (`high2`) | 0680-0718 | 32 | 1-8 |

Only the listed **173 runs** are production inputs; the ranges are not
continuous run lists. The Ge gain directory retains 184 tables for provenance
and other uses, which does not expand the production list.

## Calibration and timing conventions

Ge uses the retained base polynomial and its run-dependent gain. Si/S3 uses
the adopted source-calibration package: the selected new S3 coefficients and
accepted Si channel updates, with previous values retained where a replacement
was not accepted. The channel table above records these choices.

S3 source calibration and A/B/q form one matched package. For a matched ring
and sector, production build evaluates, in keV:

```text
E0 = hardware_gain * (p0 + p1*ADC + p2*ADC^2)
Er = q * (Ar + Br*E0r)
Es = q * (As + Bs*E0s)
```

Low beam uses the matched pixel's q; high beam uses the matched ring's q.
The pixel factor replaces the ring factor. Each face receives q once, and
the complete A/B/q correction replaces the older attenuation/front/back chain.
It is derived on the adopted source scale; those tables must be kept together.

`raw2hist` and `build_s3` retain their documented diagnostic energy stages.
Single-channel spectra cannot apply a factor that needs the matching face or
ring identity. Their S3 energy stages should therefore not be interpreted as
the complete matched-pair production correction.

The complete timing function is applied **once** to the native decoded time:

```text
t_corrected = t_native + F(run, channel, CFD state, energy/ADC)
dt = t_particle_corrected - t_Ge_corrected
```

Normal Ge CFD uses the native fine-time information; failed CFD uses coarse
time. The package already contains adopted offsets and energy-dependent terms
for all detector families. Do not add historical offsets or time-walk terms
again. Timing alignment precedes the S3 front/back pairing used for calibration
and production selection. A corrected peak position does not imply that broad
tails have disappeared.

## Build and check production events

The maintained C++ entry is `build_all/main.cpp`, with processing in
`build_all/build.cpp` and `build_all/build.h`. It reads decoded C2 data,
calibrates energies, applies the complete timing correction, forms pairs, and
applies the narrow and wide windows separately in one pass.

Current selection uses 40 <= Ge energy < 2000 keV, 15 < Si energy < 100 MeV,
and 25 < S3 energy < 100 MeV. S3 requires a unique ring-sector match with
|Er-Es| < 2.5 MeV and |tr-ts| < 80 ns. High-beam S3 rings above 8 are excluded
from production. Normal and failed CFD are retained; pileup and out-of-range
hits are rejected. This is a single-Ge selection with no Clover addback.

Run commands from the project root. Use a **new output directory** for a trial:

```bash
make -C build_all
build_all/build 0617 --inspect-parameters

build_all/build 0617 \
  --input /mnt/e/Sm150/rootfile/data_C2_0617_0617.root \
  --outdir rootfile/example_0617
```

The inspection command reads parameters only. The subsequent build reads the
whole run unless an explicit entry limit is supplied. Its outputs are:

- `data0617_build_narrow.root` and `data0617_build_wide.root`.
- `data0617_build_cpp_complete.json`, the completion and identity receipt.

`tr_event` and `tr_bg` retain prompt and random events separately, including
corrected timestamps, laboratory energies, source-hit identities, random-side
labels, and each pair's window-width normalization `alpha`. Build performs no
Doppler correction or random subtraction. Existing outputs are not overwritten.
Production files use the same naming convention in the top-level `rootfile/`.

For timing checks, compile with `make -C check_ts_align` and pass a text file
containing one build ROOT path per line to
`check_ts_align/check_build --input-list FILE --output NEW.root`.
Inputs must have the same period, narrow/wide selection, timing version, and
gate definitions, and must retain their matching completion receipts. The
checker uses stored corrected times; it does not apply the correction again.
Its spectra are already gated, so they cannot reconstruct discarded tails or
provide an independent ungated Si/S3 timing measurement.

S3 two-particle analysis is a separate branch: `s3_coin_ana` needs the complete
`s3_events` tree from `build_s3`. Production Ge-gated files retain at most one
S3 pair per entry and cannot supply the missing two-particle events.

## Doppler correction and spectra

The current Doppler model reads precomputed final **beta = v/c**, recoil
directions, and associated ring information from `include/parameter/doppler.json`.
It contains 96 entries: three periods, each with 8 Si rings and 24 S3 rings.
The adopted scaling is already included in beta; event processing applies no
additional beta or recoil-energy multiplier. Ge and particle angles are retained.
Geometry and kinematics identities are checked against the saved table sources.

For the example build above:

```bash
make -C doppler
doppler/doppler --project "$PWD" \
  --input rootfile/example_0617/data0617_build_narrow.root \
  --output rootfile/example_0617/data0617_doppler_narrow.root
```

The input's completion receipt must remain alongside it. The new output stores
one row per accepted Ge-Si or Ge-S3 pair in `pairs`, including laboratory energy
`e_raw`, recoil-corrected energy `e_dc_r` (both keV), and the original timing
and random weights. Calibration, timing, and coincidence gates are not reapplied.
Changing only Doppler parameters requires reprocessing the accepted pairs,
not rebuilding them from decoded C2 hits.

Window normalization is specific to each Ge, particle family, period, and
narrow/wide selection:

```text
alpha = prompt width / total random width
prompt weight = +1
random weight = -alpha * c
```

The current additional factor is **c = 1 for every group**. `alpha` still
accounts for the actual random-window widths, including both sides where
defined. Accumulate weighted spectra and their sums of squared weights before
combining detector groups. A universal unscaled `P-R` is valid only where the
random counts have already received the appropriate width normalization.
Preserve negative net bins and the separate prompt/random counts.

`make -C doppler test` checks the direct-beta model without reading ROOT events.
The exporter and historical comparison modes are documented in
[doppler/readme.md](doppler/readme.md). Parameter scans in `doppler_optimize`
are explicit studies and do not activate their candidates automatically.

The retained original **all-Ge** efficiency coefficients, function, and
normalization are recorded in
[ge_eff/parameters/all_ge_reference.json](ge_eff/parameters/all_ge_reference.json).
Recent source-efficiency refits remain candidates. Efficiency correction belongs
to spectrum/yield analysis; it is not an extra event-build correction.

## Dependencies and local data

The maintained native programs use ROOT with `root-config`, a C++17 compiler,
nlohmann JSON headers, and OpenSSL development libraries. Fitting targets may
also require ROOT Spectrum/Minuit2. Python analysis and plotting tools use
NumPy, SciPy, pandas, uproot, and Matplotlib as documented in their directories;
reading or running notebooks additionally requires Jupyter.

The current local decoded-data location is `/mnt/e/Sm150/rootfile/`. Production
build files are kept in the project's `rootfile/`. Decoded event files, local
runtime configuration, geometry/kinematics inputs, and output receipts are not
supplied by cloning the source repository. Provide the matching `pars/`
environment and data before running the examples. To inspect that environment
without processing events:

```bash
python3 pars/manage.py list
python3 pars/manage.py verify
```

## Local research material

The following workspace directories are excluded from Git:

| Directory | Contents |
| --- | --- |
| `codex/` | Analysis scripts, validation records, navigation, and handoff notes |
| `jupyter/` | Reading notebooks and their history |
| `gosia/` | GOSIA programs, calculation cases, and recovery results |
| `sigma/` | SIGMA shape calculations and studies |
| `fit_yield/` | Peak fitting, efficiency-corrected yields, and exports |
| `archives/` | Original recovery archives |
| `level_scheme/` | Level diagrams and historical versions |
| `tmp/` | Temporary working files |

`AGENTS.md` also remains local. Local navigation starts at
[codex/README.md](codex/README.md) and [jupyter/README.md](jupyter/README.md);
these links and the local parameter index require the populated workspace.
Production parameters and source code have their own publication scope;
generated ROOT data, notebooks, and research history are not included in source
commits.
