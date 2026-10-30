# S3 front/back normalization

## Maintained review and final export workflow

The active targeted review, confirmation and complete parameter export now
live here, alongside the C++ fitters:

| File | Responsibility |
| --- | --- |
| `review.py` | Signed spectra, counting-support diagnostics and unchanged 0230 acceptance guards |
| `calibration_workflow.py` | Select a whole-run model, freeze its choice, confirm it, or retain the whole baseline |
| `export_parameters.py` | Export the complete A/B/q package and channel, position, pixel and review tables |

The numerical solver remains `joint_calibrate.cpp`. Its fixed-offset mode
below produces the `fixed_fit` JSON supplied in the review manifest. Candidate
fitting is explicit; selection and confirmation never trigger a new fit.
Position-factor assembly and the earlier candidate comparison are maintained
in [`s3_energy_attenuation`](../s3_energy_attenuation/readme.md).

From the project root:

```bash
# Select on saved fold-1 pairs and freeze the choice before confirmation.
python3 s3_fb_cor/calibration_workflow.py select --config s3_fb_cor/workflow_0230.json --runs 0500 0553 0623 --jobs 1 --output-dir /path/to/new/selection

# Confirm those frozen choices on the configured disjoint source segments.
python3 s3_fb_cor/calibration_workflow.py confirm --config s3_fb_cor/workflow_0230.json --selection-freeze /path/to/new/selection/selection_frozen.json --jobs 1 --output-dir /path/to/new/confirmation

# Re-export all 173 runs from the existing 31-run confirmed review.
python3 s3_fb_cor/calibration_workflow.py export --config s3_fb_cor/workflow_0230.json --output-dir /path/to/new/complete_candidate
```

The last command uses the archived complete decisions selected by the config;
it does **not** consume the preceding three-run example. To export a newly
completed full review, provide both `--selection-freeze NEW_FREEZE.json` and
`--confirmation-dir NEW_CONFIRMATION/results`. The frozen and confirmed run
sets must exactly equal the targeted manifest. Other runs inherit the complete
previous candidate, with explicit status. A partial review cannot masquerade
as a completed full review.

Every invocation requires a new output directory. `--runs` selects a subset
for local checks; omit it to process the complete configured selection or
confirmation list. `--jobs` is explicit (1 by default, at most 8). Importing
modules and using `--help` perform no analysis. No polling, automatic E-drive
acquisition, production installation or build submission occurs.

`workflow_0230.json` supplies all input paths; they resolve relative to the
configuration file. Copy it for new studies. The selection manifest contains
`run`, `pairs`, `spectra`, `conditional`, `fixed_fit` and their SHA256 fields;
the reference and confirmation manifests identify their own inputs. Baseline,
previous candidate, protocol, selection freeze, review rows and confirmation
directory are explicit. The protocol is retained as provenance; acceptance
thresholds are implemented in `review.py`, not dynamically read from it.
Changing scientific thresholds requires separate validation.

Frozen files, model vectors and baseline identities are checked. Confirmation
also checks ROOT metadata and actual source-entry intervals against the
selection and reference samples. A failed confirmation returns to the coherent
baseline; it does not choose another model or splice individual strips.
Repeating these commands on previously inspected samples is reproduction,
not new independent confirmation.

Final `coefficients.json` uses `sm150_s3_targeted_review_candidate_v1`, accepted
by `build_all/build --s3-energy-file`. It contains all A/B/q coefficients and
the matched-ring requirement. `channel_coefficients.csv`,
`ring_position_factors.csv`, `pixel_coefficients.csv`, `run_decisions.csv` and
`review_decisions.csv` remain inspection tables; build reads the complete JSON.
Dependencies: Python 3, NumPy, pandas, uproot, plus the adjacent
`s3_energy_attenuation` Python modules.

The archived `codex/analysis/` scripts remain frozen provenance only. The new
entry point does not import or execute them. Historical freeze receipts may
still verify their source hashes, while newly frozen receipts record the
maintained files in these two native directories.

## Prepared complete replacement package (2026-10-30)

[joint_energy_20261030](parameters/joint_energy_20261030/README.md) contains the
frozen 173-run A/B/q candidate, full channel tables, provenance and limitations.
The associated strip/position factors are exported in
`s3_energy_attenuation/parameters/joint_energy_20261030/`.
Use the complete JSON through the existing native frontend:

```bash
./build_all/build 0634 --s3-energy-file s3_fb_cor/parameters/joint_energy_20261030/coefficients.json --outdir rootfile/s3_joint_energy_20261030
```

Run this example from the project root. The candidate replaces the old
attenuation/front-back energy chain and applies the matched-ring factor once
to both faces. It is not a drop-in pair of legacy text tables. The default
production configuration and existing ROOT files have not been changed.
See [0232](../jupyter/codex_0232_150sm_s3_replacement_parameters.ipynb) for every
run/channel value, the three-period maps, and the retained 0634 tradeoff.
The following sections describe the fitting tools and earlier staged workflows.

## Current compiled workflow

Use `build_s3/build` first to extract unmerged particle-only S3 pairs from decoded C2 data with the current complete time correction. Production Ge-coincidence build files already have selections and are not equivalent calibration inputs.

From the project root:

```sh
make -C build_s3
make -C s3_fb_cor
build_s3/build 0617 --input /mnt/e/Sm150/rootfile/data_C2_0617_0617.root --output /path/to/new/s3_0617.root
s3_fb_cor/calibrate --input /path/to/new/s3_0617.root --outdir /path/to/new/candidate_0617
```

The fitter refuses an existing output directory. It never modifies `pars/run_fb_cor`. Candidate text tables contain all 24 rings and 32 sectors; consult `results.json` for every fit, fixed reference, inherited value, uncertainty, and validation result before considering adoption.

The calibrated energy is `E_final = p0 + p1 * E_pre_fb`, with energies and p0 in keV and dimensionless p1. `E_pre_fb` already includes the existing source calibration, gain and attenuation. Sector 1 keeps its existing production p0/p1 to anchor the energy scale. All rings are fitted against sector 1; sectors 2–32 are fitted against newly corrected ring 11. This preserves the direction of the historical coefficients despite the reversed labels in old plots. A relative front/back calibration does not independently establish the absolute energy scale.

Default selection: raw multiplicity one on each face, clean quality flags, nonboundary cluster, absolute corrected time difference <=80 ns, and both pre-normalization energies in [15,100) MeV. There is no ring/sector energy-agreement cut. The source acquisition is partitioned into 40 equal entry-count blocks: even blocks fit, odd blocks validate. These are interleaved internal validation blocks, not independent experiments or exactly equal wall-time intervals.

The affine fit starts from medians in 2.5 MeV bins and the median of separated-bin slopes, then uses Huber IRLS (k=1.5, MAD scale floor 5 keV, at most 50 iterations). At least 100 training pairs and 5 MeV between energy quantiles 10% and 90% are required. Failed or unsupported channels inherit the previous coefficients with explicit flags. Failure of the ring-11 reference prevents unsupported sector updates. Two sequential local pair-tree passes avoid holding every pair in memory; the original E-drive data are not reread by this program.

`results.json` records exact medians, central-68% half-widths and the fraction beyond 2.5 MeV for all validation pairs, by channel, reference pair and acquisition block. The same held-out events are used for unnormalized, previous and candidate comparisons. A coarse no-large-degradation flag is only a diagnostic, not an adoption decision or a significance test. All-event and reference-only residuals both matter; a good reference fit can still broaden another pairing.

Quoted errors are conditional weighted regression errors. The ring-11 statistical reference is propagated into sector errors, but absolute energy calibration, model choice, time variation, x-axis errors and full cross-channel covariance are not included. Do not treat these as complete uncertainty budgets.

`diagnostics.root` provides training correlations, previous/candidate residual-versus-energy plots, complete channel histograms, and combined plots. One-dimensional residual bins are 20 keV over +/-100 MeV; two-dimensional residual panels cover +/-20 MeV, with ROOT flow bins retained. `results.json` quantiles use all selected events, including histogram overflows.

An explicit batch uses:

```sh
bash s3_fb_cor/batch_new.sh runs.txt /path/to/s3_sources /path/to/new/candidates 1
```

The input naming is `s3_RUN.root`. The worker count is explicit (default 1); choose it according to available memory. No production adoption or downstream event rebuild is automatic.

## Historical macros

`tree2gr*.cpp`, `correction_*_new.cpp`, `macro_new.cpp`, and `draw_sample*.cpp` remain for old `rootfile/si/...tr_event` files. They do not consume the new `s3_pairs` schema and are not part of the compiled workflow. Legacy `fast` graphs use the first 1/8 of a run; some checks use the first 1/32. They are prefix samples, not full-run validation. Legacy graph extraction also uses repeated tree scans and large histograms; use the compiled workflow for new analysis.

Historical fit macros can overwrite local `par/` files and lack the new fit/coverage guards. Their axis labels have been corrected without inverting coefficients. Legacy comparison macros now apply the same stated cuts, but separate old/new source files do not by themselves prove identical event membership. The current fitter evaluates all states from the same row and avoids that ambiguity. Drawing a graph no longer adds an artificial (100000,100000) point.
# Global joint-calibration candidate

Build `joint_calibrate` with `make`. Example:

```bash
./joint_calibrate --input S3_SAMPLE.root --output NEW_FIT.json --min-pairs 20
```

The input is the compact sample from `build_s3/build --sample-entries`.
Even source segments determine robust orthogonal pixel lines, followed by a
connected global strip gain/offset solution. The solver retains failed pixels,
connectivity, alternating-segment fits, and exact total coefficients. Odd
segments are excluded from fitting. Model selection and confirmation are
performed separately by the recorded analysis workflow.

High-intensity production constraints use rings 1--8 and all sectors. Outer
rings are diagnostic fits conditional on the sectors; low intensity uses all
24 rings. The old ring02 energy defines the temporary common gauge. This
constraint does not determine an absolute energy calibration or common losses.

The total response is `E = A + B*E0`, where `E0` excludes old attenuation and
front/back normalization. A compatible two-file decomposition is
`p1 = B_reference`, `k = B_run/B_reference`, `p0 = A_run`. This is a convention;
`p1` and `k` cannot be independently identified from the same equality constraint.
Replace both parameter sets together only after their candidate is adopted.
Do not apply them again to energies that already contain the old corrections.

## Temporal spectra in the joint solution

An optional mode reuses frozen even-block pixel lines and supported per-channel
temporal estimates from the fixed-pixel attenuation program:

```bash
./joint_calibrate --edge-input even_block_fit.json --temporal-input conditional_fit.json --output new_temporal_joint.json
```

This mode releases the fixed ring02 gauge and uses temporal constraints at two
energies per supported channel. Gain/affine models and weight multipliers
0.1,1,10,100 are exported. The weights use descriptive split/window/profile
stability; they are not calibrated prior uncertainties. Selection and held-out
confirmation remain external, explicit steps. At least three supported
production channels are required. The original `--input` mode remains unchanged
and was checked against its earlier complete synthetic edge/model output.

The solver still describes strip responses as `A+B*E0`. If an analysis also
uses a ring-dependent position factor, it must be applied to both faces of
that same pixel and recorded separately. Such a candidate is not automatically
representable by legacy per-strip text tables. The complete targeted-review
JSON is supported by the explicit `build_all/build --s3-energy-file` option.

## Fixed total offsets with joint slope fitting

To reduce offset/slope compensation, an explicit constrained mode fixes each
total intercept `A` at a supplied baseline and solves the connected strip
slopes together:

```bash
./joint_calibrate --edge-input even_block_fit.json --temporal-input conditional_fit.json --fixed-offset-baseline baseline.json --output new_fixed_offset_fit.json
```

The baseline contains `runs[RUN].period` and all 56 `channels`, each with an
`index`, `A_keV`, and `B`. This fixes `A` in `E=A+B*E0`, not the relative offset
of an already corrected energy. The latter would also change the total
intercept when multiplied. The four descriptive temporal weights remain
0.1, 1, 10, and 100; they are not statistical prior precisions.

Production ring/sector slopes remain coupled. Disconnected strips and high-beam
outer rings keep the supplied baseline coefficients, with explicit diagnostic
status. Replacing one strip independently can break its relations to other
strips; evaluate the complete connected solution. Choose the model before
examining a separate confirmation sample, and retain a coherent baseline if
the chosen model fails. No parameter installation or event rebuild is automatic.

