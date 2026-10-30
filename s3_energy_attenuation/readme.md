# S3 relative energy-scale correction

## Maintained candidate workflow

The active position-factor assembly, candidate comparison and selection code
is now maintained in this directory:

| File | Responsibility |
| --- | --- |
| `analysis.py` | Assemble strip/position models; same-pixel spectra and same-event front/back diagnostics |
| `candidate_workflow.py` | Compare frozen models on fold 1; confirm on fold 3 or retain the coherent baseline |
| `export_parameters.py` | Export saved per-run choices, all 56 channels, 24 position factors and 768 pixel responses |
| `workflow_io.py` | Shared explicit input, SHA256, protected output and run-record handling |

The C++ fitters remain `conditioned_extract`, `conditioned_calibrate` and
`s3_fb_cor/joint_calibrate`. Python coordinates comparison and export; it does
not replace the native event extraction or fitting programs.

From the project root, use a **new output directory** for each invocation:

```bash
# Compare three existing local samples; no raw C2 or E-drive read.
python3 s3_energy_attenuation/candidate_workflow.py evaluate --config s3_energy_attenuation/workflow_0229.json --runs 0500 0623 0701 --jobs 1 --output-dir /path/to/new/comparison

# Export all 173 runs from their saved decisions; no event reading or fitting.
python3 s3_energy_attenuation/candidate_workflow.py export --config s3_energy_attenuation/workflow_0229.json --output-dir /path/to/new/conditioned_candidate
```

`--evaluation-dir /path/to/comparison/results` explicitly selects newly
evaluated results. Export requires every run in the supplied manifest; a
three-run comparison cannot silently replace a 173-run result. Omitting
`--runs` from **evaluate** explicitly processes every manifest run.

`workflow_0229.json` is a reproducible input configuration for the archived
study, not a hardcoded path inside the algorithms. Paths are relative to the
configuration file; manifest data paths are relative to that manifest.
For new inputs, copy the configuration and supply a matching baseline,
reference manifest and candidate manifest. Each candidate row identifies
`run`, `pairs`, `spectra`, `conditional`, `coupled` and their SHA256 fields.
Run IDs are four-digit strings. Dependencies: Python 3, NumPy, pandas, uproot.

Selection thresholds and energy/time gates retain the 0229 definitions.
`0228` denotes the supplied baseline, not the current production setting.
The conditioned export is an intermediate candidate; final targeted review
and build-compatible export use
[`s3_fb_cor/calibration_workflow.py`](../s3_fb_cor/calibration_workflow.py).
That review separately reports statistical support and spectral distance.

Saved input ROOT files and research receipts remain in their task directories.
Archived Python files under `codex/analysis/` preserve historical provenance;
they are no longer the maintenance source for these stages. The active code
does not import or execute those scripts. No production files are installed
by either workflow.

## Prepared complete replacement package (2026-10-30)

[joint_energy_20261030](parameters/joint_energy_20261030/README.md) provides all
173 runs of the strip multiplier k and matched-ring multiplier q. These are
exact exports of the single complete JSON maintained under
`s3_fb_cor/parameters/joint_energy_20261030/coefficients.json`; do not edit or
apply them independently of its A/B coefficients. The existing C++ build
accepts that JSON via `--s3-energy-file` in raw-input `windows` mode.

No production parameters or ROOT files were replaced during packaging.
[0232](../jupyter/codex_0232_150sm_s3_replacement_parameters.ipynb) lists the full
parameters and validation scope. The workflows below describe calibration
studies; their outputs are not automatically selected for production.

Native C++ extracts S3 singles from decoded C2 ROOT and compares every run/channel with a fixed same-period reference.

## Build and use

~~~bash
make
./extract 472 --project .. --input /mnt/e/Sm150/rootfile/data_C2_0472_0472.root \
  --output /path/to/spectra/spectra_0472.root --sample-entries 2000000 --segments 40
./calibrate --input-dir /path/to/spectra --output-dir /path/to/new_results --run-list ../pars/run.txt
~~~

Reference spectra 0472 (high I), 0616 (low), and 0680 (high II) must be present. Existing output files/directories are protected.

- `batch_tree2hist.sh OUTPUT_DIR [RUN_LIST] [JOBS=2] [SAMPLE_ENTRIES=2000000] [INPUT_DIR]` schedules 1–8 native readers.
- `batch_fit_xy_new2.sh SPECTRA_DIR NEW_RESULT_DIR [RUN_LIST]` runs registration.
- `batch_draw_sample.sh SPECTRA_ROOT NEW_FIGURE_DIR` draws all 56 clean channel spectra.

## Source selection

Input tree `tree` uses `cid, sid, ch, evte, pileup, outofr`. Energy is the ADC calibration polynomial times hardware gain, excluding previous attenuation and front/back normalization. Channel mapping follows `include/sm150_event_calibration.h`.

Inclusive spectra retain all mapped S3 singles; clean spectra exclude pileup, out-of-range and ADC 65535. No Ge, time-coincidence, multiplicity or front/back energy-agreement gate is applied. This selection differs explicitly from historical coincidence-gated inputs.

The default sample is two million entries in 40 centered segments distributed over equal source-entry intervals. Even segments fit, odd segments check. These are not necessarily equal time intervals. `--sample-entries 0` reads all entries.

ROOT output contains `clean_ringNN/clean_sectorNN` (energy × segment), `inclusive_*`, and JSON provenance in a TNamed object plus sidecar. Bins: 50 keV, 0–150 MeV, with underflow/overflow. Source identity, exact intervals, counts, old coefficients and parameter/program hashes are recorded.

## Relative-scale estimator

Reference and target spectra are normalized independently. Search k=0.5–5 minimizes squared Hellinger distance:

`H2(k) = 1 - sum(sqrt(p_i * q_i(k)))`

Target counts are integrated over comparison-bin edges divided by k. Fit window: 20–90 MeV in reference energy; comparison bins: 250 keV. A fixed [0.25,0.5,0.25] kernel on 50-keV source bins smooths interpolation and conserves counts. Logarithmic grid search is followed by local refinement.

Odd segments validate the fixed fit; a separate fit-only 30–80 MeV window diagnoses shape sensitivity. Flags cover low counts, search boundaries, weak scale identification, split disagreement, window sensitivity and residual shape mismatch. All numerical thresholds are exported.

Distances and profile support widths are descriptive, not chi-square values or calibrated confidence intervals. Window membership changes with k. A scale-only model cannot repair arbitrary spectrum-shape changes.

## Output and adoption

`results.json` and `parameters.csv` retain old, fitted candidate and proposed k, flags, counts and validation metrics. `profiles.root` holds fit-only distance curves. The 56 files in `candidate_parameters/` use run/k columns and explicitly inherit old coefficients for flagged or unimproved branches.

Candidates are never automatically copied to `pars/run_s3_att/`. Inspect spectra and compatibility with front/back normalization before adoption. Relative registration does not establish an absolute energy calibration or a physical cause of drift.

Original macros, synthetic checks and analysis provenance from 2026-10-29 are retained locally under `codex/analysis/s3_attenuation_modernization_20261029/`. The old reference-threshold count dependence and append-only parameter output have been removed; this does not imply every historical coefficient was wrong.



# Coupling to front/back normalization

Single-channel spectrum matching and same-event front/back alignment test
different assumptions. Sector spectra mix ring angles, and a common particle
spectrum change is not uniquely identifiable as gain drift. The joint solver
in `s3_fb_cor/joint_calibrate.cpp` keeps these checks separate and exports total
`A + B*E0` coefficients before factorization. A new attenuation file must be
used with its matching front/back file; compare total responses, not isolated
`k` values defined with different normalization slopes. Existing spectrum-only
outputs remain candidates and are not implicitly promoted by the joint analysis.

## Fixed-pixel attenuation checks

`make conditioned_extract conditioned_calibrate` builds a separate, explicit
paired-event analysis. It preserves the singles workflow above.

```bash
./conditioned_extract --input paired_sample.root --output new_pixel_spectra.root
./conditioned_calibrate --input new_pixel_spectra.root --reference reference_pixel_spectra.root --baseline previous_coefficients.json --output new_conditional.json
```

The extractor reads the compact `sm150_s3_joint_sample_v1` format. It retains
four source-block folds, both faces and every pixel, with prompt/random counts
and flow bins. Clean single-hit pairs have no added energy-agreement gate.
Prompt is [-80,80) ns; random is [-880,-720) and [720,880) ns, with alpha=0.5.
The allowed source-block range is read from the source receipt's `segments`
array. This supports both historical 40-segment samples and separate
80-segment checks; the histogram split remains `block % 4`. Legacy receipts
without `segments` retain the previous 40-block bound.

The estimator compares the same ring-sector pixels across runs. It fits a
per-channel residual multiplier on the explicitly supplied total A/B baseline.
Each pixel contributes equally after its own signed P-0.5R normalization.
Negative net bins remain signed. The metric is the RMS difference between
cumulative net spectra over20-90MeV, sampled every1MeV. It is a descriptive
shape criterion, not a likelihood or confidence interval. Per-pixel estimates,
actual search curves, alternate-window and alternating-fold fits are retained.
Only folds0/2 determine coefficients; folds1/3 report selection/confirmation.

Different channels may have different coefficients. If one sector needs
different multipliers at different rings, its strip-only model may be
insufficient, or physical population/selection changes may be present. A shared
position factor applied to both faces is a separate hypothesis. It requires
the paired ring context and cannot be represented by silently replacing a
single sector attenuation table. The analysis package records any such factor
explicitly; no formal parameters or build output are changed automatically.

Separate spectral distance from statistical support when reviewing these
outputs. A guard value used for an unsupported search point is not a measured
shape discrepancy. The targeted review records the actual signed cumulative
distance, low-count/invalid support flags, and an approximate counting-noise
scale separately. Counting-noise estimates do not include fitted-parameter or
model uncertainty. Keep the same pixel set for every compared model, retain
negative net bins, and do not silently drop a failed candidate pixel.

