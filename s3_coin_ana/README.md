# S3 two-particle coincidence diagnostics

`ana` reads S3-only events, resolves the two possible ring–sector assignments,
and saves multiplicity, time, energy and angular diagnostics. Event processing
is maintained in `main.cpp`, `analysis.cpp` and `s3_tree.cpp` in this directory.

## Input

The preferred input is `build_s3` output containing the complete `s3_events`
vector tree and `sm150_s3_calibration_v1` metadata. Use explicit file paths.
The current Ge-gated `rootfile/dataRUN_build_narrow.root` and `wide.root`
retain at most one S3 pair per entry and cannot recover two-particle events.
The compact calibration `pairs`/`s3_pairs` trees are also insufficient.

Historical S3-only `tr_event` input requires explicit `--format legacy`.
Variable-length legacy arrays are read dynamically; they are not truncated to
two hits. Legacy files lack the modern quality, boundary and original-hit
identity metadata, so those controls cannot be reconstructed for them.

## Compile and run

Requirements: ROOT with `root-config`, C++17, OpenSSL, and the shared project
header `include/sm150_event_calibration.h` (including nlohmann JSON).

From the project root:

```bash
make -C s3_coin_ana
s3_coin_ana/ana 0617 --project . \
  --input codex/analysis/s3_fb_cor_modernization_20261029/data/s3_0617.root \
  --output rootfile/s3_coin_example_0617.root \
  --sample-entries 2000000 --segments 20
```

The output path, its `.json` receipt and `.partial` working file must not exist.
`--sample-entries 0` explicitly selects the complete input tree; a positive
number selects centered, nonoverlapping segments distributed over its entries.
Sampling fractions refer to stored two-sided time clusters, not beam exposure.

`--help` lists all energy, time and angle settings. For a bounded batch:

```bash
bash s3_coin_ana/batch.sh RUN_LIST INPUT_DIR NEW_OUTPUT_DIR JOBS SAMPLE_ENTRIES
```

Inputs are named `s3_RRRR.root`. The script validates the explicit list, waits
for all submitted jobs, stops dispatching after a failed job and returns a
nonzero status on failure. It does not submit the project run list implicitly.

## Calibration and event selection

The default uses stored `energy_previous_keV` and corrected `time_ns`. It
does not reapply timing offsets, time-walk or energy calibration. For an
explicit complete energy candidate, supply `--s3-energy-file FILE` using the
`sm150_s3_targeted_review_candidate_v1` schema. The ADC calibration snapshot in
the input metadata reconstructs E0, then applies `q_ring*(A+B*E0)` to each
face once, replacing the previous attenuation/front-back chain. The ring
context is recomputed for each proposed assignment. Candidate thresholds
cannot recover hits excluded by the input's earlier energy threshold.

The default diagnostic gates are strict:

- Each face: 2.5–100 MeV; |ring−sector time| < 80 ns and |Er−Es| < 2.5 MeV.
- Two particles: |difference of mean pair times| < 80 ns.
- Absolute azimuth difference and hypothesized CM angle sum: 160–200 degrees.

Boundary clusters and clusters containing any quality flags are rejected as
whole clusters. For 2+2 topology, both bijections are evaluated; exactly one
must pass. Ambiguous, unmatched and higher-multiplicity cases are counted.
Adjacent strips are flagged, not automatically merged or vetoed. Adjacent
1+2/2+1 cases are diagnostic categories, not established charge sharing.

Low/high labels mean the mean energy of the two faces. For compatibility with
the original diagnostic, the CM transformation hypothesizes low = recoil and
high = projectile with A=150/32, beam=102 MeV and deficit=2 MeV. These settings
are configurable. The deficit is a model parameter, not an independently
measured target loss, and the labels are not particle identification. Invalid
arcsine arguments are flagged rather than clipped. Geometry is read from the
complete original ring/sector angle tables and retained in the receipt.

All 24 rings are available for these diagnostics, including high-beam rings
outside the production Ge–S3 selection. There is no Ge requirement, addback,
Doppler correction, prompt/random subtraction or yield extraction here.
The input has already imposed a connected-time-cluster selection, so the
remaining time distribution is not an unbiased wide random-window spectrum.

## Outputs and plotting

The small output ROOT contains histograms, a `candidates` tree of all unique
2+2 assignments before subsequent particle-time/angular cuts, and `metadata`.
Per-candidate flags identify the later cuts. `cm_sum_deg=0` is a placeholder
when `cm_valid=false`; it is not a measured zero angle. Source hit IDs are -1
for legacy input. Accepted 1+1 events populate their own ring-energy histogram.
Histograms ending in `before` already pass the unique assignment and pair-time
gates; `after` adds both angular gates. The two energy–angle entries per
candidate and the four possible front/back combinations are not independent
reaction counts. Underflow and overflow remain in ROOT.

The JSON receipt records input UUID/size, selected entry ranges, a digest of
the read fields, input calibration provenance, geometry, gates and cutflow.
It does not claim a full byte hash of the large input file. Repeated outputs
are refused; unsuccessful partial files remain available for diagnosis.

For saved histogram plots only:

```bash
root -l -b -q 's3_coin_ana/draw_sample.cpp("ANALYSIS.root","NEW_FIGURE_DIR")'
bash s3_coin_ana/batch_draw_sample.sh NEW_FIGURE_DIR ANALYSIS.root [...]
```

The drawing macro checks input objects and uses log-y for 1D, log-z for 2D.
The project's local Chinese reading notebook is
`jupyter/codex_0234_150sm_s3_two_particle_coincidence.ipynb`; its plotting code
is embedded and requires only three small analysis ROOT outputs. The three
pilots use 2 million clusters each, not full runs or an all-173-run validation.
