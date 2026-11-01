# S3 calibration-event extraction

The current C++ entry is `build_s3/build`. It reads decoded C2 events, applies the complete particle time correction from `pars/time_alignment/active.txt`, and creates particle-only S3 calibration pairs. `set.h` is historical and no longer controls this executable.

From the project root:

```sh
make -C build_s3
build_s3/build 0617 --input /mnt/e/Sm150/rootfile/data_C2_0617_0617.root --output /path/to/new/s3_0617.root
```

Optional settings are shown by `build_s3/build --help`. The default energy threshold is 2500 keV before front/back normalization. Timing uses calibrated ADC and the complete current offset once; no random jitter is added. Native clock units are converted according to the recorded sampling rate. Current gain and attenuation coefficients are retained, and both pre-normalization and previous normalized energy are stored.

Hits are sorted with a 20000 ns buffer. Consecutive corrected-time gaps of at most 100 ns define connected clusters; a chain may span more than 100 ns. The program records all above-threshold S3 hits in two-sided clusters, without a ring/sector energy-agreement cut or adjacent-strip addback. Multiplicity means the original number of above-threshold hits, including quality-flagged hits. Exactly one hit on each face also produces a `s3_pairs` row. First/last clusters and event-quality flags remain explicit for downstream rejection. No Ge coincidence is required.

Outputs contain `s3_events`, `s3_pairs`, and JSON `metadata`, plus a JSON receipt. Source entry numbers, ADC, native/corrected timestamps, complete calibration coefficients, parameter hashes, source ROOT UUID/stat, counts, and selection settings permit readback. The original C2 file is not modified. Source stat and frozen parameter hashes are checked at completion; this is not a new whole-file input SHA256 check.

An existing output, receipt, or partial file is refused. A `.partial` ROOT is published only after completion checks. `--entries N` deliberately reads only a prefix and is recorded as partial coverage; omit it for a whole run.

For an explicit batch:

```sh
bash build_s3/batch.sh runs.txt /mnt/e/Sm150/rootfile /path/to/new/s3_sources 1
```

The output directory must be new. The final argument is the worker count (default 1); there is no implicit 120-worker launch or automatic full-run list. One raw-disk reader is a sensible starting point. Fit these outputs with the compiled program in `s3_fb_cor/`; historical `tr_event` macros use a different schema.
# Distributed joint-calibration samples

`./build RUN --sample-entries 2000000 --output NEW.root` selects centered
contiguous samples in 40 equal source-entry intervals. `--sample-entries 0`
reads all intervals; interval-edge guards still apply. This option is separate
from the legacy full-cluster mode and cannot be combined with `--entries`.

The compact `pairs` tree records the leading hit of each 100 ns face group,
all ring/sector pairs within +/-1000 ns, uncorrected calibrated energy `E0`,
previous production energy, native source identities, quality flags,
multiplicity, minor-hit energy fraction, and previous same-channel intervals.
There is no addback or front/back energy-difference cut. The 2.5 MeV extraction
threshold uses the frozen old pre-front/back energy. The previous-hit interval
includes all recorded S3 hits and is censored at the sampled interval start.
The 20 us boundary guard protects local grouping and random-sideband analysis.

The metadata stores the exact sampled intervals, a hash of the inspected decoded
columns, complete calibration/timing identities, counts, and program hash.
Existing outputs and receipts are never overwritten. This mode does not modify
production parameters or build ROOT files.

