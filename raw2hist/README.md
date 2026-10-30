# Channel singles spectra

`histogram` converts decoded ROOT hits into channel ADC and calibrated-energy
histograms. It is a calibration/monitoring diagnostic. It does not build events,
apply timing or coincidence gates, subtract random coincidences, fit peaks or
correct Doppler shifts. Normal and failed CFD hits are combined.

## Build and run

```bash
make -C raw2hist
./raw2hist/histogram --run 0617 \
  --input /mnt/e/Sm150/rootfile/data_C2_0617_0617.root \
  --output rootfile/raw2hist_check/data0617_hist.root \
  --sample-entries 10000000 --segments 20
```

The output and its JSON receipt must not already exist. Their `.partial` names
are also reserved. Existing files are never overwritten. A failed write retains
partial files for inspection. Output publication uses local filesystem hard
links; use the project filesystem for results (DrvFS volumes may not support
this operation).

Omitting `--input` selects the current E-drive C2 path above for the specified
run. An explicit decoded C1 path is also accepted. The project is found from the
compiled executable; `--project DIR` overrides it. All output paths are explicit.

`--sample-entries 0` (the default) reads the whole input. A positive number selects
centered, disjoint entry ranges across the file. Their exact boundaries and
fraction are recorded. Entry sampling is not necessarily uniform in live time
or integrated beam intensity. Source UUID, size, modification time and a SHA256
of the inspected branch values identify the sample; this digest is **not** a
full byte hash of the source ROOT.

Other options:

- `--mode both|raw|calibrated`: save ADC, energy, or both (default).
- `--detector all|ge|particle|si|s3`: select electronic crate or physical detector.
- `--format decoded`: `tree` with scalar `cid/sid/ch` (`Short_t`), `evte`
  (`UShort_t`), `pileup/outofr` (`Bool_t`). Missing or incompatible fields fail.
- `--format calibrated --mode calibrated`: explicitly accept historical
  `Double_t evte` in keV. Fill it directly without another calibration; quality
  flags and its original calibration version are unspecified. Current parameter
  files supply only physical channel mappings in this mode.

Decoded `tr_event`/`tr_bg` build files cannot reconstruct ungated singles and are
rejected. Floating energy is never silently treated as ADC.

## Formal calibration and detector mapping

The shared `include/sm150_event_calibration.h` supplies the same physical mapping
and formal energy coefficients used by the current build diagnostics: 37 Ge,
96 Si array, 24 S3 ring and 32 S3 sector channels. Periods are high I 0472--0605,
low 0616--0674 and high II 0680--0718. Raw-only `all/ge/particle` mode can inspect
source runs outside these intervals without an energy calibration.

```text
Ebase  = hardware_gain * (c0 + c1*ADC + c2*ADC^2)
EpreFB = (c0 + c1*ADC + c2*ADC^2) * hardware_gain * kGe * attenuation
Efinal = fb0 + fb1*EpreFB
```

Energies and offsets are in keV. Ge uses `pars/cali_ge.txt` and
`pars/run_k/k_RUN.txt`; particle calibration uses the period-specific
`cali_si_r460.txt`, `cali_si_r457_low_current.txt` or
`cali_si_r457_high_current.txt`. S3 also reads all 56
`pars/run_s3_att/{ring,sector}NN.txt` files and
`pars/run_fb_cor/correction_{ring,sector}_RUN.txt`. Hardware gain factors are
applied once. No random dither is added. Parameters, complete source text and
SHA256 are embedded in the result and verified unchanged after reading.

The `s3_fb_cor/parameters/joint_energy_20261030/` A/B/q package remains an explicit
candidate, not the formal default. Its sector correction may depend on a matched
ring. Singles do not supply that identity, so this program does not apply that
candidate or label a partial A/B correction as its complete result. Use matched
events and `build_all --s3-energy-file` for that comparison.

## Histograms and selection

Each `chNNNNN/` directory uses `key = 10000*cid + 100*sid + ch`:

| Histogram | Contents | Binning |
| --- | --- | --- |
| `adc_all` | Every inspected hit at this address, including quality flags | 2 ADC channels/bin, [0,65536) |
| `adc_clean` | Exclude pileup, out-of-range or ADC saturation (65535) | Same |
| `energy_base` | Clean hits after basic calibration and hardware gain | Ge: 1 keV/bin to 4096 keV; particles: 50 keV/bin to 100 MeV |
| `energy_pre_fb` | S3 clean hits after attenuation, before front/back normalization | 50 keV/bin |
| `energy_final` | Clean hits with the full formal energy chain | Same as `energy_base` |

Quality rejection is a union: a hit with multiple flags is counted once in
`rejected`; individual flag counters can overlap. No additional energy cut is
applied. Underflow/overflow counts are retained and reported. Histogram entries
and all counts include these flows; displayed plots can have narrower ranges.
All 189 physical channels are created, including zero-statistics channels.
Valid unused electronic addresses retain ADC histograms and an unmapped count;
they receive no invented energy calibration. Invalid electronic addresses are
counted separately. All 24 S3 rings are diagnostics, including high-beam outer
rings that are not used in the production Ge--S3 selection.

The output contains `metadata` plus an identical `.root.json` receipt with run,
sampling, counts, parameter records, input identity and program hashes. It holds
histograms only, not another large event tree. The accompanying notebook uses
these small outputs and keeps its figures embedded for reading without E.

## Batch and compatibility macros

```bash
./raw2hist/batch.sh RUN_LIST /mnt/e/Sm150/rootfile NEW_OUTPUT_DIR 2 10000000
```

All five arguments are required. The list is explicit, duplicate runs and
missing inputs fail before launch, and the output directory must be new. The
script stops submitting after a detected failure and waits for started jobs.
It no longer launches 120 plotting processes or implicitly processes `run.txt`.

Old `get_hist*.cpp` and `macro_*.cpp` names remain as thin wrappers around this
executable. `macro_raw`, `macro_cali`, `macro_cali_ge`, `macro_cali_si` use E C2 by
default and optionally take an input filename; the last selects the entire
particle crate (Si array and S3). The Ge function name now matches its filename.
`get_hist_cali(filename,run)` and `get_hist(run)` explicitly retain the historical
floating-energy format. Wrapper calls read the **whole** selected input and
refuse their existing output under `rootfile/raw2hist/`; use the CLI for sampling.

Redraw an existing small output from the project root:

```bash
root -l -b -q 'raw2hist/draw_sample.cpp("RESULT.root","NEW_FIGURE_DIR")'
```

This exports one log-scale ADC/energy figure per address, without reading E.
The Chinese reading notebook additionally groups spectra by physical detector,
ring and sector and lists the actual coefficients.

Validation and the three-period example are recorded locally in
`codex/analysis/raw2hist_modernization_20261030/`; the historical macros were
preserved there before this update. Formal parameters and existing event ROOT
files are unchanged.
