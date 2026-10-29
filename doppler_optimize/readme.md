# Current parameter comparison

Use the shared C++ implementation in `doppler/angle_info.cpp` and the current
`doppler/doppler` executable. A scan JSON is an array of parameter objects with
the same schema as `doppler/parameters/angle_table_baseline.json`. An input list
contains one absolute path per line to a current `pairs` ROOT file.

```bash
./doppler/doppler --project "$PWD" \
  --scan-parameters /path/to/scan_parameters.json \
  --scan-list /path/to/train_pairs.list \
  --output /tmp/train_doppler_scan.root
```

Each input is read once for all candidates. All variants use identical pairs and
weights. Histograms contain signed net counts and squared-weight variances in
0.5-keV bins over 280–800 keV, including underflow and overflow. They are grouped
by period, particle family and Ge channel. The baseline pair files retain the
full energy range and the unsubtracted event information.

Supply the scan JSON and pair list for the intended comparison; the example
paths are placeholders. `evaluate.py` requires NumPy, pandas, SciPy and uproot;
`plot_pilot.py` also requires Matplotlib. The latter reproduces the historical
six-run pilot from a supplied task directory and its saved fit/scan receipts;
those experimental files are local research products, not repository inputs.

The retained processing choice is `doppler/parameters/angle_table_baseline.json`.
The tested beta scales and five-/seven-Ge angle variants remain research
candidates. The scan and fitting tools do not install or activate them.

`evaluate.py` uses a bin-integrated Gaussian core plus linear continuum to report
centroid, core FWHM, fit quality and a diagnostic score. A poor fit is flagged;
its formal core error is not a complete peak-shape or yield uncertainty. The
script selects pilot beta candidates from training runs only, retaining separate
held-out runs. No F3 or peak yield is exported. See the task protocol and reading
notebook for scope, limits and actual acceptance decisions.

# Historical scans

The directories below retain the original scan implementations and historical
notes. Their copied models and selections differ, and their batch scripts use
old build names and launch uncontrolled per-run jobs. Do not use them for current
production or combine their parameter choices as one consistent optimization.
The S3-beta comparison macro now correctly reads `e_dc_r` rather than `e_raw`.

# optimize doppler correction

## optimize_spider_phi

### including files

+ `src_raw` folder: change for loop different spider phi angles.
+ `try` folder: save different loop files.
+ `batch_clean.sh`: clean try files.
+ `run.sh, batch_run.sh`: main run batch file.
+ `d2hist.sh, batch_d2hist.sh`: from tree event to get histograms.
+ `batch_hadd_hist.sh`: add histograms.
+ `compare_draw.cpp`: compare different results.

## op_spider_theta

0 deg

## op_spider_cut

15 MeV

## op_spider_beta

e_out*1.0

## op_s3_phi

0 deg

## op_s3_theta

0 deg

## op_s3_cut

+ 25 MeV for exp1
+ 5 MeV for exp2

## op_s3_beta

e_out*0.5

## op_ge_phi

0 deg

## op_ge_theta

0 deg

## op_coin_win

100 ns

