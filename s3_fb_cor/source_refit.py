"""Re-express source calibration, refit A/B, and preserve held-out decisions.

All outputs go to an explicit new task directory. This does not install runtime
parameters or rebuild production events. Native C++ programs do event I/O/fits.
"""
from pathlib import Path
import argparse
import concurrent.futures as futures
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from s3_energy_attenuation.analysis import differences, metrics

REFERENCES = {"high1": "0472", "low": "0616", "high2": "0680"}


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    partial.replace(path)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def native(task, name, command, output):
    """Resume only outputs whose task-owned command receipt is consistent."""
    receipt = task / "logs" / (name + ".json")
    command = [str(v) for v in command]
    signature = {"command": command, "program_sha256": sha(command[0])}
    if output.exists():
        saved = read(receipt)
        assert all(saved[k] == v for k, v in signature.items())
        assert saved["returncode"] == 0 and saved["output_sha256"] == sha(output)
        return
    receipt.parent.mkdir(exist_ok=True)
    start = time.time()
    with (task / "logs" / (name + ".log")).open("w") as log:
        result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    write(receipt, {**signature, "returncode": result.returncode,
                    "elapsed_seconds": time.time() - start,
                    "output_sha256": sha(output) if output.exists() else None})
    if result.returncode:
        raise RuntimeError(f"{name} failed; see {task / 'logs' / (name + '.log')}")


def prepare_one(task, row):
    run = row["run"]
    source_run = "0460" if int(run) <= 605 else "0457"
    output = task / "pairs" / f"pairs_{run}.root"
    native(task, f"prepare_{run}", [ROOT / "s3_fb_cor/recalibrate_source",
        "--input", row["source"], "--output", output,
        "--source-calibration", ROOT / f"cali_si/parameters/source_{source_run}_peak_centers_20261102/candidate.json",
        "--baseline", ROOT / "pars/s3_energy/joint_energy_20261030/coefficients.json",
        "--expected-sha256", row["source_sha256"]], output)
    hist = task / "histograms" / f"spectra_{run}.root"
    native(task, f"histogram_{run}", [ROOT / "s3_energy_attenuation/conditioned_extract",
        "--input", output, "--output", hist], hist)
    return run


def ab_one(task, row):
    run = row["run"]
    output = task / "ab_fits" / f"fit_{run}.json"
    native(task, f"ab_{run}", [ROOT / "s3_fb_cor/joint_calibrate",
        "--input", task / "pairs" / f"pairs_{run}.root", "--output", output,
        "--min-pairs", "20"], output)
    return run


def model_map(channels, q):
    a, b = np.zeros(56), np.ones(56)
    for c in channels:
        a[c["index"]], b[c["index"]] = c["A_keV"], c["B"]
    return {"a": a, "b": b, "radial": np.array(q)}


def summary(v, model, baseline, fold):
    raw, normalized, _, _ = differences(v, model, baseline)
    mask = ((v["dt_ns"] >= -80) & (v["dt_ns"] < 80)
            & (v["r_old_keV"] >= 15000) & (v["r_old_keV"] < 100000)
            & (v["s_old_keV"] >= 15000) & (v["s_old_keV"] < 100000))
    if fold >= 0:
        mask &= v["block"] % 4 == fold
    channels = []
    for i in range(56):
        use = mask & (v["r"] == i if i < 24 else v["s"] == i)
        channels.append({"index": i, **metrics(normalized[use])})
    return {"normalized": metrics(normalized[mask]), "raw": metrics(raw[mask]), "channels": channels}


def ab_guard(current, base, confirmation=False):
    x, y = current["normalized"], base["normalized"]
    if x["n"] < 100:
        return ["fewer_than_100_global_pairs"]
    reasons = []
    if x["width68_keV"] > (1.03 if confirmation else 1.02) * y["width68_keV"]:
        reasons.append("width_worse")
    if x["score_keV"] > (1.0 if confirmation else .99) * y["score_keV"]:
        reasons.append("score_not_improved")
    if x["tail2MeV"] > y["tail2MeV"] + .01:
        reasons.append("tail_increase")
    for a, b in zip(current["channels"], base["channels"]):
        if a["n"] >= 100 and a["score_keV"] > 1.1 * b["score_keV"] + 20:
            reasons.append(f"channel_{a['index']}_worse")
    return reasons


def evaluate_ab_one(task, row):
    run = row["run"]
    fit = read(task / "ab_fits" / f"fit_{run}.json")
    meta = fit["sample_metadata"]
    source = meta["source_recalibration"]
    base = model_map(source["channels"], source["baseline_q"])
    maps = {"equivalent": base}
    for name, model in fit["models"]["fit"].items():
        maps[name] = model_map(model["coefficients"], source["baseline_q"])
    branches = ["block", "r", "s", "dt_ns", "r_e0_keV", "s_e0_keV", "r_old_keV", "s_old_keV"]
    with uproot.open(task / "pairs" / f"pairs_{run}.root",
                     handler=uproot.source.file.MemmapSource) as f:
        v = f["pairs"].arrays(branches, library="np")
    checks = {str(fold): {name: summary(v, model, base, fold) for name, model in maps.items()}
              for fold in [1, 3, -1]}
    preferred = "gain"
    one = checks["1"]
    if (one["affine"]["normalized"]["n"] >= 100
            and one["affine"]["normalized"]["score_keV"] < .99 * one["gain"]["normalized"]["score_keV"]
            and one["affine"]["normalized"]["width68_keV"] <= 1.02 * one["gain"]["normalized"]["width68_keV"]):
        preferred = "affine"
    model = fit["models"]["fit"][preferred]
    connected = all(c["status"] in ["fixed_energy_reference", "global_fit"]
                    for c in model["coefficients"] if c["index"] < fit["active_rings"] or c["index"] >= 24)
    reasons = ab_guard(one[preferred], one["equivalent"])
    if not model["ok"] or not connected:
        reasons.append("incomplete_or_invalid_fit")
    confirmation = ab_guard(checks["3"][preferred], checks["3"]["equivalent"], True) if not reasons else []
    chosen = preferred if not reasons and not confirmation else "equivalent"
    # Reference energy maps are fixed before comparing any other run.
    if run == REFERENCES[fit["period"]]:
        chosen = "equivalent"
        reasons.append("fixed_period_reference")
    stability = []
    for name in ["gain", "affine"]:
        first = model_map(fit["models"]["fit_early_alternating"][name]["coefficients"], source["baseline_q"])
        second = model_map(fit["models"]["fit_late_alternating"][name]["coefficients"], source["baseline_q"])
        for i in range(56):
            stability.append({"index": i, "model": name, "max_prediction_difference_keV":
                float(np.max(np.abs(first["a"][i] - second["a"][i] + (first["b"][i] - second["b"][i]) * np.array([25000, 50000, 75000]))))})
    result = {"run": run, "period": fit["period"], "active_rings": fit["active_rings"],
        "reference_run": REFERENCES[fit["period"]], "preferred": preferred, "chosen": chosen,
        "selection_flags": reasons, "confirmation_flags": confirmation, "checks": checks,
        "stability": stability, "maps": {k: {n: v.tolist() for n, v in m.items()} for k, m in maps.items()},
        "fit_sha256": sha(task / "ab_fits" / f"fit_{run}.json")}
    write(task / "ab_evaluation" / f"evaluation_{run}.json", result)
    return run


def assemble_ab(task, rows):
    package = {"schema": "sm150_s3_source_rebased_ab_selection_v1", "runs": {}}
    for row in rows:
        run = row["run"]
        e = read(task / "ab_evaluation" / f"evaluation_{run}.json")
        model, old = e["maps"][e["chosen"]], e["maps"]["equivalent"]
        channels = [{"index": i, "channel": f"ring{i+1:02d}" if i < 24 else f"sector{i-23:02d}",
                    "A_keV": model["a"][i], "B": model["b"][i],
                    "old_A_keV": old["a"][i], "old_B": old["b"][i]} for i in range(56)]
        package["runs"][run] = {"period": e["period"], "reference_run": e["reference_run"],
            "active_rings": e["active_rings"], "chosen_ab": e["chosen"], "channels": channels,
            "baseline_q": old["radial"]}
    write(task / "ab_selected.json", package)


def q_one(task, row):
    run = row["run"]
    package = read(task / "ab_selected.json")
    ref = package["runs"][run]["reference_run"]
    output = task / "q_fits" / f"fit_{run}.json"
    native(task, f"q_{run}", [ROOT / "s3_energy_attenuation/conditioned_calibrate",
        "--input", task / "histograms" / f"spectra_{run}.root",
        "--reference", task / "histograms" / f"spectra_{ref}.root",
        "--baseline", task / "ab_selected.json", "--output", output], output)
    return run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["prepare", "ab", "evaluate-ab", "q"])
    p.add_argument("--task", type=Path, required=True)
    p.add_argument("--runs", nargs="*")
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    assert 1 <= args.jobs <= 8
    task = args.task.resolve()
    all_rows = read(task / "source_manifest.json")
    rows = [r for r in all_rows if not args.runs or r["run"] in args.runs]
    if args.runs:
        assert len(rows) == len(args.runs)
    jobs = min(args.jobs, 2) if args.stage == "prepare" else args.jobs
    fn = {"prepare": prepare_one, "ab": ab_one, "evaluate-ab": evaluate_ab_one, "q": q_one}[args.stage]
    start = time.time()
    with futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        pending = {pool.submit(fn, task, row): row["run"] for row in rows}
        for k, done in enumerate(futures.as_completed(pending), 1):
            run = done.result()
            print(f"{args.stage} {k}/{len(rows)} {run} elapsed={time.time()-start:.1f}s", flush=True)
    if args.stage == "evaluate-ab" and len(rows) == len(all_rows):
        assemble_ab(task, all_rows)
    write(task / f"stage_{args.stage.replace('-', '_')}_{len(rows)}.json",
          {"stage": args.stage, "runs": [r["run"] for r in rows], "elapsed_seconds": time.time()-start,
           "jobs": jobs, "program_sha256": sha(__file__)})


if __name__ == "__main__":
    main()
