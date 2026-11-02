"""Run joint-face pixel gain fits and export a complete, explicitly selected candidate."""
from pathlib import Path
import argparse
import concurrent.futures as futures
import hashlib
import json
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def fit_one(task, run, first, last, label):
    manifest = read(task / "input_manifest.json")
    baseline = Path(manifest["baseline_task"])
    assert sha(manifest["baseline_package"]) == manifest["baseline_sha256"]
    target = task / "fits" / label
    target.mkdir(parents=True, exist_ok=True)
    output = target / f"fit_{run}.json"
    assert not output.exists(), output
    command = [str(ROOT / "s3_energy_attenuation/conditioned_calibrate"), "--pixel",
               "--input", str(baseline / "histograms" / f"spectra_{run}.root"),
               "--reference", str(baseline / "histograms/spectra_0616.root"),
               "--baseline", manifest["baseline_package"], "--output", str(output),
               "--first-ring", str(first), "--last-ring", str(last)]
    start = time.time()
    with (task / "logs" / f"fit_{label}_{run}.log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
    receipt = {"run": run, "command": command, "returncode": result.returncode,
               "elapsed_seconds": time.time() - start, "output": str(output)}
    if result.returncode == 0:
        fit = read(output)
        receipt.update(accepted_pixels=fit["accepted_pixels"], pixels=len(fit["pixels"]), sha256=sha(output))
    write(task / "logs" / f"fit_{label}_{run}.json", receipt)
    if result.returncode:
        raise RuntimeError(f"Fit failed: {run}, see log")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--runs", nargs="*")
    parser.add_argument("--first-ring", type=int, default=15)
    parser.add_argument("--last-ring", type=int, default=24)
    parser.add_argument("--label", default="all_low")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    assert 1 <= args.jobs <= 8
    task = args.task.resolve()
    runs = args.runs or read(task / "input_manifest.json")["runs"]
    receipts = []
    with futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = [pool.submit(fit_one, task, run, args.first_ring, args.last_ring, args.label) for run in runs]
        for item in futures.as_completed(pending):
            row = item.result()
            receipts.append(row)
            print(row["run"], f'{row["accepted_pixels"]}/{row["pixels"]}',
                  f'{row["elapsed_seconds"]:.1f}s', f'{len(receipts)}/{len(runs)}', flush=True)
    write(task / f"stage_{args.label}.json", {"status": "complete", "results": sorted(receipts, key=lambda x: x["run"]),
          "driver_sha256": sha(__file__), "program_sha256": sha(ROOT / "s3_energy_attenuation/conditioned_calibrate")})


if __name__ == "__main__":
    main()
