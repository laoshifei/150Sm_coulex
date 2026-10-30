"""Explicit inputs, provenance and protected outputs for the S3 workflows."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import hashlib
import json
import os

PATH_KEYS = {
    "pairs": "source_sha256", "spectra": "spectra_sha256",
    "conditional": "conditional_sha256", "coupled": "coupled_sha256",
    "fixed_fit": "fixed_sha256", "edge_fit": "edge_sha256",
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def resolve(path, parent):
    path = Path(path)
    return path.resolve() if path.is_absolute() else (Path(parent) / path).resolve()


def load_config(path):
    path = Path(path).resolve()
    cfg = read_json(path)
    keys = ["baseline", "previous", "manifest", "reference_manifest", "protocol",
            "confirmation_manifest", "selection_freeze", "confirmation_dir",
            "review_rows", "evaluation_dir"]
    for key in keys:
        if key in cfg:
            cfg[key] = resolve(cfg[key], path.parent)
    cfg["config_path"] = path
    return cfg


def load_manifest(path):
    path = Path(path).resolve()
    data = read_json(path)
    rows = data if isinstance(data, list) else data["runs"]
    result = []
    seen = set()
    for source in rows:
        row = dict(source)
        run = row["run"]
        if not isinstance(run, str) or len(run) != 4 or not run.isdigit() or run in seen:
            raise ValueError(f"Expected unique four-digit run strings: {run!r}")
        seen.add(run)
        for key in PATH_KEYS:
            if key in row:
                row[key] = str(resolve(row[key], path.parent))
        result.append(row)
    return result


def choose_runs(rows, runs=None):
    if not rows:
        raise ValueError("Empty run list")
    if runs is None:
        return rows
    requested = {str(r).zfill(4) for r in runs}
    missing = requested - {r["run"] for r in rows}
    if missing or not requested:
        raise ValueError(f"Unknown or empty run selection: {sorted(missing)}")
    return [r for r in rows if r["run"] in requested]


def verify_row(row, keys):
    for key in keys:
        expected = row.get(PATH_KEYS[key])
        if not expected or digest(row[key]) != expected:
            raise ValueError(f"Missing or mismatched SHA256: run {row['run']} {key}")


def validate_baseline(document):
    import numpy as np
    for run, row in document["runs"].items():
        channels = row["channels"]
        if [c["index"] for c in channels] != list(range(56)):
            raise ValueError(f"Expected ordered 56-channel baseline for {run}")
        if row["reference_run"] not in document["runs"]:
            raise ValueError(f"Missing reference for {run}")
        for c in channels:
            if not np.isfinite([c["A_keV"], c["B"]]).all() or c["B"] <= 0:
                raise ValueError(f"Invalid baseline coefficient for {run}")


def validate_model(model):
    import numpy as np
    for key, size in [("a", 56), ("b", 56), ("radial", 24)]:
        a = np.asarray(model[key])
        if a.shape != (size,) or not np.isfinite(a).all():
            raise ValueError(f"Invalid {key} shape or coefficient")
        if key != "a" and not (a > 0).all():
            raise ValueError(f"Nonpositive {key}")


def validate_conditional(row, baseline):
    conditional = read_json(row["conditional"])
    base = baseline["runs"][row["run"]]
    if (conditional["run"] != row["run"] or conditional["period"] != base["period"]
            or conditional["reference_run"] != base["reference_run"]
            or [c["index"] for c in conditional["channels"]] != list(range(56))):
        raise ValueError(f"Mismatched conditional run/reference/channel map: {row['run']}")
    return conditional


def new_output(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=False)
    return path


def code_hashes():
    root = Path(__file__).resolve().parents[1]
    files = [p for folder in ["s3_energy_attenuation", "s3_fb_cor"]
             for p in (root / folder).glob("*.py")]
    return {str(p): digest(p) for p in sorted(files)}


def input_hashes(cfg):
    return {str(v): digest(v) for v in cfg.values()
            if isinstance(v, Path) and v.is_file()}


def run_batch(rows, worker, output, jobs=1):
    if not 1 <= jobs <= 8:
        raise ValueError("jobs must be between 1 and 8")
    results = output / "results"
    results.mkdir()
    state = {"status": "running", "pid": os.getpid(), "runs": [],
             "started": datetime.now().astimezone().isoformat()}
    write_json(output / "state.json", state)

    def one(row):
        value = worker(row)
        write_json(results / f"{row['run']}.json", value)
        return {k: value[k] for k in ["run", "status", "preferred", "selected", "final"]
                if k in value}

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(one, row): row["run"] for row in rows}
        for future in as_completed(futures):
            try:
                result = future.result()
            except Exception as error:
                result = {"run": futures[future], "status": "failed", "error": repr(error)}
            state["runs"].append(result)
            write_json(output / "state.json", state)
            print(json.dumps(result), flush=True)
    state["status"] = "complete" if all(r["status"] == "complete" for r in state["runs"]) else "failed"
    state["finished"] = datetime.now().astimezone().isoformat()
    write_json(output / "state.json", state)
    if state["status"] != "complete":
        raise RuntimeError(f"Batch incomplete; inspect {output / 'state.json'}")
    return state


def write_receipt(output, cfg, command, extra=None):
    write_json(output / "provenance.json", {
        "command": command, "inputs_sha256": input_hashes(cfg),
        "code_sha256": code_hashes(), "new_raw_event_reads": 0,
        "production_adopted": False, **(extra or {}),
    })


def write_manifest(output):
    write_json(output / "manifest.json", {
        str(p.relative_to(output)): digest(p) for p in sorted(output.rglob("*"))
        if p.is_file() and p != output / "manifest.json"
    })


def load_frozen(path, protocol):
    """Read a new freeze receipt or a preserved 0230 receipt."""
    path = Path(path).resolve()
    frozen = read_json(path)
    if digest(protocol) != frozen["protocol_sha256"]:
        raise ValueError("Frozen protocol changed")
    for group in ["files", "code_sha256", "inputs_sha256"]:
        for name, sha in frozen.get(group, {}).items():
            if digest(resolve(name, path.parent)) != sha:
                raise ValueError(f"Frozen source changed: {name}")
    choices = {r["run"]: r["preferred"] for r in frozen["choices"]}
    if len(choices) != len(frozen["choices"]):
        raise ValueError("Duplicate frozen runs")
    selections = {}
    for name in frozen["files"]:
        filename = resolve(name, path.parent)
        row = read_json(filename)
        run = row["run"]
        if row["status"] != "complete" or row["preferred"] != choices[run] or run in selections:
            raise ValueError(f"Invalid frozen selection: {run}")
        selections[run] = (filename, row)
    if set(selections) != set(choices):
        raise ValueError("Incomplete frozen selections")
    return selections


def verify_disjoint(training_path, confirmation_path):
    """Check actual entry intervals; this does not assert experiment independence."""
    old = read_json(str(training_path) + ".json")
    new = read_json(str(confirmation_path) + ".json")
    import uproot
    for path, receipt in [(training_path, old), (confirmation_path, new)]:
        with uproot.open(path, handler=uproot.source.file.MemmapSource) as root:
            metadata = json.loads(str(root["metadata"]))
        for key in ["run", "input_uuid", "input_entries", "segments"]:
            if metadata[key] != receipt[key]:
                raise ValueError(f"Pair receipt differs from ROOT metadata: {path} {key}")
    for key in ["run", "input_uuid", "input_entries"]:
        if old[key] != new[key]:
            raise ValueError(f"Different source identity: {key}")
    for a in old["segments"]:
        lo, hi = a["read_interval"]
        for b in new["segments"]:
            start, stop = b["read_interval"]
            if max(lo, start) < min(hi, stop):
                raise ValueError(f"Overlapping selection/confirmation intervals: {old['run']}")
