"""Freeze a pixel candidate and evaluate its spectra without changing fitted values."""
from pathlib import Path
import argparse
import copy
import datetime
import json
import shutil
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from s3_energy_attenuation.refit_pixels import read, sha, write
from s3_energy_attenuation.analysis import spectrum_file, interpolate_rows

GRID = np.arange(20000., 90000.1, 1000.)
DISPLAY_EDGES = np.arange(0., 120000.1, 500.)


def export(task, label):
    manifest = read(task / "input_manifest.json")
    original = read(ROOT / "s3_fb_cor/parameters/source_recalibrated_20261102/coefficients.json")
    baseline = read(manifest["baseline_package"])
    assert original["runs"] == baseline["runs"]
    package = copy.deepcopy(original)
    package.update(schema="sm150_s3_pixel_energy_v2", created=datetime.datetime.now().astimezone().isoformat(),
                   formal_replacement=False, scope="Low beam rings15-24 jointly fitted per ring-sector pixel; other settings inherited",
                   formula="E=q_pixel[matched_ring][matched_sector]*(A+B*E0); q_position is a fallback, never an extra multiplier",
                   baseline_sha256=manifest["baseline_sha256"])
    rows, fit_hashes = [], {}
    for run in manifest["runs"]:
        file = task / "fits" / label / f"fit_{run}.json"
        fit = read(file)
        assert fit["baseline_sha256"] == manifest["baseline_sha256"]
        fit_hashes[run] = sha(file)
        model = package["runs"][run]
        matrix = np.repeat(np.asarray(model["q_position"], float)[:, None], 32, axis=1)
        for r in fit["pixels"]:
            matrix[r["ring"]-1, r["sector"]-1] = r["q_selected"]
            row = {k: r.get(k) for k in ["run", "ring", "sector", "baseline_q", "q_fitted", "q_selected", "status",
                                         "q_fit_fold0", "q_fit_fold2", "q_alternate_window", "split_log_difference", "window_log_difference", "face_log_difference"]}
            row["flags"] = ";".join(r["flags"])
            for fold in ["training", "selection", "confirmation"]:
                for key in ["baseline_D", "candidate_D", "supported"]:
                    row[fold + "_" + key] = r.get(fold, {}).get(key)
            rows.append(row)
        model.update(q_pixel=matrix.tolist(), requires_pixel_context=True,
                     requires_ring_context=bool(np.any(matrix != 1.)),
                     pixel_source=str(file.relative_to(ROOT)), pixel_source_sha256=sha(file))
        for c in model["channels"]:
            c["q_application"] = "Use the common matched ring-sector q_pixel; q column is inherited ring-only baseline"
    dest = task / "candidate"
    dest.mkdir(exist_ok=False)
    shutil.copytree(Path(manifest["baseline_task"]) / "source_calibration", dest / "source_calibration")
    write(dest / "coefficients.json", package)
    pd.DataFrame(rows).to_csv(dest / "pixel_factors.csv", index=False, float_format="%.17g")
    write(task / "candidate_frozen.json", {"frozen_at": package["created"], "fit_label": label,
          "coefficients_sha256": sha(dest / "coefficients.json"), "pixel_factors_sha256": sha(dest / "pixel_factors.csv"),
          "fit_sha256": fit_hashes, "fresh_spectra_used_in_selection": False})
    print("Frozen", len(rows), "pixel decisions; accepted", sum(r["status"] == "pixel_refit_accepted" for r in rows), flush=True)


def fold_cdf(values, face, region, fold):
    a = values[(face, region)]
    return a.sum(axis=0) if fold < 0 else a[fold]


def pixel_q(model, pixels):
    matrix = np.asarray(model.get("q_pixel", np.repeat(np.asarray(model["q_position"])[:, None], 32, axis=1)), float)
    return matrix[pixels//32, pixels%32]


def projected(values, model, face, fold, pixels, grid):
    index = pixels//32 if face == "ring" else 24+pixels%32
    q = pixel_q(model, pixels)
    a = q*np.asarray(model["A_keV"])[index]
    b = q*np.asarray(model["B"])[index]
    raw = (grid[None, :]-a[:, None])/b[:, None]
    pc = fold_cdf(values, face, "prompt", fold)[pixels]
    rc = fold_cdf(values, face, "random", fold)[pixels]
    p, r = interpolate_rows(pc, raw), interpolate_rows(rc, raw)
    total = pc[:, -1]-.5*rc[:, -1]
    net = p-.5*r
    den = np.where(total > 0, total, 1.)
    return {"F": net/den[:, None], "total": total, "P": p, "R": r,
            "supported": (total >= 15) & (b > 0), "a": a, "b": b}


def energy_quantiles(values, model, face, fold, pixels):
    pc = fold_cdf(values, face, "prompt", fold)[pixels]
    rc = fold_cdf(values, face, "random", fold)[pixels]
    net = pc-.5*rc
    total = net[:, -1]
    rawq = np.full((len(pixels), 5), np.nan)
    for i, line in enumerate(net):
        if total[i] < 15:
            continue
        for j, quantile in enumerate([.10, .16, .50, .84, .90]):
            target = total[i]*quantile
            crossings = np.flatnonzero((line[:-1] <= target) & (line[1:] > target))
            if len(crossings) == 1:
                k = crossings[0]
                rawq[i, j] = 250.*(k+(target-line[k])/(line[k+1]-line[k]))
    idx = pixels//32 if face == "ring" else 24+pixels%32
    q = pixel_q(model, pixels)
    return q[:, None]*(np.asarray(model["A_keV"])[idx, None]+np.asarray(model["B"])[idx, None]*rawq)


def evaluate_set(task, origin, package, rows, curves):
    manifest = read(task / "input_manifest.json")
    baseline = read(manifest["baseline_package"])
    source = Path(manifest["baseline_task"]) / "histograms" if origin == "historical" else task / "fresh"
    runs = manifest["runs"] if origin == "historical" else manifest["pilot"]
    fold = 3 if origin == "historical" else -1
    _, reference = spectrum_file(source / "spectra_0616.root")
    refmodel = baseline["runs"]["0616"]
    for run in runs:
        _, values = spectrum_file(source / f"spectra_{run}.root")
        oldmodel, newmodel = baseline["runs"][run], package["runs"][run]
        for ring in range(1, 25):
            pixels = np.arange((ring-1)*32, ring*32)
            for face in ["ring", "sector"]:
                models = {"reference": (reference, refmodel), "baseline": (values, oldmodel), "pixel": (values, newmodel)}
                mapped = {name: projected(v, m, face, fold, pixels, GRID) for name, (v, m) in models.items()}
                use = mapped["reference"]["supported"] & mapped["baseline"]["supported"] & mapped["pixel"]["supported"]
                row = {"origin": origin, "run": run, "ring": ring, "face": face, "usable_pixels": int(use.sum()), "missing_pixels": int((~use).sum())}
                for name in ["baseline", "pixel"]:
                    delta = mapped[name]["F"]-mapped["reference"]["F"]
                    row[name+"_D"] = float(np.sqrt(np.mean(delta[use, 1:-1]**2))) if use.any() else None
                features = {name: energy_quantiles(v, m, face, fold, pixels) for name, (v, m) in models.items()}
                for name in ["baseline", "pixel"]:
                    for j, q in enumerate([10,16,50,84,90]):
                        delta = features[name][use, j]-features["reference"][use, j]
                        row[f"{name}_q{q}_median_abs_keV"] = float(np.nanmedian(np.abs(delta))) if np.isfinite(delta).any() else None
                    width = features[name][:,3]-features[name][:,1]
                    refwidth = features["reference"][:,3]-features["reference"][:,1]
                    delta = width[use]-refwidth[use]
                    row[name+"_central68_width_median_abs_keV"] = float(np.nanmedian(np.abs(delta))) if np.isfinite(delta).any() else None
                for name, (v, m) in models.items():
                    c = projected(v, m, face, fold, pixels, DISPLAY_EDGES)
                    dp, dr = np.diff(c["P"], axis=1), np.diff(c["R"], axis=1)
                    denominator = np.where(c["total"] > 0, c["total"], 1.)
                    density = (dp-.5*dr)/denominator[:,None]/.5
                    key = f"{origin}_{run}_{ring:02d}_{face}_{name}"
                    curves[key+"_equal_pixels"] = density[use].mean(axis=0) if use.any() else np.zeros(len(DISPLAY_EDGES)-1)
                    curves[key+"_pooled"] = (dp-.5*dr).sum(axis=0)/max(c["total"].sum(), 1.)/.5
                    curves[key+"_P"] = dp.sum(axis=0)
                    curves[key+"_R"] = dr.sum(axis=0)
                    row[name+"_net_total"] = float(c["total"].sum())
                    # The display interval is intentionally finite; retain both flow terms.
                    row[name+"_net_below_display"] = float((c["P"][:,0]-.5*c["R"][:,0]).sum())
                    row[name+"_net_above_display"] = float((c["total"]-c["P"][:,-1]+.5*c["R"][:,-1]).sum())
                    assert np.isclose((dp-.5*dr).sum()+row[name+"_net_below_display"]+row[name+"_net_above_display"], c["total"].sum(), atol=1e-6)
                rows.append(row)
        print("Evaluated", origin, run, flush=True)


def evaluate(task):
    frozen = read(task / "candidate_frozen.json")
    assert sha(task / "candidate/coefficients.json") == frozen["coefficients_sha256"]
    package = read(task / "candidate/coefficients.json")
    rows, curves = [], {"edges_keV": DISPLAY_EDGES}
    evaluate_set(task, "historical", package, rows, curves)
    evaluate_set(task, "fresh", package, rows, curves)
    pd.DataFrame(rows).to_csv(task / "spectral_metrics.csv", index=False, float_format="%.17g")
    np.savez_compressed(task / "display_spectra.npz", **curves)
    df = pd.DataFrame(rows)
    chosen = df[(df.ring >= 15) & (df.run != "0616")]
    summary = []
    for (origin, run), group in chosen.groupby(["origin", "run"]):
        summary.append({"origin": origin, "run": run,
                        "baseline_D": float(np.sqrt(np.mean(group.baseline_D**2))),
                        "pixel_D": float(np.sqrt(np.mean(group.pixel_D**2))),
                        "improved_face_ring": int((group.pixel_D < group.baseline_D).sum()), "face_ring_count": len(group)})
    pd.DataFrame(summary).to_csv(task / "run_metrics.csv", index=False, float_format="%.17g")
    write(task / "validation/spectra.json", {"status": "complete", "rows": len(rows),
          "candidate_sha256": frozen["coefficients_sha256"], "counts_conserved_including_flow": True,
          "fresh_source": read(task / "fresh/completion.json"), "metrics_sha256": sha(task / "spectral_metrics.csv"),
          "coefficient_changes_during_validation": False})


def package_candidate(task, tag):
    """Publish frozen numerical files; preserve separately maintained package notes."""
    assert tag and all(c in 'abcdefghijklmnopqrstuvwxyz0123456789_' for c in tag)
    manifest = read(task / 'input_manifest.json')
    frozen = read(task / 'candidate_frozen.json')
    assert sha(task / 'candidate/coefficients.json') == frozen['coefficients_sha256']
    assert read(task / 'validation/build_real.json')['status'] == 'passed'
    fb = ROOT / 's3_fb_cor/parameters' / tag
    att = ROOT / 's3_energy_attenuation/parameters' / tag
    fb.mkdir(parents=True, exist_ok=True)
    att.mkdir(parents=True, exist_ok=True)

    def copy_frozen(source, target):
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            assert sha(source) == sha(target), f'Existing package file differs: {target}'
        else:
            shutil.copy2(source, target)

    copy_frozen(task / 'candidate/coefficients.json', fb / 'coefficients.json')
    for source in (task / 'candidate/source_calibration').glob('*.txt'):
        copy_frozen(source, fb / 'source_calibration' / source.name)
    copy_frozen(task / 'candidate/pixel_factors.csv', att / 'pixel_factors.csv')
    copy_frozen(task / 'run_metrics.csv', att / 'run_metrics.csv')
    frame = pd.read_csv(Path(manifest['baseline_task']) / 'coefficients.csv', dtype={'run': str})
    frame = frame.rename(columns={'q': 'q_0244_ring_baseline'})
    frame['q_application'] = 'Read matched q_pixel when present, otherwise q_position; never multiply both'
    content = frame.to_csv(index=False, float_format='%.17g')
    destination = fb / 'coefficients.csv'
    if destination.exists():
        assert destination.read_text() == content, 'Existing A/B table differs'
    else:
        destination.write_text(content)
    pointer = {'schema': 'sm150_s3_pixel_candidate_pointer_v1',
               'path': f'../../../s3_fb_cor/parameters/{tag}/coefficients.json',
               'sha256': sha(fb / 'coefficients.json'), 'formal_replacement': False}
    if (att / 'complete_package.json').exists():
        assert read(att / 'complete_package.json') == pointer
    else:
        write(att / 'complete_package.json', pointer)
    for directory in [fb, att]:
        if not (directory / 'README.md').exists():
            (directory / 'README.md').write_text(
                '# S3 pixel candidate\n\nFrozen candidate; not activated. '
                'Use q_pixel instead of q_position, once on both matched faces. '
                f'Inputs, checks and limitations: {task.relative_to(ROOT)}.\n')
    files = list(fb.rglob('*')) + list(att.rglob('*'))
    write(task / 'package_manifest.json', {'schema': 'sm150_pixel_candidate_package_v1',
          'formal_activation': False, 'files': {str(p.relative_to(ROOT)): sha(p) for p in files if p.is_file()}})
    print('Frozen package verified/published:', fb, att)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["export", "evaluate", "package"])
    p.add_argument("--task", type=Path, required=True)
    p.add_argument("--fit-label", default="fixed_support")
    p.add_argument("--package-tag", default="low_pixel_20261102")
    args = p.parse_args()
    if args.mode == "export":
        export(args.task.resolve(), args.fit_label)
    elif args.mode == 'evaluate':
        evaluate(args.task.resolve())
    else:
        package_candidate(args.task.resolve(), args.package_tag)


if __name__ == "__main__":
    main()
