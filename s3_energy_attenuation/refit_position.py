"""Validate refitted ring factors on both faces and export a coupled candidate.

The fit is performed by conditioned_calibrate. This module applies frozen
selection/confirmation rules; it never tunes a model on confirmation segments.
"""
from pathlib import Path
import argparse
import concurrent.futures as futures
import json
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from s3_fb_cor.source_refit import read, write, sha, model_map
from s3_fb_cor.review import cumulative_net
from s3_energy_attenuation.analysis import spectrum_file


def mapped_cdf(values, face, fold, pixels, model):
    indices = pixels // 32 if face == "ring" else 24 + pixels % 32
    q = model["radial"][pixels // 32]
    return cumulative_net(values, face, fold, pixels,
                          q * model["a"][indices], q * model["b"][indices])


def compare_ring(values, reference, refmap, baseline, candidate, ring, eligible, fold):
    pixels = np.asarray(eligible, int)
    if not len(pixels):
        return {"fold": fold, "ring": ring + 1, "pixels": [], "status": "insufficient_common_pixels", "baseline_D": None, "candidate_D": None}
    valid = np.ones(len(pixels), bool)
    reference_cdfs, baseline_cdfs = {}, {}
    for face in ["ring", "sector"]:
        rr = mapped_cdf(reference, face, fold, pixels, refmap)
        bb = mapped_cdf(values, face, fold, pixels, baseline)
        reference_cdfs[face], baseline_cdfs[face] = rr, bb
        valid &= rr["reliable"] & bb["reliable"]
    pixels = pixels[valid]
    result = {"fold": fold, "ring": ring + 1, "pixels": pixels.tolist(),
              "baseline_D": None, "candidate_D": None, "faces": {}}
    if len(pixels) < 3:
        return {**result, "status": "insufficient_common_pixels"}
    valid_candidate = True
    for face in ["ring", "sector"]:
        rr, bb = reference_cdfs[face], baseline_cdfs[face]
        cc = mapped_cdf(values, face, fold, pixels, candidate)
        valid_candidate &= bool(np.all(cc["reliable"]))
        db = float(np.sqrt(np.mean((bb["F"][valid, 1:-1] - rr["F"][valid, 1:-1]) ** 2)))
        dc = float(np.sqrt(np.mean((cc["F"][:, 1:-1] - rr["F"][valid, 1:-1]) ** 2)))
        noise = float(np.sqrt(np.mean(cc["variance_F"][:, 1:-1] + rr["variance_F"][valid, 1:-1])))
        result["faces"][face] = {"baseline_D": db, "candidate_D": dc,
            "candidate_noise_floor": noise, "baseline_net": bb["net"][valid].tolist(),
            "candidate_net": cc["net"].tolist(), "reference_net": rr["net"][valid].tolist()}
    result["baseline_D"] = float(np.sqrt(np.mean([c["baseline_D"] ** 2 for c in result["faces"].values()])))
    result["candidate_D"] = float(np.sqrt(np.mean([c["candidate_D"] ** 2 for c in result["faces"].values()])))
    result["status"] = "comparable" if valid_candidate else "candidate_insufficient_support"
    return result


def q_guard(result, confirmation=False):
    if result["status"] != "comparable":
        return [result["status"]]
    old, new = result["baseline_D"], result["candidate_D"]
    reasons = []
    if confirmation:
        if new > old + .002:
            reasons.append("confirmation_distance_worse")
    elif new > .95 * old or old - new < .001:
        reasons.append("selection_improvement_below_threshold")
    for face, row in result["faces"].items():
        if row["candidate_D"] > row["baseline_D"] + .002:
            reasons.append(face + "_distance_worse")
    return reasons


def evaluate_one(task, row):
    run = row["run"]
    ab = read(task / "ab_evaluation" / f"evaluation_{run}.json")
    qfit = read(task / "q_fits" / f"fit_{run}.json")
    ref = read(task / "ab_evaluation" / f"evaluation_{ab['reference_run']}.json")
    base = {k: np.array(v) for k, v in ab["maps"]["equivalent"].items()}
    chosen = {k: np.array(v) for k, v in ab["maps"][ab["chosen"]].items()}
    refmap = {k: np.array(v) for k, v in ref["maps"]["equivalent"].items()}
    _, values = spectrum_file(task / "histograms" / f"spectra_{run}.root")
    _, reference = spectrum_file(task / "histograms" / f"spectra_{ab['reference_run']}.root")
    rows = []
    channels = {r["index"]: r for r in qfit["channels"]}
    for i in range(24):
        c = channels[i]
        candidate = {k: v.copy() for k, v in chosen.items()}
        candidate["radial"][i] = c["q"]
        selection = compare_ring(values, reference, refmap, chosen, candidate, i, c["eligible_pixels"], 1)
        confirmation = compare_ring(values, reference, refmap, chosen, candidate, i, c["eligible_pixels"], 3)
        flags = list(c["training_flags"])
        if run == ab["reference_run"]:
            flags = ["fixed_reference"]
        elif i >= ab["active_rings"]:
            flags = ["excluded_high_beam_outer_ring"]
        else:
            flags += q_guard(selection)
        confirm_flags = q_guard(confirmation, True) if not flags else []
        accepted = not flags and not confirm_flags
        if accepted:
            chosen["radial"][i] = c["q"]
        rows.append({"run": run, "period": ab["period"], "ring": i+1, "production": i < ab["active_rings"],
            "q_old": float(base["radial"][i]), "q_fitted": c["q"], "q_selected": float(chosen["radial"][i]),
            "q_status": "refit_accepted" if accepted else ("excluded_high_beam_outer_ring" if i >= ab["active_rings"] else "equivalent_inherited"),
            "training_and_selection_flags": flags, "confirmation_flags": confirm_flags,
            "selection": selection, "confirmation": confirmation,
            "fit_fold0": c["q_fit_fold0"], "fit_fold2": c["q_fit_fold2"], "alternate_window": c["q_alternate_window"],
            "profile_support": c["q_profile_support"]})
    # A/B changes must also pass cross-run spectral checks on both faces.
    final_checks = []
    global_flags = []
    for fold in [1, 3]:
        for i in range(ab["active_rings"]):
            check = compare_ring(values, reference, refmap, base, chosen, i, channels[i]["eligible_pixels"], fold)
            final_checks.append(check)
            if ab["chosen"] != "equivalent":
                if check["status"] != "comparable":
                    global_flags.append(f"fold{fold}_ring{i+1}_unsupported")
                elif check["candidate_D"] > check["baseline_D"] + .02:
                    global_flags.append(f"fold{fold}_ring{i+1}_spectra_worse")
        complete = [c for c in final_checks if c["fold"] == fold and c["status"] == "comparable"]
        if complete and ab["chosen"] != "equivalent":
            d0 = np.sqrt(np.mean([c["baseline_D"] ** 2 for c in complete]))
            d1 = np.sqrt(np.mean([c["candidate_D"] ** 2 for c in complete]))
            if d1 > d0 + .002:
                global_flags.append(f"fold{fold}_global_spectra_worse")
    if global_flags:
        chosen = {k: v.copy() for k, v in base.items()}
        for c in rows:
            if c["production"]:
                c["q_status"] = "equivalent_inherited_after_ab_spectral_guard"
                c["q_selected"] = c["q_old"]
    result = {"run": run, "period": ab["period"], "reference_run": ab["reference_run"],
        "active_rings": ab["active_rings"], "ab_proposed": ab["chosen"],
        "ab_selected": "equivalent" if global_flags else ab["chosen"],
        "ab_spectral_flags": global_flags, "rings": rows, "final_proposal_spectral_checks": final_checks,
        "map": {k: v.tolist() for k, v in chosen.items()},
        "baseline_map": {k: v.tolist() for k, v in base.items()},
        "q_fit_sha256": sha(task / "q_fits" / f"fit_{run}.json"),
        "ab_evaluation_sha256": sha(task / "ab_evaluation" / f"evaluation_{run}.json")}
    write(task / "final_evaluation" / f"evaluation_{run}.json", result)
    return run


def export(task, source_rows):
    package = {"schema": "sm150_s3_targeted_review_candidate_v1", "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "formula": "Er=q[matched_ring]*(Ar+Br*E0r); Es=q[matched_ring]*(As+Bs*E0s)",
        "scope": "173-run source-recalibrated candidate; high-beam rings1-8, low-beam rings1-24",
        "reference": "Current physical scales of 0472/0616/0680; all source and A/B coefficients must be used together",
        "formal_replacement": False, "runs": {}}
    channels, ring_rows, summary_rows = [], [], []
    for source in source_rows:
        run = source["run"]
        e = read(task / "final_evaluation" / f"evaluation_{run}.json")
        meta = read(task / "pairs" / f"pairs_{run}.root.json")["source_recalibration"]
        source_channels = {c["index"]: c for c in meta["channels"]}
        m = e["map"]
        cr = []
        for i in range(56):
            old = source_channels[i]
            production = i < e["active_rings"] or i >= 24
            row = {"run": run, "period": e["period"], "reference_run": e["reference_run"],
                "index": i, "channel": old["channel"], "hardware_key": old["key"], "production": production,
                "A_keV": m["a"][i], "B": m["b"][i], "old_A_keV": old["old_A_keV"], "old_B": old["old_B"],
                "equivalent_A_keV": old["A_keV"], "equivalent_B": old["B"],
                "new_c0_keV": old["source_coefficients"][0], "new_c1_keV_per_ADC": old["source_coefficients"][1],
                "hardware_gain": old["hardware_gain"],
                "status": e["ab_selected"] if production else "excluded_high_beam_outer_ring",
                "q": m["radial"][i] if i < 24 else None,
                "q_application": "same matched-ring q on both faces"}
            cr.append(row)
            channels.append(row)
        package["runs"][run] = {"period": e["period"], "reference_run": e["reference_run"], "active_rings": e["active_rings"],
            "A_keV": m["a"], "B": m["b"], "q_position": m["radial"], "channels": cr,
            "requires_ring_context": any(q != 1 for q in m["radial"]),
            "decision": e["ab_selected"], "source": "new_source_ADC_refit_with_heldout_guards"}
        for c in e["rings"]:
            rr = {k: v for k, v in c.items() if k not in ["selection", "confirmation"]}
            for fold in ["selection", "confirmation"]:
                for field in ["baseline_D", "candidate_D", "status"]:
                    rr[fold + "_" + field] = c[fold][field]
            rr["training_and_selection_flags"] = ";".join(rr["training_and_selection_flags"])
            rr["confirmation_flags"] = ";".join(rr["confirmation_flags"])
            ring_rows.append(rr)
        summary_rows.append({"run": run, "period": e["period"], "reference_run": e["reference_run"],
            "active_rings": e["active_rings"], "ab_proposed": e["ab_proposed"], "ab_selected": e["ab_selected"],
            "ab_spectral_flags": ";".join(e["ab_spectral_flags"]),
            "q_updates": sum(c["q_status"] == "refit_accepted" for c in e["rings"]),
            "sample_entries": source["sample_entries"], "source_entries": source["source_entries"], "sample_fraction": source["sample_fraction"],
            "clean_prompt_pairs": meta["prompt"], "clean_random_pairs": meta["random"],
            "source_input_sha256": source["source_sha256"]})
    write(task / "coefficients.json", package)
    pd.DataFrame(channels).to_csv(task / "coefficients.csv", index=False, float_format="%.17g")
    pd.DataFrame(ring_rows).to_csv(task / "position_factors.csv", index=False, float_format="%.17g")
    pd.DataFrame(summary_rows).to_csv(task / "run_summary.csv", index=False, float_format="%.17g")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--task", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=4)
    args = p.parse_args()
    assert 1 <= args.jobs <= 8
    task = args.task.resolve()
    rows = read(task / "source_manifest.json")
    start = time.time()
    with futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {pool.submit(evaluate_one, task, row): row["run"] for row in rows}
        for k, done in enumerate(futures.as_completed(pending), 1):
            print(f"position-check {k}/{len(rows)} {done.result()} elapsed={time.time()-start:.1f}s", flush=True)
    export(task, rows)
    write(task / "stage_final_evaluation.json", {"runs": len(rows), "elapsed_seconds": time.time()-start,
          "program_sha256": sha(__file__), "jobs": args.jobs})


if __name__ == "__main__":
    main()
