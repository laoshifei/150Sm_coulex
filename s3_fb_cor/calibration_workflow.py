"""Explicit targeted selection, disjoint confirmation and final export."""
from pathlib import Path
from datetime import datetime
import argparse
import json
import sys
import numpy as np

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from s3_fb_cor.review import (
    candidates, maps_from_channels, spectrum_file, load_pairs,
    spectra, pair_metrics, guard, distance_score, model_json,
)
from s3_energy_attenuation import workflow_io as io
from s3_energy_attenuation.workflow_io import digest

def select_run(row, base, prev, reference):
    run = row["run"]
    conditional=json.loads(Path(row['conditional']).read_text());fixed=json.loads(Path(row['fixed_fit']).read_text())
    active=24 if conditional['period']=='low' else 8;ref=conditional['reference_run']
    models=candidates(base['runs'][run],prev['runs'][run],fixed)
    refmap=maps_from_channels(base['runs'][ref]['channels'])
    _,hist=spectrum_file(row['spectra']);v,meta=load_pairs(row['pairs'])
    spectral=spectra(hist,reference[ref],refmap,models,conditional,active,1)
    global_metrics,channels=pair_metrics(v,models,active,1)
    flags={name:guard(name,spectral,global_metrics,channels) for name in models if name not in ['formal','0228']}
    passing=[name for name in flags if name.startswith('fixed_offset_') and not flags[name]]
    preferred=min(passing,key=lambda n:distance_score(spectral,n)) if passing else ('0229' if not flags['0229'] else '0228')
    result={'status':'complete','run':run,'period':conditional['period'],'reference_run':ref,'active_rings':active,'preferred':preferred,
            'flags':flags,'spectral':spectral,'global_metrics':global_metrics,'channels':channels,'models':{n:model_json(m) for n,m in models.items()},
            'spectral_scores':{n:distance_score(spectral,n) for n in models},'source':row,
            'source_sample_metadata':{k:meta[k] for k in ['sample_entries','input_entries','sample_fraction']},
            'scope':'Reused fold1 model selection; disjoint confirmation samples not read.'}
    return result


def validate_confirmation(result, selection_path, selected):
    """Prevent exports from mixing models or selecting again after confirmation."""
    if (result["status"] != "complete" or result["run"] != selected["run"]
            or result["selection_result_sha256"] != digest(selection_path)
            or result["preferred"] != selected["preferred"]
            or result["models"] != selected["models"]):
        raise ValueError(f"Mismatched confirmation: {selected['run']}")
    preferred = selected["preferred"]
    flags = guard(preferred, result["spectral"], result["global_metrics"],
                  result["channels"], confirmation=True)
    final = preferred if preferred == "0228" or not flags else "0228"
    decision = ("baseline_retained_at_selection" if preferred == "0228" else
                ("candidate_confirmed" if not flags else "confirmation_failed_retained_0228"))
    if result["confirmation_flags"] != flags or result["final"] != final or result["decision"] != decision:
        raise ValueError(f"Confirmation rule mismatch: {selected['run']}")
    for model in result["models"].values():
        io.validate_model(model)


def validate_selection(selected, base, prev):
    run = selected["run"]
    baseline = base["runs"][run]
    previous = prev["runs"][run]
    if (selected["reference_run"] != baseline["reference_run"]
            or selected["period"] != baseline["period"]):
        raise ValueError(f"Selection reference/period mismatch: {run}")
    expected = candidates(baseline, previous, {"models": {}})
    for name, model in expected.items():
        if model_json(model) != selected["models"][name]:
            raise ValueError(f"Selection baseline changed: {run} {name}")
    for model in selected["models"].values():
        io.validate_model(model)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["select", "confirm", "export"])
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="New directory; existing paths are never overwritten")
    parser.add_argument("--runs", nargs="+", help="Optional subset for select/confirm")
    parser.add_argument("--jobs", type=int, choices=range(1, 9), default=1)
    parser.add_argument("--selection-freeze", type=Path, help="Override frozen selection receipt")
    parser.add_argument("--confirmation-dir", type=Path, help="Override saved confirmation results for export")
    args = parser.parse_args(argv)
    if args.command == "export" and args.runs:
        parser.error("Export requires every targeted run in the manifest")
    cfg = io.load_config(args.config)
    for key in ["selection_freeze", "confirmation_dir"]:
        value = getattr(args, key)
        if value:
            cfg[key] = value.resolve()
    base = io.read_json(cfg["baseline"])
    prev = io.read_json(cfg["previous"])
    io.validate_baseline(base)
    manifest = io.load_manifest(cfg["manifest"])
    training_refs = {r["run"]: r for r in io.load_manifest(cfg["reference_manifest"])}
    output = io.new_output(args.output_dir)
    io.write_receipt(output, cfg, args.command)

    if args.command == "select":
        rows = io.choose_runs(manifest, args.runs)
        required = {base["runs"][r["run"]]["reference_run"] for r in rows}
        reference = {}
        for run in sorted(required):
            io.verify_row(training_refs[run], ["spectra"])
            reference[run] = spectrum_file(training_refs[run]["spectra"])[1]

        def worker(row):
            io.verify_row(row, ["pairs", "spectra", "conditional", "fixed_fit"])
            io.validate_conditional(row, base)
            result = select_run(row, base, prev, reference)
            validate_selection(result, base, prev)
            return result

        state = io.run_batch(rows, worker, output, args.jobs)
        freeze = {
            "frozen": state["finished"], "choices": sorted(state["runs"], key=lambda r: r["run"]),
            "files": {str(p.relative_to(output)): digest(p) for p in sorted((output / "results").glob("*.json"))},
            "protocol_sha256": digest(cfg["protocol"]), "code_sha256": io.code_hashes(),
            "inputs_sha256": {str(cfg[k]): digest(cfg[k]) for k in
                              ["baseline", "previous", "manifest", "reference_manifest", "protocol"]},
        }
        io.write_json(output / "selection_frozen.json", freeze)
    else:
        selections = io.load_frozen(cfg["selection_freeze"], cfg["protocol"])
        for _, selected in selections.values():
            validate_selection(selected, base, prev)
        if args.command == "confirm":
            rows = io.choose_runs([{"run": r} for r in selections], args.runs)
            sources = {r["run"]: r for r in io.load_manifest(cfg["confirmation_manifest"])}
            required = {selections[r["run"]][1]["reference_run"] for r in rows}
            reference = {}
            for run in sorted(required):
                io.verify_row(sources[run], ["pairs", "spectra"])
                io.verify_disjoint(training_refs[run]["pairs"], sources[run]["pairs"])
                reference[run] = spectrum_file(sources[run]["spectra"])[1]

            def worker(row):
                run = row["run"]
                path, selected = selections[run]
                source = sources[run]
                io.verify_row(source, ["pairs", "spectra"])
                io.verify_row(selected["source"], ["conditional"])
                io.verify_disjoint(selected["source"]["pairs"], source["pairs"])
                result = confirm_run(source, selected, path, base, reference)
                validate_confirmation(result, path, selected)
                return result

            io.run_batch(rows, worker, output, args.jobs)
        else:
            from s3_fb_cor.export_parameters import export_parameters
            expected = {r["run"] for r in manifest}
            files = {p.stem: p for p in cfg["confirmation_dir"].glob("*.json")}
            if set(files) != expected or set(selections) != expected:
                raise ValueError("Confirmation/frozen run sets must exactly match the targeted manifest")
            confirmed = {}
            for run, path in files.items():
                row = io.read_json(path)
                validate_confirmation(row, *selections[run])
                confirmed[run] = row
            for run, row in prev["runs"].items():
                if run not in base["runs"]:
                    raise ValueError(f"Missing baseline run: {run}")
                io.validate_model({"a": row["A_keV"], "b": row["B"], "radial": row["q_position"]})
            result = export_parameters(base, prev, io.read_json(cfg["review_rows"]), confirmed, output)
            io.write_receipt(output, cfg, "export", {
                "confirmation_sha256": {r: digest(p) for r, p in files.items()},
                "event_files_read": 0,
            })
            print(json.dumps(result), flush=True)
    io.write_manifest(output)
    return 0


def confirm_run(row, selected, selection_path, base, reference):
    run=row['run'];ref=selected['reference_run']
    assert digest(row['pairs'])==row['source_sha256'] and digest(row['spectra'])==row['spectra_sha256']
    models={n:{k:np.asarray(v) for k,v in m.items()} for n,m in selected['models'].items()}
    conditional=json.loads(Path(selected['source']['conditional']).read_text())
    _,hist=spectrum_file(row['spectra']);v,meta=load_pairs(row['pairs'])
    refmap=maps_from_channels(base['runs'][ref]['channels']);active=selected['active_rings']
    spectral=spectra(hist,reference[ref],refmap,models,conditional,active,-1)
    global_metrics,channels=pair_metrics(v,models,active,-1)
    preferred=selected['preferred'];flags=guard(preferred,spectral,global_metrics,channels,confirmation=True)
    final=preferred if preferred=='0228' or not flags else '0228'
    status='baseline_retained_at_selection' if preferred=='0228' else ('candidate_confirmed' if not flags else 'confirmation_failed_retained_0228')
    result={'status':'complete','run':run,'period':selected['period'],'reference_run':ref,'active_rings':active,
            'preferred':preferred,'final':final,'decision':status,'confirmation_flags':flags,
            'spectral':spectral,'global_metrics':global_metrics,'channels':channels,'models':selected['models'],
            'spectral_scores':{n:distance_score(spectral,n) for n in models},'source':row,
            'source_sample_metadata':{k:meta[k] for k in ['sample_entries','input_entries','sample_fraction']},
            'selection_result_sha256':digest(selection_path),
            'scope':'All saved confirmation pairs; entry intervals checked against the supplied selection and reference samples. Same experiment; not a new independent experiment or a claim that no earlier project study inspected any such data.'}
    return result


if __name__ == "__main__":
    raise SystemExit(main())
