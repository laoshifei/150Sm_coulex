"""Position-factor assembly and frozen candidate comparison.

The numerical evaluator preserves the 0229 selection and confirmation rules.
Paths and execution are supplied explicitly through the command-line frontend.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from s3_energy_attenuation.analysis import (
    model_maps, maps_from_channels, spectrum_file, spectral_metrics,
    load_pairs, differences, pair_masks, metrics,
)
from s3_energy_attenuation import workflow_io as io

def evaluate(row,baseline,reference_paths):
 run=row['run'];conditional=json.loads(Path(row['conditional']).read_text());coupled=json.loads(Path(row['coupled']).read_text());base=baseline['runs'][run]
 active=24 if conditional['period']=='low' else 8;ref=conditional['reference_run']
 maps=model_maps(base,conditional,coupled);refmaps=maps_from_channels(baseline['runs'][ref]['channels'])
 _,hist=spectrum_file(row['spectra']);_,refhist=spectrum_file(reference_paths[ref])
 spectral=[]
 for fold in [1,3]:spectral.extend(spectral_metrics(hist,refhist,refmaps,maps,conditional,active,fold))
 v,meta=load_pairs(row['pairs']);pair_summary={};channel_rows=[]
 for name,model in maps.items():
  raw,d,_,_=differences(v,model,maps['0228']);pair_summary[name]={}
  for fold in [1,3]:
   common,prompt,random,production=pair_masks(v,active,fold);selected=common&prompt
   pair_summary[name][str(fold)]={'normalized':metrics(d[selected&production]),'raw':metrics(raw[selected&production]),'random':int(np.sum(common&random&production))}
   for c in range(56):
    channel=(v['r']==c) if c<24 else ((v['s']==c)&production)
    channel_rows.append({'index':c,'model':name,'fold':fold,**metrics(d[selected&channel])})
 def spectral_score(name,fold):
  a=[x['D'] for x in spectral if x['model']==name and x['fold']==fold and x['production']]
  return float(np.sqrt(np.mean(np.square(a)))) if a else None
 def guards(name,fold,confirmation=False):
  flags=[];before=pair_summary['0228'][str(fold)]['normalized'];after=pair_summary[name][str(fold)]['normalized']
  factor=1.05 if confirmation else 1.03
  if after['n']<100:flags.append('too_few_front_back_pairs');return flags
  for key in ['width68_keV','score_keV']:
   if after[key]>factor*before[key]+1:flags.append('front_back_'+key+'_worse')
  if after['tail2MeV']>before['tail2MeV']+.01:flags.append('front_back_tail_increase')
  old={x['index']:x for x in channel_rows if x['model']=='0228' and x['fold']==fold}
  for item in channel_rows:
   if item['model']!=name or item['fold']!=fold or not(item['index']<active or item['index']>=24) or item['n']<100:continue
   before=old[item['index']]
   if item['score_keV']>1.10*before['score_keV']+20:flags.append('channel_front_back_worse_'+str(item['index']))
  d0=spectral_score('0228',fold);d1=spectral_score(name,fold)
  if d0 is None or d1 is None:flags.append('too_few_spectral_pixels');return flags
  if confirmation:
   if d1>d0+.002:flags.append('confirmation_spectrum_worse')
  elif d1>.95*d0 or d0-d1<.001:flags.append('spectral_improvement_below_threshold')
  old={x['channel_index']:x for x in spectral if x['model']=='0228' and x['fold']==fold and x['production']}
  bad=[]
  for x in spectral:
   if x['model']==name and x['fold']==fold and x['production']:
    if x['invalid_pixels'] or x['D']>old[x['channel_index']]['D']+.02:bad.append(x['channel_index'])
  if len(bad)>max(1,int(.2*len(old))):flags.append('spectral_channels_worse_'+str(bad))
  return flags
 selection={name:guards(name,1) for name in maps if name not in ['formal','0228','direct']}
 passing=[name for name,flags in selection.items() if not flags]
 preferred=min(passing,key=lambda name:spectral_score(name,1)) if passing else '0228'
 if run==ref:preferred='0228'
 confirm=guards(preferred,3,True) if preferred!='0228' else []
 chosen=preferred if not confirm else '0228'
 selected=maps[chosen]
 return {'status':'complete','run':run,'period':conditional['period'],'reference_run':ref,'active_rings':active,'source':row,
  'source_sample_metadata':{k:meta[k] for k in ['sample_entries','sample_fraction','input_entries','sampled_columns_sha256']},
  'spectral_metrics':spectral,'pair_summary':pair_summary,'channel_pair_metrics':channel_rows,
  'spectral_scores':{name:{str(fold):spectral_score(name,fold) for fold in [1,3]} for name in maps},
  'selection_flags':selection,'preferred':preferred,'confirmation_flags':confirm,'selected':chosen,
  'selected_a':selected['a'].tolist(),'selected_b':selected['b'].tolist(),'selected_radial':selected['radial'].tolist(),
  'conditional_sha256':hashlib.sha256(Path(row['conditional']).read_bytes()).hexdigest(),'coupled_sha256':hashlib.sha256(Path(row['coupled']).read_bytes()).hexdigest(),
  'scope':'Internal held-out comparison of frozen candidates; no production adoption. Radial variant applies the same position factor to both faces.'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["evaluate", "export"])
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="New directory; existing paths are never overwritten")
    parser.add_argument("--runs", nargs="+", help="Optional run subset for evaluate")
    parser.add_argument("--jobs", type=int, choices=range(1, 9), default=1)
    parser.add_argument("--evaluation-dir", type=Path,
                        help="Saved run JSON directory for export (overrides config)")
    args = parser.parse_args(argv)
    if args.command == "export" and args.runs:
        parser.error("Export uses the complete manifest; supply a separate manifest for a subset")
    cfg = io.load_config(args.config)
    if args.evaluation_dir:
        cfg["evaluation_dir"] = args.evaluation_dir.resolve()
    baseline = io.read_json(cfg["baseline"])
    io.validate_baseline(baseline)
    rows = io.choose_runs(io.load_manifest(cfg["manifest"]), args.runs)
    output = io.new_output(args.output_dir)
    io.write_receipt(output, cfg, args.command)

    if args.command == "evaluate":
        reference_rows = {r["run"]: r for r in io.load_manifest(cfg["reference_manifest"])}
        required = {baseline["runs"][r["run"]]["reference_run"] for r in rows}
        references = {}
        for run in sorted(required):
            io.verify_row(reference_rows[run], ["spectra"])
            references[run] = reference_rows[run]["spectra"]

        def worker(row):
            io.verify_row(row, ["pairs", "spectra", "conditional", "coupled"])
            conditional = io.validate_conditional(row, baseline)
            for model in model_maps(baseline["runs"][row["run"]], conditional,
                                    io.read_json(row["coupled"])).values():
                io.validate_model(model)
            return evaluate(row, baseline, references)

        io.run_batch(rows, worker, output, args.jobs)
    else:
        from s3_energy_attenuation.export_parameters import export_parameters
        for row in rows:
            e = io.read_json(cfg["evaluation_dir"] / f"{row['run']}.json")
            if e["status"] != "complete" or e["source"] != row:
                raise ValueError(f"Stale or incomplete evaluation: {row['run']}")
            io.verify_row(row, ["conditional", "coupled"])
            conditional = io.read_json(row["conditional"])
            if [c["index"] for c in conditional["channels"]] != list(range(56)):
                raise ValueError("Conditional channel order differs from baseline")
            model = model_maps(baseline["runs"][row["run"]], conditional,
                               io.read_json(row["coupled"]))[e["selected"]]
            io.validate_model(model)
            for key, saved in [("a", "selected_a"), ("b", "selected_b"), ("radial", "selected_radial")]:
                if not np.array_equal(model[key], e[saved]):
                    raise ValueError(f"Saved selection coefficients changed: {row['run']}")
        result = export_parameters(rows, baseline, cfg["evaluation_dir"], output)
        io.write_receipt(output, cfg, "export", {
            "saved_decisions_sha256": {r["run"]: io.digest(cfg["evaluation_dir"] / f"{r['run']}.json") for r in rows},
            "event_files_read": 0,
        })
        print(json.dumps(result), flush=True)
    io.write_manifest(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
