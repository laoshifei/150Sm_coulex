"""Freeze and validate empirical-flank S3 candidates without activating them."""
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
from s3_energy_attenuation.anchor_review import TASK, SOURCE, PACKAGE_PATH, sha, write


def records(label, run):
    data = json.loads((TASK/'results'/label/f'fit_{run}.json').read_text())
    return {(r['ring'], r['sector'], r['face']): r['fits'] for r in data['records']}


def decide():
    base = json.loads(PACKAGE_PATH.read_text())
    old = json.loads((ROOT/'s3_fb_cor/parameters/source_recalibrated_20261102/coefficients.json').read_text())
    reference = records('combined', '0616')
    rows = []
    for run in ['0641', '0659', '0674']:
        current = records('combined', run)
        for ring in range(20, 25):
            for sector in range(1, 33):
                q_current = base['runs'][run]['q_pixel'][ring-1][sector-1]
                q_reference = base['runs']['0616']['q_pixel'][ring-1][sector-1]
                row = {'run': run, 'ring': ring, 'sector': sector,
                       'q_0244': old['runs'][run]['q_position'][ring-1],
                       'q_0245': q_current, 'q_trial': q_current}
                flags, ratios, errors, window_spreads, split_spreads = [], [], [], [], []
                for face in ['ring', 'sector']:
                    ref, cur = reference[ring, sector, face], current[ring, sector, face]
                    for name in ['main', 'narrow', 'wide', 'even', 'odd']:
                        for origin, fit in [('reference', ref[name]), ('current', cur[name])]:
                            if not fit['ok']:
                                flags.extend(f'{face}:{origin}:{name}:{flag}' for flag in fit['flags'])
                    if 'mu' not in ref['main'] or 'mu' not in cur['main']:
                        continue
                    rmu, cmu = ref['main']['mu'], cur['main']['mu']
                    ratio = q_reference*rmu/cmu
                    error = np.hypot(ref['main']['mu_error']/rmu, cur['main']['mu_error']/cmu)
                    row.update({f'reference_mu_{face}_MeV': q_reference*rmu,
                                f'uncorrected_mu_{face}_MeV': cmu,
                                f'q_{face}': ratio,
                                f'conditional_log_error_{face}': error,
                                f'delta_0245_{face}_MeV': q_current*cmu-q_reference*rmu,
                                f'sigma_reference_{face}_MeV': q_reference*ref['main']['sigma'],
                                f'sigma_uncorrected_{face}_MeV': cur['main']['sigma']})
                    ratios.append(ratio)
                    errors.append(error)
                    variants = {}
                    for name in ['narrow', 'wide', 'even', 'odd']:
                        if 'mu' in ref[name] and 'mu' in cur[name]:
                            variants[name] = q_reference*ref[name]['mu']/cur[name]['mu']
                    if len(variants) == 4:
                        window = max(abs(np.log(variants[name]/ratio)) for name in ['narrow', 'wide'])
                        split = abs(np.log(variants['even']/variants['odd']))
                        row[f'window_log_difference_{face}'] = window
                        row[f'split_log_difference_{face}'] = split
                        window_spreads.append(window)
                        split_spreads.append(split)
                        if window > .03:
                            flags.append(face+':window_ratio_unstable')
                        if split > .05:
                            flags.append(face+':split_ratio_unstable')
                if len(ratios) == 2:
                    common = float(np.sqrt(np.prod(ratios)))
                    face_spread = float(abs(np.log(ratios[0]/ratios[1])))
                    # Do not count two correlated readouts as independent events.
                    error = float(max(errors))
                    row.update(q_anchor=common, conditional_log_error=error,
                               face_log_difference=face_spread)
                    if face_spread > .05:
                        flags.append('two_face_ratios_disagree')
                else:
                    flags.append('missing_main_face_location')
                if flags:
                    row['decision'] = 'retain_unresolved'
                else:
                    shift = np.log(row['q_anchor']/q_current)
                    if abs(shift) <= max(.02, 2*row['conditional_log_error']):
                        row['decision'] = 'retain_consistent'
                    else:
                        row['decision'] = 'retune_up' if shift > 0 else 'retune_down'
                        row['q_trial'] = row['q_anchor']
                for face in ['ring', 'sector']:
                    if f'uncorrected_mu_{face}_MeV' in row:
                        row[f'delta_trial_{face}_MeV'] = row['q_trial']*row[f'uncorrected_mu_{face}_MeV']-row[f'reference_mu_{face}_MeV']
                        row[f'width_ratio_trial_{face}'] = row['q_trial']*row[f'sigma_uncorrected_{face}_MeV']/row[f'sigma_reference_{face}_MeV']
                row['flags'] = ';'.join(sorted(set(flags)))
                rows.append(row)
    frame = pd.DataFrame(rows)
    assert len(frame) == 480 and not frame[['run', 'ring', 'sector']].duplicated().any()
    destination = TASK/'candidate'
    destination.mkdir(exist_ok=False)
    frame.to_csv(destination/'pixel_decisions.csv', index=False, float_format='%.17g')
    package = copy.deepcopy(base)
    for row in frame.itertuples():
        package['runs'][row.run]['q_pixel'][row.ring-1][row.sector-1] = row.q_trial
    package.update(created=datetime.datetime.now().astimezone().isoformat(), formal_replacement=False,
                   scope='Three low-beam runs, rings20-24, empirical falling-flank relative-scale candidates',
                   previous_candidate_sha256=sha(PACKAGE_PATH),
                   feature_note='Location is a local empirical falling flank, not an identified elastic endpoint. Width/tails are not forced to match.')
    write(destination/'coefficients.json', package)
    shutil.copytree(SOURCE/'candidate/source_calibration', destination/'source_calibration')
    frozen = {'created': package['created'], 'package_sha256': sha(destination/'coefficients.json'),
              'decision_sha256': sha(destination/'pixel_decisions.csv'),
              'data_label': 'combined historical0244 and fresh0245; previously inspected',
              'new20M_spectra_inspected': False,
              'decision_counts': frame.decision.value_counts().to_dict(),
              'script_sha256': {name: sha(ROOT/'s3_energy_attenuation'/name)
                                for name in ['anchor_review.py', 'anchor_results.py']}}
    write(TASK/'candidate_frozen.json', frozen)
    print(json.dumps(frozen['decision_counts']), flush=True)


def validate():
    frozen = json.loads((TASK/'candidate_frozen.json').read_text())
    assert sha(TASK/'candidate/coefficients.json') == frozen['package_sha256']
    decision = pd.read_csv(TASK/'candidate/pixel_decisions.csv', dtype={'run': str})
    reference = records('validation', '0616')
    lookup = {run: records('validation', run) for run in ['0641', '0659', '0674']}
    rows = []
    for item in decision.to_dict('records'):
        run, ring, sector = item['run'], item['ring'], item['sector']
        row = dict(item)
        flags, old_deltas, new_deltas = [], [], []
        for face in ['ring', 'sector']:
            ref = reference[ring, sector, face]['main']
            cur = lookup[run][ring, sector, face]['main']
            if not ref['ok'] or not cur['ok']:
                flags.extend(f'{face}:reference:{x}' for x in ref['flags'])
                flags.extend(f'{face}:current:{x}' for x in cur['flags'])
            if 'mu' not in ref or 'mu' not in cur:
                flags.append(face+':no_validation_location')
                continue
            old_delta = item['q_0245']*cur['mu']-ref['mu']
            new_delta = item['q_trial']*cur['mu']-ref['mu']
            uncertainty = np.hypot(item['q_trial']*cur['mu_error'], ref['mu_error'])
            row.update({f'validation_delta_0245_{face}_MeV': old_delta,
                        f'validation_delta_trial_{face}_MeV': new_delta,
                        f'validation_delta_error_{face}_MeV': uncertainty,
                        f'validation_width_ratio_{face}': item['q_trial']*cur['sigma']/ref['sigma']})
            old_deltas.append(abs(old_delta))
            new_deltas.append(abs(new_delta))
            # This is a descriptive compatibility band, not an absolute energy certificate.
            if abs(new_delta) > max(1.5, 3*uncertainty):
                flags.append(face+':remaining_landmark_offset')
            if abs(new_delta) > abs(old_delta)+.5:
                flags.append(face+':landmark_worsened')
        proposed = item['decision'] in ['retune_up', 'retune_down']
        if proposed and len(old_deltas) == 2 and np.mean(new_deltas) >= np.mean(old_deltas):
            flags.append('no_average_landmark_improvement')
        row['validation_flags'] = ';'.join(sorted(set(flags)))
        row['validation_status'] = ('confirmed_change' if not flags else 'change_not_confirmed') if proposed else 'no_change_proposed'
        row['q_final_candidate'] = item['q_trial'] if row['validation_status'] == 'confirmed_change' else item['q_0245']
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(TASK/'validated_decisions.csv', index=False, float_format='%.17g')
    summary = {'status': 'complete', 'frozen_package_sha256': frozen['package_sha256'],
               'proposed_changes': int(decision.decision.isin(['retune_up', 'retune_down']).sum()),
               'validation_counts': out.validation_status.value_counts().to_dict(),
               'counts_by_run': out.groupby(['run', 'validation_status']).size().reset_index(name='count').to_dict('records'),
               'criteria': 'Both empirical flank fits supported; per-face remaining offset <=max(1.5MeV,3conditional_sigma); neither face worsens by>0.5MeV; mean absolute offset improves.',
               'formal_activation': False, 'production_rebuild': False}
    good = out[out.validation_status == 'confirmed_change']
    all_proposed = out[out.decision.isin(['retune_up', 'retune_down'])]
    summary['all_proposed_before_validation_selection'] = {}
    for face in ['ring', 'sector']:
        summary[face] = {'median_abs_old_MeV': float(good[f'validation_delta_0245_{face}_MeV'].abs().median()) if len(good) else None,
                         'median_abs_new_MeV': float(good[f'validation_delta_trial_{face}_MeV'].abs().median()) if len(good) else None}
        summary['all_proposed_before_validation_selection'][face] = {
            'median_abs_old_MeV': float(all_proposed[f'validation_delta_0245_{face}_MeV'].abs().median()),
            'median_abs_new_MeV': float(all_proposed[f'validation_delta_trial_{face}_MeV'].abs().median()),
            'available': int(all_proposed[f'validation_delta_trial_{face}_MeV'].notna().sum())}
    write(TASK/'validation/independent.json', summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


def package():
    frame = pd.read_csv(TASK/'validated_decisions.csv', dtype={'run': str})
    original = json.loads(PACKAGE_PATH.read_text())
    model = copy.deepcopy(original)
    for row in frame.itertuples():
        model['runs'][row.run]['q_pixel'][row.ring-1][row.sector-1] = row.q_final_candidate
    selected = frame[frame.validation_status == 'confirmed_change']
    model.update(created=datetime.datetime.now().astimezone().isoformat(), formal_replacement=False,
                 scope='Empirical relative falling-flank updates at 166 validated pixels in0641/0659/0674, rings20-24; all other values inherit0245',
                 previous_candidate_sha256=sha(PACKAGE_PATH),
                 feature_note='A local spectral landmark, not an identified physical endpoint; width/tails not forced to match. Held-out data were used to accept/reject proposals, with no refit.')
    for run in ['0641', '0659', '0674']:
        model['runs'][run]['pixel_source'] = str((TASK/'validated_decisions.csv').relative_to(ROOT))
        model['runs'][run]['pixel_source_sha256'] = sha(TASK/'validated_decisions.csv')
    destination = TASK/'final_candidate'
    destination.mkdir(exist_ok=False)
    write(destination/'coefficients.json', model)
    shutil.copytree(SOURCE/'candidate/source_calibration', destination/'source_calibration')
    fb = ROOT/'s3_fb_cor/parameters/feature_anchor_20261102'
    att = ROOT/'s3_energy_attenuation/parameters/feature_anchor_20261102'
    fb.mkdir(exist_ok=False)
    att.mkdir(exist_ok=False)
    shutil.copy2(destination/'coefficients.json', fb/'coefficients.json')
    shutil.copytree(destination/'source_calibration', fb/'source_calibration')
    shutil.copy2(TASK/'validated_decisions.csv', att/'pixel_decisions.csv')
    selected.to_csv(att/'confirmed_changes.csv', index=False, float_format='%.17g')
    write(att/'complete_package.json', {'schema': 'sm150_s3_pixel_candidate_pointer_v1',
          'path': '../../../s3_fb_cor/parameters/feature_anchor_20261102/coefficients.json',
          'sha256': sha(fb/'coefficients.json'), 'formal_replacement': False})
    description = ('# S3 falling-flank candidate\n\n'
        'This package is not activated. It updates only confirmed positions in low-beam '
        'runs0641,0659,0674, rings20-24; all other numbers inherit0245. '
        'Use the matched q_pixel once on both faces, replacing q_position. '
        'The empirical falling flank is not an identified elastic endpoint or an absolute-energy calibration. '
        'Widths and low-energy tails are not forced to match.\n\n'
        'The complete package retains173runs and requires the supplied source calibration. '
        'See codex/analysis/s3_energy_feature_anchor_20261102 for protocols, '
        'all480decisions, new20M-C2 validation, uncertainty limitations and plots.\n')
    (fb/'README.md').write_text(description)
    (att/'README.md').write_text(description)
    files = list(fb.rglob('*'))+list(att.rglob('*'))
    write(TASK/'package_manifest.json', {'files': {str(p.relative_to(ROOT)): sha(p) for p in files if p.is_file()},
          'confirmed_changes': len(selected), 'up': int((selected.q_final_candidate > selected.q_0245).sum()),
          'down': int((selected.q_final_candidate < selected.q_0245).sum()), 'formal_activation': False})
    print('Published candidate; changes', len(selected), 'up', (selected.q_final_candidate > selected.q_0245).sum(),
          'down', (selected.q_final_candidate < selected.q_0245).sum(), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['freeze', 'validate', 'package'])
    args = parser.parse_args()
    {'freeze': decide, 'validate': validate, 'package': package}[args.mode]()
