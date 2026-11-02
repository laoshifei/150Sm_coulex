"""Review relative S3 energy scales using local falling flanks, not whole spectra.

The fitted location is an empirical landmark, not an identified physical endpoint.
No runtime parameters or production ROOT files are changed by this program.
"""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import multiprocessing
import sys
import time

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
from scipy.special import ndtr, xlogy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from s3_energy_attenuation.analysis import spectrum_file

SOURCE = ROOT / 'codex/analysis/s3_low_pixel_energy_20261102'
TASK = ROOT / 'codex/analysis/s3_energy_feature_anchor_20261102'
PACKAGE_PATH = SOURCE / 'candidate/coefficients.json'
WINDOWS = {'main': (-6., 10.), 'narrow': (-4., 8.), 'wide': (-8., 12.)}
ALPHA = .5
WORKER_VALUES = None
WORKER_MODEL = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def spectrum(values, model, face, pixel, folds):
    """Return unweighted integer P/R bins in the A/B energy coordinate (q=1)."""
    channel = pixel // 32 if face == 'ring' else 24 + pixel % 32
    a = model['A_keV'][channel] / 1000.
    b = model['B'][channel]
    edges = a + b * np.arange(801) * .25
    p = np.diff(values[(face, 'prompt')][folds, pixel].sum(axis=0))
    r = np.diff(values[(face, 'random')][folds, pixel].sum(axis=0))
    assert b > 0 and np.all(p >= 0) and np.all(r >= 0)
    return edges, p, r


def bin_model(edges, params):
    """Exact integral of H*Phi((mu-E)/sigma)+C, with H,C in counts/MeV."""
    h, mu, sigma, c = params
    z = (np.asarray(edges) - mu) / sigma
    primitive = z * ndtr(-z) - np.exp(-.5*z*z)/np.sqrt(2*np.pi)
    return np.maximum(h*sigma*np.diff(primitive) + c*np.diff(edges), 0.)


def profiled_random(signal, p, r):
    """Analytic nonnegative random expectation in joint P/R Poisson likelihood."""
    aa = ALPHA*(1+ALPHA)
    bb = (1+ALPHA)*signal - ALPHA*(p+r)
    disc = np.sqrt(bb*bb + 4*aa*r*signal)
    # The rational form avoids cancellation when the fitted signal is large.
    positive = bb >= 0
    result = np.empty_like(signal)
    result[positive] = np.divide(2*r[positive]*signal[positive],
                                  disc[positive]+bb[positive],
                                  out=np.zeros(np.count_nonzero(positive)),
                                  where=(disc[positive]+bb[positive]) > 0)
    result[~positive] = (disc[~positive]-bb[~positive])/(2*aa)
    return np.maximum(result, 0.)


def deviance_residual(observed, expected):
    expected = np.maximum(expected, 1.e-12)
    dev = 2*(expected-observed+xlogy(observed, observed/expected))
    return np.sign(observed-expected)*np.sqrt(np.maximum(dev, 0.))


def residual(params, edges, p, r):
    signal = bin_model(edges, params)
    background = profiled_random(signal, p, r)
    return np.r_[deviance_residual(p, signal+ALPHA*background),
                 deviance_residual(r, background)]


def initial_location(edges, p, r):
    x = .5*(edges[1:]+edges[:-1])
    density = (p-ALPHA*r)/np.diff(edges)
    smooth = gaussian_filter1d(density, 1./np.median(np.diff(edges)))
    allowed = (x >= 20) & (x <= 110)
    if not allowed.any():
        return None
    indices = np.flatnonzero(allowed)
    peak = indices[np.argmax(smooth[indices])]
    if smooth[peak] <= 0 or peak == indices[0]:
        return None
    crossings = np.flatnonzero((np.arange(len(x)) > peak) & (smooth <= .5*smooth[peak]))
    if not len(crossings):
        return None
    k = crossings[0]
    if x[k] > 115 or k == 0:
        return None
    fraction = (.5*smooth[peak]-smooth[k-1])/(smooth[k]-smooth[k-1])
    return {'seed': float(x[k-1]+fraction*(x[k]-x[k-1])),
            'peak': float(x[peak]), 'height': float(smooth[peak])}


def fit_flank(edges, p, r, seed_info, window):
    if seed_info is None:
        return {'ok': False, 'flags': ['no_resolved_high_side_seed']}
    seed = seed_info['seed']
    lo, hi = seed+window[0], seed+window[1]
    use = np.flatnonzero((edges[:-1] >= lo) & (edges[1:] <= hi))
    if len(use) < 12:
        return {'ok': False, 'flags': ['too_few_fit_bins']}
    e = edges[use[0]:use[-1]+2]
    pp, rr = p[use], r[use]
    net = float(np.sum(pp-ALPHA*rr))
    record = {'window': [float(e[0]), float(e[-1])], 'seed': seed,
              'peak_seed': seed_info['peak'], 'net': net, 'prompt': float(pp.sum()),
              'random': float(rr.sum()), 'bin_count': len(use)}
    if net < 80:
        return dict(record, ok=False, flags=['insufficient_local_net_counts'])
    height = max(seed_info['height'], 1.)
    bounds = ([0., seed-5., .3, 0.], [height*15., seed+5., 12., height*5.])
    starts, results = [], []
    for sigma in [1.5, 4.]:
        x0 = [height*1.5, seed, sigma, max(.005*height, .01)]
        fit = least_squares(residual, x0, bounds=bounds, args=(e, pp, rr),
                            max_nfev=350, ftol=1.e-9, xtol=1.e-9, gtol=1.e-8,
                            x_scale='jac')
        starts.append({'params': fit.x.tolist(), 'deviance': float(2*fit.cost),
                       'success': bool(fit.success), 'nfev': fit.nfev})
        results.append(fit)
    fit = min(results, key=lambda f: f.cost)
    rank = int(np.linalg.matrix_rank(fit.jac))
    cov = np.linalg.pinv(fit.jac.T @ fit.jac)
    error = float(np.sqrt(max(cov[1, 1], 0)))
    h, mu, sigma, c = map(float, fit.x)
    dev = float(2*fit.cost)
    df = len(pp)-4  # One profiled random nuisance per P/R bin.
    flags = []
    if not fit.success:
        flags.append('optimizer_not_converged')
    if rank < 4 or error > 1.5 or error/max(mu, 1.) > .03:
        flags.append('location_uncertain')
    if min(mu-bounds[0][1], bounds[1][1]-mu) < .1:
        flags.append('location_at_search_boundary')
    if sigma <= .31 or sigma >= 11.9:
        flags.append('width_at_search_boundary')
    if mu-e[0] < .5 or e[-1]-mu < 2*sigma:
        flags.append('insufficient_flank_support')
    if h/(h+c) < .5:
        flags.append('weak_drop_contrast')
    if dev/df > 3.:
        flags.append('local_model_mismatch')
    if abs(starts[0]['params'][1]-starts[1]['params'][1]) > .1:
        flags.append('start_dependence')
    record.update(ok=not flags, flags=flags, params=[h, mu, sigma, c],
                  mu=mu, sigma=sigma, mu_error=error, covariance=cov.tolist(),
                  deviance=dev, dof=df, deviance_per_dof=dev/df,
                  rank=rank, starts=starts)
    return record


def fit_one(argument):
    run, face, pixel, spectra, model = argument
    edges, p, r = spectrum(spectra, model, face, pixel, [0, 1, 2, 3])
    seed = initial_location(edges, p, r)
    fits = {name: fit_flank(edges, p, r, seed, window)
            for name, window in WINDOWS.items()}
    for label, folds in [('even', [0, 2]), ('odd', [1, 3])]:
        e, pp, rr = spectrum(spectra, model, face, pixel, folds)
        fits[label] = fit_flank(e, pp, rr, seed, WINDOWS['main'])
    return {'run': run, 'ring': pixel//32+1, 'sector': pixel%32+1,
            'face': face, 'fits': fits}


def fit_worker(argument):
    return fit_one((*argument, WORKER_VALUES, WORKER_MODEL))


def run_fits(args):
    global WORKER_VALUES, WORKER_MODEL
    package = json.loads(PACKAGE_PATH.read_text())
    protocol = json.loads((TASK/'protocol.json').read_text())
    assert sha(PACKAGE_PATH) == protocol['inputs']['s3_fb_cor/parameters/low_pixel_20261102/coefficients.json']
    directory = TASK/'results'/args.label
    directory.mkdir(exist_ok=False)
    started = time.time()
    sectors = [1, 8, 16, 24, 32] if args.pilot else list(range(1, 33))
    rings = [24] if args.pilot else list(range(20, 25))
    if args.validation:
        assert not args.combine_historical
        frozen = json.loads((TASK/'candidate_frozen.json').read_text())
        assert sha(TASK/'candidate/coefficients.json') == frozen['package_sha256']
        assert (TASK/'fresh/completion.json').exists()
    for run in ['0616', '0641', '0659', '0674']:
        path = (TASK if args.validation else SOURCE)/'fresh'/f'spectra_{run}.root'
        if args.validation:
            receipt = json.loads((TASK/'fresh'/f'receipt_{run}.json').read_text())
            assert sha(path) == receipt['histogram_sha256']
        else:
            assert sha(path) == protocol['inputs'][str(path.relative_to(ROOT))]
        metadata, WORKER_VALUES = spectrum_file(path)
        inputs = {str(path): sha(path)}
        if args.combine_historical:
            other = ROOT/'codex/analysis/s3_source_calibration_refit_20261102/histograms'/f'spectra_{run}.root'
            other_meta, other_values = spectrum_file(other)
            assert metadata['alpha'] == other_meta['alpha'] == ALPHA
            assert metadata['input_metadata']['input_uuid'] == other_meta['input_metadata']['input_uuid']
            left = [x['read_interval'] for x in metadata['input_metadata']['segments']]
            right = [x['read_interval'] for x in other_meta['input_metadata']['segments']]
            assert sum(max(0, min(b, d)-max(a, c)) for a, b in left for c, d in right) == 0
            WORKER_VALUES = {key: value + other_values[key] for key, value in WORKER_VALUES.items()}
            inputs[str(other)] = sha(other)
        WORKER_MODEL = package['runs'][run]
        arguments = [(run, face, (ring-1)*32+sector-1)
                     for ring in rings for sector in sectors for face in ['ring', 'sector']]
        # Fork shares immutable input arrays; each worker returns only fit records.
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs,
                mp_context=multiprocessing.get_context('fork')) as pool:
            records = list(pool.map(fit_worker, arguments, chunksize=2))
        write(directory/f'fit_{run}.json', {'run': run, 'input': str(path),
              'input_sha256': sha(path), 'package_sha256': sha(PACKAGE_PATH),
              'combined_inputs': inputs, 'records': records})
        print(run, len(records), 'faces', sum(x['fits']['main']['ok'] for x in records),
              'main supported', f'{time.time()-started:.1f}s', flush=True)
    write(directory/'completion.json', {'status': 'complete', 'elapsed_seconds': time.time()-started,
          'script_sha256': sha(__file__), 'pilot': args.pilot, 'jobs': args.jobs})


def self_test():
    from scipy.integrate import quad
    from scipy.optimize import minimize_scalar
    e = np.linspace(20, 90, 281)
    params = np.array([100., 58., 2.1, .8])
    b = bin_model(e, params)
    direct = np.array([quad(lambda x: params[0]*ndtr((params[1]-x)/params[2])+params[3], a, c)[0]
                       for a, c in zip(e[:-1], e[1:])])
    assert np.max(np.abs(b-direct)) < 1.e-9
    for signal, p, r in [(10., 15., 4.), (0., 7., 6.), (20., 0., 0.), (1.e5, 1.e5, 1.)]:
        fitted = profiled_random(np.array([signal]), np.array([p]), np.array([r]))[0]
        def obj(x):
            return signal+(1+ALPHA)*x-xlogy(p, max(signal+ALPHA*x, 1e-15))-xlogy(r, max(x, 1e-15))
        numerical = minimize_scalar(obj, bounds=(0., max(p+r+signal, 20.)), method='bounded')
        assert obj(fitted) <= numerical.fun+1.e-7
    trials = []
    for gain, width in [(1., 1.5), (.7, 3.), (.9, 5.)]:
        truth = np.array([1000., 60.*gain, width, .2])
        rr = np.ones(len(e)-1)*.2
        pp = bin_model(e, truth)+ALPHA*rr
        fit = fit_flank(e, pp, rr, {'seed': truth[1], 'peak': truth[1]-5., 'height': 1000.}, (-8., 15.))
        assert abs(fit['mu']-truth[1]) < .005
        trials.append({'gain': gain, 'width': width, 'true_mu': truth[1], 'fit_mu': fit['mu']})
    write(TASK/'validation/mathematics.json', {'status': 'passed', 'bin_integrals': 'numerical quadrature agrees',
          'profile_random': 'independent scalar minimization agrees', 'known_location_recovery': trials,
          'limitation': 'Analytic-model checks do not validate the model on damaged-detector spectra.'})
    print('Mathematical checks passed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['self-test', 'fit'])
    parser.add_argument('--label', default='main')
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--combine-historical', action='store_true')
    parser.add_argument('--validation', action='store_true')
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if args.mode == 'self-test':
        self_test()
    else:
        run_fits(args)


if __name__ == '__main__':
    main()
