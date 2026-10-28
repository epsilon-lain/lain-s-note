"""A reproducible neural-regression experiment for the beta-coordinate model.

Labels and finite-tail certificates are generated with rational arithmetic.
Training itself uses ordinary floating-point scikit-learn MLPRegressor.
The experiment is finite; it makes no claim about pi being badly approximable.

Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 for reproducible small runs.
"""

from fractions import Fraction as F
from pathlib import Path
import json
import time
import warnings

import numpy as np
import sklearn
from sklearn.linear_model import LinearRegression
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler

from rational_base_gamma_checks import pi_bounds
from base_regression_checks import fit_interval, interval_ratio
from same_base_curve_line_checks import exact_digits


OUTPUT = Path(__file__).with_name('base_ml_multiscale_output')
POST_DIGITS = 12
SEED = 2848


def exact_multiscale_checks():
    import random
    rng = random.Random(SEED)
    for _ in range(500):
        den = rng.randrange(1, 40)
        base = 8 + F(rng.randrange(den + 1), den)
        line = 1 / (base - 1)
        prefix = rng.randrange(1, 13)
        while True:
            phase = F(rng.randrange(10000), 10000)
            if (base * phase).__floor__() != 1:
                break
        z = line + (phase - line) / base ** prefix
        digits, rem = exact_digits(z, base, prefix + POST_DIGITS)
        assert digits[:prefix] == [1] * prefix
        assert digits[prefix] != 1 and rem[prefix] == phase
        assert abs(phase - line) >= F(1, 72)
        global_D = sum(((d - 1) * base ** (-2*j)
                        for j, d in enumerate(digits, 1)), F())
        tail_D = sum(((d - 1) * base ** (-2*j)
                      for j, d in enumerate(digits[prefix:], 1)), F())
        assert global_D == tail_D / base ** (2 * prefix)
        global_Q = sum(((d-1)**2 * base ** (-j) / (base-1)
                        for j, d in enumerate(digits, 1)), F())
        tail_Q = sum(((d-1)**2 * base ** (-j) / (base-1)
                      for j, d in enumerate(digits[prefix:], 1)), F())
        assert global_Q == tail_Q / base ** prefix
        # Slopes of the finite curve completed by an all-one tail.
        epsilon0 = sum(((d-1)*base**(-j) for j,d in enumerate(digits,1)), F())
        epsilon1 = sum(((d-1)*base**(-j) for j,d in enumerate(digits[1:],1)), F())
        slope0 = -6*epsilon0 + 6*(base+1)*global_D
        slope1 = -6*epsilon1 + 6*(base+1)*sum(
            ((d-1)*base**(-2*j) for j,d in enumerate(digits[1:],1)), F())
        assert base**2*slope0 == slope1 + 6*digits[0] - 6*(base-1)*(line+epsilon1)
    return 500


def midpoint(interval):
    return (interval[0] + interval[1]) / 2


def sample(base, xlo, xhi):
    fit = fit_interval(base, xlo, xhi, length=52)
    line = 1 / (base - 1)
    error = xlo - line, xhi - line
    assert error[0] * error[1] > 0
    error_abs = min(abs(error[0]), abs(error[1])), max(abs(error[0]), abs(error[1]))
    prefix = next(j for j, d in enumerate(fit['digits']) if d != 1)
    assert prefix + POST_DIGITS < len(fit['digits'])
    tail_lower = sum((F(d) * base ** (-j)
                      for j, d in enumerate(fit['digits'][prefix:prefix+POST_DIGITS], 1)), F())
    tail_upper = tail_lower + base ** (-POST_DIGITS)
    tail_mid = (tail_lower + tail_upper) / 2
    # Fit the 12-layer local curve, filling the unobserved tail by height 1/2.
    # Its mean is tail_mid; the two fitted parameters preserve this mean.
    local_compressed = sum((F(d) * base ** (-2*j)
                            for j,d in enumerate(fit['digits'][prefix:prefix+POST_DIGITS],1)), F())
    local_slope = 6*(base+1)*local_compressed + 3*base**(-2*POST_DIGITS) - 6*tail_mid
    local_intercept = tail_mid - local_slope/2
    assert local_intercept + local_slope/2 == tail_mid
    local_phase_lo = line + base ** prefix * error[0]
    local_phase_hi = line + base ** prefix * error[1]
    assert tail_lower <= local_phase_lo <= local_phase_hi <= tail_upper
    assert min(abs(local_phase_lo - line), abs(local_phase_hi - line)) >= F(1, 72)
    finite_error = base ** (-prefix) * (tail_mid - line)
    actual_error = midpoint(error)
    assert finite_error * actual_error > 0
    relative_bound = F(36) * base ** (-POST_DIGITS)
    assert abs(finite_error - actual_error) <= relative_bound * error_abs[0]
    sq_error = error_abs[0] ** 2, error_abs[1] ** 2
    H_interval = interval_ratio(fit['quadratic_correction'], sq_error)
    H = midpoint(H_interval)
    slope = midpoint(fit['slope'])
    mse = midpoint(fit['mse'])
    # V = Q/|epsilon|, Q = MSE + epsilon^2 + slope^2/12.
    V = (mse + actual_error ** 2 + slope ** 2 / 12) / abs(actual_error)
    assert H * actual_error > 0
    assert abs(H) <= 8160
    assert F(1, 8) <= V <= 17
    q = base.numerator - base.denominator
    p = 3 * base.numerator - 2 * base.denominator
    return {
        'base': str(base), 'b': float(base), 'prefix': prefix,
        'tail_mid': float(tail_mid), 'finite_error': float(finite_error),
        'local_slope': float(local_slope), 'local_intercept': float(local_intercept),
        'epsilon': float(actual_error), 'slope': float(slope), 'mse': float(mse),
        'H': float(H), 'V': float(V),
        'H_label_interval_width': float(H_interval[1] - H_interval[0]),
        'slope_label_interval_width': float(fit['slope'][1] - fit['slope'][0]),
        'mse_label_interval_width': float(fit['mse'][1] - fit['mse'][0]),
        'finite_phase_relative_error': float(abs(finite_error / actual_error - 1)),
        'finite_phase_relative_error_bound': float(relative_bound),
        'q': q, 'pi_fraction': f'{p}/{q}',
    }


def fit_network(inputs, targets):
    scaler_x = StandardScaler().fit(inputs)
    scaler_y = StandardScaler().fit(targets)
    model = MLPRegressor(hidden_layer_sizes=(32, 32), activation='tanh',
                         solver='lbfgs', alpha=1e-6, random_state=SEED,
                         max_iter=10000, max_fun=150000, tol=1e-7)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model.fit(scaler_x.transform(inputs), scaler_y.transform(targets))
    messages = [str(w.message) for w in caught]
    return model, scaler_x, scaler_y, messages


def predict(network, inputs):
    model, sx, sy, _ = network
    return sy.inverse_transform(model.predict(sx.transform(inputs)))


def phase_prediction(network, rows):
    inputs = np.array([[r['b'], r['tail_mid']] for r in rows])
    raw = predict(network, inputs)
    H_abs = np.clip(np.exp(np.clip(raw[:, 0], -30, 30)), 0, 8160)
    V = np.clip(np.exp(np.clip(raw[:, 1], -30, 30)), 1/8, 17)
    eps = np.array([r['finite_error'] for r in rows])
    H = np.sign(eps) * H_abs
    slope = -6 * eps + eps ** 2 * H
    mse = np.maximum(0, np.abs(eps) * V - eps ** 2 - slope ** 2 / 12)
    return np.column_stack((slope, mse)), np.column_stack((H, V))


def metrics(truth, pred):
    error = pred - truth
    return {'slope_mae': float(np.mean(np.abs(error[:, 0]))),
            'slope_max_absolute_error': float(np.max(np.abs(error[:, 0]))),
            'mse_mae': float(np.mean(np.abs(error[:, 1]))),
            'mse_max_absolute_error': float(np.max(np.abs(error[:, 1])))}


def export_network(network):
    model, sx, sy, messages = network
    return {'activation': 'tanh', 'output_activation': model.out_activation_,
            'input_mean': sx.mean_.tolist(), 'input_scale': sx.scale_.tolist(),
            'target_mean': sy.mean_.tolist(), 'target_scale': sy.scale_.tolist(),
            'coefs': [array.tolist() for array in model.coefs_],
            'intercepts': [array.tolist() for array in model.intercepts_],
            'iterations': int(model.n_iter_), 'training_loss': float(model.loss_),
            'warnings': messages}


def predict_exported(model, inputs):
    value = (inputs - np.asarray(model['input_mean'])) / np.asarray(model['input_scale'])
    for i, (weights, bias) in enumerate(zip(model['coefs'], model['intercepts'])):
        value = value @ np.asarray(weights) + np.asarray(bias)
        if i < len(model['coefs']) - 1:
            value = np.tanh(value)
    return value * np.asarray(model['target_scale']) + np.asarray(model['target_mean'])


def plot_results(data, test_indices, raw_pred_all, phase_pred_all, critical, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)
    b = np.array([r['b'] for r in data])
    slope = np.array([r['slope'] for r in data])
    axes[0].plot(b, slope, color='#142d4e', linewidth=1.8, label='Exact coordinate fit')
    axes[0].plot(b, raw_pred_all[:, 0], color='#d2692a', linewidth=1.1, label='Raw-base neural net')
    axes[0].plot(b, phase_pred_all[:, 0], color='#188268', linewidth=1.0, linestyle='--', label='Local-phase neural net')
    axes[0].set(xlabel='Base n', ylabel='Coordinate-fit slope s(n)', title='The coarse curve can look accurate')
    axes[0].legend(frameon=False, fontsize=8)
    labels = ['355/113', '103993/33102', '104348/33215', '245850922/78256779']
    rows = [r for r in critical if r['pi_fraction'] in labels]
    q = [r['q'] for r in rows]
    for key, label, color, marker in [
        ('true_slope_score', 'Exact slope score', '#142d4e', 'o'),
        ('raw_network_score', 'Raw-base neural net', '#d2692a', 'x'),
        ('local_phase_network_score', 'Local-phase neural net', '#188268', '+')]:
        axes[1].plot(q, [r[key] for r in rows], marker=marker, linestyle='none',
                     color=color, markersize=7, label=label)
    axes[1].set(xscale='log', yscale='log', xlabel='Reduced denominator q',
                ylabel='q^2 |predicted slope| / 6', title='Fine rational scales expose the error')
    axes[1].legend(frameon=False, fontsize=8)
    idx = np.array(test_indices)
    exact_H = np.array([r['H'] for r in data])[idx]
    predicted_H = np.array(result['phase_test_predictions'])[:, 0]
    axes[2].scatter(exact_H, predicted_H, s=12, alpha=0.65, color='#188268')
    span = [min(exact_H.min(), predicted_H.min()), max(exact_H.max(), predicted_H.max())]
    axes[2].plot(span, span, color='#142d4e', linewidth=1)
    axes[2].set(xlabel='Exact scaled correction H', ylabel='Neural estimate of H',
                title='Learn the exit phase after rescaling')
    fig.suptitle('Pi digit-coordinate regression: global fit versus local scale', fontsize=13)
    fig.savefig(OUTPUT / 'base_ml_multiscale.png', dpi=170)
    fig.savefig(OUTPUT / 'base_ml_multiscale.svg')
    plt.close(fig)


def plot_summary(critical):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    wanted = ['355/113', '103993/33102', '104348/33215', '245850922/78256779']
    rows = [next(r for r in critical if r['pi_fraction'] == name) for name in wanted]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    x = np.arange(len(rows))
    for key, label, color, marker, size in [
        ('true_slope_score', 'Exact coordinate-fit score', '#142d4e', 'o', 7),
        ('raw_network_score', 'Base-only neural network', '#d2692a', 'x', 8),
        ('local_phase_network_score', 'Multiscale neural network', '#188268', '+', 10)]:
        axis.plot(x, [r[key] for r in rows], linestyle='none', marker=marker,
                  color=color, markersize=size, label=label)
    axis.set_xticks(x, ['355/113', '103993/33102', '104348/33215', '245850922/\n78256779'])
    axis.set(yscale='log', ylabel='Normalized slope score: q² |s| / 6',
             xlabel='Held-out rational approximation p/q',
             title='Retaining the fine rational scale in a learned fit')
    axis.grid(axis='y', alpha=0.18)
    axis.set_xlim(-0.35, 3.35)
    axis.legend(frameon=False, loc='upper left', fontsize=9)
    fig.savefig(OUTPUT/'critical_scale_summary.png', dpi=170)
    fig.savefig(OUTPUT/'critical_scale_summary.svg')
    plt.close(fig)


def main():
    started = time.time()
    OUTPUT.mkdir(exist_ok=True)
    plo, phi = pi_bounds()
    xlo, xhi = plo - 3, phi - 3
    data_path = OUTPUT / 'dataset.json'
    if data_path.exists():
        data = json.loads(data_path.read_text())
        if 'local_slope' not in data[0]:
            data = [sample(8 + F(i, 1000), xlo, xhi) for i in range(1001)]
            data_path.write_text(json.dumps(data, separators=(',', ':')))
    else:
        data = [sample(8 + F(i, 1000), xlo, xhi) for i in range(1001)]
        data_path.write_text(json.dumps(data, separators=(',', ':')))
    train_idx = [i for i in range(len(data)) if i % 5 in (1, 2, 3)]
    validation_idx = [i for i in range(len(data)) if i % 5 == 4]
    test_idx = [i for i in range(len(data)) if i % 5 == 0]
    X_raw = np.array([[r['b']] for r in data])
    Y_raw = np.array([[r['slope'], r['mse']] for r in data])
    X_phase = np.array([[r['b'], r['tail_mid']] for r in data])
    Y_phase = np.array([[np.log(abs(r['H'])), np.log(r['V'])] for r in data])
    linear = LinearRegression().fit(X_raw[train_idx], Y_raw[train_idx])
    raw = fit_network(X_raw[train_idx], Y_raw[train_idx])
    phase = fit_network(X_phase[train_idx], Y_phase[train_idx])
    raw_all = predict(raw, X_raw)
    phase_all, phase_features_all = phase_prediction(phase, data)
    result = {'configuration': {
        'grid': '8+i/1000, i=0,...,1000', 'train': len(train_idx),
        'validation': len(validation_idx), 'test': len(test_idx),
        'split_rule': 'train i%5 in {1,2,3}; validation i%5=4; test i%5=0',
        'hyperparameter_selection': 'Architecture fixed; optimizer limits/tolerance adjusted after convergence warnings; exploratory fixed-split evaluation',
        'network': '32,32 tanh; lbfgs; alpha=1e-6; max_iter=10000; tol=1e-7; seed=2848',
        'raw_inputs': ['base'], 'raw_outputs': ['slope', 'absolute MSE'],
        'phase_inputs': ['base', 'midpoint from 12 post-exit digits'],
        'phase_outputs': ['log |H|', 'log V'],
        'structured_reconstruction': 'epsilon_hat=n^-N*(tail_mid-1/(n-1)); s=-6*epsilon_hat+epsilon_hat^2*H_hat',
        'post_exit_digits': POST_DIGITS, 'numpy': np.__version__, 'sklearn': sklearn.__version__},
        'exact_multiscale_checks': exact_multiscale_checks(),
        'label_checks': {
            'max_slope_interval_width': max(r['slope_label_interval_width'] for r in data),
            'max_mse_interval_width': max(r['mse_label_interval_width'] for r in data),
            'max_H_interval_width': max(r['H_label_interval_width'] for r in data),
            'uniform_phase_relative_error_bound': float(F(36) / 8 ** POST_DIGITS)},
        'validation_metrics': {}, 'test_metrics': {}}
    for name, idx in [('validation_metrics', validation_idx), ('test_metrics', test_idx)]:
        result[name] = {
            'linear_base_regression': metrics(Y_raw[idx], linear.predict(X_raw[idx])),
            'raw_base_neural_network': metrics(Y_raw[idx], raw_all[idx]),
            'local_phase_neural_network': metrics(Y_raw[idx], phase_all[idx])}
    critical_bases = [F(121, 15), F(129, 16), F(37789, 4687), F(37918, 4703),
                      F(89337364, 11080585)]
    critical_rows = [sample(base, xlo, xhi) for base in critical_bases]
    critical_raw = predict(raw, np.array([[r['b']] for r in critical_rows]))
    critical_phase, critical_HV = phase_prediction(phase, critical_rows)
    critical_report = []
    for r, raw_prediction, phase_prediction_row, learned_features in zip(
            critical_rows, critical_raw, critical_phase, critical_HV):
        q = r['q']
        row = {'base': r['base'], 'pi_fraction': r['pi_fraction'], 'q': q,
               'leading_ones': r['prefix'], 'true_slope': r['slope'],
               'true_slope_score': q*q*abs(r['slope'])/6,
               'true_kappa': q*q*abs(r['epsilon']),
               'raw_network_score': q*q*abs(float(raw_prediction[0]))/6,
               'local_phase_network_score': q*q*abs(float(phase_prediction_row[0]))/6,
               'local_analytic_score': q*q*abs(r['finite_error']),
               'finite_phase_relative_error': r['finite_phase_relative_error'],
               'exact_H': r['H'], 'predicted_H': float(learned_features[0]),
               'exact_V': r['V'], 'predicted_V': float(learned_features[1])}
        row['raw_score_relative_error'] = abs(row['raw_network_score']/row['true_slope_score']-1)
        row['phase_score_relative_error'] = abs(row['local_phase_network_score']/row['true_slope_score']-1)
        critical_report.append(row)
    result['critical_holdouts'] = critical_report
    result['network_training'] = {'raw': {'iterations': raw[0].n_iter_, 'warnings': raw[3]},
                                  'phase': {'iterations': phase[0].n_iter_, 'warnings': phase[3]}}
    B = float(1 + 1 / midpoint((xlo, xhi)))
    result['raw_prediction_at_special_base'] = predict(raw, np.array([[B]]))[0].tolist()
    root_lo, root_hi = 8.0, 9.0
    assert predict(raw, np.array([[root_lo]]))[0,0] > 0
    assert predict(raw, np.array([[root_hi]]))[0,0] < 0
    for _ in range(70):
        mid = (root_lo + root_hi) / 2
        if mid == root_lo or mid == root_hi:
            break
        if predict(raw, np.array([[mid]]))[0,0] > 0:
            root_lo = mid
        else:
            root_hi = mid
    result['raw_network_zero_estimate'] = (root_lo + root_hi) / 2
    result['raw_network_zero_error_approx'] = result['raw_network_zero_estimate'] - B
    result['phase_test_predictions'] = phase_features_all[test_idx].tolist()
    result['phase_test_correction_mae'] = float(np.mean(np.abs(
        phase_features_all[test_idx, 0] - np.array([r['H'] for r in data])[test_idx])))
    models = {'configuration': result['configuration'], 'raw_model': export_network(raw),
              'phase_model': export_network(phase)}
    raw_export_error = np.max(np.abs(predict_exported(models['raw_model'], X_raw[test_idx])-raw_all[test_idx]))
    phase_export_error = np.max(np.abs(predict_exported(models['phase_model'], X_phase[test_idx])-predict(phase, X_phase[test_idx])))
    assert raw_export_error < 1e-12 and phase_export_error < 1e-12
    result['export_checks_max_absolute_error'] = {'raw': float(raw_export_error), 'phase': float(phase_export_error)}
    (OUTPUT/'models.json').write_text(json.dumps(models, separators=(',', ':')))
    plot_results(data, test_idx, raw_all, phase_all, critical_report, result)
    plot_summary(critical_report)
    result.pop('phase_test_predictions')
    result['elapsed_seconds'] = time.time()-started
    (OUTPUT/'results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2))
    print(json.dumps({key:result[key] for key in ['test_metrics','critical_holdouts','network_training']},
                     ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
