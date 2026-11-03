#!/usr/bin/env python3
"""Fit descriptive curves to the reconstructed pi CF prefix, with temporal checks.

Requires numpy, scipy, matplotlib. First reconstruct prefix.u32le.gz with
../windows/analyze_windows.py (gmpy2). No claim of a new law for fixed pi.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = Path(__file__).resolve().parent
PREFIX_SHA = 'ce0fd0c028276918ac51d04e122486dc628d283c0c6d2972cf27e2d4eae21739'
LN2 = math.log(2)
LEVY10 = math.pi**2 / (12 * LN2 * math.log(10))


def states(a, terminal):
    """Forward denominator ratios/logs; backward inverse complete quotients."""
    beta = np.zeros(len(a))
    logq = np.zeros(len(a))
    for n in range(1, len(a)):
        factor = float(a[n]) + beta[n-1]
        beta[n] = 1 / factor
        logq[n] = logq[n-1] + math.log10(factor)
    inverse_tail = np.empty(len(a) + 1)
    inverse_tail[-1] = terminal
    for j in range(len(a)-1, 0, -1):
        inverse_tail[j] = 1 / (float(a[j]) + inverse_tail[j+1])
    u = inverse_tail[1:len(a)]
    quality = u / (1 + beta[:-1] * u)
    return logq, quality


def linear_fit(x, y):
    xc, yc = x.mean(), y.mean()
    slope = float(np.dot(x-xc, y-yc) / np.dot(x-xc, x-xc))
    intercept = float(yc - slope*xc)
    return {'slope': slope, 'intercept': intercept}


def errors(actual, prediction):
    residual = actual-prediction
    return {'rmse': float(np.sqrt(np.mean(residual**2))),
            'mae': float(np.mean(abs(residual))),
            'maximum_absolute_error': float(np.max(abs(residual))),
            'last_error': float(residual[-1])}


def quality_cdf(c):
    c = np.asarray(c)
    return np.where(c <= .5, c/LN2,
                    (1-c+np.log(2*np.maximum(c, 1e-300)))/LN2)


def empirical_cdf(data, thresholds):
    return np.searchsorted(np.sort(data), thresholds, side='right')/len(data)


def tail_fit(data, kmin=10):
    z = data[data >= kmin].astype(float)
    # Discrete Pareto P(A=k | A>=kmin)=(kmin/k)^alpha-(kmin/(k+1))^alpha.
    def nll(alpha, observations=z):
        return float(np.sum(alpha*np.log(observations/kmin)
            - np.log(-np.expm1(-alpha*np.log1p(1/observations)))))
    fit = minimize_scalar(nll, bounds=(.4, 1.6), method='bounded',
                          options={'xatol': 1e-12})
    return {'minimum_a': kmin, 'tail_samples': len(z), 'alpha': float(fit.x),
            'tail_mass': len(z)/len(data), 'training_negative_log_likelihood': float(fit.fun)}


def tail_nll(data, alpha, kmin=10):
    z = data[data >= kmin].astype(float)
    return float(np.sum(alpha*np.log(z/kmin)
        -np.log(-np.expm1(-alpha*np.log1p(1/z)))))


def gauss_conditional_nll(data, kmin=10):
    z = data[data >= kmin].astype(float)
    # P(A=k)=ln(1+1/(k(k+2)))/ln2, normalized to A>=kmin.
    logp = np.log(np.log1p(1/(z*(z+2)))) - math.log(math.log1p(1/kmin))
    return float(-np.sum(logp))


def make_plot(a, ix, split, logq, c, fit, result, outdir):
    font = Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none'})
    blue, orange, gray = '#185f9e', '#c35524', '#525866'
    fig, axs = plt.subplots(2, 2, figsize=(13.6, 9.5), layout='constrained')
    fig.suptitle('π 连分数：整段曲线拟合与后半段检验', fontsize=20, weight='bold')
    fig.get_layout_engine().set(hspace=.12, wspace=.07)
    ax = axs[0, 0]
    thin = ix[::100]
    ax.plot(thin/10000, logq[thin]/10000, color=blue, label='逐项递推所得分母')
    ax.plot(thin/10000, (fit['slope']*thin+fit['intercept'])/10000,
            color=orange, ls='--', label='仅用前半段拟合')
    ax.axvline(split/10000, ls=':', color=gray, lw=1)
    ax.set(title='① 分母的对数近似线性增长', xlabel='收敛分数下标（万）',
           ylabel='log10 q（万）')
    ax.legend(loc='upper left', fontsize=9)
    ax.text(.98, .04, f"训练斜率 {fit['slope']:.6f}\n验证独立拟合斜率 "
            f"{result['denominator']['test_fit']['slope']:.6f}", transform=ax.transAxes,
            ha='right', va='bottom', fontsize=10)
    ax = axs[0, 1]
    residual = logq[thin]-(fit['slope']*thin+fit['intercept'])
    ax.plot(thin/10000, residual, color=blue, lw=1.2)
    ax.axhline(0, lw=.8, color=gray)
    ax.axvline(split/10000, ls=':', color=gray, lw=1)
    ax.axvspan(split/10000, ix[-1]/10000, color=orange, alpha=.06)
    ax.set(title='② 放大残差：直线无法逐项准确外推', xlabel='收敛分数下标（万）',
           ylabel='实际 log10 q − 训练直线')
    ax.text(.04, .05, f"后半段均方根误差：{result['denominator']['test_prediction_errors']['rmse']:.1f} 个十进制位",
            transform=ax.transAxes, fontsize=10)
    ax = axs[1, 0]
    k = np.unique(np.rint(np.geomspace(1, 1e4, 180)).astype(int))
    traina, testa = a[ix[ix < split]], a[ix[ix >= split]]
    survival = np.array([(testa >= t).mean() for t in k])
    alpha = result['partial_quotient_tail']['train_fit']['alpha']
    mass = result['partial_quotient_tail']['train_fit']['tail_mass']
    ax.loglog(k, survival, color=blue, label='后半段实际频率')
    ax.loglog(k, np.log1p(1/k)/LN2, ls='--', color=gray, label='Gauss 通用基线')
    kfit = k[k >= 10]
    ax.loglog(kfit, mass*(10/kfit)**alpha, ls=':', color=orange,
              label=f'前半段尾部拟合：指数 {alpha:.4f}')
    ax.set(title='③ 大部分商出现频率接近 1/K', xlabel='大项门槛 K', ylabel='a ≥ K 的比例')
    ax.legend(fontsize=9)
    ax = axs[1, 1]
    thresholds = np.geomspace(1e-4, .5, 200)
    testc = c[ix[ix >= split]]
    ax.loglog(thresholds, empirical_cdf(testc, thresholds), color=blue,
              label='后半段实际频率')
    ax.loglog(thresholds, thresholds/LN2, ls='--', color=gray,
              label='Lenstra 通用基线：c / ln 2')
    slope = result['quality_distribution']['train_linear_cdf_slope']
    ax.loglog(thresholds, slope*thresholds, ls=':', color=orange,
              label=f'前半段频率拟合：{slope:.4f}c')
    ax.set(title='④ 更小的归一化误差，出现频率近似正比下降',
           xlabel='质量门槛 c（越小越好）', ylabel='C ≤ c 的比例')
    ax.legend(fontsize=9)
    for ax in axs.flat:
        ax.grid(alpha=.17, which='major')
        ax.tick_params(labelsize=9)
    fig.supxlabel('实际使用下标 64–453228，共 453,165 项；样本来自纪录终点前缀，结果为探索分析。',
                  fontsize=10, color=gray)
    fig.savefig(outdir/'curve_fits.png', dpi=160)
    fig.savefig(outdir/'curve_fits.svg')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, default=HERE.parent/'windows/prefix.u32le.gz')
    parser.add_argument('--outdir', type=Path, default=HERE)
    args = parser.parse_args()
    args.outdir.mkdir(exist_ok=True, parents=True)
    raw = gzip.decompress(args.prefix.read_bytes())
    assert hashlib.sha256(raw).hexdigest() == PREFIX_SHA
    a = np.frombuffer(raw, dtype='<u4')
    assert len(a) == 453294 and a[0] == 3 and a[-1] == 12996958
    logq, c = states(a, 0.)
    _, c_alternate = states(a, 1.)
    # Omit the final 64 scorable convergents and the final unscorable coefficient.
    # At least 65 known coefficients enter every retained inverse forward tail.
    start, end = 64, len(a)-65
    split = (start+end)//2
    ix = np.arange(start, end)
    train, test = ix[ix < split], ix[ix >= split]
    fit = linear_fit(train.astype(float), logq[train])
    prediction = fit['slope']*test+fit['intercept']
    # A quadratic is a comparator, with no parameter selection on the test half.
    ztrain, ztest = train/len(a), test/len(a)
    quadratic = np.polynomial.polynomial.polyfit(ztrain, logq[train], 2)
    prediction_quadratic = np.polynomial.polynomial.polyval(ztest, quadratic)
    quality_thresholds = np.array([1e-4, 1e-3, .01, .05, .1, .2, .5])
    fit_grid = np.linspace(.01, .5, 50)
    ftrain = empirical_cdf(c[train], fit_grid)
    slope_c = float(np.dot(fit_grid, ftrain)/np.dot(fit_grid, fit_grid))
    test_grid = np.linspace(.001, .5, 1000)
    ftest = empirical_cdf(c[test], test_grid)
    fit_tail = tail_fit(a[train])
    test_tail = tail_fit(a[test])
    thresholds_a = [10, 100, 1000, 10000]
    result = {
        'input_sha256_u32le': PREFIX_SHA,
        'prefix_terms': len(a), 'used_positions': len(ix),
        'train_indices_inclusive': [int(train[0]), int(train[-1])],
        'test_indices_inclusive': [int(test[0]), int(test[-1])],
        'denominator': {
            'definition': 'log10(q_n) = slope*n + intercept',
            'train_fit': fit, 'test_fit': linear_fit(test.astype(float), logq[test]),
            'all_fit': linear_fit(ix.astype(float), logq[ix]),
            'generic_levy_slope': LEVY10,
            'test_prediction_errors': errors(logq[test], prediction),
            'quadratic_coefficients_for_n_divided_by_prefix_length': quadratic.tolist(),
            'quadratic_test_prediction_errors': errors(logq[test], prediction_quadratic),
            'endpoint_log10_q_453292': float(logq[453292]),
        },
        'partial_quotient_tail': {
            'method': 'Discrete Pareto maximum marginal likelihood above fixed kmin=10; dependence ignored for fitting, no IID standard errors claimed.',
            'train_fit': fit_tail, 'test_fit_for_shape_comparison_only': test_tail,
            'test_tail_negative_log_likelihood_fitted_train_alpha': tail_nll(a[test], fit_tail['alpha']),
            'test_tail_negative_log_likelihood_alpha_one': tail_nll(a[test], 1),
            'test_tail_negative_log_likelihood_conditional_gauss': gauss_conditional_nll(a[test]),
            'thresholds': [{'minimum_a': k, 'train_count': int((a[train] >= k).sum()),
                'test_count': int((a[test] >= k).sum()),
                'test_rate': float((a[test] >= k).mean()),
                'train_fitted_test_rate': fit_tail['tail_mass']*(10/k)**fit_tail['alpha'],
                'gauss_rate': math.log1p(1/k)/LN2} for k in thresholds_a]
        },
        'quality_distribution': {
            'definition': 'C_n=q_n^2*abs(pi-p_n/q_n)=u_(n+1)/(1+beta_n*u_(n+1))',
            'train_linear_cdf_slope': slope_c,
            'generic_cdf_slope': 1/LN2,
            'cdf_fit_grid': {'minimum': .01, 'maximum': .5, 'count': 50},
            'held_out_cdf_rmse_fit': float(np.sqrt(np.mean((ftest-slope_c*test_grid)**2))),
            'held_out_cdf_rmse_generic': float(np.sqrt(np.mean((ftest-test_grid/LN2)**2))),
            'observed_terminal_zero_vs_one_max_difference': float(np.max(abs(c[ix]-c_alternate[ix]))),
            'thresholds': [{'maximum_C': float(t), 'train_count': int((c[train] <= t).sum()),
                'test_count': int((c[test] <= t).sum()),
                'test_rate': float((c[test] <= t).mean()),
                'train_fitted_test_rate': slope_c*float(t), 'generic_rate': float(t/LN2)}
                for t in quality_thresholds],
            'full_cdf_max_absolute_difference_test': float(np.max(abs(empirical_cdf(c[test], np.linspace(0,1,1001))-quality_cdf(np.linspace(0,1,1001)))))
        },
        'cautions': [
            'The prefix is selected by a record endpoint. Removing the endpoint does not create a random sample.',
            'The temporal test is exploratory and reuses a prefix already examined in an earlier analysis. It is not a pristine validation sample.',
            'Gauss, Levy, and Lenstra laws apply to almost every real number; they are comparison models, not theorems asserting these limits for pi.',
            'Numerical curves use float64 recurrences; exact candidate verification is in the parent analysis. They are not new rigorous rational certificates.',
            'Tail endpoint uncertainty is suppressed by 65 omitted terminal positions; binary floating-point rounding remains.',
            'Cumulative fits can have excellent R2 while yielding large prediction errors. CDF residuals do not imply calibrated significance.',
            'No fitted curve identifies the position of the next record or supports extrapolating beyond the observed range.'
        ]
    }
    (args.outdir/'fit_results.json').write_text(json.dumps(result, indent=2)+'\n')
    make_plot(a, ix, split, logq, c, fit, result, args.outdir)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
