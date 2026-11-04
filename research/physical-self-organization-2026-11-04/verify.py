#!/usr/bin/env python3
"""Small independent checks; no pi approximation search or annealing claim."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parent
mp.mp.dps = 65


def screened_width(b, xi):
    """Sum Eq. (1) with an explicit positive geometric tail bound."""
    s = 1 / mp.mpf(xi)
    z = mp.exp(-s * b)
    upper_b = 2 * mp.exp(s) + mp.exp(-s)
    count = max(1, int(mp.ceil(mp.log(mp.mpf('1e-48') * (1-z) / upper_b)
                              / mp.log(z))))
    total = mp.mpf(0)
    for n in range(1, count + 1):
        r = mp.mpf(n * b)
        # Stable identity for r times the second difference; g=1.
        bracket = (2 * (mp.cosh(s)-1) + 2 * mp.sinh(s)/r
                   + 2 * (mp.cosh(s)+mp.sinh(s)/r)/(r*r-1))
        total += mp.exp(-s*r) * bracket
    tail = upper_b * z**(count+1) / (1-z)
    assert tail < mp.mpf('1e-46')
    return total, tail


def balanced_energy(rho, xi, neighbor_count=800):
    total = mp.mpf(0)
    for k in range(1, neighbor_count+1):
        v = k/rho
        d = int(mp.floor(v))
        t = v-d
        potential = lambda r: mp.exp(-mp.mpf(r)/xi)/r
        total += rho*((1-t)*potential(d)+t*potential(d+1))
    return total


def phi(u):
    u = mp.mpf(u)
    if not u:
        return mp.pi**2/3
    z = mp.exp(-u)
    return u*u/mp.expm1(u)-2*u*mp.log1p(-z)+2*mp.polylog(2,z)


def run_checks():
    out = {}
    root_checks = 0
    for q in range(2, 24):
        for p in range(1, q):
            if math.gcd(p, q) != 1:
                continue
            assert sorted((p*j) % q for j in range(q)) == list(range(q))
            # Actual complex Fourier average, independently of divisibility sum.
            theta = mp.mpf('0.173')
            indices = [1, 2, 3, 5, q, 2*q]
            direct = sum(sum(mp.cos(2*mp.pi*k*(theta+mp.mpf(j*p)/q))/k
                             for k in indices) for j in range(q))/q
            filtered = sum(mp.cos(2*mp.pi*k*theta)/k
                           for k in indices if k % q == 0)
            assert abs(direct-filtered) < mp.mpf('1e-58')
            root_checks += 1
    out['fourier_filter'] = {'coprime_pairs': root_checks, 'passed': True}

    for r in range(2, 101):
        diff = r*(Fraction(1,r-1)+Fraction(1,r+1)-Fraction(2,r))
        assert diff == Fraction(2,r*r-1)
    coulomb = []
    nmax = 10000
    for b in [2, 3, 7, 22, 113, 355]:
        partial = mp.fsum(2/(mp.mpf(n*b)**2-1) for n in range(1,nmax+1))
        exact = 1-(mp.pi/b)*mp.cot(mp.pi/b)
        tail_bound = mp.mpf(2)/((b*b-1)*nmax)
        assert 0 < exact-partial < tail_bound
        coulomb.append({'density_denominator': b, 'closed_width': str(exact),
                        'scaled_b2_width': str(b*b*exact),
                        'partial_sum_error': str(exact-partial),
                        'analytic_tail_bound': str(tail_bound)})
    out['coulomb_identity'] = {'passed': True, 'rows': coulomb}

    # Independent check via slopes of the balanced-state energy.
    slope_rows = []
    xi = mp.mpf(5)
    eps = mp.mpf('1e-15')
    for a,b in [(1,3),(2,5),(3,8),(7,22)]:
        rho = mp.mpf(a)/b
        # One-sided intervals avoid evaluating floor exactly at a rational kink.
        left = (balanced_energy(rho-eps,xi)-balanced_energy(rho-2*eps,xi))/eps
        right = (balanced_energy(rho+2*eps,xi)-balanced_energy(rho+eps,xi))/eps
        width,tail = screened_width(b,xi)
        assert abs(right-left-width) < mp.mpf('1e-38')
        slope_rows.append({'density': f'{a}/{b}', 'slope_jump': str(right-left),
                           'defect_sum': str(width)})
    out['independent_energy_slopes'] = {'passed': True, 'rows': slope_rows}

    fixed_xi = []
    xi = mp.mpf(5)
    prefactor = 2*(mp.cosh(1/xi)-1)
    for b in [100, 300, 1000]:
        width,tail = screened_width(b,xi)
        # For tiny widths the absolute tail target needs tightening by rescaling.
        s = 1/xi
        def bracket(r):
            return 2*(mp.cosh(s)-1)+2*mp.sinh(s)/r+2*(mp.cosh(s)+mp.sinh(s)/r)/(r*r-1)
        scaled = mp.fsum(mp.exp(-s*(n-1)*b)*bracket(mp.mpf(n*b))
                         for n in range(1,10))/prefactor
        fixed_xi.append({'b': b, 'xi': str(xi), 'ratio_to_fixed_xi_asymptotic': str(scaled)})
    assert abs(mp.mpf(fixed_xi[-1]['ratio_to_fixed_xi_asymptotic'])-1) < mp.mpf('0.011')
    out['fixed_screening_asymptotic'] = {'passed': True, 'rows': fixed_xi}

    joint = []
    for u in [mp.mpf('0.25'),mp.mpf(1),mp.mpf(3)]:
        errors = []
        for b in [20,80,320]:
            width,tail = screened_width(b,mp.mpf(b)/u)
            relative = abs(b*b*width/phi(u)-1)
            errors.append(relative)
            joint.append({'u': str(u),'b': b,'scaled_width': str(b*b*width),
                          'predicted_limit': str(phi(u)),'relative_error': str(relative)})
        assert errors[2] < errors[1] < errors[0]
        assert errors[2] < mp.mpf('0.00003')
    out['joint_limit'] = {'passed': True, 'rows': joint}

    # All occupation states of a 12-site periodic ring. The pair interaction is
    # a periodic-image sum; omitted self-image energies are constant at fixed N.
    length = 12
    xi_float = 2.0
    pair = {}
    for i in range(length):
        for j in range(i+1,length):
            d = j-i
            pair[i,j] = sum(math.exp(-abs(d+m*length)/xi_float)/abs(d+m*length)
                            for m in range(-12,13))
    minimum = [math.inf]*(length+1)
    for mask in range(1 << length):
        sites = [i for i in range(length) if mask >> i & 1]
        energy = sum(pair[i,j] for i,j in itertools.combinations(sites,2))
        minimum[len(sites)] = min(minimum[len(sites)],energy)
    maximum_difference = 0.0
    for count in range(1,length):
        sites = [i for i in range(length)
                 if ((i+1)*count)//length-(i*count)//length]
        energy = sum(pair[i,j] for i,j in itertools.combinations(sites,2))
        maximum_difference = max(maximum_difference, abs(energy-minimum[count]))
        assert abs(energy-minimum[count]) < 1e-12
    out['finite_ground_states'] = {'passed': True,'sites': length,
                                   'states_checked': 1 << length,
                                   'max_energy_difference': maximum_difference,
                                   'not_a_dynamical_convergence_test': True}
    return out


def plot_scales():
    import numpy as np
    from scipy.special import k0
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    b_values = np.unique(np.geomspace(2,600,180).astype(int))
    x = b_values.astype(float)
    fig, axes = plt.subplots(1,2,figsize=(11.6,4.4),layout='constrained')
    reference = x*x*(1-(np.pi/x)/np.tan(np.pi/x))
    axes[0].plot(x,reference,color='#26394c',lw=2,label='Local Coulomb limit')
    for xi in [5,20,80]:
        values = []
        s = 1.0/xi
        for b in b_values:
            ns = np.arange(1,max(2,int(math.ceil(48/(s*b))))+1,dtype=float)
            r = ns*b
            bracket = (4*np.sinh(s/2)**2+2*np.sinh(s)/r
                       +2*(np.cosh(s)+np.sinh(s)/r)/(r*r-1))
            values.append(b*b*np.sum(np.exp(-s*r)*bracket))
        axes[0].plot(x,values,lw=1.8,label=fr'Screening length $\xi={xi}$')
    axes[0].axhline(np.pi**2/3,color='#718096',ls=':',lw=1)
    axes[0].set(xscale='log',xlabel='Denominator b of particle density a/b',
                ylabel=r'$b^2 W_\xi(b)/g$',title='Collective stability to chemical potential')
    axes[0].legend(frameon=False,fontsize=8)

    q = np.arange(2,151,dtype=float)
    for h in [.02,.05,.1]:
        z = 2*np.pi*h*q
        # Leading harmonic peak-to-trough, not a claimed exact barrier.
        axes[1].semilogy(q,8*k0(z),lw=1.8,label=fr'Row separation $h={h}$')
    axes[1].set(xlabel='Denominator q of rigid-chain spacing p/q',
                ylabel=r'Leading harmonic peak-to-trough / $|g|$',
                title='Coulomb substrate: high harmonics are filtered')
    axes[1].legend(frameon=False,fontsize=8)
    for ax in axes:
        ax.grid(alpha=.18)
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Physical locking scales — neither panel is a pi approximation error',fontsize=12)
    fig.savefig(ROOT/'scales.png',dpi=180)
    fig.savefig(ROOT/'scales.svg')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot',action='store_true')
    args = parser.parse_args()
    results = run_checks()
    (ROOT/'verification.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
    if args.plot:
        plot_scales()
    print(json.dumps({key: {'passed': value['passed']} for key,value in results.items()},indent=2))
