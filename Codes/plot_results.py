# -*- coding: utf-8 -*-
"""
plot_results.py - plot the CSV files written by extract_results.py.
Run in an ordinary Python with matplotlib:   python plot_results.py PREFIX
(PREFIX is the name used by extract_results.py, e.g. PS_v5_C6_dqs_m025)
Writes PREFIX_stress_strain.png, PREFIX_energy_checks.png, PREFIX_lateral.png
"""
import csv
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load(path):
    with open(path) as f:
        rd = csv.DictReader(f)
        cols = dict((k, []) for k in rd.fieldnames)
        for row in rd:
            for k in rd.fieldnames:
                try:
                    cols[k].append(float(row[k]))
                except ValueError:
                    cols[k].append(float('nan'))
    return cols


def main():
    prefix = sys.argv[1]
    h = load(prefix + '_history.csv')
    eps = [100.0 * e for e in h['strain']]

    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(eps, h['stress_MPa'], lw=1.6)
    ax.set_xlabel('Nominal strain (%)')
    ax.set_ylabel('Nominal stress (MPa)')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(prefix + '_stress_strain.png', dpi=200)

    fig, axs = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    axs[0].plot(eps, h['ALLIE_mJ'], label='ALLIE')
    axs[0].plot(eps, h['ALLKE_mJ'], label='ALLKE')
    axs[0].plot(eps, h['ALLAE_mJ'], label='ALLAE')
    axs[0].set_ylabel('Energy (mJ)')
    axs[0].legend()
    axs[0].grid(alpha=0.3)
    axs[1].plot(eps, h['KE_over_IE'], label='ALLKE / ALLIE')
    axs[1].plot(eps, h['AE_over_IE'], label='ALLAE / ALLIE')
    axs[1].axhline(0.10, color='k', ls='--', lw=0.8)
    axs[1].set_ylabel('Energy ratio')
    axs[1].set_xlabel('Nominal strain (%)')
    axs[1].legend()
    axs[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(prefix + '_energy_checks.png', dpi=200)

    try:
        lat = load(prefix + '_lateral.csv')
        fig, ax = plt.subplots(figsize=(6, 4.2))
        ax.plot([100.0 * e for e in lat['strain_load']], lat['poisson_ratio_eff'], lw=1.6)
        ax.set_xlabel('Nominal strain (%)')
        ax.set_ylabel('Effective Poisson ratio')
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(prefix + '_lateral.png', dpi=200)
    except IOError:
        pass
    print('plots written with prefix', prefix)


main()
