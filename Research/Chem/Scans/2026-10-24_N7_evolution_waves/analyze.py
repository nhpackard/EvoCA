"""CHEM-N7 analysis: results.csv -> per-condition table (printed and
summary.csv) and evolution.png (one panel per measure, seeds as dots).

Regime: a run is 'waves' in a window when site_period finds a period there.
Demography: alive density, turnover (births per alive cell per step) and
generation time 1 / turnover. Activity slopes are per step; populations
turn over at different rates, so they are also given per generation
(slope x generation time)."""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = ['W1', 'W2', 'S1', 'T1', 'C0']
LABEL = {'W1': 'W1 waves\none-way', 'W2': 'W2 waves\ntwo-way', 'S1': 'S1 speckle\n(wdiff 0.08)',
         'T1': 'T1 one-way\ntox 1.0', 'C0': 'C0 waste\noff'}
LATE = ('period_t10000', 'period_t15000', 'period_t20000')
MEASURES = [
    ('wave_windows_late', 'late windows with a period (of 3)'),
    ('alive_density_mean', 'alive density'),
    ('gen_time', 'generation time (steps)'),
    ('dyn_excess_pc_slope', 'dyn excess slope (per step)'),
    ('dyn_excess_per_gen', 'dyn excess per generation'),
    ('eg_excess_pc_slope', 'eg excess slope (per step)'),
    ('eg_excess_per_gen', 'eg excess per generation'),
    ('a_new_final', 'A_new (final)'),
    ('n_distinct_genomes_mean', 'distinct genomes (mean)'),
    ('unique_top_genomes', 'dominant-genome changes'),
    ('negene_final', 'egenes per cell (final)'),
    ('cog_specificity_final', 'egene specificity (final)'),
]


def load():
    d = pd.read_csv(os.path.join(HERE, 'results.csv')).copy()
    gen = 1.0 / d['births_per_alive_step']
    return pd.concat([d, pd.DataFrame({
        'wave_windows_late': d[list(LATE)].notna().sum(axis=1),
        'gen_time': gen,
        'dyn_excess_per_gen': d['dyn_excess_pc_slope'] * gen,
        'eg_excess_per_gen': d['eg_excess_pc_slope'] * gen})], axis=1)


def main():
    d = load()
    cols = [m for m, _ in MEASURES] + ['period_t20000', 'births_per_alive_step', 'excess_pc_slope', 'a_new_slope',
                                       'mean_intake_mean', 'repro_interval_mean']
    g = d.groupby('condition')[cols]
    summ = pd.concat({'mean': g.mean(), 'sd': g.std()}, axis=1).reindex(ORDER)
    summ.to_csv(os.path.join(HERE, 'summary.csv'))
    pd.set_option('display.width', 250)
    print(g.mean().reindex(ORDER).T.round(4).to_string())
    print('\nsd across seeds:')
    print(g.std().reindex(ORDER).T.round(4).to_string())
    print('\nperiod by window (seed mean, nan = no period in any seed):')
    pw = [c for c in d.columns if c.startswith('period_t')]
    print(d.groupby('condition')[pw].mean().reindex(ORDER).round(0).to_string())
    figure(d)


def figure(d):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    nc = 4
    nr = -(-len(MEASURES) // nc)
    fig, axes = plt.subplots(nr, nc, figsize=(3.4 * nc, 2.7 * nr), squeeze=False)
    x = np.arange(len(ORDER))
    for ax, (m, title) in zip(axes.flat, MEASURES):
        for i, c in enumerate(ORDER):
            v = d.loc[d.condition == c, m].astype(float).values
            ax.scatter(np.full(len(v), i) + np.linspace(-0.12, 0.12, len(v)), v, s=14, color='#2a78d6', zorder=3)
            if len(v):
                ax.hlines(np.nanmean(v), i - 0.25, i + 0.25, color='#1c5cab', lw=2, zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels([LABEL[c] for c in ORDER], fontsize=6.5)
        ax.set_title(title, fontsize=8.5)
        ax.tick_params(axis='y', labelsize=7)
        ax.grid(axis='y', alpha=0.25, lw=0.5)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    for ax in axes.flat[len(MEASURES):]:
        ax.axis('off')
    fig.suptitle('CHEM-N7: evolution in vs out of the wave regime (N=512, 20000 steps, 4 seeds; bar = mean)', fontsize=9.5)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'evolution.png'), dpi=110)
    plt.close(fig)


if __name__ == '__main__':
    main()
