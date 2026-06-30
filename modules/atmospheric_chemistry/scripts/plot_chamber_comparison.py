#!/usr/bin/env python3
"""
plot_chamber_comparison.py -- Plot MOOSE vs F0AM gold comparison figures.

Generates 2x2 panel figures (C5H8, OH, NOx, O3) for each chamber scenario
(S1, S2, S3, S2b), comparing MOOSE CSV output against F0AM gold reference.

Usage:
    # Plot all 4 scenarios using default paths
    python3 scripts/plot_chamber_comparison.py

    # Check files only (no plotting)
    python3 scripts/plot_chamber_comparison.py --check

    # Custom paths for a single scenario
    python3 scripts/plot_chamber_comparison.py \\
        --moose /path/to/moose.csv \\
        --gold /path/to/gold.csv \\
        --save /path/to/output.pdf

Dependencies: numpy, matplotlib
"""

import argparse, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Molecules/cm^3 -> ppb conversion at 298K, 1013mbar
M = 2.46e19

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 9,
    'axes.labelsize': 10, 'axes.titlesize': 11,
    'legend.fontsize': 7, 'figure.dpi': 150,
    'savefig.dpi': 200, 'savefig.bbox': 'tight',
    'lines.linewidth': 1.3, 'axes.linewidth': 0.7,
    'grid.alpha': 0.3, 'grid.linestyle': ':',
})

# Default scenarios (relative paths from script location)
SCENARIOS = [
    ('S1 (NO2=0.1ppb)',  'vs_F0AM_chamber_S1_box.csv',   'vs_F0AM_chamber_S1_box.csv'),
    ('S2 (NO2=1ppb)',    'vs_F0AM_chamber_S2_box.csv',   'vs_F0AM_chamber_S2_box.csv'),
    ('S3 (NO2=10ppb)',   'vs_F0AM_chamber_S3_box.csv',   'vs_F0AM_chamber_S3_box.csv'),
    ('S2b (jcorr=10)',   'vs_F0AM_chamber_S2b_box.csv',  'vs_F0AM_chamber_S2b_box.csv'),
]


def load_csv(path):
    """Load a MOOSE or gold CSV; return {name: array} with species in ppb."""
    with open(path) as f:
        r = csv.reader(f)
        h = next(r)
        rows = list(r)
    data = {name: [] for name in h}
    for row in rows:
        for i, name in enumerate(h):
            data[name].append(float(row[i]))
    data = {k: np.array(v) for k, v in data.items()}
    # Convert molec/cm^3 to ppb for species columns (skip 'time')
    for name in h[1:]:
        data[name] = data[name] * 1.0e9 / M
    return data


def plot_panel(ax, gold, moose, sp, ylabel, title, scale=1.0):
    """Plot a single species panel on the given axes."""
    g = gold.get(sp, np.zeros(1)) * scale
    m = moose.get(sp, np.zeros(1)) * scale
    ax.plot(gold['time'] / 3600, g, '-', color='#1f77b4', label='F0AM gold', lw=1.2)
    ax.plot(moose['time'] / 3600, m, 'o--', color='#d62728',
            label='MOOSE', ms=3, lw=1.0)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='best', framealpha=0.8)
    ax.grid(True)


def plot_scenario(gold, moose, label, save_path):
    """Create a 2x2 figure comparing MOOSE vs F0AM for one scenario."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle(f'MOOSE vs F0AM — {label}', fontsize=11, fontweight='bold')

    panels = [
        (axes[0, 0], 'C5H8', 'C$_5$H$_8$ (ppb)', 'Isoprene', 1.0),
        (axes[0, 1], 'OH',   'OH (ppt)',          'OH Radical', 1.0e6),
        (axes[1, 0], None,   'NO$_x$ (ppb)',      'NO + NO$_2$', None),
        (axes[1, 1], 'O3',   'O$_3$ (ppb)',       'Ozone', 1.0),
    ]
    for ax, sp, yl, title, scale in panels:
        if sp is None:
            # NO + NO2 sum
            g = gold.get('NO', np.zeros(1)) + gold.get('NO2', np.zeros(1))
            m = moose.get('NO', np.zeros(1)) + moose.get('NO2', np.zeros(1))
            ax.plot(gold['time'] / 3600, g, '-', color='#1f77b4', label='F0AM gold', lw=1.2)
            ax.plot(moose['time'] / 3600, m, 'o--', color='#d62728',
                    label='MOOSE', ms=3, lw=1.0)
        else:
            plot_panel(ax, gold, moose, sp, yl, title, scale)
        ax.set_xlabel('Time (hours)')

    plt.tight_layout()
    plt.savefig(save_path)
    print(f'  Saved: {save_path}')
    plt.close(fig)


def check_files(actions_dir, gold_dir):
    """Check existence and size of all scenario files."""
    all_ok = True
    for label, gold_f, moose_f in SCENARIOS:
        for desc, d, f in [('gold', gold_dir, gold_f), ('moose', actions_dir, moose_f)]:
            path = os.path.join(d, f)
            if os.path.exists(path):
                sz = os.path.getsize(path)
                print(f'  [ok] {desc:>5}  {f:40s}  {sz//1024:4d} KB')
            else:
                print(f'  [!!] {desc:>5}  {f:40s}  MISSING')
                all_ok = False
    return all_ok


def main():
    parser = argparse.ArgumentParser(
        description='Plot MOOSE vs F0AM chamber comparison figures')
    parser.add_argument('--moose', help='MOOSE CSV path (single-scenario mode)')
    parser.add_argument('--gold', help='F0AM gold CSV path (single-scenario mode)')
    parser.add_argument('--save', help='Output PDF path')
    parser.add_argument('--check', action='store_true',
                        help='Only verify files exist, no plotting')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    actions_dir = os.path.join(script_dir, '../test/tests/actions')
    gold_dir = os.path.join(actions_dir, 'gold')

    # Single-scenario mode
    if args.moose and args.gold:
        if not os.path.exists(args.moose) or not os.path.exists(args.gold):
            print('ERROR: specified --moose or --gold file not found')
            sys.exit(1)
        gold = load_csv(args.gold)
        moose = load_csv(args.moose)
        save = args.save or 'chamber_comparison.pdf'
        plot_scenario(gold, moose, os.path.basename(args.gold), save)
        return

    # Check mode
    if args.check:
        ok = check_files(actions_dir, gold_dir)
        sys.exit(0 if ok else 1)

    # Default: plot all 4 scenarios
    print('Loading gold and MOOSE data...')
    for label, gold_f, moose_f in SCENARIOS:
        gpath = os.path.join(gold_dir, gold_f)
        mpath = os.path.join(actions_dir, moose_f)
        if not os.path.exists(gpath) or not os.path.exists(mpath):
            print(f'SKIP {label}: missing file (gold={os.path.exists(gpath)}, moose={os.path.exists(mpath)})')
            continue
        gold = load_csv(gpath)
        moose = load_csv(mpath)
        print(f'  {label}: gold {len(gold["time"])} tp, MOOSE {len(moose["time"])} tp')

        safe = label.replace(' ', '_').replace('(', '').replace(')', '').replace('=', '')
        save = args.save or os.path.join(actions_dir, f'chamber_{safe}.pdf')
        plot_scenario(gold, moose, label, save)

    print('Done.')


if __name__ == '__main__':
    main()
