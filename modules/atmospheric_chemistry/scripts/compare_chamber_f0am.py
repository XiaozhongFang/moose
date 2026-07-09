#!/usr/bin/env python3
"""
compare_chamber_f0am.py -- Compare MOOSE CSV output against F0AM gold CSV.

Computes per-species ratios at each common time point, categorizes results
(good/ok/low/high/moose_zero), and prints a summary report.

Usage:
    # Compare checked-in chamber scenarios using default paths
    python3 scripts/compare_chamber_f0am.py

    # Compare a single scenario
    python3 scripts/compare_chamber_f0am.py \\
        --moose /path/to/vs_F0AM_chamber_box.csv \\
        --gold /path/to/gold/vs_F0AM_chamber_box.csv
"""

import argparse, csv, os, sys
from collections import OrderedDict


def load_csv(path):
    """Load a CSV file; return (header, {time: {species: value}})."""
    with open(path) as f:
        r = csv.reader(f)
        header = next(r)
        data = {}
        for row in r:
            t = float(row[0])
            data[t] = {header[i].strip(): float(row[i]) for i in range(1, len(header))}
    return header, data


def compare(gold_path, moose_path, rel_err=0.5, abs_zero=1e-20):
    """Compare MOOSE vs gold, print categorized summary. Return (n_total, n_bad)."""
    _, gold = load_csv(gold_path)
    _, moose = load_csv(moose_path)

    common_t = sorted(set(gold.keys()) & set(moose.keys()))
    if not common_t:
        print('ERROR: No common time points between gold and MOOSE CSV')
        return 0, 0

    common_sp = sorted(set(gold[common_t[0]].keys()) & set(moose[common_t[0]].keys()))
    tf = common_t[-1]

    results = OrderedDict()
    for sp in common_sp:
        vg = gold[tf][sp]
        vm = moose[tf][sp]

        if abs(vg) < 1e-30 and abs(vm) < 1e-30:
            results[sp] = ('both_zero', 1.0)
        elif abs(vm) < abs_zero and abs(vg) > 1e-10:
            results[sp] = ('moose_zero', 0.0)
        elif abs(vg) < 1e-30:
            results[sp] = ('gold_zero', float('inf'))
        else:
            ratio = vm / vg
            if 0.8 <= ratio <= 1.2:
                results[sp] = ('good', ratio)
            elif 0.5 <= ratio <= 2.0:
                results[sp] = ('ok', ratio)
            elif ratio < 0.5:
                results[sp] = ('low', ratio)
            else:
                results[sp] = ('high', ratio)

    # Summary
    print(f'{"="*70}')
    print(f'F0AM vs MOOSE at t={tf:.0f}s  ({os.path.basename(gold_path)})')
    print(f'{"="*70}')
    for cat, label in [('moose_zero', 'ZERO in MOOSE (MOOSE=0, F0AM>0)'),
                       ('low', 'LOW (ratio < 0.5)'),
                       ('high', 'HIGH (ratio > 2.0)')]:
        bad = [(sp, r) for sp, (c, r) in results.items() if c == cat]
        if bad:
            bad.sort(key=lambda x: x[1])
            print(f'\n--- {label}: {len(bad)} species ---')
            for sp, ratio in bad[:30]:
                vg = gold[tf].get(sp, 0)
                vm = moose[tf].get(sp, 0)
                print(f'  {sp:>20}: F0AM={vg:12.4e}  MOOSE={vm:12.4e}  ratio={ratio:8.4f}')

    counts = {}
    for c, _ in results.values():
        counts[c] = counts.get(c, 0) + 1

    print(f'\n{"="*70}')
    print(f'SUMMARY: {len(common_sp)} species, {len(common_t)} time points')
    for cat in ['good', 'ok', 'low', 'high', 'moose_zero', 'both_zero', 'gold_zero']:
        if cat in counts:
            print(f'  {cat:15}: {counts[cat]:4d}')

    n_bad = counts.get('low', 0) + counts.get('high', 0) + counts.get('moose_zero', 0)
    return len(common_sp), n_bad


def main():
    parser = argparse.ArgumentParser(
        description='Compare MOOSE CSV output against F0AM gold CSV')
    parser.add_argument('--moose', help='MOOSE output CSV path')
    parser.add_argument('--gold', help='F0AM gold CSV path')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Relative error threshold (default: 0.5 = 50 percent)')
    args = parser.parse_args()

    if args.moose and args.gold:
        compare(args.gold, args.moose)
        return

    # Default: compare chamber scenarios with checked-in gold.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    chamber_dir = os.path.join(script_dir, '../test/tests/chamber')
    gold_dir = os.path.join(chamber_dir, 'gold')

    scenarios = [
        ('vs_F0AM_chamber_S1_box.csv',   'vs_F0AM_chamber_S1_box.csv'),
        ('vs_F0AM_chamber_S2_box.csv',   'vs_F0AM_chamber_S2_box.csv'),
        ('vs_F0AM_chamber_S2b_box.csv',  'vs_F0AM_chamber_S2b_box.csv'),
        ('vs_F0AM_chamber_S3_box.csv',   'vs_F0AM_chamber_S3_box.csv'),
    ]
    total_bad = 0
    for gold_f, moose_f in scenarios:
        gpath = os.path.join(gold_dir, gold_f)
        mpath = os.path.join(chamber_dir, moose_f)
        if not os.path.exists(gpath) or not os.path.exists(mpath):
            print(f'SKIP {moose_f}: missing file')
            continue
        n, nb = compare(gpath, mpath)
        total_bad += nb
        print()

    print(f'Total bad species across all scenarios: {total_bad}')


if __name__ == '__main__':
    sys.exit(main())
