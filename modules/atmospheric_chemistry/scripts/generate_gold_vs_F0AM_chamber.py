#!/usr/bin/env python3
"""
generate_gold_vs_F0AM_chamber_box.py -- Generate gold CSV files from F0AM .mat output.

Reads ChamberExampleOutput.mat (S1/S2/S3) and ChamberExampleHighLightsOutput.mat (S2b),
extracts species concentrations, removes F0AM-internal species (ONE, RO2, CH3ONO),
reorders columns to match MOOSE output, and writes to the MOOSE gold directory.

Dependencies: scipy, numpy (no MATLAB required)

Usage:
    # Default paths (for this repo structure)
    python3 scripts/generate_gold_vs_F0AM_chamber_box.py

    # Custom paths
    python3 scripts/generate_gold_vs_F0AM_chamber_box.py \\
        --mat /path/to/ChamberExampleOutput.mat \\
        --mat-s2b /path/to/ChamberExampleHighLightsOutput.mat \\
        --moose-csv /path/to/vs_F0AM_chamber_box.csv \\
        --gold-dir /path/to/gold
"""

import argparse, csv, os
import numpy as np
import scipy.io as sio

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def unwrap(v):
    """Unnest scipy-loaded MATLAB scalars/arrays."""
    while isinstance(v, np.ndarray) and v.shape == (1, 1):
        v = v.item()
    return v


def read_air_den(S_arr):
    """Read air number density Met.M (molec/cm³) from the S structure."""
    Met = S_arr['Met'].flatten()[0]
    M_field = Met['M']
    while hasattr(M_field, 'item') and hasattr(M_field, 'shape') and M_field.shape == (1, 1):
        M_field = M_field.item()
    return float(np.asarray(M_field).flatten()[0])


def load_conc(conc_struct, sp):
    """Extract a single species' concentration array from the Conc struct."""
    fd = conc_struct[0, 0][sp]
    while isinstance(fd, np.ndarray) and fd.shape == (1, 1):
        fd = fd.item()
    return np.asarray(fd).flatten().astype(float)


def extract_scenario(S_arr, step, target_times, air_den):
    """Extract species data at target times for a given step (or all if step=None)."""
    cnames = [str(c[0]) for c in unwrap(S_arr['Cnames']).flatten()]
    time_all = unwrap(S_arr['Time']).flatten()
    step_idx = unwrap(S_arr['StepIndex']).flatten().astype(int)
    conc_s = S_arr['Conc']

    mask = np.ones(len(time_all), dtype=bool) if step is None else (step_idx == step)
    t_filt = time_all[mask]

    result = np.zeros((len(target_times), len(cnames)))
    for ti, tt in enumerate(target_times):
        idx = np.argmin(np.abs(t_filt - tt))
        if idx >= len(t_filt):
            continue
        for si, sp in enumerate(cnames):
            vals = load_conc(conc_s, sp)
            result[ti, si] = vals[mask][idx] * air_den / 1.0e9
    return cnames, result


def clean_reorder(cnames_raw, conc_raw, moose_sp):
    """Drop F0AM-internal pseudo-species and reorder to MOOSE column order.
    ONE is a FACSIMILE placeholder (conc=1), CH3ONO is an external add-on
    not in the MOOSE mechanism. RO2 is kept as a derived diagnostic."""
    drop = {'ONE', 'CH3ONO'}
    keep = [i for i, sp in enumerate(cnames_raw) if sp not in drop]
    cnames_clean = [cnames_raw[i] for i in keep]
    conc_clean = conc_raw[:, keep]

    order = [-1] * len(moose_sp)
    for i, msp in enumerate(moose_sp):
        for j, fsp in enumerate(cnames_clean):
            if msp == fsp:
                order[i] = j
                break

    conc_out = np.zeros((conc_clean.shape[0], len(moose_sp)))
    for i, j in enumerate(order):
        conc_out[:, i] = conc_clean[:, j] if j >= 0 else 0.0
    return conc_out


def write_gold(path, header, target_times, conc):
    """Write gold CSV with the given header and concentration matrix."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for ti, t in enumerate(target_times):
            row = [str(int(t))]
            for si in range(len(header) - 1):
                row.append('{:.12e}'.format(conc[ti, si]))
            w.writerow(row)
    sz = os.path.getsize(path)
    print(f'  {os.path.basename(path)}: {len(header)-1} sp x {len(target_times)} tp ({sz//1024} KB)')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate gold CSV files from F0AM .mat output')
    parser.add_argument('--mat', default='/mnt/d/BaiduSyncdisk/Code/F0AM/Runs/ChamberExampleOutput.mat',
                        help='Path to ChamberExampleOutput.mat (default: F0AM Runs dir)')
    parser.add_argument('--mat-s2b',
                        default='/mnt/d/BaiduSyncdisk/Code/F0AM/Runs/ChamberExampleHighLightsOutput.mat',
                        help='Path to ChamberExampleHighLightsOutput.mat (default: F0AM Runs dir)')
    parser.add_argument('--moose-csv',
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             '../test/tests/actions/vs_F0AM_chamber_S1_box.csv'),
                        help='MOOSE CSV to use as column-order reference')
    parser.add_argument('--gold-dir',
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             '../test/tests/actions/gold'),
                        help='Output directory for gold CSVs')
    args = parser.parse_args()

    target_times = np.arange(1000, 11000, 1000)
    s2b_targets  = np.arange(1000, 4000, 1000)

    print('=== Generating F0AM gold CSVs ===')

    # 1. Load .mat files
    d1 = sio.loadmat(args.mat)
    d2 = sio.loadmat(args.mat_s2b)

    air_den = read_air_den(d1['S'])
    print(f'Air density: {air_den:.4e} molec/cm^3 (from Met.M)')

    # 2. Read MOOSE column order and ensure RO2 is included
    with open(args.moose_csv) as f:
        moose_header = next(csv.reader(f))
    moose_sp = moose_header[1:]  # species list (skip 'time')
    # RO2 is now always output as a ScalarVariable by AtmosphericChemistryAction box mode.
    # If the reference MOOSE CSV predates this change, append RO2 explicitly.
    if 'RO2' not in moose_sp:
        moose_sp = moose_sp + ['RO2']
        moose_header = ['time'] + moose_sp
        print('  Added RO2 column (MOOSE CSV predates RO2 output)')

    # 3. Export S1, S2, S3
    scenarios = [
        ('vs_F0AM_chamber_S1_box.csv',  1, target_times, 'S1 (NO2=0.1ppb)'),
        ('vs_F0AM_chamber_S2_box.csv',  2, target_times, 'S2 (NO2=1ppb)'),
        ('vs_F0AM_chamber_S3_box.csv',  3, target_times, 'S3 (NO2=10ppb)'),
    ]
    for fname, step, times, label in scenarios:
        cnames, conc = extract_scenario(d1['S'], step, times, air_den)
        conc_out = clean_reorder(cnames, conc, moose_header[1:])
        write_gold(os.path.join(args.gold_dir, fname), moose_header, times, conc_out)
        print(f'  {label}')

    # 4. Export S2b (restart from S2, jcorr=10)
    print('\nS2b (NO2=1ppb, jcorr=10, restart):')
    cnames_s2b, conc_s2b = extract_scenario(d2['S'], None, s2b_targets, air_den)
    conc_s2b_out = clean_reorder(cnames_s2b, conc_s2b, moose_header[1:])
    write_gold(os.path.join(args.gold_dir, 'vs_F0AM_chamber_S2b_box.csv'),
               moose_header, s2b_targets, conc_s2b_out)

    print(f'\nDone. Gold files in {args.gold_dir}')


if __name__ == '__main__':
    main()
