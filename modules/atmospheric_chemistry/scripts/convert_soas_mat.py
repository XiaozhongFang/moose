#!/usr/bin/env python3
"""
convert_soas_mat.py — Convert F0AM SOAS observation .mat files to
MOOSE-readable CSV format.

Usage:
    python3 convert_soas_mat.py <mat_file> [--output <csv_file>]

Arguments:
    mat_file        Path to SOAS observation .mat file (required)
    --output, -o    Output CSV file path (default: SOAS_DielCycle.csv in
                    the same directory as the input file)
    --hour, -H      Output only the initial condition for the given hour
                    (default: all 24 hours)

Examples:
    # Convert all 24 hours of data
    python3 convert_soas_mat.py Obs_SOAS_CampaignAvg_60min.mat

    # Specify output path
    python3 convert_soas_mat.py Obs_SOAS_CampaignAvg_60min.mat -o ../database/SOAS_DielCycle.csv

    # Extract only hour 0 initial condition
    python3 convert_soas_mat.py Obs_SOAS_CampaignAvg_60min.mat --hour 0

Output CSV columns:
    Time_h, T(K), P(mbar), RH(%), SZA(deg), M(molec/cm3), BLheight(m),
    NO, NO2, O3, OH, CO, H2O2, PAN, C2H5NO3, IC3H7NO3,
    C5H8, APINENE, BPINENE, LIMONENE, C2H4, C2H6, C3H8, ...
    H2 (constant 550 ppb), CH4 (constant 1770 ppb)

Dependencies: scipy, numpy
"""

import scipy.io as sio
import numpy as np
import argparse
import os
import sys

# Observed species fields in the SOAS struct (units: ppb)
SPECIES_FIELDS = [
    'NO', 'NO2', 'O3', 'OH', 'CO', 'H2O2', 'PAN', 'C2H5NO3', 'IC3H7NO3',
    'C5H8', 'APINENE', 'BPINENE', 'LIMONENE',
    'C2H4', 'C2H6', 'C3H8', 'IC4H10', 'IC5H12', 'NC5H12', 'NC6H14', 'NC10H22',
    'BENZENE', 'TOLUENE', 'EBENZ', 'TM124B', 'TM135B', 'MXYL', 'OXYL', 'PXYL', 'BENZAL',
    'CH3CHO', 'C2H5CHO', 'C3H7CHO', 'HOCH2CHO', 'GLYOX', 'CH3OH', 'C2H5OH',
    'ACETOL', 'BIACET', 'MACR', 'MVK', 'HCHO', 'CH3COCH3',
    'C2H2', 'C3H6', 'NC4H10', 'DMS', 'HNO3', 'HO2', 'IEPOX', 'ISOPOOH', 'MEK', 'MPAN',
]

# Environmental variable fields (raw units)
ENV_FIELDS = ['T', 'P', 'RH', 'SZA', 'M', 'BLheight']

# Constant species (not in SOAS observations, set manually in F0AM scripts)
CONSTANT_SPECIES = {'H2': 550.0, 'CH4': 1770.0}  # ppb


def get_field_1d(soas, name):
    """Extract a 1-D numpy array from the SOAS struct."""
    val = soas[name]
    if val.ndim == 0:
        val = val[()]
    return np.asarray(val, dtype=float).ravel()


def main():
    parser = argparse.ArgumentParser(
        description='Convert F0AM SOAS observation .mat files to MOOSE CSV format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s Obs_SOAS_CampaignAvg_60min.mat
  %(prog)s Obs_SOAS_CampaignAvg_60min.mat -o SOAS_DielCycle.csv
  %(prog)s Obs_SOAS_CampaignAvg_60min.mat --hour 0
        """
    )
    parser.add_argument('mat_file', help='Path to SOAS observation .mat file')
    parser.add_argument('--output', '-o', default=None,
                        help='Output CSV file path (default: SOAS_DielCycle.csv in input dir)')
    parser.add_argument('--hour', '-H', type=int, default=None,
                        help='Output only the initial condition for the given hour')
    args = parser.parse_args()

    # Check input file
    if not os.path.isfile(args.mat_file):
        print(f"Error: file not found '{args.mat_file}'", file=sys.stderr)
        sys.exit(1)

    # Default output path
    out_path = args.output
    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(args.mat_file)) or '.'
        out_path = os.path.join(out_dir, 'SOAS_DielCycle.csv')

    # Read .mat file
    mat = sio.loadmat(args.mat_file, squeeze_me=True)
    if 'SOAS' not in mat:
        print("Error: 'SOAS' struct not found in .mat file", file=sys.stderr)
        sys.exit(1)
    soas = mat['SOAS']

    time_arr = get_field_1d(soas, 'Time')
    n = len(time_arr)
    print(f"SOAS: {n} time points, {len(soas.dtype.names)} fields")

    # Determine output row range
    if args.hour is not None:
        hour_indices = [args.hour]
    else:
        hour_indices = range(n)

    # Collect column names
    all_cols = ['Time_h'] + ENV_FIELDS + SPECIES_FIELDS + list(CONSTANT_SPECIES.keys())

    # Build data rows
    rows = []
    for i in hour_indices:
        row = {'Time_h': int(time_arr[i])}
        for f in ENV_FIELDS:
            vals = get_field_1d(soas, f)
            row[f] = float(vals[i]) if i < len(vals) else 0.0
        for f in SPECIES_FIELDS:
            if f in soas.dtype.names:
                vals = get_field_1d(soas, f)
                row[f] = float(vals[i]) if i < len(vals) else 0.0
            else:
                row[f] = 0.0
        for name, val in CONSTANT_SPECIES.items():
            row[name] = val
        rows.append(row)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)

    # Write CSV
    with open(out_path, 'w') as f:
        f.write(','.join(all_cols) + '\n')
        for row in rows:
            f.write(','.join(str(row[c]) for c in all_cols) + '\n')

    print(f"Written: {len(rows)} rows x {len(all_cols)} cols -> {out_path}")

    # Print initial condition summary
    if rows:
        print(f"\nInitial conditions (hour {rows[0]['Time_h']}):")
        for f in ENV_FIELDS:
            print(f"  {f}: {rows[0][f]:.6g}")
        print(f"  Constant H2: {CONSTANT_SPECIES['H2']} ppb")
        print(f"  Constant CH4: {CONSTANT_SPECIES['CH4']} ppb")
        for f in SPECIES_FIELDS[:8]:
            if f in rows[0]:
                print(f"  {f}: {rows[0][f]:.6g} ppb")
        print("  ...")


if __name__ == '__main__':
    main()
