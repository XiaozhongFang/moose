#!/usr/bin/env python3
"""
Generate gold CSV for vs_AtChem2_inorg_box from AtChem2 reference output.

Column order matches AtChem2's output files exactly:
  Photolysis: J1,J2,J3,J4,J5,J6,J7,J8,J11,J12,...,J61
  Env: M,TEMP,PRESS,RH,H2O,DEC,BLHEIGHT,DILUTE,JFAC,ROOF,ASA,RO2
  Species: CH3NO3,CH3O,...,SO3 (29 species in mechanism order)

Usage:
    python3 generate_atchem2_gold.py \
        --atchem2 .reasonix/docs/AtChem2/model/output \
        --output modules/atmospheric_chemistry/test/tests/actions/gold/vs_AtChem2_inorg_box.csv
"""

import argparse, sys
import numpy as np
from pathlib import Path


def load_atchem2(path):
    var = np.genfromtxt(path, max_rows=1, dtype=str)
    df = np.genfromtxt(path, skip_header=1)
    if df.ndim == 1:
        df = df.reshape(-1, 1)
    return var, df


def find_time_indices(t_ref, dt=900.0, end_time=43200.0):
    t_moose = np.arange(dt, end_time + dt / 2, dt)
    indices = [np.argmin(np.abs(t_ref - tm)) for tm in t_moose]
    return indices, t_moose


def generate_gold(atchem2_dir, output_path):
    atchem2 = Path(atchem2_dir)

    sp_var, sp_df = load_atchem2(atchem2 / "speciesConcentrations.output")
    ph_var, ph_df = load_atchem2(atchem2 / "photolysisRates.output")
    ev_var, ev_df = load_atchem2(atchem2 / "environmentVariables.output")

    t_ref = sp_df[:, 0]
    indices, t_moose = find_time_indices(t_ref)

    # Build lookup dicts
    sp_data = {v: sp_df[:, i] for i, v in enumerate(sp_var)}
    ph_data = {v: ph_df[:, i] for i, v in enumerate(ph_var)}
    ev_data = {v: ev_df[:, i] for i, v in enumerate(ev_var)}

    # ── Column order: AtChem2 photolysis + AtChem2 env + MOOSE species ──
    # Photolysis (AtChem2 order from photolysisRates.output)
    j_cols = [
        "J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8",
        "J11", "J12", "J13", "J14", "J15", "J16", "J17", "J18",
        "J19", "J20", "J21", "J22", "J23", "J24",
        "J31", "J32", "J33", "J34", "J35", "J41",
        "J51", "J52", "J53", "J54", "J55", "J56", "J61",
    ]

    # Environment (AtChem2 order; MOOSE uses block names → map here)
    # AtChem2: M, TEMP, PRESS, RH, H2O, DEC, BLHEIGHT, DILUTE, JFAC, ROOF, ASA, RO2
    # MOOSE block names: M_env, TEMP, PRESS, RH, H2O, DEC, BLHEIGHT, DILUTE, JFAC, ROOF, ASA, RO2_sum
    env_cols = [
        "M", "TEMP", "PRESS", "RH", "H2O", "DEC",
        "BLHEIGHT", "DILUTE", "JFAC", "ROOF", "ASA", "RO2",
    ]
    # AtChem2→MOOSE name mapping for env columns
    env_moose_names = [
        "M_env", "TEMP", "PRESS", "RH", "H2O", "DEC",
        "BLHEIGHT", "DILUTE", "JFAC", "ROOF", "ASA", "RO2_sum",
    ]

    # Species (MOOSE mechanism order)
    species_cols = [
        "CH3NO3", "CH3O", "CH3O2", "CH3O2NO2", "CH3OH", "CH3OOH",
        "CH4", "CL", "CO", "H2", "H2O2", "HCHO",
        "HNO3", "HO2", "HO2NO2", "HONO", "HSO3", "N2O5",
        "NA", "NO", "NO2", "NO3", "O", "O1D", "O3", "OH",
        "SA", "SO2", "SO3",
    ]

    # Solar params (MCMSolarPostprocessor columns).
    #
    # NOTE: AtChem2 computes these identically (Madronich 1993, zenith_data_mod /
    # solar_functions_mod) and PLOTS them (see plot-atchem2-numpy.py), but does
    # NOT output them to any text file.  The gold values below are computed from
    # the same Madronich 1993 formula.  This is NOT a true independent validation
    # (both codes use identical equations) — it only verifies the implementation
    # is bug-free, not that the formula itself is correct.
    solar_cols = ["cosx", "secx", "lha", "sinld", "cosld", "eqtime", "lat", "lon"]

    all_cols = j_cols + env_moose_names + solar_cols + species_cols
    header = ["time"] + all_cols

    # Pre-compute solar parameters for each time step (Madronich 1993)
    # lat=51.51, lon=0.13, day=21, month=6, year=2010
    pi = 3.14159265358979323846
    lat_deg, lon_deg = 51.51, 0.13
    day, month, year = 21, 6, 2010
    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (year % 4 == 0 and year % 100 != 0) or year % 400 == 0:
        days_in_months[1] = 29
    doy = sum(days_in_months[:month - 1]) + day
    days_in_year = 366 if days_in_months[1] == 29 else 365
    theta_day = 2.0 * pi * doy / days_in_year
    dec = (0.006918 - 0.399912 * np.cos(theta_day) + 0.070257 * np.sin(theta_day)
           - 0.006758 * np.cos(2 * theta_day) + 0.000907 * np.sin(2 * theta_day)
           - 0.002697 * np.cos(3 * theta_day) + 0.001480 * np.sin(3 * theta_day))
    eqt = (0.000075 + 0.001868 * np.cos(theta_day) - 0.032077 * np.sin(theta_day)
           - 0.014615 * np.cos(2 * theta_day) - 0.040849 * np.sin(2 * theta_day))
    lat_rad = lat_deg * pi / 180.0
    sinld = np.sin(lat_rad) * np.sin(dec)
    cosld = np.cos(lat_rad) * np.cos(dec)

    def solar_params(t_sec):
        frac_hour = (t_sec / 3600.0) % 24.0
        lha = pi * (frac_hour / 12.0 - (1.0 + lon_deg / 180.0)) + eqt
        cosx = np.cos(lha) * cosld + sinld
        secx = 1.0 / cosx if cosx > 1e-10 else 1e2
        if cosx <= 0.0:
            cosx = 0.0
            secx = 1e2
        return cosx, secx, lha, sinld, cosld, eqt

    # Build rows
    rows = []
    for step, aidx in enumerate(indices):
        row = [t_moose[step]]

        # Photolysis
        for col in j_cols:
            src = ph_data.get(col)
            row.append(src[aidx] if src is not None and aidx < len(src) else 0.0)

        # Environment
        for col, _ in zip(env_cols, env_moose_names):
            src = ev_data.get(col)
            row.append(src[aidx] if src is not None and aidx < len(src) else 0.0)

        # Solar params (computed from Madronich 1993 — same as MCMBoxModel)
        cosx, secx, lha, _sinld, _cosld, _eqt = solar_params(t_moose[step])
        row.extend([cosx, secx, lha, sinld, cosld, eqt, lat_deg, lon_deg])

        # Species
        for col in species_cols:
            src = sp_data.get(col)
            row.append(src[aidx] if src is not None and aidx < len(src) else 0.0)

        rows.append(row)

    with open(output_path, "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(f"{v:.15e}" for v in row) + "\n")

    n_ph = sum(1 for c in j_cols if c in ph_data)
    n_ev = sum(1 for c in env_cols if c in ev_data)
    n_sp = sum(1 for c in species_cols if c in sp_data)
    print(f"Gold: {len(t_moose)} steps × {len(all_cols)} columns")
    print(f"  J: {n_ph}/{len(j_cols)}  Env: {n_ev}/{len(env_cols)}  Sp: {n_sp}/{len(species_cols)}")


def main():
    parser = argparse.ArgumentParser(description="Generate AtChem2 gold CSV")
    parser.add_argument("--atchem2", required=True, help="AtChem2 output directory")
    parser.add_argument("--output", required=True, help="Output gold CSV path")
    args = parser.parse_args()
    atchem2 = Path(args.atchem2)
    if not atchem2.is_dir():
        print(f"Error: not a directory: {atchem2}")
        sys.exit(1)
    generate_gold(atchem2, args.output)


if __name__ == "__main__":
    main()
