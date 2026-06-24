#!/usr/bin/env python3
"""
Plot MOOSE vs AtChem2 results — matching plot-atchem2-numpy.py format.

- 3×3 grid per page, each variable its own subplot
- Time in seconds, scientific notation on Y axis
- MOOSE data (black solid line) overlaid with AtChem2 reference (red dashed)
- Column order matches AtChem2 output files:
  Species: CH3NO3 ... SO3 (29 species)
  Photolysis: J1, J2, J3, ... (AtChem2 order)
  Environment: M, TEMP, PRESS, RH, H2O, DEC, BLHEIGHT, DILUTE, JFAC, ROOF, ASA, RO2

Usage:
    python3 plot_atchem2.py [--moose CSV] [--atchem2 DIR] [--output PDF]
"""

import argparse, csv, sys
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent  # modules/atmospheric_chemistry
MOOSE_DIR = MODULE_DIR.parents[1]  # moose root
ACTIONS_DIR = MODULE_DIR / "test/tests/actions"
ATCHEM2_DIR = MOOSE_DIR / ".reasonix/docs/AtChem2/model/output"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight",
    "lines.linewidth": 1.0, "axes.linewidth": 0.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7,
})

# AtChem2 column order (must match generate_atchem2_gold.py)
J_ORDER = [
    "J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8",
    "J11", "J12", "J13", "J14", "J15", "J16", "J17", "J18",
    "J19", "J20", "J21", "J22", "J23", "J24",
    "J31", "J32", "J33", "J34", "J35", "J41",
    "J51", "J52", "J53", "J54", "J55", "J56", "J61",
]
ENV_ATCHEM2_ORDER = [
    "M", "TEMP", "PRESS", "RH", "H2O", "DEC",
    "BLHEIGHT", "DILUTE", "JFAC", "ROOF", "ASA", "RO2",
]
# MOOSE env block names → AtChem2 names
ENV_MOOSE_TO_ATCHEM2 = {
    "M_env": "M", "TEMP": "TEMP", "PRESS": "PRESS", "RH": "RH",
    "H2O": "H2O", "DEC": "DEC", "BLHEIGHT": "BLHEIGHT",
    "DILUTE": "DILUTE", "JFAC": "JFAC", "ROOF": "ROOF",
    "ASA": "ASA", "RO2_sum": "RO2",
}


def load_moose_csv(path):
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [[float(v) for v in r] for r in reader if r]
    arr = np.array(rows)
    return arr[:, 0], header, arr[:, 1:]


def load_atchem2(path):
    var = np.genfromtxt(path, max_rows=1, dtype=str)
    df = np.genfromtxt(path, skip_header=1)
    if df.ndim == 1:
        df = df.reshape(-1, 1)
    return var, df


def plot_grid(t_moose, moose_data, moose_cols, atchem_var, atchem_df,
              title_prefix, pdf):
    """Plot 3×3 grid with MOOSE (black) vs AtChem2 (red dashed) overlay."""
    nc = len(moose_cols)
    if nc == 0:
        return
    # Build AtChem2 lookup: name → data column
    atchem_lookup = {}
    if atchem_var is not None and atchem_df is not None:
        for i, v in enumerate(atchem_var):
            atchem_lookup[v] = atchem_df[:, i]

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(11, 7))
    axs = axs.ravel()
    j = 0
    for i in range(nc):
        ax = axs[j]
        name = moose_cols[i]
        ax.plot(t_moose, moose_data[:, i], linestyle="-", color="black",
                linewidth=0.8, label="MOOSE")
        if name in atchem_lookup:
            ax.plot(atchem_df[:, 0], atchem_lookup[name],
                    linestyle="--", color="red", linewidth=0.6,
                    label="AtChem2")
        ax.legend(fontsize=6, frameon=False, loc="upper right")
        ax.set(title=name, xlabel="seconds", ylabel="")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "%.1e" % x))
        if j == 8:
            fig.tight_layout()
            pdf.savefig(fig)
            fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(11, 7))
            axs = axs.ravel()
            j = 0
        else:
            j = j + 1
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot MOOSE vs AtChem2")
    parser.add_argument("--moose", "-m",
                        default=str(ACTIONS_DIR / "vs_AtChem2_inorg_box.csv"))
    parser.add_argument("--atchem2", "-a", default=str(ATCHEM2_DIR))
    parser.add_argument("--output", "-o", default="moose_vs_atchem2.pdf")
    args = parser.parse_args()

    t, cols, data = load_moose_csv(Path(args.moose))
    col_name = cols[1:]  # skip "time"

    # Build column name → index lookup
    col_idx = {c: i for i, c in enumerate(col_name)}

    # Load AtChem2 reference
    atchem2 = Path(args.atchem2)
    sp_var, sp_df = None, None
    ph_var, ph_df = None, None
    ev_var, ev_df = None, None
    if atchem2.is_dir():
        if (atchem2 / "speciesConcentrations.output").exists():
            sp_var, sp_df = load_atchem2(atchem2 / "speciesConcentrations.output")
        if (atchem2 / "photolysisRates.output").exists():
            ph_var, ph_df = load_atchem2(atchem2 / "photolysisRates.output")
        if (atchem2 / "environmentVariables.output").exists():
            ev_var, ev_df = load_atchem2(atchem2 / "environmentVariables.output")

    # Species columns (MOOSE mechanism order)
    species_cols = [c for c in col_name
                    if c not in J_ORDER
                    and c not in ENV_MOOSE_TO_ATCHEM2
                    and c not in ENV_ATCHEM2_ORDER]
    sp_indices = [col_idx[c] for c in species_cols if c in col_idx]
    sp_names = [c for c in species_cols if c in col_idx]

    with PdfPages(args.output) as pdf:
        # ── Page 1: Species ──
        if sp_names:
            sp_data = data[:, sp_indices]
            plot_grid(t, sp_data, sp_names, sp_var, sp_df, "Species", pdf)

        # ── Page 2: Photolysis (AtChem2 order) ──
        j_names = [j for j in J_ORDER if j in col_idx]
        if j_names:
            j_indices = [col_idx[j] for j in j_names]
            j_data = data[:, j_indices]
            plot_grid(t, j_data, j_names, ph_var, ph_df, "Photolysis", pdf)

        # ── Page 3: Environment (AtChem2 order) ──
        env_display = []
        env_indices = []
        for atchem_name in ENV_ATCHEM2_ORDER:
            # Find MOOSE column that maps to this AtChem2 name
            moose_name = None
            for mk, mv in ENV_MOOSE_TO_ATCHEM2.items():
                if mv == atchem_name and mk in col_idx:
                    moose_name = mk
                    break
            if moose_name:
                env_display.append(atchem_name)
                env_indices.append(col_idx[moose_name])
        if env_display:
            env_data = data[:, env_indices]
            plot_grid(t, env_data, env_display, ev_var, ev_df, "Environment", pdf)

    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
