#!/usr/bin/env python3
"""
Plot MOOSE vs AtChem2 atmospheric chemistry results.

Matches the layout of AtChem2's plot-atchem2-numpy.py:
  - 3×3 grid per page, each variable its own subplot
  - Time in seconds on X axis, scientific notation on Y axis
  - MOOSE data (black solid) overlaid with AtChem2 reference (red dashed)

Output: multi-page PDF with 4 sections:
  Page 1 — Species concentrations (29 species)
  Page 2 — Photolysis rates (35 J values: J1-J8, J11-J24, J31-J35, J41, J51-J56, J61)
  Page 3 — Environment variables (M, TEMP, PRESS, RH, H2O, DEC, BLHEIGHT, DILUTE, JFAC, ROOF, ASA, RO2)
  Page 4 — Solar parameters (cosx, secx, lha, sinld, cosld, eqtime, lat, lon)

Arguments:
  --moose, -m    Path to MOOSE CSV output (vs_AtChem2_inorg_box.csv).
                 Default: modules/atmospheric_chemistry/test/tests/actions/vs_AtChem2_inorg_box.csv
  --atchem2, -a  Directory containing AtChem2 output files:
                   speciesConcentrations.output
                   photolysisRates.output
                   environmentVariables.output
                 Also tries <dir>/model/output/ subdirectory.
                 If files not found, plots MOOSE data only (no overlay).
                 Default: .reasonix/docs/AtChem2/model/output/
  --output, -o   Output PDF path.  Default: moose_vs_atchem2.pdf

Usage:
    python3 plot_atchem2.py
    python3 plot_atchem2.py --moose my_run.csv --output my_plot.pdf
    python3 plot_atchem2.py --atchem2 ~/AtChem2/model/output
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
        rows = []
        for lineno, r in enumerate(reader, start=2):
            if not r:
                continue
            try:
                rows.append([float(v) for v in r])
            except ValueError as exc:
                raise ValueError(
                    f"Non-numeric data in {path} at line {lineno}: {r}") from exc
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
    step = max(1, len(t_moose) // 8)
    for i in range(nc):
        ax = axs[j]
        name = moose_cols[i]

        # MOOSE: black solid line + hollow ○ markers
        ax.plot(t_moose, moose_data[:, i], linestyle="-", color="black",
                linewidth=1.0, label="MOOSE")
        ax.plot(t_moose[::step], moose_data[::step, i],
                linestyle="none", marker="o", color="black",
                markersize=3.5, markerfacecolor="none", markeredgewidth=0.8)

        # AtChem2: red dashed line + hollow □ markers
        if name in atchem_lookup:
            atchem_t = atchem_df[:, 0]
            atchem_y = atchem_lookup[name]
            ax.plot(atchem_t, atchem_y, linestyle="--", color="#d62728",
                    linewidth=0.5, alpha=0.8, label="AtChem2")
            ax.plot(atchem_t[::step], atchem_y[::step],
                    linestyle="none", marker="s", color="#d62728",
                    markersize=3.5, markerfacecolor="none", markeredgewidth=0.8)

        ax.legend(fontsize=6, frameon=True, loc="best",
                  fancybox=False, edgecolor="gray", facecolor="white",
                  markerscale=0.8)
        ax.set(title=name, xlabel="seconds", ylabel="")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "%.1e" % x))
        j = j + 1
        if j == 9:
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(11, 7))
            axs = axs.ravel()
            j = 0
    if j > 0:
        # Hide unused subplots on the last page
        for k in range(j, 9):
            axs[k].set_visible(False)
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

    # Load AtChem2 reference — try the given directory, and if files are not
    # found there, look in the standard model output location.
    atchem2 = Path(args.atchem2)
    sp_var, sp_df = None, None
    ph_var, ph_df = None, None
    ev_var, ev_df = None, None

    def _try_load(base_dir):
        """Try loading AtChem2 output files from base_dir; return True if any found."""
        nonlocal sp_var, sp_df, ph_var, ph_df, ev_var, ev_df
        sp_file = base_dir / "speciesConcentrations.output"
        ph_file = base_dir / "photolysisRates.output"
        ev_file = base_dir / "environmentVariables.output"
        found = False
        if sp_file.exists():
            sp_var, sp_df = load_atchem2(str(sp_file))
            found = True
        if ph_file.exists():
            ph_var, ph_df = load_atchem2(str(ph_file))
            found = True
        if ev_file.exists():
            ev_var, ev_df = load_atchem2(str(ev_file))
            found = True
        return found

    if atchem2.is_dir():
        if not _try_load(atchem2):
            # Also try model/output/ subdirectory (common AtChem2 layout)
            alt = atchem2 / "model" / "output"
            if alt.is_dir():
                _try_load(alt)

    if sp_var is None and ph_var is None and ev_var is None:
        print("Note: AtChem2 reference data not found — plotting MOOSE data only.\n"
              "      Expected files: speciesConcentrations.output, photolysisRates.output,\n"
              "      environmentVariables.output in a model/output/ directory.\n"
              "      Run with --atchem2 <path-to-model-output> to add AtChem2 overlay.",
              file=sys.stderr)

    # Solar columns (separate from species)
    SOLAR_ORDER = ["cosx", "secx", "lha", "sinld", "cosld", "eqtime", "lat", "lon"]

    # Species columns (MOOSE mechanism order)
    species_cols = [c for c in col_name
                    if c not in J_ORDER
                    and c not in ENV_MOOSE_TO_ATCHEM2
                    and c not in ENV_ATCHEM2_ORDER
                    and c not in SOLAR_ORDER]
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
        # Reverse mapping: AtChem2 name → MOOSE column name (built once)
        atchem2_to_moose = {v: k for k, v in ENV_MOOSE_TO_ATCHEM2.items()}
        env_display = []
        env_indices = []
        for atchem_name in ENV_ATCHEM2_ORDER:
            moose_name = atchem2_to_moose.get(atchem_name)
            if moose_name and moose_name in col_idx:
                env_display.append(atchem_name)
                env_indices.append(col_idx[moose_name])
        if env_display:
            env_data = data[:, env_indices]
            plot_grid(t, env_data, env_display, ev_var, ev_df, "Environment", pdf)

        # ── Page 4: Solar parameters ──
        sol_names = [s for s in SOLAR_ORDER if s in col_idx]
        if sol_names:
            sol_indices = [col_idx[s] for s in sol_names]
            sol_data = data[:, sol_indices]
            plot_grid(t, sol_data, sol_names, None, None, "Solar", pdf)

    n_types = (1 if sp_names else 0) + (1 if j_names else 0) + \
              (1 if env_display else 0) + (1 if sol_names else 0)
    print(f"Saved: {args.output}  ({len(sp_names)} species, {len(j_names)} photolysis, "
          f"{len(env_display)} env, {len(sol_names)} solar → {n_types} pages)")


if __name__ == "__main__":
    main()
