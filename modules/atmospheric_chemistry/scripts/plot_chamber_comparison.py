#!/usr/bin/env python3
"""
Compare MOOSE Chamber box results against F0AM ExampleSetup_Chamber reference.

Plots matching F0AM ExampleSetup_Chamber.m figures (S1 = low NO2 0.1 ppb):
  1. C5H8 (isoprene) time series
  2. OH radical time series (ppt)
  3. NO + NO2 time series
  4. O3 time series

Usage:
    python3 plot_chamber_comparison.py [--save PATH]
"""

import argparse, sys, os
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────
ACTIONS_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test/tests/actions")
F0AM_MAT = os.path.expanduser(
    "~/git_repo/moose/.reasonix/docs/BaiduSyncdisk/Code/F0AM/Runs/ChamberExampleOutput.mat")
MOOSE_CSV = os.path.join(ACTIONS_DIR, "vs_F0AM_chamber_box.csv")

M = 2.46e19  # molecules/cm³ at 298K, 1013mbar

# ── Plot style ─────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.labelsize": 10, "axes.titlesize": 11,
    "legend.fontsize": 7, "figure.dpi": 150,
    "savefig.dpi": 200, "savefig.bbox": "tight",
    "lines.linewidth": 1.3, "axes.linewidth": 0.7,
    "grid.alpha": 0.3, "grid.linestyle": ":",
})

F0AM_C = "#1f77b4"
MOOSE_C = "#d62728"


def load_f0am(mat_path):
    """Load F0AM Chamber S1 (low NO2)."""
    mat = sio.loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    S = mat["S"]
    t = np.asarray(S.Time, dtype=float).ravel()
    conc = {}
    for name in S.Conc._fieldnames:
        vals = np.asarray(getattr(S.Conc, name), dtype=float)
        if vals.ndim > 1 and vals.shape[1] > 1:
            vals = vals[:, 0]  # take first column if 2D
        conc[name] = vals.ravel()
    return t, conc


def load_moose(csv_path):
    """Load MOOSE CSV, return {name: array} with species in ppb."""
    import csv as c
    with open(csv_path) as f:
        r = c.reader(f)
        h = next(r)
        rows = list(r)
    data = {}
    for i, name in enumerate(h):
        data[name] = np.array([float(row[i]) for row in rows])
    # Convert molecules/cm³ → ppb for species columns
    ppb = {}
    for k, v in data.items():
        if k == "time":
            ppb[k] = v
        else:
            ppb[k] = v * 1e9 / M
    return ppb


def main():
    parser = argparse.ArgumentParser(
        description="Compare MOOSE Chamber vs F0AM reference")
    parser.add_argument("--save", help="Save figure (PDF)")
    args = parser.parse_args()

    if not os.path.exists(F0AM_MAT):
        print(f"F0AM output: {F0AM_MAT}")
        print("  -> NOT FOUND. Run ExampleSetup_Chamber.m first.")
        sys.exit(1)
    if not os.path.exists(MOOSE_CSV):
        print(f"MOOSE CSV: {MOOSE_CSV}")
        sys.exit(1)

    f0am_t, f0am = load_f0am(F0AM_MAT)
    moose = load_moose(MOOSE_CSV)

    print(f"F0AM: {len(f0am_t)} points, C5H8 range [{f0am['C5H8'][0]:.2f}, {f0am['C5H8'][-1]:.2f}] ppb")
    print(f"MOOSE: {len(moose['time'])} points, C5H8 range [{moose['C5H8'][0]:.2f}, {moose['C5H8'][-1]:.2f}] ppb")

    # ── Figure ──
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    def plot_pair(ax, sp, ylabel, title, scale=1.0, unit=""):
        f0am_y = f0am.get(sp, np.zeros(1)) * scale
        moose_y = moose.get(sp, np.zeros(1)) * scale
        ax.plot(f0am_t / 3600, f0am_y, "-", color=F0AM_C, label="F0AM", lw=1.2)
        ax.plot(moose["time"] / 3600, moose_y, "o--", color=MOOSE_C,
                label="MOOSE", ms=3, lw=1.0)
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc="best", framealpha=0.8)
        ax.grid(True)

    plot_pair(axes[0, 0], "C5H8", "C$_5$H$_8$ (ppb)", "Isoprene")
    plot_pair(axes[0, 1], "OH", "OH (ppt)", "OH Radical", scale=1e6)
    plot_pair(axes[1, 0], "NO", "NO$_x$ (ppb)", "NO + NO$_2$")
    # Add NO2 to NO
    if "NO2" in f0am and "NO2" in moose:
        f0am_y2 = f0am["NO"] + f0am["NO2"]
        moose_y2 = moose["NO"] + moose["NO2"]
        axes[1, 0].clear()
        axes[1, 0].plot(f0am_t / 3600, f0am_y2, "-", color=F0AM_C, label="F0AM", lw=1.2)
        axes[1, 0].plot(moose["time"] / 3600, moose_y2, "o--", color=MOOSE_C,
                        label="MOOSE", ms=3, lw=1.0)
        axes[1, 0].set_xlabel("Time (hours)")
        axes[1, 0].set_ylabel("NO + NO$_2$ (ppb)")
        axes[1, 0].set_title("NOx")
        axes[1, 0].legend(loc="best", framealpha=0.8)
        axes[1, 0].grid(True)

    plot_pair(axes[1, 1], "O3", "O$_3$ (ppb)", "Ozone")

    fig.suptitle("MOOSE vs F0AM — Chamber Box (610sp, NO$_2$=0.1 ppb, BOTTOMUP photolysis)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()

    out = args.save or os.path.join(ACTIONS_DIR, "chamber_comparison.pdf")
    plt.savefig(out)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
