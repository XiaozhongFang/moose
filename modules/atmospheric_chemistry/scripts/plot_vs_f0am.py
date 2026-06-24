#!/usr/bin/env python3
"""
F0AM-comparison plot for MOOSE atmospheric_chemistry gold CSVs.

Matches the layout of F0AM LearnF0AM_ODE.mlx figure_2:
  subplot 1 — concentration time series (A,C on left Y; B on right Y)
  subplot 2 — B concentration log-scale (with Analytical Solution comparison)
  subplot 3 — reaction rates
  subplot 4 — total mass

Usage:
    python3 plot_vs_f0am.py [--save PATH] [--no-show]
"""

import argparse, sys, csv
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent  # modules/atmospheric_chemistry
MOOSE_DIR = MODULE_DIR.parents[1]  # moose root
ACTIONS_DIR = MODULE_DIR / "test/tests/actions"
GOLD_DIR = ACTIONS_DIR / "gold"

# ── Publication style ─────────────────────────────────────────
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.fontsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "lines.linewidth": 1.3,
    "axes.linewidth": 0.7,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
})

# Reaction labels for tutorial_5sp
RXN_NAMES = [
    "A + B → C + B",
    "B → loss",
    "A → loss",
    "B + B → loss",
    "→ A",
    "C → loss",
]
# Rate constants
K = [1e-3, 1e-2, 1e-4, 1e-1, 5e-1, 1e-4]
# Reactant index pairs per reaction (species order: ONE, RO2, A, B, C)
IG = [(2,3), (3,), (2,), (3,3), (), (4,)]


def load_csv(path):
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [[float(v) for v in r] for r in reader if r]
    arr = np.array(rows)
    return arr[:, 0], header[1:], arr[:, 1:]


def compute_rates(t, y, cols):
    """Compute reaction rates from F0AM-style ODE system."""
    n_steps = len(t)
    n_rx = len(K)
    rates = np.zeros((n_steps, n_rx))
    C = np.zeros((n_steps, 5))
    for name, idx in [("ONE",0), ("RO2",1), ("A",2), ("B",3), ("C",4)]:
        if name in cols:
            j = cols.index(name)
            C[:, idx] = y[:, j]
        else:
            C[:, idx] = 1.0 if name == "ONE" else 0.0

    for r in range(n_rx):
        ig = IG[r]
        if len(ig) == 2:
            G = C[:, ig[0]] * C[:, ig[1]]
        elif len(ig) == 1:
            G = C[:, ig[0]]
        else:
            G = np.ones(n_steps)  # zero-order
        rates[:, r] = K[r] * G
    return rates


def plot_vs_f0am(path, save=None, show=False):
    t, cols, y = load_csv(path)

    # ── Build F0AM-compatible concentration arrays ──
    # Column names may be "A"/"B"/"C" (box mode) or "A_avg"/"B_avg"/"C_avg" (coupled mode).
    # Build a lookup that matches both conventions.
    col_lookup = {}
    for i, col in enumerate(cols):
        base = col.replace("_avg", "").replace("_val", "")
        col_lookup[base] = i

    c = {}
    for name in ["ONE", "RO2", "A", "B", "C"]:
        if name in col_lookup:
            c[name] = y[:, col_lookup[name]]
        else:
            c[name] = None

    # ── Figure ──
    fig = plt.figure(figsize=(10, 7.5), facecolor="white")
    c_A = "#1f77b4"   # blue
    c_C = "#d62728"   # red
    c_B = "#2ca02c"   # green

    # ── 计算物种 B 的解析解 ──
    if c["B"] is not None:
        k2 = K[1]  # 0.01
        k4 = K[3]  # 0.1
        b0 = c["B"][0]  # 自动提取初始浓度 [B]_0
        exp_k2t = np.exp(-k2 * t)
        b_analytical = (k2 * b0 * exp_k2t) / (k2 + 2.0 * k4 * b0 * (1.0 - exp_k2t))

    # ── Subplot 1: concentration time series (dual Y) ──
    ax1 = plt.subplot(2, 2, 1)
    ax1b = ax1.twinx()
    if c["A"] is not None and c["C"] is not None:
        ax1.plot(t, c["A"], color=c_A, lw=2.5, label="A (reactant)")
        ax1.plot(t, c["C"], color=c_C, lw=2.5, label="C (product)")
    if c["B"] is not None:
        ax1b.plot(t, c["B"], color=c_B, ls="-", lw=2, label="B (Simulation)")
        ax1b.plot(t, b_analytical, color="black", ls=":", lw=1.5, label="B (Analytical)")
        
    ax1.set_ylabel("Major species (A, C) [ppb]", fontweight="bold")
    ax1b.set_ylabel("Radical (B) [ppb]", fontweight="bold", color=c_B)
    ax1b.tick_params(axis="y", colors=c_B)
    ax1.set_xlabel("Time (s)", fontweight="bold")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", frameon=False)
    ax1.set_title("Species concentration evolution", fontweight="bold")
    ax1.grid(True)

    # ── Subplot 2: B log scale (含解析解对比) ──
    ax2 = plt.subplot(2, 2, 2)
    if c["B"] is not None:
        ax2.semilogy(t, c["B"], color=c_B, lw=2.5, label="Simulation")
        ax2.semilogy(t, b_analytical, color="black", ls=":", lw=1.5, label="Analytical")
        ax2.legend(frameon=False, loc="best")
    ax2.set_xlabel("Time (s)", fontweight="bold")
    ax2.set_ylabel("[B] (ppb)", fontweight="bold")
    ax2.set_title("Radical B (log scale)", fontweight="bold")
    ax2.grid(True)

    # ── Subplot 3: reaction rates ──
    ax3 = plt.subplot(2, 2, 3)
    rates = compute_rates(t, y, cols)
    rate_colors = ["#030303", "#2ca02c", "#d62728", "#9467bd", "#ff7f0e", "#8c564b"]
    for r in range(len(K)):
        ax3.plot(t, rates[:, r], color=rate_colors[r], lw=2, label=RXN_NAMES[r])
    ax3.set_xlabel("Time (s)", fontweight="bold")
    ax3.set_ylabel("Reaction rate (ppb/s)", fontweight="bold")
    ax3.legend(fontsize=7, frameon=False, loc="best")
    ax3.set_title("Reaction rates", fontweight="bold")
    ax3.grid(True)

    # ── Subplot 4: total mass ──
    ax4 = plt.subplot(2, 2, 4)
    if c["A"] is not None and c["B"] is not None and c["C"] is not None:
        total = c["A"] + c["B"] + c["C"]
        ax4.plot(t, total, color="#444444", lw=2.5)
    ax4.set_xlabel("Time (s)", fontweight="bold")
    ax4.set_ylabel("Σ(A + B + C) (ppb)", fontweight="bold")
    ax4.set_title("Total mass (source + sink)", fontweight="bold")
    ax4.grid(True)

    fig.suptitle("F0AM / MOOSE BDF2 — kinetic simulation results",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()

    if save:
        fig.savefig(save)
        print(f"Saved: {save}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="F0AM figure_2 style plot")
    parser.add_argument("--file", "-f",
                        default=str(ACTIONS_DIR / "vs_F0AM_tutorial5_box.csv"),
                        help="CSV file to plot")
    parser.add_argument("--save", "-s", help="Save figure (.pdf/.png)")
    parser.add_argument("--no-show", "-n", action="store_true", help="Don't display")
    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        path = GOLD_DIR / args.file
    if not path.exists():
        print(f"File not found: {args.file}")
        sys.exit(1)

    plot_vs_f0am(path, save=args.save, show=not args.no_show)


if __name__ == "__main__":
    main()