#!/usr/bin/env python3
"""Plot F0AM chamber comparison figures for test and benchmark outputs.

The default mode reproduces the checked chamber CSVDiff figure set for S1, S2,
S2b, and S3.  The benchmark mode reads CSVs produced by
benchmark_chamber_solvers.py and writes the same figure layout for each
solver/scenario pair, so ordinary chamber tests and KPP/PETSc benchmark runs
use one plotting path.

Examples
--------
Plot checked-in chamber outputs:

    python3 scripts/plot_chamber_comparison.py

Check the checked-in chamber outputs without plotting:

    python3 scripts/plot_chamber_comparison.py --check

Plot one CSV against one F0AM gold CSV:

    python3 scripts/plot_chamber_comparison.py \
        --moose test/tests/chamber/vs_F0AM_chamber_S1_box.csv \
        --gold test/tests/chamber/gold/vs_F0AM_chamber_S1_box.csv \
        --save /tmp/chamber_S1.png

Plot all solver CSVs in a benchmark run directory:

    python3 scripts/plot_chamber_comparison.py \
        --run-dir test/tests/chamber/kpp_chamber/solver_runs/f0am_accuracy_timing \
        --solvers kpp_rosenbrock,petsc_bdf \
        --output-dir test/tests/chamber/kpp_chamber/solver_runs/f0am_accuracy_timing/figures

Write early-time atmospheric-chemistry diagnostics:

    python3 scripts/plot_chamber_comparison.py \
        --diagnostics-csv test/tests/chamber/chamber_diagnostics.csv
"""

import argparse
import csv
import os
import re
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
CHAMBER_DIR = MODULE_DIR / "test" / "tests" / "chamber"
GOLD_DIR = CHAMBER_DIR / "gold"
MECHANISM_FILE = MODULE_DIR / "doc" / "content" / "modules" / "atmospheric_chemistry" / "database" / "MCMv331_Inorg_Isoprene.fac"

# Molecules/cm^3 -> ppb conversion at the chamber test conditions.
AIR_DENSITY = 2.462203318410e19
EARLY_WINDOW_SECONDS = 1010.0

SCENARIOS = [
    ("S1", "S1 (NO2=0.1ppb)", "vs_F0AM_chamber_S1_box.csv"),
    ("S2", "S2 (NO2=1ppb)", "vs_F0AM_chamber_S2_box.csv"),
    ("S2b", "S2b restart (NO2=1ppb, jcorr=10)", "vs_F0AM_chamber_S2b_box.csv"),
    ("S3", "S3 (NO2=10ppb)", "vs_F0AM_chamber_S3_box.csv"),
]
SCENARIO_BY_NAME = {name: (label, gold) for name, label, gold in SCENARIOS}

# S2b CSVDiff uses restart-relative seconds; figures place it after S2.
RESTART_TIME_OFFSETS = {"S2b": 10800.0}
DIAGNOSTIC_COLUMNS = ("OH", "RO2", "HO2", "CH3O2", "C5H8", "NO", "NO2", "O3",
                      "MVK", "MACR", "HCHO")

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "lines.linewidth": 1.3,
    "axes.linewidth": 0.7,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
})


def resolve_path(path):
    """Resolve a user path relative to cwd, module root, then repo-style test dir."""
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    for base in (Path.cwd(), MODULE_DIR, CHAMBER_DIR):
        resolved = (base / candidate).resolve()
        if resolved.exists():
            return resolved
    return (Path.cwd() / candidate).resolve()


def parse_solvers(value):
    return [item for item in value.split(",") if item.strip()]


def read_ro2_species(path=MECHANISM_FILE):
    if not path.exists():
        return []
    text = path.read_text()
    match = re.search(r"\bRO2\s*=\s*(.*?)\s*;", text, flags=re.DOTALL)
    if not match:
        return []
    names = []
    for item in match.group(1).replace("\n", " ").split("+"):
        name = re.sub(r"[^A-Za-z0-9_]", "", item.strip())
        if name:
            names.append(name)
    return names


RO2_SPECIES = read_ro2_species()


def load_csv(path):
    """Load a MOOSE or F0AM gold CSV; return {column: numpy array}.

    Species concentrations are converted from molec/cm^3 to ppb; the time
    column remains in seconds.
    """
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = [[float(value) for value in row] for row in reader if row]
    if not rows:
        raise ValueError(f"CSV has no data rows: {path}")

    arr = np.array(rows)
    data = {}
    for index, name in enumerate(header):
        values = arr[:, index]
        if index > 0:
            values = values * 1.0e9 / AIR_DENSITY
        data[name] = values
    if "RO2" not in data and RO2_SPECIES:
        ro2_values = np.zeros_like(data["time"])
        found = False
        for name in RO2_SPECIES:
            if name in data:
                ro2_values += data[name]
                found = True
        if found:
            data["RO2"] = ro2_values
    return data


def species(data, name):
    return data.get(name, np.zeros_like(data["time"]))


def with_time_offset(data, offset):
    shifted = dict(data)
    shifted["time"] = data["time"] + offset
    return shifted


def apply_plot_time_offset(scenario, gold, moose):
    offset = RESTART_TIME_OFFSETS.get(scenario, 0.0)
    if offset == 0.0:
        return gold, moose
    return with_time_offset(gold, offset), with_time_offset(moose, offset)


def infer_scenario(*paths):
    joined = " ".join(str(path) for path in paths)
    for scenario, _label, _gold_file in sorted(SCENARIOS, key=lambda item: len(item[0]), reverse=True):
        if scenario in joined:
            return scenario
    return ""


def plot_species(ax, gold, moose, column, ylabel, title, scale=1.0, legend_label="MOOSE"):
    gold_values = species(gold, column) * scale
    moose_values = species(moose, column) * scale
    ax.plot(gold["time"] / 3600.0, gold_values, "-", color="#1f77b4",
            label="F0AM gold", lw=1.2)
    ax.plot(moose["time"] / 3600.0, moose_values, "o--", color="#d62728",
            label=legend_label, ms=3, lw=1.0)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def plot_scenario(gold, moose, label, save_path, legend_label="MOOSE"):
    """Create the chamber test figure for one scenario."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle(f"MOOSE vs F0AM - {label}", fontsize=11, fontweight="bold")

    plot_species(axes[0, 0], gold, moose, "C5H8", "C$_5$H$_8$ (ppb)",
                 "Isoprene", legend_label=legend_label)
    plot_species(axes[0, 1], gold, moose, "OH", "OH (ppt)",
                 "OH Radical", scale=1.0e6, legend_label=legend_label)

    gold_nox = species(gold, "NO") + species(gold, "NO2")
    moose_nox = species(moose, "NO") + species(moose, "NO2")
    axes[1, 0].plot(gold["time"] / 3600.0, gold_nox, "-", color="#1f77b4",
                    label="F0AM gold", lw=1.2)
    axes[1, 0].plot(moose["time"] / 3600.0, moose_nox, "o--", color="#d62728",
                    label=legend_label, ms=3, lw=1.0)
    axes[1, 0].set_ylabel("NO$_x$ (ppb)")
    axes[1, 0].set_title("NO + NO$_2$")

    plot_species(axes[1, 1], gold, moose, "O3", "O$_3$ (ppb)",
                 "Ozone", legend_label=legend_label)

    for ax in axes.flat:
        ax.set_xlabel("Time (hours)")
        ax.legend(loc="best", framealpha=0.8)
        ax.grid(True)

    fig.tight_layout()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path)
    print(f"  Saved: {save_path}")
    plt.close(fig)


def plot_radical_diagnostics(gold, moose, label, save_path, legend_label="MOOSE"):
    """Create an early-time radical-family figure for chamber analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(f"Early Radical Diagnostics - {label}", fontsize=11, fontweight="bold")
    specs = [
        ("OH", "OH (ppt)", 1.0e6),
        ("RO2", "RO$_2$ (ppt)", 1.0e6),
        ("HO2", "HO$_2$ (ppt)", 1.0e6),
        ("CH3O2", "CH$_3$O$_2$ (ppt)", 1.0e6),
    ]
    for ax, (column, ylabel, scale) in zip(axes.flat, specs):
        plot_species(ax, gold, moose, column, ylabel, column, scale=scale,
                     legend_label=legend_label)
        ax.set_xlim((min(gold["time"].min(), moose["time"].min()) / 3600.0,
                     min(EARLY_WINDOW_SECONDS + min(gold["time"].min(), moose["time"].min()),
                         max(gold["time"].max(), moose["time"].max())) / 3600.0))
        ax.set_xlabel("Time (hours)")
        ax.legend(loc="best", framealpha=0.8)
        ax.grid(True)
    fig.tight_layout()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path)
    print(f"  Saved: {save_path}")
    plt.close(fig)


def diagnostic_rows(gold, moose, scenario, label, legend_label):
    rows = []
    times = gold["time"]
    early = times <= (times.min() + EARLY_WINDOW_SECONDS)
    for column in DIAGNOSTIC_COLUMNS:
        if column not in gold or column not in moose:
            continue
        gold_values = gold[column][early]
        moose_values = moose[column][early]
        time_values = times[early]
        if gold_values.size == 0:
            continue
        denom = np.maximum(np.maximum(np.abs(gold_values), np.abs(moose_values)), 1.0e-300)
        rel = np.abs(moose_values - gold_values) / denom
        index = int(np.argmax(rel))
        rows.append({
            "scenario": scenario,
            "label": label,
            "series": legend_label,
            "species": column,
            "time": f"{time_values[index]:.0f}",
            "max_rel_error": f"{rel[index]:.6e}",
            "max_abs_error": f"{abs(moose_values[index] - gold_values[index]):.6e}",
        })
    return rows


def write_diagnostics_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["scenario", "label", "series", "species", "time",
                  "max_rel_error", "max_abs_error"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {path}")


def check_default_files(chamber_dir, gold_dir):
    all_ok = True
    for _scenario, _label, gold_file in SCENARIOS:
        for desc, directory in (("gold", gold_dir), ("moose", chamber_dir)):
            path = directory / gold_file
            if path.exists():
                print(f"  [ok] {desc:>5}  {gold_file:40s}  {path.stat().st_size // 1024:4d} KB")
            else:
                print(f"  [!!] {desc:>5}  {gold_file:40s}  MISSING")
                all_ok = False
    return all_ok


def plot_default(output_dir, diagnostics_csv=None):
    print("Loading chamber test and F0AM gold data...")
    output_dir = output_dir or CHAMBER_DIR
    diagnostics = []
    for scenario, label, gold_file in SCENARIOS:
        gold_path = GOLD_DIR / gold_file
        moose_path = CHAMBER_DIR / gold_file
        if not gold_path.exists() or not moose_path.exists():
            print(f"SKIP {scenario}: missing gold or MOOSE CSV")
            continue
        gold = load_csv(gold_path)
        moose = load_csv(moose_path)
        gold, moose = apply_plot_time_offset(scenario, gold, moose)
        save_path = output_dir / f"chamber_{scenario}_comparison.png"
        plot_scenario(gold, moose, label, save_path)
        plot_radical_diagnostics(gold, moose, label,
                                 output_dir / f"chamber_{scenario}_radicals.png")
        diagnostics.extend(diagnostic_rows(gold, moose, scenario, label, "MOOSE"))
    if diagnostics_csv:
        write_diagnostics_csv(diagnostics_csv, diagnostics)


def plot_benchmark_run(run_dir, solvers, output_dir, legend_label, diagnostics_csv=None):
    print(f"Loading benchmark run: {run_dir}")
    output_dir = output_dir or (run_dir / "figures")
    diagnostics = []
    for solver_name in solvers:
        for scenario, (label, gold_file) in SCENARIO_BY_NAME.items():
            moose_path = run_dir / f"{solver_name}_{scenario}.csv"
            gold_path = GOLD_DIR / gold_file
            if not moose_path.exists():
                print(f"SKIP {solver_name} {scenario}: missing {moose_path.name}")
                continue
            if not gold_path.exists():
                print(f"SKIP {solver_name} {scenario}: missing gold {gold_path.name}")
                continue
            gold = load_csv(gold_path)
            moose = load_csv(moose_path)
            gold, moose = apply_plot_time_offset(scenario, gold, moose)
            save_path = output_dir / f"{solver_name}_{scenario}_comparison.png"
            plot_scenario(gold, moose, label, save_path, legend_label=legend_label)
            plot_radical_diagnostics(gold, moose, label,
                                     output_dir / f"{solver_name}_{scenario}_radicals.png",
                                     legend_label=legend_label)
            diagnostics.extend(diagnostic_rows(gold, moose, scenario, label, solver_name))
    if diagnostics_csv:
        write_diagnostics_csv(diagnostics_csv, diagnostics)


def main():
    parser = argparse.ArgumentParser(
        description="Plot F0AM chamber comparison figures for test and benchmark CSV outputs")
    parser.add_argument("--moose", help="MOOSE/solver CSV path for single-scenario mode")
    parser.add_argument("--gold", help="F0AM gold CSV path for single-scenario mode")
    parser.add_argument("--save", help="Output path for single-scenario mode")
    parser.add_argument("--run-dir", help="Benchmark solver run directory")
    parser.add_argument("--solvers", default="kpp_rosenbrock,petsc_bdf",
                        help="Comma-separated solver names for --run-dir mode")
    parser.add_argument("--output-dir", help="Directory for generated figures")
    parser.add_argument("--legend-label", default="MOOSE",
                        help="Label used for the non-F0AM curve")
    parser.add_argument("--diagnostics-csv",
                        help="Write early-time diagnostic error summary CSV")
    parser.add_argument("--check", action="store_true",
                        help="Verify default chamber files exist without plotting")
    args = parser.parse_args()

    output_dir = resolve_path(args.output_dir) if args.output_dir else None
    diagnostics_csv = resolve_path(args.diagnostics_csv) if args.diagnostics_csv else None

    if args.check:
        sys.exit(0 if check_default_files(CHAMBER_DIR, GOLD_DIR) else 1)

    if args.moose or args.gold:
        if not args.moose or not args.gold:
            raise SystemExit("ERROR: --moose and --gold must be supplied together")
        moose_path = resolve_path(args.moose)
        gold_path = resolve_path(args.gold)
        if not moose_path.exists() or not gold_path.exists():
            raise SystemExit("ERROR: specified --moose or --gold file not found")
        save_path = resolve_path(args.save) if args.save else Path("chamber_comparison.png").resolve()
        scenario = infer_scenario(gold_path, moose_path)
        gold = load_csv(gold_path)
        moose = load_csv(moose_path)
        gold, moose = apply_plot_time_offset(scenario, gold, moose)
        plot_scenario(gold, moose, gold_path.stem, save_path)
        if diagnostics_csv:
            write_diagnostics_csv(
                diagnostics_csv,
                diagnostic_rows(gold, moose, scenario, gold_path.stem, args.legend_label))
        return

    if args.run_dir:
        run_dir = resolve_path(args.run_dir)
        if not run_dir.exists():
            raise SystemExit(f"ERROR: benchmark run directory not found: {run_dir}")
        plot_benchmark_run(run_dir, parse_solvers(args.solvers), output_dir,
                           args.legend_label, diagnostics_csv)
        return

    plot_default(output_dir, diagnostics_csv)


if __name__ == "__main__":
    main()
