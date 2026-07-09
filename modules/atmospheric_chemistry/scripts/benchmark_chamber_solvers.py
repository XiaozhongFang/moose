#!/usr/bin/env python3
"""Benchmark chamber runtime and F0AM accuracy across box solver backends.

The chamber timing unit is the current four-test workflow:
S1, S2, S2b restart from S2, and S3. The F0AM baseline supplied with
--f0am-seconds is interpreted as the total compute time for those same four
work items, not a per-scenario time.

Examples
--------
List all solver selections:

    python3 scripts/benchmark_chamber_solvers.py --list-solvers

Generate solver-specific chamber inputs without running them:

    python3 scripts/benchmark_chamber_solvers.py \
        --solvers kpp_rosenbrock,petsc_bdf \
        --write-inputs-only \
        --output-dir kpp_chamber/solver_runs/input_check

Run the accuracy/timing benchmark against the F0AM gold CSV files:

    python3 scripts/benchmark_chamber_solvers.py \
        --solvers kpp_rosenbrock,petsc_bdf \
        --rtol 5e-4 --atol 5e-8 \
        --compare-gold --plot-comparison \
        --f0am-seconds 7 \
        --output-dir kpp_chamber/solver_runs/f0am_accuracy_timing
"""

import argparse
import csv
import datetime as _dt
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
REPO_ROOT = MODULE_DIR.parent.parent
CHAMBER_DIR = MODULE_DIR / "test" / "tests" / "chamber"

SCENARIOS = [
    ("S1", "vs_F0AM_chamber_S1_box.i"),
    ("S2", "vs_F0AM_chamber_S2_box.i"),
    ("S2b", "vs_F0AM_chamber_S2b_box.i"),
    ("S3", "vs_F0AM_chamber_S3_box.i"),
]
GOLD_SCENARIOS = {
    "S1": "vs_F0AM_chamber_S1_box.csv",
    "S2": "vs_F0AM_chamber_S2_box.csv",
    "S2b": "vs_F0AM_chamber_S2b_box.csv",
    "S3": "vs_F0AM_chamber_S3_box.csv",
}

DIRECT_SOLVERS = {
    "petsc_bdf": ("petsc_ts", "bdf"),
    "petsc_arkimex": ("petsc_ts", "arkimex"),
    "petsc_eimex": ("petsc_ts", "eimex"),
    "petsc_rosw": ("petsc_ts", "rosw"),
    "petsc_mimex": ("petsc_ts", "mimex"),
    "petsc_beuler": ("petsc_ts", "beuler"),
    "petsc_cn": ("petsc_ts", "cn"),
    "petsc_rk": ("petsc_ts", "rk"),
    "petsc_theta": ("petsc_ts", "theta"),
    "petsc_ssp": ("petsc_ts", "ssp"),
    "petsc_sundials": ("petsc_ts", "sundials"),
    "sundials": ("sundials", "bdf"),
    "moose_implicit": ("moose_implicit", "bdf"),
}

KPP_SOLVERS = ("kpp_rosenbrock", "kpp_sdirk", "kpp_runge_kutta")
KPP_INTEGRATORS = {
    "kpp_rosenbrock": "rosenbrock",
    "kpp_sdirk": "sdirk",
    "kpp_runge_kutta": "runge_kutta",
}
KPP_FAC = MODULE_DIR / "doc" / "content" / "modules" / "atmospheric_chemistry" / "database" / "MCMv331_Inorg_Isoprene.fac"
KPP_BASE_DIR = CHAMBER_DIR / "kpp_chamber" / "generated_mechanisms"
KPP_BUILD_MAKEFILE = MODULE_DIR / "kpp" / "build" / "Makefile"
FAC_TO_KPP = SCRIPT_DIR / "fac_to_kpp.py"
PLOT_CHAMBER = SCRIPT_DIR / "plot_chamber_comparison.py"
DEFAULT_GOLD_DIR = CHAMBER_DIR / "gold"
DEFAULT_IGNORE_COLUMNS = ("ONE", "RO2", "CH3ONO")
SCENARIO_TOLERANCES = {
    "S1": (1.0e-1, 1.0e4),
    "S2": (5.0e-2, 1.0),
    "S2b": (5.0e-2, 1.0e-2),
    "S3": (1.0e-1, 1.8e4),
}
DIAGNOSTIC_COLUMNS = ("OH", "RO2", "C5H8", "NO", "NO2", "O3", "HO2", "CH3O2",
                      "MVK", "MACR", "HCHO")
EARLY_WINDOW_SECONDS = 1010.0

TIMING_FIELDS = [
    "solver", "repeat", "scenario", "status", "seconds", "returncode",
    "accuracy_status", "max_rel_error", "max_abs_error", "compared_columns",
    "accuracy_message", "gold", "input", "log",
]
SUMMARY_FIELDS = ["solver", "repeat", "status", "total_seconds",
                  "f0am_seconds", "moose_over_f0am", "accuracy_status",
                  "max_rel_error", "max_abs_error", "compared_scenarios",
                  "accuracy_failed_scenarios", "failed_scenarios",
                  "analysis"]


def find_app(path_arg):
    if path_arg:
        path = Path(path_arg).resolve()
        if not os.access(path, os.X_OK):
            raise SystemExit(f"ERROR: app is not executable: {path}")
        return path

    for path in (MODULE_DIR / "atmospheric_chemistry-opt",
                 REPO_ROOT / "atmospheric_chemistry-opt"):
        if os.access(path, os.X_OK):
            return path
    raise SystemExit("ERROR: atmospheric_chemistry-opt not found. Build the module first.")


def block_bounds(lines, block_name):
    start = None
    for i, line in enumerate(lines):
        if line.strip() == f"[{block_name}]":
            start = i
            break
    if start is None:
        raise ValueError(f"block [{block_name}] not found")

    for i in range(start + 1, len(lines)):
        if lines[i].strip() == "[]":
            return start, i
    raise ValueError(f"block [{block_name}] has no closing []")


def set_param(lines, block_name, param, value):
    start, end = block_bounds(lines, block_name)
    pattern = re.compile(rf"^(\s*){re.escape(param)}\s*=")
    for i in range(start + 1, end):
        match = pattern.match(lines[i])
        if match:
            lines[i] = f"{match.group(1)}{param} = {value}"
            return

    block_indent = re.match(r"^(\s*)", lines[start]).group(1)
    lines.insert(end, f"{block_indent}  {param} = {value}")


def input_relpath(path):
    return str(path.resolve().relative_to(CHAMBER_DIR.resolve()))


def display_path(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(CHAMBER_DIR.resolve()))
    except ValueError:
        return str(path)


def resolve_existing_path(path):
    if path.is_absolute():
        return path
    for base in (Path.cwd(), MODULE_DIR, REPO_ROOT):
        candidate = (base / path).resolve()
        if candidate.exists():
            return candidate
    return (Path.cwd() / path).resolve()


def output_base(run_dir, solver_name, scenario):
    return run_dir / f"{solver_name}_{scenario}"


def ensure_kpp_mechanism(solver_name, build):
    integrator = KPP_INTEGRATORS[solver_name]
    model = f"chamber_mcm_{integrator}"
    out_dir = KPP_BASE_DIR / model
    mech = out_dir / f"{model}.kpp"

    subprocess.run(
        [sys.executable, str(FAC_TO_KPP), str(KPP_FAC),
         "--output-dir", str(out_dir), "--model", model, "--integrator", integrator],
        cwd=REPO_ROOT,
        check=True,
    )
    if build:
        subprocess.run(
            ["make", "-f", str(KPP_BUILD_MAKEFILE), f"MECH={mech}"],
            cwd=REPO_ROOT,
            check=True,
        )
    return mech


def generate_input(template_path, run_dir, solver_name, solver, solver_type,
                   scenario, s2_restart_base, keep_csv, rtol, atol,
                   mechanism_file=None):
    text = template_path.read_text()
    lines = text.splitlines()

    if mechanism_file is not None:
        set_param(lines, "Box", "mechanism_file", f"'{input_relpath(mechanism_file)}'")
    set_param(lines, "Box", "chem_solver", solver)
    set_param(lines, "Box", "chem_solver_type", solver_type)
    set_param(lines, "Box", "chem_solver_rtol", rtol)
    set_param(lines, "Box", "chem_solver_atol", atol)
    if solver_name in KPP_SOLVERS:
        set_param(lines, "Box", "output_ro2_sum", "false")
    set_param(lines, "Outputs", "console", "false")
    set_param(lines, "Outputs", "csv", "true" if keep_csv else "false")
    set_param(lines, "Outputs", "checkpoint", "true" if scenario == "S2" else "false")
    # Outputs/file_base is resolved relative to the chamber run directory used
    # by this script. Problem/restart_file_base is resolved relative to the
    # generated input file, so S2b uses the local S2 checkpoint name below.
    set_param(lines, "Outputs", "file_base",
              f"'{input_relpath(output_base(run_dir, solver_name, scenario))}'")

    if scenario == "S2b":
        if not s2_restart_base:
            raise ValueError("S2b requires an S2 restart base")
        set_param(lines, "Problem", "restart_file_base", f"'{s2_restart_base}_cp/LATEST'")

    generated = run_dir / f"{solver_name}_{scenario}.i"
    generated.write_text("\n".join(lines) + "\n")
    return generated


def run_case(app, input_path, log_path, timeout):
    start = time.perf_counter()
    with open(log_path, "w") as log:
        proc = subprocess.run(
            [str(app), "-i", input_relpath(input_path)],
            cwd=CHAMBER_DIR,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            check=False,
        )
    elapsed = time.perf_counter() - start
    return proc.returncode, elapsed


def parse_solvers(value):
    if value == "all":
        return list(DIRECT_SOLVERS) + list(KPP_SOLVERS)

    solvers = []
    for item in value.split(","):
        name = item.strip()
        if not name:
            continue
        if name not in DIRECT_SOLVERS and name not in KPP_SOLVERS:
            valid = ", ".join(list(DIRECT_SOLVERS) + list(KPP_SOLVERS) + ["all"])
            raise SystemExit(f"ERROR: unknown solver '{name}'. Valid values: {valid}")
        solvers.append(name)
    return solvers


def parse_column_list(value):
    return [item for item in re.split(r"[,\s]+", value.strip()) if item]


def write_rows(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames,
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_numeric_csv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    return fieldnames, rows


def read_ro2_species(path=KPP_FAC):
    if not path.exists():
        return []
    match = re.search(r"\bRO2\s*=\s*(.*?)\s*;", path.read_text(), flags=re.DOTALL)
    if not match:
        return []
    species = []
    for item in match.group(1).replace("\n", " ").split("+"):
        name = re.sub(r"[^A-Za-z0-9_]", "", item.strip())
        if name:
            species.append(name)
    return species


RO2_SPECIES = read_ro2_species()


def row_value(row, column, header):
    if column in row:
        return float(row[column])
    if column == "RO2" and RO2_SPECIES:
        total = 0.0
        found = False
        for species in RO2_SPECIES:
            if species in header:
                total += float(row[species])
                found = True
        if found:
            return total
    raise KeyError(column)


def compare_csv_to_gold(moose_csv, gold_csv, rel_err, abs_zero, ignore_columns):
    """Compare one solver CSV against a F0AM gold CSV using CSVDiff-like tolerances."""
    if not moose_csv.exists():
        return {
            "status": "FAIL",
            "max_rel_error": "",
            "max_abs_error": "",
            "compared_columns": "0",
            "message": f"missing output CSV: {display_path(moose_csv)}",
        }
    if not gold_csv.exists():
        return {
            "status": "NO_GOLD",
            "max_rel_error": "",
            "max_abs_error": "",
            "compared_columns": "0",
            "message": f"missing gold CSV: {display_path(gold_csv)}",
        }

    ignored = set(ignore_columns)
    gold_header, gold_rows = read_numeric_csv(gold_csv)
    moose_header, moose_rows = read_numeric_csv(moose_csv)
    gold_columns = [name for name in gold_header if name not in ignored]
    moose_columns = set(name for name in moose_header if name not in ignored)
    common_columns = [name for name in gold_columns if name in moose_columns]
    missing_columns = [name for name in gold_columns if name not in moose_columns]

    messages = []
    status = "OK"
    if missing_columns:
        status = "FAIL"
        shown = ",".join(missing_columns[:8])
        suffix = "..." if len(missing_columns) > 8 else ""
        messages.append(f"missing columns: {shown}{suffix}")
    if len(gold_rows) != len(moose_rows):
        status = "FAIL"
        messages.append(f"row count differs: moose={len(moose_rows)} gold={len(gold_rows)}")
    if not common_columns:
        status = "FAIL"
        messages.append("no comparable columns")

    n_rows = min(len(gold_rows), len(moose_rows))
    max_rel = 0.0
    max_abs = 0.0
    fail_count = 0
    fail_examples = []
    for row_index in range(n_rows):
        gold_row = gold_rows[row_index]
        moose_row = moose_rows[row_index]
        for column in common_columns:
            try:
                gold_value = float(gold_row[column])
                moose_value = float(moose_row[column])
            except (KeyError, ValueError):
                status = "FAIL"
                fail_count += 1
                if len(fail_examples) < 5:
                    fail_examples.append(f"{column}@row{row_index}: nonnumeric")
                continue
            if not math.isfinite(gold_value) or not math.isfinite(moose_value):
                status = "FAIL"
                fail_count += 1
                if len(fail_examples) < 5:
                    fail_examples.append(f"{column}@row{row_index}: nonfinite")
                continue

            gold_check = 0.0 if abs(gold_value) < abs_zero else gold_value
            moose_check = 0.0 if abs(moose_value) < abs_zero else moose_value
            diff = abs(moose_check - gold_check)
            if gold_check == 0.0 and moose_check == 0.0:
                rel = 0.0
            else:
                rel = diff / max(abs(gold_check), abs(moose_check))
            max_rel = max(max_rel, rel)
            max_abs = max(max_abs, diff)
            if rel > rel_err:
                status = "FAIL"
                fail_count += 1
                if len(fail_examples) < 5:
                    fail_examples.append(
                        f"{column}@row{row_index}: rel={rel:.3e} abs={diff:.3e}")

    if fail_count:
        messages.append(f"{fail_count} values exceed tolerances")
        messages.extend(fail_examples)

    return {
        "status": status,
        "max_rel_error": f"{max_rel:.6e}" if common_columns else "",
        "max_abs_error": f"{max_abs:.6e}" if common_columns else "",
        "compared_columns": str(len(common_columns)),
        "message": "; ".join(messages),
    }


def scenario_tolerances(scenario, args):
    if args.rel_err is not None or args.abs_zero is not None:
        rel = args.rel_err if args.rel_err is not None else 5.0e-2
        zero = args.abs_zero if args.abs_zero is not None else 2.0e-6
        return rel, zero
    return SCENARIO_TOLERANCES.get(scenario, (5.0e-2, 2.0e-6))


def diagnostic_errors(moose_csv, gold_csv):
    """Return early-time relative errors for atmospheric-chemistry diagnostics."""
    if not moose_csv.exists() or not gold_csv.exists():
        return []

    gold_header, gold_rows = read_numeric_csv(gold_csv)
    moose_header, moose_rows = read_numeric_csv(moose_csv)
    moose_columns = set(moose_header)
    rows = []
    n_rows = min(len(gold_rows), len(moose_rows))
    for column in DIAGNOSTIC_COLUMNS:
        if column not in gold_header:
            continue
        if column not in moose_columns and not (column == "RO2" and RO2_SPECIES):
            continue
        max_rel = 0.0
        max_abs = 0.0
        max_time = ""
        for row_index in range(n_rows):
            gold_row = gold_rows[row_index]
            try:
                time_value = float(gold_row["time"])
            except (KeyError, ValueError):
                continue
            if time_value > EARLY_WINDOW_SECONDS:
                continue
            try:
                gold_value = row_value(gold_row, column, gold_header)
                moose_value = row_value(moose_rows[row_index], column, moose_header)
            except (KeyError, ValueError):
                continue
            diff = abs(moose_value - gold_value)
            rel = diff / max(abs(gold_value), abs(moose_value), 1.0e-300)
            if rel > max_rel:
                max_rel = rel
                max_abs = diff
                max_time = f"{time_value:g}"
        rows.append({
            "column": column,
            "max_rel": f"{max_rel:.6e}",
            "max_abs": f"{max_abs:.6e}",
            "time": max_time,
        })
    return rows


def write_analysis(path, timing_rows, summary_rows, diagnostic_rows, f0am_seconds, compare_gold,
                   figure_dir=None):
    """Write a compact Markdown archive of solver timing and accuracy results."""
    lines = [
        "# Chamber Solver Timing Analysis",
        "",
        f"Generated: {_dt.datetime.now().isoformat(timespec='seconds')}",
        f"F0AM baseline scope: S1 + S2 + S2b restart + S3 ({f0am_seconds:.3f} s)",
    ]
    if compare_gold:
        lines.append("Accuracy gold: F0AM CSVs for S1, S2, S2b, and S3.")
    else:
        lines.append("Accuracy gold: not checked in this run.")
    if figure_dir is not None:
        lines.append(f"Comparison figures: `{input_relpath(figure_dir)}`")

    lines.extend([
        "",
        "## Solver Totals",
        "",
        "| Solver | Repeat | Status | Accuracy | Total seconds | vs F0AM | Max rel error | Failed scenarios |",
        "|--------|--------|--------|----------|---------------|---------|---------------|------------------|",
    ])
    for row in summary_rows:
        failed_text = row.get("failed_scenarios", "")
        if row.get("accuracy_failed_scenarios"):
            failed_text = (failed_text + ";" if failed_text else "") + \
                "accuracy=" + row.get("accuracy_failed_scenarios")
        lines.append(
            "| {solver} | {repeat} | {status} | {accuracy_status} | {total_seconds} | "
            "{moose_over_f0am} | {max_rel_error} | {failed_text} |".format(
                **{key: row.get(key, "") for key in SUMMARY_FIELDS}, failed_text=failed_text,
            )
        )

    lines.extend([
        "",
        "## Scenario Timings",
        "",
        "| Solver | Scenario | Status | Seconds | Accuracy | Max rel error | Compared columns |",
        "|--------|----------|--------|---------|----------|---------------|------------------|",
    ])
    for row in timing_rows:
        lines.append(
            "| {solver} | {scenario} | {status} | {seconds} | {accuracy_status} | "
            "{max_rel_error} | {compared_columns} |".format(
                **{key: row.get(key, "") for key in TIMING_FIELDS})
        )

    if diagnostic_rows:
        lines.extend([
            "",
            "## Early-Time Diagnostic Errors",
            "",
            "Maximum relative errors for selected atmospheric-chemistry diagnostics "
            "at output times <= 1000 s. Concentrations are compared in molec/cm^3.",
            "",
            "| Solver | Scenario | Species | Time (s) | Max rel error | Max abs error |",
            "|--------|----------|---------|----------|---------------|---------------|",
        ])
        for row in diagnostic_rows:
            lines.append(
                "| {solver} | {scenario} | {column} | {time} | {max_rel} | {max_abs} |".format(
                    **row)
            )

    path.write_text("\n".join(lines) + "\n")


def write_comparison_figures(run_dir, solvers):
    """Archive chamber comparison figures using the same script as chamber tests."""
    figure_dir = run_dir / "figures"
    diagnostics_csv = run_dir / "chamber_diagnostic_errors.csv"
    subprocess.run(
        [sys.executable, str(PLOT_CHAMBER),
         "--run-dir", str(run_dir),
         "--solvers", ",".join(solvers),
         "--output-dir", str(figure_dir),
         "--diagnostics-csv", str(diagnostics_csv)],
        cwd=MODULE_DIR,
        check=True,
    )
    return figure_dir


def main():
    parser = argparse.ArgumentParser(description="Benchmark chamber solver runtime.")
    parser.add_argument("--app", help="Path to atmospheric_chemistry-opt")
    parser.add_argument("--solvers", default="all",
                        help="Comma-separated solver list, or 'all' to cover direct FAC and generated KPP solvers")
    parser.add_argument("--repeat", type=int, default=1,
                        help="Number of repeats per solver")
    parser.add_argument("--f0am-seconds", type=float, default=6.0,
                        help="F0AM total seconds for S1+S2+S2b+S3")
    parser.add_argument("--rtol", default="5e-4")
    parser.add_argument("--atol", default="5e-8")
    parser.add_argument("--timeout", type=float, default=3600.0,
                        help="Per-scenario timeout in seconds")
    parser.add_argument("--keep-csv", action="store_true",
                        help="Keep MOOSE CSV output for each generated input")
    parser.add_argument("--compare-gold", action="store_true",
                        help="Compare S1/S2/S2b/S3 output CSV files against F0AM gold CSV files")
    parser.add_argument("--gold-dir", type=Path, default=DEFAULT_GOLD_DIR,
                        help="Directory containing F0AM chamber gold CSV files")
    parser.add_argument("--rel-err", type=float,
                        help="Override the scenario-specific relative tolerance used by --compare-gold")
    parser.add_argument("--abs-zero", type=float,
                        help="Override the scenario-specific absolute zero threshold used by --compare-gold")
    parser.add_argument("--ignore-columns", default=" ".join(DEFAULT_IGNORE_COLUMNS),
                        help="Columns ignored by --compare-gold")
    parser.add_argument("--analysis-file", default="chamber_solver_analysis.md",
                        help="Markdown timing analysis path relative to --output-dir")
    parser.add_argument("--plot-comparison", action="store_true",
                        help="Write chamber comparison figures for solver CSV outputs")
    parser.add_argument("--output-dir",
                        help="Benchmark output directory relative to chamber test dir")
    parser.add_argument("--write-inputs-only", action="store_true",
                        help="Generate solver/scenario input files without running")
    parser.add_argument("--list-solvers", action="store_true",
                        help="List solver names and exit")
    args = parser.parse_args()

    if args.list_solvers:
        print("Direct FAC chamber solvers:")
        for name, (solver, solver_type) in DIRECT_SOLVERS.items():
            print(f"  {name:16s} chem_solver={solver} chem_solver_type={solver_type}")
        print("Generated KPP chamber solvers:")
        for name in KPP_SOLVERS:
            print(f"  {name:16s} integrator={KPP_INTEGRATORS[name]}")
        return 0

    if args.write_inputs_only and args.compare_gold:
        raise SystemExit("ERROR: --compare-gold cannot be used with --write-inputs-only")
    if args.write_inputs_only and args.plot_comparison:
        raise SystemExit("ERROR: --plot-comparison cannot be used with --write-inputs-only")
    if args.compare_gold or args.plot_comparison:
        args.keep_csv = True

    solvers = parse_solvers(args.solvers)
    ignore_columns = parse_column_list(args.ignore_columns)
    gold_dir = resolve_existing_path(args.gold_dir)
    if not gold_dir.exists() and not args.write_inputs_only:
        raise SystemExit(f"ERROR: gold directory does not exist: {gold_dir}")
    kpp_mechanisms = {}
    for solver_name in solvers:
        if solver_name in KPP_SOLVERS:
            kpp_mechanisms[solver_name] = ensure_kpp_mechanism(
                solver_name, build=not args.write_inputs_only)

    app = None if args.write_inputs_only else find_app(args.app)

    timestamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.output_dir:
        run_dir = Path(args.output_dir)
        if not run_dir.is_absolute():
            run_dir = CHAMBER_DIR / run_dir
    else:
        run_dir = CHAMBER_DIR / "kpp_chamber" / "solver_runs" / timestamp
    run_dir = run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    try:
        run_dir.relative_to(CHAMBER_DIR.resolve())
    except ValueError as exc:
        raise SystemExit("ERROR: --output-dir must be inside test/tests/chamber") from exc

    timing_rows = []
    summary_rows = []
    diagnostic_rows = []
    timing_csv = run_dir / "chamber_solver_timing.csv"
    summary_csv = run_dir / "chamber_solver_summary.csv"
    analysis_path = Path(args.analysis_file)
    if not analysis_path.is_absolute():
        analysis_path = run_dir / analysis_path

    interrupted = False
    try:
        for solver_name in solvers:
            if solver_name in KPP_SOLVERS:
                solver, solver_type = solver_name, "bdf"
                mechanism_file = kpp_mechanisms[solver_name]
            else:
                solver, solver_type = DIRECT_SOLVERS[solver_name]
                mechanism_file = None
            for repeat in range(1, args.repeat + 1):
                total = 0.0
                failed = []
                accuracy_failed = []
                compared_scenarios = []
                repeat_max_rel = 0.0
                repeat_max_abs = 0.0
                s2_base = f"{solver_name}_S2"

                for scenario, template in SCENARIOS:
                    template_path = CHAMBER_DIR / template
                    generated = generate_input(template_path, run_dir, solver_name, solver,
                                               solver_type, scenario, s2_base,
                                               args.keep_csv, args.rtol, args.atol,
                                               mechanism_file=mechanism_file)
                    log_path = run_dir / f"{solver_name}_{scenario}.log"

                    if args.write_inputs_only:
                        status = "WRITE_ONLY"
                        elapsed = 0.0
                        returncode = ""
                    elif scenario == "S2b" and "S2" in failed:
                        status = "SKIP"
                        elapsed = 0.0
                        returncode = ""
                    else:
                        print(f"Running {solver_name} {scenario} ...", flush=True)
                        try:
                            returncode, elapsed = run_case(app, generated, log_path, args.timeout)
                        except subprocess.TimeoutExpired:
                            returncode = "timeout"
                            elapsed = args.timeout
                        status = "OK" if returncode == 0 else "FAIL"
                        if status != "OK":
                            failed.append(scenario)
                        total += elapsed

                    accuracy = {
                        "status": "WRITE_ONLY" if args.write_inputs_only else "NOT_CHECKED",
                        "max_rel_error": "",
                        "max_abs_error": "",
                        "compared_columns": "0",
                        "message": "",
                        "gold": "",
                    }
                    if args.compare_gold and not args.write_inputs_only:
                        if status != "OK":
                            accuracy.update({
                                "status": "SKIP",
                                "message": "solver run did not complete",
                                "gold": display_path(gold_dir / GOLD_SCENARIOS[scenario]),
                            })
                        else:
                            gold_csv = gold_dir / GOLD_SCENARIOS[scenario]
                            moose_csv = output_base(run_dir, solver_name, scenario).with_suffix(".csv")
                            rel_err, abs_zero = scenario_tolerances(scenario, args)
                            accuracy = compare_csv_to_gold(
                                moose_csv, gold_csv, rel_err,
                                abs_zero, ignore_columns)
                            accuracy["gold"] = display_path(gold_csv)
                            compared_scenarios.append(scenario)
                            for diag in diagnostic_errors(moose_csv, gold_csv):
                                diag["solver"] = solver_name
                                diag["scenario"] = scenario
                                diagnostic_rows.append(diag)
                            if accuracy["status"] != "OK":
                                accuracy_failed.append(scenario)
                            if accuracy["max_rel_error"]:
                                repeat_max_rel = max(repeat_max_rel,
                                                     float(accuracy["max_rel_error"]))
                            if accuracy["max_abs_error"]:
                                repeat_max_abs = max(repeat_max_abs,
                                                     float(accuracy["max_abs_error"]))

                    timing_rows.append({
                        "solver": solver_name,
                        "repeat": repeat,
                        "scenario": scenario,
                        "status": status,
                        "seconds": f"{elapsed:.6f}",
                        "returncode": returncode,
                        "accuracy_status": accuracy["status"],
                        "max_rel_error": accuracy["max_rel_error"],
                        "max_abs_error": accuracy["max_abs_error"],
                        "compared_columns": accuracy["compared_columns"],
                        "accuracy_message": accuracy["message"],
                        "gold": accuracy["gold"],
                        "input": input_relpath(generated),
                        "log": input_relpath(log_path),
                    })
                    write_rows(timing_csv, TIMING_FIELDS, timing_rows)

                if args.write_inputs_only:
                    summary_status = "WRITE_ONLY"
                    accuracy_status = "WRITE_ONLY"
                    ratio = ""
                elif failed:
                    summary_status = "FAIL"
                    accuracy_status = "FAIL" if accuracy_failed else "SKIP"
                    ratio = f"{total / args.f0am_seconds:.3f}" if args.f0am_seconds > 0 else ""
                elif accuracy_failed:
                    summary_status = "FAIL"
                    accuracy_status = "FAIL"
                    ratio = f"{total / args.f0am_seconds:.3f}" if args.f0am_seconds > 0 else ""
                else:
                    summary_status = "OK"
                    accuracy_status = "OK" if args.compare_gold else "NOT_CHECKED"
                    ratio = f"{total / args.f0am_seconds:.3f}" if args.f0am_seconds > 0 else ""

                summary_rows.append({
                    "solver": solver_name,
                    "repeat": repeat,
                    "status": summary_status,
                    "total_seconds": f"{total:.6f}" if not args.write_inputs_only else "",
                    "f0am_seconds": args.f0am_seconds,
                    "moose_over_f0am": ratio,
                    "accuracy_status": accuracy_status,
                    "max_rel_error": f"{repeat_max_rel:.6e}" if compared_scenarios else "",
                    "max_abs_error": f"{repeat_max_abs:.6e}" if compared_scenarios else "",
                    "compared_scenarios": ";".join(compared_scenarios),
                    "accuracy_failed_scenarios": ";".join(accuracy_failed),
                    "failed_scenarios": ";".join(failed),
                    "analysis": input_relpath(analysis_path),
                })
                write_rows(summary_csv, SUMMARY_FIELDS, summary_rows)
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrupted; partial benchmark results were preserved.", file=sys.stderr)

    figure_dir = None
    if args.plot_comparison and not interrupted:
        figure_dir = write_comparison_figures(run_dir, solvers)

    write_rows(timing_csv, TIMING_FIELDS, timing_rows)
    write_rows(summary_csv, SUMMARY_FIELDS, summary_rows)
    write_analysis(analysis_path, timing_rows, summary_rows, diagnostic_rows,
                   args.f0am_seconds, args.compare_gold, figure_dir)

    print(f"Wrote {input_relpath(timing_csv)}")
    print(f"Wrote {input_relpath(summary_csv)}")
    print(f"Wrote {input_relpath(analysis_path)}")
    if figure_dir is not None:
        print(f"Wrote {input_relpath(figure_dir)}")
    print("F0AM baseline scope: S1 + S2 + S2b restart + S3")
    if interrupted:
        return 130
    if any(row.get("status") == "FAIL" for row in summary_rows):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
