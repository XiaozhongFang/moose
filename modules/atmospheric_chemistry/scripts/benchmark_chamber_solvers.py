#!/usr/bin/env python3
"""Benchmark chamber runtime across box chemical solver backends.

The chamber timing unit is the current four-test workflow:
S1, S2, S2b restart from S2, and S3. The F0AM baseline supplied with
--f0am-seconds is interpreted as the total compute time for those same four
work items, not a per-scenario time.
"""

import argparse
import csv
import datetime as _dt
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
KPP_BASE_DIR = CHAMBER_DIR / "kpp_chamber"
KPP_BUILD_MAKEFILE = MODULE_DIR / "kpp" / "build" / "Makefile"
FAC_TO_KPP = SCRIPT_DIR / "fac_to_kpp.py"

TIMING_FIELDS = ["solver", "repeat", "scenario", "status", "seconds",
                 "returncode", "input", "log"]
SUMMARY_FIELDS = ["solver", "repeat", "status", "total_seconds",
                  "f0am_seconds", "moose_over_f0am", "failed_scenarios"]


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


def write_rows(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


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

    solvers = parse_solvers(args.solvers)
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
        run_dir = CHAMBER_DIR / "benchmark_runs" / timestamp
    run_dir = run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    try:
        run_dir.relative_to(CHAMBER_DIR.resolve())
    except ValueError as exc:
        raise SystemExit("ERROR: --output-dir must be inside test/tests/chamber") from exc

    timing_rows = []
    summary_rows = []
    timing_csv = run_dir / "chamber_solver_timing.csv"
    summary_csv = run_dir / "chamber_solver_summary.csv"

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

                    timing_rows.append({
                        "solver": solver_name,
                        "repeat": repeat,
                        "scenario": scenario,
                        "status": status,
                        "seconds": f"{elapsed:.6f}",
                        "returncode": returncode,
                        "input": input_relpath(generated),
                        "log": input_relpath(log_path),
                    })
                    write_rows(timing_csv, TIMING_FIELDS, timing_rows)

                if args.write_inputs_only:
                    summary_status = "WRITE_ONLY"
                    ratio = ""
                elif failed:
                    summary_status = "FAIL"
                    ratio = ""
                else:
                    summary_status = "OK"
                    ratio = f"{total / args.f0am_seconds:.3f}" if args.f0am_seconds > 0 else ""

                summary_rows.append({
                    "solver": solver_name,
                    "repeat": repeat,
                    "status": summary_status,
                    "total_seconds": f"{total:.6f}" if not args.write_inputs_only else "",
                    "f0am_seconds": args.f0am_seconds,
                    "moose_over_f0am": ratio,
                    "failed_scenarios": ";".join(failed),
                })
                write_rows(summary_csv, SUMMARY_FIELDS, summary_rows)
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrupted; partial benchmark results were preserved.", file=sys.stderr)

    write_rows(timing_csv, TIMING_FIELDS, timing_rows)
    write_rows(summary_csv, SUMMARY_FIELDS, summary_rows)

    print(f"Wrote {input_relpath(timing_csv)}")
    print(f"Wrote {input_relpath(summary_csv)}")
    print("F0AM baseline scope: S1 + S2 + S2b restart + S3")
    return 130 if interrupted else 0


if __name__ == "__main__":
    sys.exit(main())
