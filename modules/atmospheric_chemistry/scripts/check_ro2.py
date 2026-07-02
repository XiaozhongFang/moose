#!/usr/bin/env python3
"""Run the atmospheric_chemistry module to extract RO2 species names from a
.fac mechanism file.

Usage:
    cd <work_dir>
    python3 check_ro2.py <mechanism.fac> [-o ro2_species.txt]

The script runs from the current working directory (or the output file's
directory if -o is given), generates a minimal MOOSE input, runs
atmospheric_chemistry-opt (expected at ../atmospheric_chemistry-opt
relative to the script), and parses the RO2_SPECIES(N): name1,name2,...
line from the console output.

Examples:
    python3 ~/git_repo/moose/modules/atmospheric_chemistry/scripts/check_ro2.py \
        ~/git_repo/moose/modules/atmospheric_chemistry/doc/content/modules/atmospheric_chemistry/database/MCMv331_Inorg_Isoprene.fac \
        -o ro2_secies.txt
"""
import argparse, subprocess, sys, re, os, tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BINARY = str(SCRIPT_DIR.parent / "atmospheric_chemistry-opt")


def run_module_get_ro2(fac_file, binary, workdir):
    """Run atmospheric_chemistry-opt in workdir and extract RO2 species names."""
    fac_path = str(Path(os.path.expanduser(fac_file)).resolve())
    bin_path = str(Path(binary).resolve())

    if not os.path.exists(bin_path):
        raise FileNotFoundError(f"binary not found: {bin_path}")
    if not os.path.exists(fac_path):
        raise FileNotFoundError(f"mechanism file not found: {fac_path}")
    if not os.path.isdir(workdir):
        os.makedirs(workdir, exist_ok=True)

    # mechanism_file must be relative to workdir
    rel_path = os.path.relpath(fac_path, workdir)
    input_content = (
        "[Mesh]\n"
        "  [gen]\n"
        "    type = GeneratedMeshGenerator\n"
        "    dim = 1\n"
        "    nx = 1\n"
        "  []\n"
        "[]\n"
        "\n"
        "[AtmosphericChemistry]\n"
        "  mode = box\n"
        "  mechanism_file = '" + rel_path + "'\n"
        "  temperature = 298.0\n"
        "  air_density = 2.46e19\n"
        "  photolysis_scheme = MCM_SZA\n"
        "  jfac = 1.0\n"
        "[]\n"
        "\n"
        "[VectorPostprocessors]\n"
        "  [ro2_list]\n"
        "    type = MCMRO2ListPostprocessor\n"
        "    box_model = box_model\n"
        "  []\n"
        "[]\n"
        "\n"
        "[Executioner]\n"
        "  type = Transient\n"
        "  dt = 1\n"
        "  end_time = 1\n"
        "  nl_max_its = 1\n"
        "[]\n"
        "\n"
        "[Outputs]\n"
        "  console = true\n"
        "  csv = false\n"
        "[]\n"
    )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.i', dir=workdir,
                                     delete=False) as f:
        f.write(input_content)
        input_path = f.name

    try:
        result = subprocess.run(
            [bin_path, '-i', os.path.basename(input_path)],
            cwd=workdir,
            capture_output=True, text=True, timeout=300
        )

        # Try to find RO2_SPECIES(N): name1,name2,... in stdout
        m = re.search(r'RO2_SPECIES\((\d+)\):\s*(\S+)', result.stdout)
        if m:
            return int(m.group(1)), m.group(2).split(',')

        # Fallback: search the full output (stdout + stderr) for species names
        full_output = result.stdout + '\n' + result.stderr
        # Look for "Detected N RO2 species" as a fallback indicator
        m2 = re.search(r'Detected\s+(\d+)\s+RO2\s+species', full_output)
        if m2:
            print(f"[parser found {m2.group(1)} RO2 species, but RO2_SPECIES "
                  "line not emitted — check MCMRO2ListPostprocessor setup]",
                  file=sys.stderr)
        else:
            print("[no RO2 detection output found — module may have failed "
                  "before parsing]", file=sys.stderr)
            print("=== STDOUT (last 1000 chars) ===", file=sys.stderr)
            print(result.stdout[-1000:], file=sys.stderr)
            print("=== STDERR (last 1000 chars) ===", file=sys.stderr)
            print(result.stderr[-1000:], file=sys.stderr)
        return 0, []
    finally:
        if os.path.exists(input_path):
            os.unlink(input_path)


def main():
    parser = argparse.ArgumentParser(
        description="Extract RO2 species from a .fac mechanism via the module")
    parser.add_argument("mechanism",
                        help="Path to .fac mechanism file")
    parser.add_argument("-o", "--output", default=None,
                        help="Output file (default: print to stdout, "
                        "workdir = current directory)")
    args = parser.parse_args()

    if args.output:
        workdir = os.path.dirname(os.path.abspath(args.output))
    else:
        workdir = os.getcwd()

    count, names = run_module_get_ro2(args.mechanism, DEFAULT_BINARY, workdir)

    # Sort alphabetically for reproducible output
    names = sorted(names)

    if args.output:
        with open(args.output, 'w') as f:
            for s in names:
                f.write(s + '\n')
        print(f"Wrote {count} RO2 species to {args.output}")
    else:
        # Aligned column output
        max_len = max((len(s) for s in names), default=0)
        cols = max(1, 80 // (max_len + 2))
        print(f"RO2 species ({count}):")
        for i, s in enumerate(names, 1):
            print(f"  {s:{max_len}}", end="")
            if i % cols == 0:
                print()
        if len(names) % cols:
            print()


if __name__ == "__main__":
    main()
