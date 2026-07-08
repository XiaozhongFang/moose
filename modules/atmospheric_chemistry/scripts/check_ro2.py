#!/usr/bin/env python3
"""Extract and validate RO2 species for atmospheric_chemistry mechanisms.

For FACSIMILE mechanisms with an explicit

    RO2 = A + B + C ;

declaration, this script reads that declaration directly from the mechanism
file. This is the preferred source for regenerating RO2 CSVDiff gold files.

Use --run-app to also run atmospheric_chemistry-opt and compare the parser's
detected species names against the explicit declaration. For simplified
mechanisms without an explicit RO2 declaration, --run-app reports the parser's
fallback result.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BINARY = SCRIPT_DIR.parent / "atmospheric_chemistry-opt"


def mechanism_path(path):
  return Path(os.path.expanduser(path)).resolve()


def strip_fac_comments(text):
  lines = []
  for line in text.splitlines():
    if line.lstrip().startswith("*"):
      continue
    lines.append(line)
  return "\n".join(lines)


def extract_explicit_ro2_species(fac_file):
  """Return RO2 species from the explicit FACSIMILE RO2 declaration."""
  path = mechanism_path(fac_file)
  text = strip_fac_comments(path.read_text())
  match = re.search(r"(?ims)^\s*RO2\s*=\s*(.*?)\s*;", text)
  if not match:
    return []

  names = []
  seen = set()
  duplicates = []
  for token in match.group(1).replace("&", " ").split("+"):
    name = token.strip()
    if not name:
      continue
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name):
      raise ValueError(f"invalid RO2 species token in {path}: {name!r}")
    if name in seen:
      duplicates.append(name)
    seen.add(name)
    names.append(name)

  if duplicates:
    dup_list = ", ".join(sorted(set(duplicates)))
    raise ValueError(f"duplicate RO2 species in {path}: {dup_list}")
  return names


def build_input(fac_file, workdir):
  rel_path = os.path.relpath(str(fac_file), str(workdir))
  return f"""[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
    mechanism_file = '{rel_path}'
    temperature = 298.0
    air_density = 2.46e19
    photolysis_scheme = MCM_SZA
    jfac = 1.0
  []
[]

[VectorPostprocessors]
  [ro2_list]
    type = MCMRO2ListPostprocessor
    box_model = box_model
  []
[]

[Executioner]
  type = Transient
  dt = 1
  end_time = 1
  nl_max_its = 1
[]

[Outputs]
  console = true
  csv = false
[]
"""


def run_module_get_ro2(fac_file, binary, workdir):
  """Run atmospheric_chemistry-opt and return detected RO2 species names."""
  fac_path = mechanism_path(fac_file)
  bin_path = mechanism_path(binary)
  workdir = Path(workdir).resolve()

  if not bin_path.exists():
    raise FileNotFoundError(f"binary not found: {bin_path}")
  if not fac_path.exists():
    raise FileNotFoundError(f"mechanism file not found: {fac_path}")
  workdir.mkdir(parents=True, exist_ok=True)

  with tempfile.NamedTemporaryFile(mode="w", suffix=".i", prefix="check_ro2_",
                                   dir=workdir, delete=False) as handle:
    handle.write(build_input(fac_path, workdir))
    input_path = Path(handle.name)

  try:
    result = subprocess.run(
        [str(bin_path), "-i", input_path.name],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if result.returncode != 0:
      raise RuntimeError(
          f"{bin_path.name} failed with exit code {result.returncode}\n"
          f"=== STDOUT ===\n{result.stdout[-2000:]}\n"
          f"=== STDERR ===\n{result.stderr[-2000:]}"
      )

    match = re.search(r"RO2_SPECIES\((\d+)\):\s*([^\r\n]*)", result.stdout)
    if not match:
      raise RuntimeError(
          "RO2_SPECIES line not found in application output\n"
          f"=== STDOUT ===\n{result.stdout[-2000:]}\n"
          f"=== STDERR ===\n{result.stderr[-2000:]}"
      )

    count = int(match.group(1))
    names_text = match.group(2).strip()
    names = [] if not names_text else [item.strip() for item in names_text.split(",") if item.strip()]
    if count != len(names):
      raise RuntimeError(f"RO2_SPECIES count {count} does not match {len(names)} parsed names")
    return names
  finally:
    input_path.unlink(missing_ok=True)


def compare_species(expected, actual):
  expected_set = set(expected)
  actual_set = set(actual)
  missing = sorted(expected_set - actual_set)
  extra = sorted(actual_set - expected_set)

  if not missing and not extra and len(expected) == len(actual):
    return True

  if missing:
    print("Missing in module output:", ", ".join(missing), file=sys.stderr)
  if extra:
    print("Extra in module output:", ", ".join(extra), file=sys.stderr)
  if len(expected) != len(actual):
    print(f"Count mismatch: expected {len(expected)}, got {len(actual)}", file=sys.stderr)
  return False


def write_species_list(path, names):
  with open(path, "w") as handle:
    for name in names:
      handle.write(name + "\n")


def write_gold_csv(path, names):
  with open(path, "w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["ro2_count", *names])
    writer.writerow([len(names), *([1] * len(names))])


def print_species(names):
  print(f"RO2 species ({len(names)}):")
  for name in names:
    print(name)


def main():
  parser = argparse.ArgumentParser(
      description="Extract RO2 species names from .fac mechanisms.")
  parser.add_argument("mechanism", help="Path to .fac mechanism file")
  parser.add_argument("-o", "--output", help="Write one RO2 species name per line")
  parser.add_argument("--gold-csv", help="Write CSVDiff gold CSV from the explicit RO2 declaration")
  parser.add_argument("--run-app", action="store_true",
                      help="Run atmospheric_chemistry-opt and compare/print detected RO2 species")
  parser.add_argument("--binary", default=str(DEFAULT_BINARY),
                      help=f"Path to atmospheric_chemistry-opt (default: {DEFAULT_BINARY})")
  parser.add_argument("--workdir", default=None,
                      help="Working directory for --run-app (default: output directory or current directory)")
  args = parser.parse_args()

  explicit_names = extract_explicit_ro2_species(args.mechanism)
  module_names = None

  if args.run_app:
    if args.workdir:
      workdir = args.workdir
    elif args.output:
      workdir = os.path.dirname(os.path.abspath(args.output))
    elif args.gold_csv:
      workdir = os.path.dirname(os.path.abspath(args.gold_csv))
    else:
      workdir = os.getcwd()
    module_names = run_module_get_ro2(args.mechanism, args.binary, workdir)
    if explicit_names and not compare_species(explicit_names, module_names):
      return 1

  if args.gold_csv:
    if not explicit_names:
      print("No explicit RO2 declaration found; refusing to write gold CSV", file=sys.stderr)
      return 1
    write_gold_csv(args.gold_csv, explicit_names)

  output_names = explicit_names or module_names
  if not output_names:
    print("No explicit RO2 declaration found. Use --run-app to inspect fallback detection.",
          file=sys.stderr)
    return 1

  if args.output:
    write_species_list(args.output, output_names)
  elif not args.gold_csv:
    print_species(output_names)

  if args.gold_csv:
    print(f"Wrote {len(explicit_names)} RO2 species to {args.gold_csv}")
  if args.output:
    print(f"Wrote {len(output_names)} RO2 species to {args.output}")
  if args.run_app and explicit_names:
    print(f"Module RO2 detection matched {len(explicit_names)} explicit species")

  return 0


if __name__ == "__main__":
  sys.exit(main())
