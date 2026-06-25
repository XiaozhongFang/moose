#!/usr/bin/env python3
"""Compare RO2 species detected from a .fac mechanism with the MCM reference list.

Usage:
    # Compare against reference list using pre-extracted RO2 file:
    python3 check_ro2.py <ref_file> [detected_file]

    # Extract RO2 species directly from a .fac mechanism file:
    python3 check_ro2.py <ref_file> --fac <mechanism.fac> [-o ro2_detected.txt]

    ref_file:     path to peroxy-radicals reference (e.g., peroxy-radicals_v3.3.1.dat)
    detected_file: path to ro2_detected.txt (default: ro2_detected.txt)

The RO2 detection from .fac extracts all species in the VARIABLE block whose
names end with 'O2' (MCM peroxy radical naming convention).
"""
import argparse, sys, re
from pathlib import Path


def extract_ro2_from_fac(fac_file):
    """Extract RO2 species from the VARIABLE block of a .fac file."""
    with open(fac_file) as f:
        content = f.read()

    # Find VARIABLE block: starts with VARIABLE, ends at first '%' or lone ';'
    m = re.search(r'VARIABLE\s*\n(.*?)(?=\n\s*[%;]|\n\s*$)', content, re.DOTALL | re.IGNORECASE)
    if not m:
        print(f"ERROR: no VARIABLE block found in {fac_file}", file=sys.stderr)
        sys.exit(1)

    var_block = m.group(1)
    # Use same detection logic as MCMFacsimileParser (name-based heuristic):
    # - Species ending in "O2" (excluding known false positives)
    # - Species containing "RO2"
    non_ro2 = {'HO2', 'NO2', 'SO2', 'H2O2', 'O2', 'N2O2',
               'NO3', 'HNO3', 'CO2', 'CLO2', 'CL2O2', 'BRO2'}
    all_species = re.findall(r'\b(\w+)\b', var_block)
    ro2_set = set()
    for sp in all_species:
        if sp in non_ro2:
            continue
        if (len(sp) >= 3 and sp.endswith('O2')) or 'RO2' in sp:
            ro2_set.add(sp)
    return sorted(ro2_set)


def main():
    parser = argparse.ArgumentParser(
        description="Compare RO2 species lists (from .fac or pre-extracted file)")
    parser.add_argument("ref_file", help="Path to peroxy-radicals reference file")
    parser.add_argument("detected_file", nargs="?", default=None,
                        help="Path to detected RO2 file (default: ro2_detected.txt)")
    parser.add_argument("--fac", help="Extract RO2 directly from .fac mechanism file")
    parser.add_argument("-o", "--output", help="Save detected RO2 to file")
    args = parser.parse_args()

    ref_path = Path(args.ref_file)
    if not ref_path.exists():
        print(f"RO2 reference file not found: {ref_path}")
        sys.exit(1)

    with open(ref_path) as f:
        ref_set = set(l.strip() for l in f if l.strip())

    # Determine detected set
    if args.fac:
        det_set = set(extract_ro2_from_fac(args.fac))
        if args.output:
            with open(args.output, 'w') as f:
                for s in sorted(det_set):
                    f.write(s + '\n')
            print(f"Wrote {len(det_set)} RO2 species to {args.output}")
    elif args.detected_file:
        det_path = Path(args.detected_file)
        if not det_path.exists():
            print(f"Detected file not found: {det_path}")
            sys.exit(1)
        with open(det_path) as f:
            det_set = set(l.strip() for l in f if l.strip())
    else:
        det_path = Path("ro2_detected.txt")
        if not det_path.exists():
            print(f"Detected file not found: {det_path}")
            print("Use --fac <mechanism.fac> to extract RO2 from a .fac file.")
            sys.exit(1)
        with open(det_path) as f:
            det_set = set(l.strip() for l in f if l.strip())

    missing = det_set - ref_set
    extra = ref_set - det_set
    print(f"Detected: {len(det_set)}, Reference: {len(ref_set)}")
    print(f"Missing from reference: {len(missing)}")
    if missing:
        print("\n".join(sorted(missing)[:20]))
    print(f"In reference but not detected: {len(extra)}")


if __name__ == "__main__":
    main()
