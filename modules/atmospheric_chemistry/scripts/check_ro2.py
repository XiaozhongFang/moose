#!/usr/bin/env python3
"""Compare detected RO2 species with MCM reference list.

Usage:
    python3 check_ro2.py <ref_file> [detected_file]

    ref_file:     path to peroxy-radicals reference (e.g., peroxy-radicals_v3.3.1.dat)
    detected_file: path to ro2_detected.txt (default: ro2_detected.txt)
"""
import argparse, sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Compare RO2 species lists")
    parser.add_argument("ref_file", help="Path to peroxy-radicals reference file")
    parser.add_argument("detected_file", nargs="?", default="ro2_detected.txt",
                        help="Path to detected RO2 file (default: ro2_detected.txt)")
    args = parser.parse_args()

    ref_path = Path(args.ref_file)
    det_path = Path(args.detected_file)

    if not ref_path.exists():
        print(f"RO2 reference file not found: {ref_path}")
        sys.exit(1)

    with open(ref_path) as f:
        ref_set = set(l.strip() for l in f if l.strip())

    if not det_path.exists():
        print(f"Detected file not found: {det_path}")
        sys.exit(1)

    with open(det_path) as f:
        det_set = set(l.strip() for l in f if l.strip())

    missing = det_set - ref_set
    extra = ref_set - det_set
    print(f"Detected: {len(det_set)}, Reference: {len(ref_set)}")
    print(f"Missing from reference: {len(missing)}")
    if missing:
        print("\n".join(sorted(missing)[:20]))
    print(f"\nIn reference but not detected: {len(extra)}")


if __name__ == "__main__":
    main()
