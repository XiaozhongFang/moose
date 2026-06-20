#!/usr/bin/env python3
"""Compare detected RO2 species with MCM reference list."""
import sys, os

ref = sys.argv[1] if len(sys.argv) > 1 else '/home/fangxiaozhong/git_repo/AtChem2/mcm/peroxy-radicals_v3.3.1'
detected = sys.argv[2] if len(sys.argv) > 2 else 'ro2_detected.txt'

if os.path.exists(ref):
    with open(ref) as f: ref_set = set(l.strip() for l in f if l.strip())
else:
    print(f'RO2 reference file not found: {ref}')
    sys.exit(1)

if os.path.exists(detected):
    with open(detected) as f: det_set = set(l.strip() for l in f if l.strip())
    missing = det_set - ref_set
    extra = ref_set - det_set
    print(f'Detected: {len(det_set)}, Reference: {len(ref_set)}')
    print(f'Missing from reference: {len(missing)}')
    if missing: print('\n'.join(sorted(missing)[:20]))
    print(f'\nIn reference but not detected: {len(extra)}')
else:
    print(f'Detected file not found: {detected}')
