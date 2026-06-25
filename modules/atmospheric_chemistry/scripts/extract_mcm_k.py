#!/usr/bin/env python3
"""
Extract MCM standard rate constant definitions from F0AM MCMv331_K.m
and inject them into a .fac mechanism file.

Usage:
  python3 extract_mcm_k.py <MCMv331_K.m> <mechanism.fac> [--output combined.fac]

Workflow:
  1. Parse all rate constant definition blocks from MCMv331_K.m
  2. Extract intermediate variables + final expression for each constant
  3. Convert to fparser-compatible format (.* -> *, .^ -> ^, T -> TEMP)
  4. Scan mechanism.fac for actually referenced K names
  5. Inject only the referenced K definitions

Examples:
  # Inject K constants in-place (overwrites mechanism.fac)
  python3 extract_mcm_k.py MCMv331_K.m my_mechanism.fac

  # Write to a separate output file, keeping the original unchanged
  python3 extract_mcm_k.py MCMv331_K.m my_mechanism.fac -o combined.fac
"""

import re, sys, argparse
from collections import OrderedDict

def parse_mcm_k(k_file):
    """Parse MCMv331_K.m, return dict of {Kname: [definition_lines]}."""
    
    with open(k_file) as f:
        lines = f.readlines()
    
    result = OrderedDict()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        # Detect K definition block start: "Knames{i} = 'KXXX';"
        m = re.match(r"Knames\{i\}\s*=\s*'(\w+)'\s*;", line)
        if m:
            kname = m.group(1)
            body_lines = []
            i += 1
            # Collect intermediate variable definitions until krx(:,i) = ...
            while i < len(lines):
                l = lines[i].strip()
                if re.match(r'krx\(:,i\)\s*=', l):
                    # This is the final expression line
                    final_match = re.search(r'krx\(:,i\)\s*=\s*(.+?)\s*;', l)
                    if final_match:
                        final_expr = final_match.group(1).strip()
                        # Convert MATLAB syntax
                        final_expr = convert_matlab_expr(final_expr)
                        body_lines.append(f"{kname} = {final_expr} ;")
                    break
                # Skip blank lines and pure comment lines
                if l and not l.startswith('%') and not l.startswith('%%'):
                    # Intermediate variable definition: "VAR = EXPR ;"
                    vm = re.match(r'(\w+)\s*=\s*(.+?)\s*;', l)
                    if vm:
                        vname = vm.group(1)
                        vexpr = vm.group(2).strip()
                        vexpr = convert_matlab_expr(vexpr)
                        body_lines.append(f"{vname} = {vexpr} ;")
                i += 1
            if body_lines:
                result[kname] = body_lines
        i += 1
    
    return result


def convert_matlab_expr(expr):
    """Convert MATLAB expression to fparser-compatible format."""
    # .* -> * , ./ -> / , .^ -> ^
    expr = expr.replace('.*', '*').replace('./', '/').replace('.^', '^')
    # T -> TEMP (MCM .fac convention uses TEMP for temperature)
    # Use word-boundary replacement to avoid matching T inside other names
    expr = re.sub(r'\bT\b', 'TEMP', expr)
    # M, H2O stay as-is (already in _func_params)
    return expr


def find_referenced_k(fac_file, all_k_names):
    """Detect K constants referenced in .fac, or return all if --all flag set.
    
    Args:
        fac_file: Path to .fac mechanism file
        all_k_names: If non-empty, return this set directly (--all mode)
    """
    if all_k_names:
        return all_k_names

    referenced = set()
    with open(fac_file) as f:
        content = f.read()

    # Reaction lines: % KXXX : or % EXPR*KXXX*... :
    for m in re.finditer(r'%\s+([^:]+)\s*:', content):
        for token in re.findall(r'\b(K[A-Za-z0-9]+)\b', m.group(1)):
            referenced.add(token)

    # Coefficient definitions and reaction expressions
    for m in re.finditer(r'=\s*([^;]+)\s*;', content):
        for token in re.findall(r'\b(K[A-Za-z0-9]+)\b', m.group(1)):
            referenced.add(token)

    return referenced


def main():
    parser = argparse.ArgumentParser(
        description='Extract K constant definitions from MCMv331_K.m and inject into .fac file'
    )
    parser.add_argument('k_file', help='Path to MCMv331_K.m')
    parser.add_argument('fac_file', help='Target .fac mechanism file')
    parser.add_argument('--all', action='store_true',
                        help='Inject ALL 33 MCM K constants (recommended for complete .fac files)')
    parser.add_argument('--output', '-o', help='Output file (default: overwrite the .fac file)')
    args = parser.parse_args()
    
    k_defs = parse_mcm_k(args.k_file)
    print(f"Parsed {len(k_defs)} K constants from {args.k_file}")
    
    all_k = set(k_defs.keys()) if args.all else set()
    referenced = find_referenced_k(args.fac_file, all_k)
    
    matched = referenced & set(k_defs.keys())
    missing = referenced - set(k_defs.keys())
    unused = set(k_defs.keys()) - referenced
    
    if args.all:
        print(f"Injecting ALL {len(matched)} K constants (--all)")
    else:
        print(f"Mechanism references: {sorted(referenced)}")
        print(f"Matched: {len(matched)}, Unreferenced: {len(unused)}")
    if missing:
        print(f"WARNING: not found in MCMv331_K.m: {sorted(missing)}")

    if not matched:
        print("No K definitions to inject")
        return
    
    # Build the K definition block
    k_block_lines = ["", "* MCM v3.3.1 standard rate constants (from MCMv331_K.m) ;",
                     "* Auto-generated by extract_mcm_k.py ;", ""]
    
    for kname in sorted(matched):
        for line in k_defs[kname]:
            k_block_lines.append(line)
        k_block_lines.append("")  # blank line between constants
    
    k_block = '\n'.join(k_block_lines)
    
    # Read the original .fac, remove any previously injected old K definition block
    with open(args.fac_file) as f:
        fac_content = f.read()

    # Remove old auto-generated K definition block
    fac_content = re.sub(
        r'\n\*\s*MCM v3\.3\.1 standard rate constants.*?(?=\n\* MCMv3|$|\Z)',
        '', fac_content, flags=re.DOTALL
    )

    # Insert before the VARIABLE line
    var_match = re.search(r'\nVARIABLE\b', fac_content)
    if var_match:
        insert_pos = var_match.start()
        out_content = fac_content[:insert_pos] + k_block + '\n' + fac_content[insert_pos:]
    else:
        out_content = fac_content.rstrip() + '\n' + k_block + '\n'
    
    out_path = args.output or args.fac_file
    with open(out_path, 'w') as f:
        f.write(out_content)
    
    print(f"Written: {out_path}")

if __name__ == '__main__':
    main()
