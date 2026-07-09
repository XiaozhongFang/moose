#!/usr/bin/env python3
"""Convert a FACSIMILE MCM mechanism to a minimal KPP C mechanism.

This converter is intentionally narrow: it targets the MCM/AtChem2-style
FACSIMILE files used by the chamber benchmark. It preserves species order,
RO2, rate definitions, and reaction equations closely enough to build a KPP
shared library for runtime benchmarking.

Examples
--------
Convert the chamber FAC mechanism to a KPP Rosenbrock model:

    python3 scripts/fac_to_kpp.py \
        doc/content/modules/atmospheric_chemistry/database/MCMv331_Inorg_Isoprene.fac \
        --output-dir test/tests/chamber/kpp_chamber/generated_mechanisms/chamber_mcm_rosenbrock \
        --model chamber_mcm_rosenbrock \
        --integrator rosenbrock

Build the generated KPP shared library:

    make -f kpp/build/Makefile \
        MECH=test/tests/chamber/kpp_chamber/generated_mechanisms/chamber_mcm_rosenbrock/chamber_mcm_rosenbrock.kpp
"""

import argparse
import re
from pathlib import Path


BASE_FIXED = ("M", "O2", "N2", "H2O")
KPP_DUMMY_SPECIES = ("ONE", "PROD", "HV")
PHYSICAL_GLOBALS = ("M", "O2", "N2", "H2O")
BASE_RATE_NAMES = {
    "TEMP",
    "M",
    "O2",
    "N2",
    "H2O",
    "RO2",
}


def is_identifier_char(ch):
    return ch.isalnum() or ch == "_"


def find_left_operand(expr, op_pos):
    i = op_pos - 1
    while i >= 0 and expr[i].isspace():
        i -= 1
    if i < 0:
        raise ValueError(f"missing left operand before ** in {expr!r}")
    if expr[i] == ")":
        depth = 1
        i -= 1
        while i >= 0:
            if expr[i] == ")":
                depth += 1
            elif expr[i] == "(":
                depth -= 1
                if depth == 0:
                    return i, op_pos
            i -= 1
        raise ValueError(f"unbalanced parentheses before ** in {expr!r}")
    while i >= 0 and (is_identifier_char(expr[i]) or expr[i] in ".+-"):
        if expr[i] in "+-" and i > 0 and expr[i - 1] not in "Ee":
            break
        i -= 1
    return i + 1, op_pos


def find_right_operand(expr, op_end):
    i = op_end
    while i < len(expr) and expr[i].isspace():
        i += 1
    if i >= len(expr):
        raise ValueError(f"missing right operand after ** in {expr!r}")
    if expr[i] == "(":
        depth = 1
        i += 1
        while i < len(expr):
            if expr[i] == "(":
                depth += 1
            elif expr[i] == ")":
                depth -= 1
                if depth == 0:
                    return op_end, i + 1
            i += 1
        raise ValueError(f"unbalanced parentheses after ** in {expr!r}")
    if expr[i] in "+-":
        i += 1
    while i < len(expr) and (is_identifier_char(expr[i]) or expr[i] == "."):
        i += 1
    return op_end, i


def c_power_to_pow(expr):
    """Convert KPP/FAC '**' powers to C pow() calls for inline C code."""
    while "**" in expr:
        op_pos = expr.find("**")
        left_start, left_end = find_left_operand(expr, op_pos)
        right_start, right_end = find_right_operand(expr, op_pos + 2)
        left = expr[left_start:left_end].strip()
        right = expr[right_start:right_end].strip()
        expr = expr[:left_start] + f"pow({left}, {right})" + expr[right_end:]
    return expr


def clean_line(line):
    line = line.rstrip()
    if not line:
        return ""
    return line


def split_statements(lines):
    fixed = []
    in_reactions = False
    for raw in lines:
        line = clean_line(raw)
        stripped = line.strip()
        if not stripped:
            continue
        if "Reaction definitions" in stripped:
            in_reactions = True
        if in_reactions and stripped[0] not in "%*" and fixed:
            fixed[-1] = fixed[-1].rstrip() + " " + stripped
        else:
            fixed.append(line)

    statements = []
    buf = ""
    for line in fixed:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("*"):
            continue
        buf = (buf + " " + stripped).strip()
        while ";" in buf:
            stmt, buf = buf.split(";", 1)
            stmt = stmt.strip()
            if stmt:
                statements.append(stmt)
            buf = buf.strip()
    if buf:
        statements.append(buf)
    return statements


def convert_expr(expr, *, c_inline=False):
    expr = expr.strip()
    expr = re.sub(r"(\d)[dD]([+-]?\d+)", r"\1E\2", expr)
    expr = expr.replace("@", "**")
    expr = expr.replace("^", "**")
    expr = re.sub(r"\*\*(?!\()([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)", r"**(\1)", expr)
    expr = re.sub(r"\bEXP\s*\(", "EXP(", expr, flags=re.IGNORECASE)
    expr = re.sub(r"\bLOG10\s*\(", "LOG10(", expr, flags=re.IGNORECASE)
    expr = re.sub(r"J<([0-9]+)>", r"J\1", expr)
    if c_inline:
        expr = expr.replace("EXP(", "exp(").replace("LOG10(", "log10(")
        expr = c_power_to_pow(expr)
    return expr


def normalize_equation(equation):
    if "=" not in equation:
        return equation
    lhs, rhs = equation.split("=", 1)
    rhs_terms = [term.strip() for term in rhs.split("+") if term.strip()]
    if not rhs_terms:
        rhs_terms = ["PROD"]
    return lhs.strip() + " = " + " + ".join(rhs_terms)


def canonical_term(term):
    term = term.strip()
    if not term:
        return None
    match = re.fullmatch(
        r"(?:(\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\s*)?([A-Za-z][A-Za-z0-9_]*)",
        term,
    )
    if not match:
        return (term, 1.0)
    coeff = float(match.group(1)) if match.group(1) else 1.0
    return (match.group(2), coeff)


def canonical_side(side):
    terms = {}
    for raw in side.split("+"):
        parsed = canonical_term(raw)
        if parsed is None:
            continue
        name, coeff = parsed
        terms[name] = terms.get(name, 0.0) + coeff
    return tuple(sorted((name, coeff) for name, coeff in terms.items() if coeff != 0.0))


def canonical_equation_key(equation):
    lhs, rhs = equation.split("=", 1)
    return (canonical_side(lhs), canonical_side(rhs))


def equation_species(equation):
    names = set()
    if "=" not in equation:
        return names
    lhs, rhs = equation.split("=", 1)
    for side in (lhs, rhs):
        for raw in side.split("+"):
            parsed = canonical_term(raw)
            if parsed is None:
                continue
            name, _ = parsed
            if name not in KPP_DUMMY_SPECIES:
                names.add(name)
    return names


def merge_duplicate_reactions(reactions):
    merged = {}
    order = []
    for rate, equation in reactions:
        key = canonical_equation_key(equation)
        if key not in merged:
            order.append(key)
            merged[key] = [equation, []]
        merged[key][1].append(rate)

    out = []
    for key in order:
        equation, rates = merged[key]
        if len(rates) == 1:
            out.append((rates[0], equation))
        else:
            out.append(("+".join(f"({rate})" for rate in rates), equation))
    return out


def is_number(expr):
    return re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?", expr.strip()) is not None


def is_identifier(expr):
    return re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", expr.strip()) is not None


def parse_fac(path):
    statements = split_statements(path.read_text().splitlines())
    species = []
    rates = []
    ro2_species = []
    reactions = []

    for stmt in statements:
        s = stmt.strip()
        if not s:
            continue
        if s.startswith("CONSTANT "):
            s = s[len("CONSTANT "):].strip()
        if s.startswith("VARIABLE"):
            species.extend(s[len("VARIABLE"):].split())
            continue
        if s.startswith("%"):
            body = s[1:].strip()
            if ":" not in body:
                continue
            rate, equation = body.split(":", 1)
            reactions.append((convert_expr(rate, c_inline=True), normalize_equation(equation.strip())))
            continue
        if "=" in s:
            name, expr = s.split("=", 1)
            name = name.strip()
            expr = expr.strip()
            if name == "RO2":
                ro2_species = [term.strip() for term in expr.split("+") if term.strip()]
            else:
                rates.append((name, convert_expr(expr, c_inline=True)))

    # KPP fixed/dummy species are not active ODE variables. KPP also drops
    # declared species that never appear in the reaction graph, so filter them
    # here to keep #DEFVAR aligned with generated NVAR/SPC_NAMES.
    ignored_species = set(BASE_FIXED) | set(KPP_DUMMY_SPECIES)
    active_species = set()
    for _, equation in reactions:
        active_species.update(equation_species(equation))
    defvar = [sp for sp in species if sp not in ignored_species and sp in active_species]
    return defvar, ro2_species, rates, merge_duplicate_reactions(reactions)


def species_atom_line(name):
    return f"{name} = IGNORE;"


def wrap_c_assignment(lhs, rhs, indent="  ", width=100):
    line = f"{indent}{lhs} = {rhs};"
    if len(line) <= width:
        return [line]
    out = []
    current = f"{indent}{lhs} = "
    for part in rhs.split(" + "):
        token = part if not current.strip().endswith("=") else part
        extra = token if current.strip().endswith("=") else " + " + token
        if len(current) + len(extra) + 1 > width:
            out.append(current + " \\")
            current = indent + "  " + token
        else:
            current += extra
    out.append(current + ";")
    return out


def write_kpp(out_dir, model, defvar, ro2_species, rates, reactions, integrator):
    out_dir.mkdir(parents=True, exist_ok=True)

    (out_dir / f"{model}.kpp").write_text(
        f"#INCLUDE    {model}.def\n"
        "#LANGUAGE   C\n"
        "#EQNTAGS    ON\n"
        "#STOICMAT   OFF\n"
        "#HESSIAN    OFF\n"
        "#MEX        OFF\n"
        f"#INTEGRATOR {integrator}\n"
    )
    (out_dir / f"{model}.def").write_text(
        f"#INCLUDE {model}.spc\n"
        f"#INCLUDE {model}.eqn\n"
        "\n"
        "#LOOKATALL\n"
    )

    spc_lines = ["#DEFVAR"]
    spc_lines.extend(species_atom_line(sp) for sp in defvar)
    spc_lines.extend(["", "#DEFFIX"])
    spc_lines.extend(species_atom_line(sp) for sp in BASE_FIXED)
    (out_dir / f"{model}.spc").write_text("\n".join(spc_lines) + "\n")

    j_refs = sorted(
        {int(j) for rate, _ in reactions for j in re.findall(r"\bJ([0-9]+)\b", rate)}
        | {int(j) for _, expr in rates for j in re.findall(r"\bJ([0-9]+)\b", expr)}
    )

    eqn_lines = [
        f"// {model} converted from FACSIMILE",
        "",
        "#INLINE C_GLOBAL",
    ]
    c_globals = ["RO2", *PHYSICAL_GLOBALS]
    for name, _ in rates:
        if name not in c_globals:
            c_globals.append(name)
    for j in j_refs:
        name = f"J{j}"
        if name not in c_globals:
            c_globals.append(name)
    reaction_rates = []
    for i, (rate, equation) in enumerate(reactions, start=1):
        if is_number(rate) or is_identifier(rate):
            reaction_rates.append((rate, equation, None, None))
        else:
            helper = f"KPP_RATE_{i}"
            reaction_rates.append((helper, equation, helper, rate))
            c_globals.append(helper)

    eqn_lines.extend(f"extern double {name};" for name in c_globals)
    eqn_lines.extend(["#ENDINLINE", "", "#INLINE C_DATA"])
    eqn_lines.extend(f"double {name};" for name in c_globals)
    eqn_lines.extend(["#ENDINLINE", "", "#INLINE C_RCONST"])

    ro2_terms = [f"VAR[ind_{sp}]" for sp in ro2_species if sp in defvar]
    eqn_lines.append("  RO2 = 0.0;")
    for term in ro2_terms:
        eqn_lines.append(f"  RO2 += {term};")
    for name, expr in rates:
        eqn_lines.extend(wrap_c_assignment(name, expr))
    for _, _, helper, expr in reaction_rates:
        if helper:
            eqn_lines.extend(wrap_c_assignment(helper, expr))
    eqn_lines.extend(["#ENDINLINE", "", "#EQUATIONS"])

    for i, (rate, equation, _, _) in enumerate(reaction_rates, start=1):
        eqn_lines.append(f"<{i}> {equation} : {rate} ;")

    (out_dir / f"{model}.eqn").write_text("\n".join(eqn_lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fac", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--model", default="chamber_mcm")
    parser.add_argument("--integrator", default="rosenbrock")
    args = parser.parse_args()

    defvar, ro2_species, rates, reactions = parse_fac(args.fac)
    write_kpp(args.output_dir, args.model, defvar, ro2_species, rates, reactions, args.integrator)
    print(f"Wrote {args.output_dir / (args.model + '.kpp')}")
    print(f"species={len(defvar)} ro2={len(ro2_species)} rates={len(rates)} reactions={len(reactions)}")


if __name__ == "__main__":
    main()
