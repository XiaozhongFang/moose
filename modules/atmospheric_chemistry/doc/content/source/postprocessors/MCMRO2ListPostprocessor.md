# MCMRO2ListPostprocessor

!syntax description /VectorPostprocessors/MCMRO2ListPostprocessor

## Description

`MCMRO2ListPostprocessor` is a `GeneralVectorPostprocessor` that outputs the RO2 (peroxy radical) species detected by the `MCMBoxModel` from the parsed mechanism file. It produces one count vector and one flag vector per detected species name:

- **`ro2_count`** (1 element): total number of detected RO2 species
- **`<species name>`** (1 element): `1` when that species was detected as RO2

These vectors are output to a separate CSV file (`<file_base>_<vpp_name>_<timestep>.csv`). The gold file from a CSVDiff test verifies the count and the detected species-name set by comparing CSV headers.

## Underlying Detection

The RO2 species list is populated by `MCMFacsimileParser::parse()`:

1. **`.fac` files**: Parses the explicit `RO2 = A + B + C + ...` declaration line.
2. **`.kpp` files**: Parses the `RO2 = & ... )` section using KPP's variable naming convention (`C(ind_NAME)`).
3. **Fallback**: If neither is present, uses a name-based heuristic (`ends with "O2"` or contains `RO2`, with exclusions for false positives like `HO2`, `NO2`, `SO2`, etc.).

The explicit declaration (method 1 or 2) is always preferred and is the recommended form for accurate RO2 detection. The fallback keeps simplified mechanisms such as `tutorial_5sp.fac`, where `RO2` is declared directly as a mechanism species, working without a separate family declaration.

## Example

```moose
[VectorPostprocessors]
  [ro2_list]
    type = MCMRO2ListPostprocessor
    box_model = box_model
  []
[]
```

## Gold CSV Format

The VPP output CSV has the following format (shown for a mechanism with 3 RO2 species):

```
ro2_count,CH3O2,NISOPO2,ISOP34O2
3,1,1,1
```

Each species column is a presence flag. This makes CSVDiff compare species names rather than internal species indices, which can change when mechanism species ordering changes.

The test at `test_ro2_detection.i` uses `MCMRO2ListPostprocessor` to validate that `MCMv331_Inorg_Isoprene.fac` (610 species, 1974 reactions) has exactly 117 RO2 species from its explicit `RO2 = ... ;` declaration. `test_ro2_detection_tutorial_5sp.i` validates the fallback path for a simplified mechanism that declares `RO2` directly as a species.

## Console Output

In addition to the VPP CSV output, `MCMRO2ListPostprocessor` prints a parseable line to the console:

```
RO2_SPECIES(117): CH3O2,NISOPO2,ISOP34O2,...
```

This is used by `check_ro2.py --run-app` to compare the parser output against the explicit mechanism declaration or inspect fallback detection for mechanisms without an explicit RO2 family.

## Validation Workflow

```bash
# Regenerate the checked-in CSVDiff gold from the explicit mechanism declaration:
python3 modules/atmospheric_chemistry/scripts/check_ro2.py \
    modules/atmospheric_chemistry/doc/content/.../MCMv331_Inorg_Isoprene.fac \
    --gold-csv modules/atmospheric_chemistry/test/tests/ro2_family/gold/test_ro2_detection_ro2_list_0001.csv

# Compare the parser output with the explicit declaration:
python3 modules/atmospheric_chemistry/scripts/check_ro2.py \
    modules/atmospheric_chemistry/doc/content/.../MCMv331_Inorg_Isoprene.fac \
    --run-app

# Inspect fallback detection for a simplified mechanism with direct RO2 species:
python3 modules/atmospheric_chemistry/scripts/check_ro2.py \
    modules/atmospheric_chemistry/doc/content/.../tutorial_5sp.fac \
    --run-app
```

!syntax parameters /VectorPostprocessors/MCMRO2ListPostprocessor

!syntax inputs /VectorPostprocessors/MCMRO2ListPostprocessor

!syntax children /VectorPostprocessors/MCMRO2ListPostprocessor
