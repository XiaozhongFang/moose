# MCMRO2ListPostprocessor

!syntax description /VectorPostprocessors/MCMRO2ListPostprocessor

## Description

`MCMRO2ListPostprocessor` is a `GeneralVectorPostprocessor` that outputs the RO2 (peroxy radical) species detected by the `MCMBoxModel` from the parsed mechanism file. It produces two vectors:

- **`ro2_count`** (1 element): total number of detected RO2 species
- **`ro2_species`** (N elements): 0-based species indices for each detected peroxy radical

These vectors are output to a separate CSV file (`<file_base>_<vpp_name>_<timestep>.csv`) in a row-per-element format. The gold file from a CSVDiff test verifies both the count and every individual species index.

## Underlying Detection

The RO2 species list is populated by `MCMFacsimileParser::parse()`:

1. **`.fac` files**: Parses the explicit `RO2 = A + B + C + ...` declaration line.
2. **`.kpp` files**: Parses the `RO2 = & ... )` section using KPP's variable naming convention (`C(ind_NAME)`).
3. **Fallback**: If neither is present, uses a name-based heuristic (`ends with "O2"` with exclusions for false positives like `HO2`, `NO2`, `SO2`, etc.).

The explicit declaration (method 1 or 2) is always preferred and is the recommended form for accurate RO2 detection.

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

The VPP output CSV has the following format (shown for a mechanism with 3 RO2 species at indices 7, 42, and 103):

```
ro2_count,ro2_species
3,7
0,42
0,103
```

The first row gives the total count, subsequent rows list each species index. The `0` in the `ro2_count` column for rows 2+ is the VPP format — the vector is repeated as one element per row.

The test at `test_ro2_detection.i` uses `MCMRO2ListPostprocessor` to validate that `MCMv331_Inorg_Isoprene.fac` (610 species, 1974 reactions) has exactly 117 RO2 species with the correct indices.

## Console Output

In addition to the VPP CSV output, `MCMRO2ListPostprocessor` prints a parseable line to the console:

```
RO2_SPECIES(117): CH3O2,NISOPO2,ISOP34O2,...
```

This is used by `check_ro2.py` to extract the actual RO2 species NAMES (not indices) from the module.

## Validation Workflow

```bash
# Gold CSV (indices) is maintained manually in test/tests/actions/gold/

# Verify the RO2 species NAMES via the module (ground truth):
python3 modules/atmospheric_chemistry/scripts/check_ro2.py \
    modules/atmospheric_chemistry/doc/content/.../MCMv331_Inorg_Isoprene.fac \
    -o ro2_species.txt
```

!syntax parameters /VectorPostprocessors/MCMRO2ListPostprocessor

!syntax inputs /VectorPostprocessors/MCMRO2ListPostprocessor

!syntax children /VectorPostprocessors/MCMRO2ListPostprocessor
