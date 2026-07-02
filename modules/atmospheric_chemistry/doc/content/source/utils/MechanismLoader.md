# MechanismLoader

!syntax description /Utils/MechanismLoader (standalone, not a MOOSE object)

`MechanismLoader` is a standalone utility that orchestrates the loading of an MCM
Facsimile-format (`.fac`) chemical mechanism. It encapsulates:

1. **`MCMFacsimileParser` invocation** — Parses the `.fac` file
2. **Photolysis file path resolution** — 3-tier search (as-is → input file directory → mechanism file directory)
3. **Full photolysis parameter set loading** — Reads all J\<N\> entries from the photolysis-rates file

The returned `MechanismData` struct provides all parsed mechanism data in a
self-contained format, replacing the ad-hoc parsing and path-resolution logic
previously embedded in `AtmosphericChemistryAction`.

## Motivation

Before `MechanismLoader`, the `AtmosphericChemistryAction` constructor performed
both path resolution and mechanism parsing inline, with identical photolysis path
resolution appearing twice (constructor + `actCoupledAddMaterial`). This led to:

- Code duplication (~45 lines of identical path resolution)
- Maintenance risk (bug fixes needed to be applied in two places)
- Constructor bloat (~160 lines for data loading + validation)

`MechanismLoader` centralizes all file I/O and parsing into a single call:

```cpp
MechanismData data = MechanismLoader::load(
    mechanism_file, photo_path, mcm_version, peroxy_path, input_file_dirs);
```

## API

### `MechanismLoader::load()`

```cpp
static MechanismData load(
    const std::string & mechanism_file,
    const std::string & photo_path,
    const std::string & mcm_version,
    const std::string & peroxy_path,
    const std::vector<std::string> & input_file_dirs);
```

### `MechanismData` struct

| Field | Type | Description |
|-------|------|-------------|
| `species` | `vector<string>` | Species names in order |
| `reactions` | `vector<Reaction>` | Each reaction has `rate_expression`, `reactants`, `products` |
| `stoichiometric_matrix` | `vector<vector<Real>>` | `stoichiometry[species_idx][reaction_idx]` |
| `rate_coefficients` | `map<string,string>` | Name → raw expression |
| `converted_coefficients` | `map<string,string>` | Name → fparser-converted expression |
| `eval_order` | `vector<string>` | Coefficient names in topological eval order |
| `reaction_rate_expressions` | `vector<string>` | Per-reaction rate expression |
| `ro2_species` | `vector<string>` | Peroxy radical species |
| `resolved_photo_path` | `string` | Resolved photolysis file path |
| `j_numbers_all` / `j_cl_values` / ... | `vector<...>` | Full photolysis parameter set |

## Usage

`MechanismLoader` is typically used by `AtmosphericChemistryAction` during construction.
It can also be used independently for mechanism analysis:

```cpp
MechanismData data = MechanismLoader::load(
    "mechanism.fac",
    "mcm_photolysis_rates_v3.3.1.dat",
    "v3.3.1",
    "mcm_peroxy_radicals_v3.3.1.dat",
    input_file_dirs);
```
