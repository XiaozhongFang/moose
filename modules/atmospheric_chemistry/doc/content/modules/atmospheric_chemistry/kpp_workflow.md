# KPP Workflow for Atmospheric Chemistry

This page documents the supported Kinetic PreProcessor (KPP) workflow for
atmospheric chemistry. It covers the bundled `small_strato` test
mechanism, switching to another KPP mechanism file, and using a user-defined
mechanism.

The module supports two KPP runtime paths:

- Box mode: `[AtmosphericChemistry/Box]` loads a mechanism-specific shared
  library and advances the chemical state with the KPP generated integrator.
  Supported values are `chem_solver = kpp_rosenbrock`, `kpp_sdirk`, and
  `kpp_runge_kutta`.
- Coupled mode: `[AtmosphericChemistry/Coupled]` loads the same generated
  shared library, evaluates KPP `Fun` and `Jac_SP` at quadrature points, and
  lets the MOOSE `[Executioner]` perform the fully coupled implicit solve.
  KPP's sparse analytical Jacobian is used for the chemical source Jacobian.

In box mode, MOOSE executioner settings advance only the outer time loop; the
chemical substep is owned by the selected KPP integrator. In coupled mode, the
MOOSE executioner owns the nonlinear and linear solve.

## Supported Scope

The KPP integration in this module is a runtime workflow for generated C
mechanisms. It is not a replacement for the full upstream KPP command-line
toolchain. The supported generated-library surface is:

- variable species metadata from `SPC_NAMES` and `NVAR`
- box integration through the adapter symbol `kpp_integrate`
- RHS evaluation through KPP `Fun`
- sparse analytical Jacobian evaluation through KPP `Jac_SP`, `LU_IROW`,
  `LU_ICOL`, and `LU_CROW`
- common physical globals when exported by the mechanism, including `TEMP`,
  `TIME`, `SUN`, `AIR`, and `H2O`
- common fixed species names `M`, `AIR`, `O2`, `N2`, and `H2O`

Per-reaction rates, separated species production/loss rates, and RO2 metadata
are not exported through a stable KPP C API. Those diagnostics remain MCM/FAC
features unless the mechanism-specific KPP metadata is added separately.

## Runtime Architecture

The box KPP path has two separate stages:

```
KPP mechanism files
  root.kpp, root.spc, root.eqn, included files
        |
        v
kpp/build/Makefile
  runs kpp, compiles generated C sources, links adapter/kpp_adapter.c
        |
        v
mechanism directory/kpp_build_<root>/libkpp_<root>.so
        |
        v
KppBoxIntegrator
  dlopen() shared library
  kpp_init()
  kpp_get_nvar()
  kpp_get_species_name(i)
  kpp_set_conc()
  kpp_integrate()
  kpp_get_conc()
        |
        v
MCMBoxModel::execute()
  reads scalar concentrations, integrates one time step, writes scalars back
```

Coupled mode uses the same build stage, then evaluates the generated RHS and
Jacobian inside MOOSE:

```
libkpp_<root>.so
        |
        v
KPPGeneratedMechanism
  dlopen() shared library
  Fun(VAR, FIX, RCONST, Ydot)
  Jac_SP(VAR, FIX, RCONST, JVS)
        |
        v
KPPMechanismMaterial
  stores kpp_rhs and kpp_jacobian_dense material properties
        |
        v
KPPChemicalSourceKernel + TimeDerivative
  assembled by the MOOSE nonlinear solve
```

Both paths get variable species names directly from the generated KPP library.
The order is KPP's active variable order, not the textual order a separate
parser happens to read from the `.spc` file. This matters for mechanisms where
KPP reorders variables.

## Mechanism File Layout

KPP uses one root input file, usually named `root.kpp`. The upstream KPP manual
calls the file stem `ROOT`; KPP-generated source names, monitor arrays, and
driver files are derived from that root.

A minimal mechanism used by this module should follow this layout:

```text
my_mechanism/
  my_mechanism.kpp
  my_mechanism.spc
  my_mechanism.eqn
  atoms.kpp or other include files
```

The root `.kpp` file should select C output and one supported integrator:

```text
#MODEL      my_mechanism
#LANGUAGE   C
#INTEGRATOR rosenbrock
```

Supported `#INTEGRATOR` values are:

| KPP `#INTEGRATOR` | MOOSE `chem_solver` | Runtime entry point |
|-------------------|---------------------|---------------------|
| `rosenbrock` | `kpp_rosenbrock` | direct `Rosenbrock(...)` call |
| `sdirk` | `kpp_sdirk` | generated `INTEGRATE(TIN, TOUT)` |
| `runge_kutta` | `kpp_runge_kutta` | generated `INTEGRATE(TIN, TOUT)` |

For coupled mode, the KPP integrator choice only affects the generated library
build. The coupled MOOSE solve uses KPP `Fun` and `Jac_SP`, not KPP's time
integrator.

The species file usually declares atoms, variable species, and fixed species:

```text
#INCLUDE atoms.kpp

#DEFVAR
O3  = O + O + O;
NO  = N + O;
NO2 = N + O + O;

#DEFFIX
O2  = O + O;
M   = IGNORE;
```

Important conventions:

- `#DEFVAR` species are the active variables advanced by KPP. These become
  MOOSE scalar variables in box mode and FE variables in coupled mode.
- `#DEFFIX` species are fixed inside the KPP mechanism. They are initialized by
  KPP code and are not emitted as ordinary MOOSE scalar output unless separate
  output objects are added.
- `#INCLUDE` files are parsed by KPP as if they were part of the same input.
  KPP can search its own directories and the current build directory. The
  module build helper copies direct `#INCLUDE` files that are next to the
  `.kpp` file; keep custom include files beside the mechanism unless the
  upstream KPP search path is deliberately used.
- Keep the `.kpp` file stem, `#MODEL` name, and `.spc/.eqn` file stems aligned
  for user mechanisms. KPP supports more complex predefined-model search paths,
  but matching names make the generated shared library and MOOSE auto-discovery
  predictable.

## Build the Atmospheric Chemistry Module

From the repository root:

```bash
cd modules/atmospheric_chemistry
make -j4
```

The KPP backend also requires the `kpp` executable to be available on `PATH`.
If KPP is installed outside `PATH`, set `KPP_BIN` when building the mechanism
library:

```bash
make -f kpp/build/Makefile \
  MECH=test/tests/kpp/kpp_small_strato/small_strato.kpp \
  KPP_BIN=/path/to/kpp
```

## Build the Bundled `small_strato` Library

From the repository root:

```bash
make -f modules/atmospheric_chemistry/kpp/build/Makefile \
  MECH=modules/atmospheric_chemistry/test/tests/kpp/kpp_small_strato/small_strato.kpp
```

Or from `modules/atmospheric_chemistry`:

```bash
make -f kpp/build/Makefile \
  MECH=test/tests/kpp/kpp_small_strato/small_strato.kpp
```

The output library is:

```text
modules/atmospheric_chemistry/test/tests/kpp/kpp_small_strato/kpp_build_small_strato/libkpp_small_strato.so
```

Verify that the adapter symbols are present:

```bash
nm -D modules/atmospheric_chemistry/test/tests/kpp/kpp_small_strato/kpp_build_small_strato/libkpp_small_strato.so \
  | grep -E 'kpp_init|kpp_get_nvar|kpp_get_species_name|kpp_integrate'
```

## Run the Bundled Test

The MOOSE input is:

```text
modules/atmospheric_chemistry/test/tests/kpp/kpp_small_strato.i
```

The box input uses:

```moose
[AtmosphericChemistry]
  [Box]
    mechanism_file = 'kpp_small_strato/small_strato.kpp'
    temperature = 270.0
    chem_solver = kpp_rosenbrock
  []
[]
```

Run through the test harness from `modules/atmospheric_chemistry`:

```bash
./run_tests --heavy -C test/tests/kpp --re kpp_small_strato
```

The test specification builds the `small_strato` shared library before running
the CSV comparison. If you run the input manually, build the shared library
first with the command in the previous section.

Or run the input directly from the test directory:

```bash
cd modules/atmospheric_chemistry/test/tests/kpp
../../../atmospheric_chemistry-opt -i kpp_small_strato.i --no-gdb-backtrace
```

The test writes:

```text
modules/atmospheric_chemistry/test/tests/kpp/kpp_small_strato.csv
```

The gold file follows the original KPP `small_strato` box example:

- start time: `12 * 3600` seconds
- end time: `12 * 3600 + 3 * 24 * 3600` seconds
- time step: `0.25 * 3600` seconds
- temperature used by the KPP mechanism: `270 K`

The gold CSV includes `time_h`, `O2`, and `N` columns from the KPP reference
driver. The MOOSE scalar CSV contains MOOSE time in seconds and the active
`#DEFVAR` species. The test therefore ignores `time_h`, `O2`, and `N` during
CSV comparison.

The coupled smoke-test input is:

```text
modules/atmospheric_chemistry/test/tests/kpp/kpp_small_strato_coupled.i
```

It uses the same shared library but runs:

```moose
[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'kpp_small_strato/small_strato.kpp'
    temperature = 270.0
    press = 1013.0
    chem_solver = kpp_rosenbrock
  []
[]
```

This path assembles `TimeDerivative + KPPChemicalSourceKernel` for each KPP
variable species and uses the generated sparse Jacobian in the Newton Jacobian.

## Operator-Split FV Transport with KPP Chamber Chemistry

The chamber split-coupling smoke test uses a single mesh-wide `TransientMultiApp`
instead of one sub-application per grid point. The parent application solves a
small distributed finite-volume street-canyon flow and transports selected
species fields. The child application clones that same mesh and runs
`[AtmosphericChemistry/Coupled]` with the generated F0AM chamber KPP Rosenbrock
mechanism, so every local grid degree of freedom owns its own 610-species
chemical state.

The diagnostic CSV averages in the parent and child inputs are only output
checks. Coupling is done by `MultiAppCopyTransfer` on the species fields at
`timestep_end`; no spatial average is used to drive chemistry.

The checked inputs are:

```text
modules/atmospheric_chemistry/test/tests/vs_F0AM_tutorial5_split_fv.i
modules/atmospheric_chemistry/test/tests/vs_F0AM_tutorial5_split_sub.i
```

The child input uses the F0AM chamber BottomUp lamp spectrum:

```moose
[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'chamber/kpp_chamber/generated_mechanisms/chamber_mcm_rosenbrock/chamber_mcm_rosenbrock.kpp'
    chem_solver = kpp_rosenbrock
    photolysis_scheme = BOTTOMUP
    lamp_flux_file = 'ExampleLightFlux.txt'
    bottomup_data_dir = '../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup'
  []
[]
```

Run the heavy split-coupling check from the module directory:

```bash
cd modules/atmospheric_chemistry
./run_tests -C test/tests --heavy --re 'build_chamber_kpp_rosenbrock_inputs|build_chamber_kpp_rosenbrock_library|vs_F0AM_chamber_split_kpp_fv'
```

For manual debugging, first generate and build the KPP chamber library, then run
the parent input directly:

```bash
cd modules/atmospheric_chemistry/test/tests
python3 ../../scripts/benchmark_chamber_solvers.py \
    --solvers kpp_rosenbrock \
    --write-inputs-only \
    --output-dir kpp_chamber/solver_runs/split_input_check
make -f ../../kpp/build/Makefile \
    MECH=chamber/kpp_chamber/generated_mechanisms/chamber_mcm_rosenbrock/chamber_mcm_rosenbrock.kpp
../../atmospheric_chemistry-opt -i vs_F0AM_tutorial5_split_fv.i Outputs/exodus=false
```

KPP generated C mechanisms store runtime state in process-global arrays. The
wrapper serializes access inside a rank when evaluating `Fun` and `Jac_SP`, while
the parent mesh and cloned chemistry sub-application can still be decomposed
across MPI ranks for HPC runs.

For the split test, both parent and child inputs use GMRES with block-Jacobi and
rank-local LU sub-solves. This matches the operator-split assumption that
chemistry is local to each grid field after transport. A global parallel direct
solve couples unrelated chemistry blocks and can route the small test through
MUMPS/METIS ordering, which is slower and less robust for this local chemistry
workload.

## Use a Different KPP Mechanism

For another KPP mechanism, repeat the mechanism-library build with the new
`.kpp` file:

```bash
make -f modules/atmospheric_chemistry/kpp/build/Makefile \
  MECH=/absolute/or/relative/path/to/my_mechanism.kpp
```

By default the build helper writes the shared library next to the mechanism:

```text
path/to/kpp_build_my_mechanism/libkpp_my_mechanism.so
```

Then point the MOOSE input at the same `.kpp` file:

```moose
[AtmosphericChemistry]
  [Box]
    mechanism_file = 'path/to/my_mechanism.kpp'
      chem_solver = kpp_rosenbrock
  []
[]
```

For coupled mode, use the same mechanism path under `[Coupled]`:

```moose
[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'path/to/my_mechanism.kpp'
    chem_solver = kpp_rosenbrock
  []
[]
```

At runtime, `KppBoxIntegrator` auto-discovers the library using this convention:

```text
<mechanism_dir>/kpp_build_<mechanism_stem>/libkpp_<mechanism_stem>.so
```

For example, if the input says:

```moose
mechanism_file = 'mechanisms/cases/urban.kpp'
```

the default library path is:

```text
mechanisms/cases/kpp_build_urban/libkpp_urban.so
```

Use `KPP_LIB` only when the library is intentionally stored somewhere else:

```bash
KPP_LIB=/scratch/kpp/libkpp_urban.so \
  modules/atmospheric_chemistry/atmospheric_chemistry-opt -i urban_box.i
```

When running multiple mechanisms, build one shared library per `.kpp` file.
Avoid setting a global `KPP_LIB` while comparing mechanisms unless every run is
intended to use that same override.

## Add a User-Defined Mechanism

Use this checklist for a new mechanism maintained outside the test directory:

1. Put all mechanism files in one mechanism directory.

   ```text
   my_case/kpp/
     my_case.kpp
     my_case.spc
     my_case.eqn
     atoms.kpp
     custom_rates.kpp
   ```

2. Use a root `.kpp` file that selects C output and one supported integrator.

   ```text
   #MODEL      my_case
   #LANGUAGE   C
   #INTEGRATOR rosenbrock
   ```

3. Declare chemically integrated species in `#DEFVAR`.

   Every `#DEFVAR` species should have an initial condition in the MOOSE input,
   unless the simulation intentionally relies on a zero default.

4. Declare fixed species in `#DEFFIX`.

   Fixed species are owned by the KPP generated code. The wrapper updates
   common background fixed species named `M`, `AIR`, `O2`, `N2`, and `H2O`.
   Other fixed species remain at the values initialized by the KPP mechanism.

5. Build the shared library.

   ```bash
   make -f modules/atmospheric_chemistry/kpp/build/Makefile \
     MECH=my_case/kpp/my_case.kpp
   ```

6. Use the same `.kpp` path in the MOOSE input.

   ```moose
   [AtmosphericChemistry]
     [Box]
       mechanism_file = 'my_case/kpp/my_case.kpp'
       chem_solver = kpp_rosenbrock
       chem_solver_rtol = 1e-6
       chem_solver_atol = 1e-10
     []
   []
   ```

7. Match the Executioner times to the reference driver when validating.

   KPP examples often report time in hours but integrate in seconds. The MOOSE
   `[Executioner]` times are seconds unless the model explicitly defines another
   unit convention.

8. Rebuild the shared library whenever a `.kpp`, `.spc`, `.eqn`, included file,
   or the module KPP adapter changes.

For `#INTEGRATOR sdirk` or `#INTEGRATOR runge_kutta`, set the matching
`chem_solver` in box mode so the input documents the generated library being
used. The adapter calls the generated `INTEGRATE(TIN, TOUT)` entry point for
those integrators.

## Validation and Gold Files

When validating a custom mechanism, compare against an independent KPP driver or
another trusted reference before updating a MOOSE gold file.

Typical differences to account for:

- KPP reference drivers may output time in hours; MOOSE CSV output uses the
  MOOSE `time` value.
- KPP reference drivers may output fixed species and monitor species that are
  not MOOSE scalar variables.
- Very small stiff-chemistry concentrations can underflow differently across
  compilers and BLAS/LAPACK implementations. Use a physically meaningful
  absolute-zero threshold in CSV comparisons for near-zero species.

For the bundled test, direct CSV comparison can be reproduced with:

```bash
python3 python/mooseutils/csvdiff.py \
  modules/atmospheric_chemistry/test/tests/kpp/gold/kpp_small_strato.csv \
  modules/atmospheric_chemistry/test/tests/kpp/kpp_small_strato.csv \
  --ignore-fields time_h O2 N \
  --relative-tolerance 5.0e-2 \
  --abs-zero 1.0e-1
```

## Troubleshooting

### `cannot load KPP shared library`

Build the shared library for the same `.kpp` path used by the input, or set
`KPP_LIB` to the intended library. Check the default path convention:

```text
<mechanism_dir>/kpp_build_<mechanism_stem>/libkpp_<mechanism_stem>.so
```

### `failed to resolve KPP symbols`

The library was probably built with an old adapter or is not a module-compatible
KPP library. Rebuild it with `modules/atmospheric_chemistry/kpp/build/Makefile`
and verify `kpp_get_nvar`, `kpp_get_species_name`, and `kpp_integrate` with
`nm -D`.

### Concentrations stay constant

Make sure the library was rebuilt after adapter changes and that the input time
interval is nonzero. For `small_strato`, the output should include 900-second
steps from 43200 to 302400 seconds.

### CSV columns do not match the gold file

Check whether the reference file includes KPP monitor or fixed-species columns.
Ignore those columns in the test specification, or add explicit MOOSE outputs
for the additional quantities if they are part of the validation requirement.

### KPP cannot find an included file

Keep direct include files beside the root `.kpp` file, or configure the upstream
KPP search path deliberately. The module build helper copies direct `#INCLUDE`
files from the mechanism directory into the build directory.
