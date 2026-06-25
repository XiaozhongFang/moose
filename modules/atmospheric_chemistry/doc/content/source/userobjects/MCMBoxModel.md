# MCMBoxModel

!syntax description /UserObjects/MCMBoxModel

`MCMBoxModel` is the central computational engine for 0-D atmospheric chemistry
box model simulations in MOOSE. It manages the chemical system matrices
(stoichiometric coefficients, reactant indices, rate constants) and provides
the core $\mathrm{d}C/\mathrm{d}t$ computation following F0AM's `dydt_eval.m`
algorithm. The stoichiometric matrix uses a pluggable `StoichMatrix` storage
backend, selectable via the [!param](/UserObjects/MCMBoxModel/stoich_format) parameter
(analogous to PETSc's `-mat_type`).

### Core Algorithm

The ODE system follows F0AM's matrix formulation:

\begin{equation}
\frac{\mathrm{d}C}{\mathrm{d}t} = \mathbf{f}^T \cdot \left( \mathbf{k} \odot
\prod \mathbf{C}[\mathbf{iG}] \right)
\end{equation}

where:
- $\mathbf{f}$ [nRx $\times$ nSp]: stoichiometric coefficient matrix (sparse, ~99.9% zeros for full MCM)
- $\mathbf{iG}$ [nRx $\times$ 3]: reactant index triplets (padded with ONE=0 for pseudo-first-order)
- $\mathbf{k}$ [nRx]: rate constant vector

### Stoichiometric Matrix Storage (`stoich_format`)

The stoichiometric matrix $\mathbf{f}$ is stored in a `StoichMatrix` backend that
supports four formats, selected via [!param](/UserObjects/MCMBoxModel/stoich_format):

| Format | Storage | Memory (full MCM) | Best For |
|--------|---------|-------------------|----------|
| `CSR` (default) | Compressed Sparse Row, net coefficients | ~3 MB | HPC large mechanisms, PETSc AIJ-compatible |
| `COO` | AtChem2-style split reactant/product vectors | ~5 MB | Loss/production rate diagnostics |
| `DENSE` | Dense 2D array `f[reaction][species]` | ~800 MB | Tiny mechanisms (< 50 species), zero indirection overhead |
| `CSC` | Compressed Sparse Column + CSR forward index | ~6 MB | Column queries ("which reactions involve species X?") |

All formats expose an identical `forEachInRow(r, fn)` iteration interface — compute
methods (`computeDCdt`, Jacobian assembly) are format-agnostic. The compiler fully
inlines the `switch`-based dispatch, so format selection has zero runtime overhead.

The CSR format mirrors PETSc's internal AIJ storage and is the recommended choice for
HPC production runs. COO mirrors AtChem2's `clhs`/`crhs` arrays.

### Photolysis Methods

Two photolysis rate calculation methods are supported:

| Method | Description | Input Required |
|--------|-------------|---------------|
| `MCM_SZA` (default) | MCM parameterization $J = C_L \cos^C_{MM}(\theta) \exp(-C_{NN}/\cos\theta) \cdot F_{JFAC}$ | Solar zenith angle $\theta$ |
| `HYBRID` | TUV-based 4D lookup table (F0AM Hybrid J) | SZA, albedo, O$_3$ column, altitude |
| `BOTTOMUP` | Cross-section × quantum-yield × lamp-flux integration (F0AM Jmethod=1) | Lamp flux file, CS/QY data |

BOTTOMUP is the standard method for laboratory chamber experiments.
Data files are generated from an F0AM installation via `scripts/generate_bottomup_jmap.py`.
See [BottomUpJIntegrator](/utils/BottomUpJIntegrator.md) for details.

### Constrained Species (AtChem2 Mode)

Species can be marked as "constrained" (fixed to observations). Constrained species
still participate in reaction rate calculations but their concentrations are not
solved by the ODE system — reducing the effective problem size from $N$ to $N-C$.

### Solar Cycle

Multi-day simulations track the solar zenith angle using Madronich (1993)
parameterization, updating photolysis rates at each time step. Convergence
mode (`checkConvergence()`) enables diurnal steady-state detection.

### Air Density Calculation

If [!param](/UserObjects/MCMBoxModel/press) is set (> 0 mbar), the air number
density $M$ is computed dynamically from the ideal gas law
(AtChem2 `calcAirDensity` equivalent):

$$ M = 10^{-6} \cdot \frac{N_A}{R} \cdot \frac{p \times 100}{T} $$

where $N_A = 6.02214129\times 10^{23}$ mol⁻¹, $R = 8.3144621$ J mol⁻¹ K⁻¹.
If `press` is 0 (default), the fixed `air_density` parameter is used directly.

### Humidity Calculation

If [!param](/UserObjects/MCMBoxModel/rh) is set (≥ 0 %), the water vapor
concentration is computed from relative humidity using the Vaisala (2013)
formula (AtChem2 `convertRHtoH2O` equivalent):

$$ p_{\text{sat}} = 6.116441 \cdot 10^{\frac{7.591386(T-273.15)}{T-273.15+240.7263}} \quad\text{(mbar)} $$
$$ \text{H}_2\text{O} = \frac{\text{RH}/100 \cdot p_{\text{sat}}}{p - \text{RH}/100 \cdot p_{\text{sat}}} \cdot M $$

If `rh` is −1 (default, AtChem2 sentinel), the fixed `water_vapor` parameter
is used directly.

### Dilution (Box Mode)

If [!param](/UserObjects/MCMBoxModel/dilute) is set (> 0 s⁻¹), first-order
dilution is applied at each timestep:

$$ \frac{\mathrm{d}C_i}{\mathrm{d}t} \mathrel{-}= k_{\text{dil}} \cdot (C_i - C_i^{\text{bkgd}}) $$

This matches AtChem2's DILUTE parameter.  Background concentrations default
to zero (clean air dilution).

### Lazy Initialization

`MCMBoxModel` supports lazy initialization: if the MOOSE framework calls
`computeDCdt()` before `GeneralUserObject::initialize()`, the mechanism
is loaded on first use.  This ensures correct chemistry from the very
first time step regardless of framework initialization ordering.

### Key API Methods

| Method | Description |
|--------|-------------|
| `loadMechanism(mech)` | Populate matrices from parsed mechanism |
| `computeDCdt(C, dC)` | Compute $\mathrm{d}C/\mathrm{d}t$ |
| `computeJacobianTriplets(C, J)` | Analytical Jacobian as triplets |
| `enableHybridPhotolysis(dir)` | Switch to TUV lookup tables |
| `loadBottomUpData(dir, flux)` | Load lamp flux + reaction map for BOTTOMUP scheme |
| `updatePhotolysisBottomUp()` | Evaluate all J-values via CS×QY×Flux integration |
| `updatePhotolysis(sza, alb, o3, alt)` | Update photolysis rates |
| `setSolarCycle(lat, lon, d, m, y)` | Configure solar cycle |
| `setConstrainedSpecies(names)` | Mark constrained species |
| `setDilution(kdil, bg)` | Configure dilution |
| `reactionRate(r, C)` | Single reaction rate |
| `allReactionRates(C, rates)` | All reaction rates |

## Example Input File Syntax

```moose
[UserObjects]
  [box_model]
    type = MCMBoxModel
    mechanism_file = 'mechanism.fac'
  []
[]
```

!syntax parameters /UserObjects/MCMBoxModel

!syntax inputs /UserObjects/MCMBoxModel

!syntax children /UserObjects/MCMBoxModel
