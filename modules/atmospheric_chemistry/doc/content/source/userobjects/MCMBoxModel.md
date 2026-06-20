# MCMBoxModel

!syntax description /UserObjects/MCMBoxModel

`MCMBoxModel` is the central computational engine for 0-D atmospheric chemistry
box model simulations in MOOSE. It manages the chemical system matrices
(stoichiometric coefficients, reactant indices, rate constants) and provides
the core $\mathrm{d}C/\mathrm{d}t$ computation following F0AM's `dydt_eval.m`
algorithm.

### Core Algorithm

The ODE system follows F0AM's matrix formulation:

\begin{equation}
\frac{\mathrm{d}C}{\mathrm{d}t} = \mathbf{f}^T \cdot \left( \mathbf{k} \odot
\prod \mathbf{C}[\mathbf{iG}] \right)
\end{equation}

where:
- $\mathbf{f}$ [nRx $\times$ nSp]: stoichiometric coefficient matrix
- $\mathbf{iG}$ [nRx $\times$ 2]: reactant index pairs (-1 for pseudo-first-order)
- $\mathbf{k}$ [nRx]: rate constant vector

### Photolysis Methods

Two photolysis rate calculation methods are supported:

| Method | Description | Input Required |
|--------|-------------|---------------|
| `MCM_SZA` (default) | MCM parameterization $J = C_L \cos^C_{MM}(\theta) \exp(-C_{NN}/\cos\theta) \cdot F_{JFAC}$ | Solar zenith angle $\theta$ |
| `HYBRID` | TUV-based 4D lookup table (F0AM Hybrid J) | SZA, albedo, O$_3$ column, altitude |

### Constrained Species (AtChem2 Mode)

Species can be marked as "constrained" (fixed to observations). Constrained species
still participate in reaction rate calculations but their concentrations are not
solved by the ODE system — reducing the effective problem size from $N$ to $N-C$.

### Solar Cycle

Multi-day simulations track the solar zenith angle using Madronich (1993)
parameterization, updating photolysis rates at each time step. Convergence
mode (`checkConvergence()`) enables diurnal steady-state detection.

### Key API Methods

| Method | Description |
|--------|-------------|
| `loadMechanism(mech)` | Populate matrices from parsed mechanism |
| `computeDCdt(C, dC)` | Compute $\mathrm{d}C/\mathrm{d}t$ |
| `computeJacobianTriplets(C, J)` | Analytical Jacobian as triplets |
| `enableHybridPhotolysis(dir)` | Switch to TUV lookup tables |
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
