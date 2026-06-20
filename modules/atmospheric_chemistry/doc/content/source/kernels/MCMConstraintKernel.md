# MCMConstraintKernel

!syntax description /Kernels/MCMConstraintKernel

`MCMConstraintKernel` constrains a chemical species concentration to a prescribed
time-dependent function. This implements the "constrained species" mode from AtChem2,
where certain species (e.g., NO, NO$_2$) are fixed to observed values rather than
solved by the chemical ODE system.

The residual is simply:
\begin{equation}
R = u - f(t, \mathbf{x})
\end{equation}
where $f(t, \mathbf{x})$ is a MOOSE [`Function`](Functions/index.md) providing
the time-dependent constrained value. The Jacobian contribution is diagonal:
\begin{equation}
\frac{\partial R}{\partial u} = 1
\end{equation}

When used alongside [`MCMBoxModel`](MCMBoxModel.md), constrained species still
participate in the chemical reaction rate calculations (they consume/produce other
species), but their own concentrations do not evolve — the ODE solver is
effectively solving a reduced system of $N - C$ equations where $C$ is the
number of constrained species.

## Example Input File Syntax

!listing test/tests/actions/mcm_constrained.i block=Kernels

```moose
[Kernels]
  [no_constraint]
    type = MCMConstraintKernel
    variable = NO
    function = no_observed_function
  []
  [no2_constraint]
    type = MCMConstraintKernel
    variable = NO2
    function = no2_observed_function
  []
[]
```

!syntax parameters /Kernels/MCMConstraintKernel

!syntax inputs /Kernels/MCMConstraintKernel

!syntax children /Kernels/MCMConstraintKernel
