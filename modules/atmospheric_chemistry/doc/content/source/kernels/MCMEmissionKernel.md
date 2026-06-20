# MCMEmissionKernel

!syntax description /Kernels/MCMEmissionKernel

`MCMEmissionKernel` adds a time-dependent emission source term for atmospheric
chemistry species. This implements the emission terms from F0AM's
`Chem/Emission/` scripts (e.g., biogenic isoprene, soil NO, anthropogenic).

The weak form is:
\begin{equation}
R = -E(t, \mathbf{x}) \, \phi_{\text{test}}
\end{equation}
where $E(t, \mathbf{x})$ [molec/cm$^3$/s] is the emission rate provided by a
MOOSE [`Function`](Functions/index.md). The diagonal Jacobian contribution is zero
(the emission rate does not depend on concentration).

## Example Input File Syntax

```moose
[Functions]
  [isoprene_emission]
    type = ParsedFunction
    expression = 'if(t < 21600, 0, 5.0e6)'  # sunrise at 6am
  []
[]

[Kernels]
  [emit_c5h8]
    type = MCMEmissionKernel
    variable = C5H8
    function = isoprene_emission
  []
[]
```

!syntax parameters /Kernels/MCMEmissionKernel

!syntax inputs /Kernels/MCMEmissionKernel

!syntax children /Kernels/MCMEmissionKernel
