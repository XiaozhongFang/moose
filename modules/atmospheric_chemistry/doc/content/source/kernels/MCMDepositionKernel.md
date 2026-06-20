# MCMDepositionKernel

!syntax description /Kernels/MCMDepositionKernel

`MCMDepositionKernel` implements first-order dry deposition loss for atmospheric
chemistry species. The deposition rate follows the standard resistance-in-series
formulation where the overall deposition velocity $v_d$ is converted to a
pseudo-first-order rate constant $k_{\text{dep}} = v_d / H$ ($H$ = boundary
layer height).

The weak form is:
\begin{equation}
R = k_{\text{dep}} \, u \, \phi_{\text{test}}
\end{equation}
with diagonal Jacobian:
\begin{equation}
\frac{\partial R}{\partial u} = k_{\text{dep}} \, \phi_{\text{trial}} \, \phi_{\text{test}}
\end{equation}

Typical deposition velocities for atmospheric species:
- O$_3$: 0.5–1.0 cm/s (over vegetation)
- HNO$_3$: 2–5 cm/s (highly soluble)
- NO$_2$: 0.1–0.5 cm/s
- SO$_2$: 0.5–1.5 cm/s
- H$_2$O$_2$: 1–2 cm/s

With a boundary layer height of ~1000 m, the corresponding $k_{\text{dep}}$
values are on the order of $10^{-5}$ to $10^{-4}$ s$^{-1}$.

## Example Input File Syntax

```moose
[Kernels]
  [o3_deposition]
    type = MCMDepositionKernel
    variable = O3
    rate = 5.0e-6   # /s (0.5 cm/s at H=1000m)
  []
  [hno3_deposition]
    type = MCMDepositionKernel
    variable = HNO3
    rate = 3.0e-5   # /s (3 cm/s at H=1000m)
  []
[]
```

!syntax parameters /Kernels/MCMDepositionKernel

!syntax inputs /Kernels/MCMDepositionKernel

!syntax children /Kernels/MCMDepositionKernel
