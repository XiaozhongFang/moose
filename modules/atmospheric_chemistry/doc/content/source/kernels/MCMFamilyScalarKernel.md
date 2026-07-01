# MCMFamilyScalarKernel

!syntax description /ScalarKernels/MCMFamilyScalarKernel

`MCMFamilyScalarKernel` implements the DAE family conservation residual for the slack variable of each chemical family. It replaces the standard `ChemistryODEKernel` for species designated as DAE slack variables.

The kernel computes:
\begin{equation}
R = \frac{\mathrm{d}u}{\mathrm{d}t} - \frac{\mathrm{d}C_{\text{slack}}}{\mathrm{d}t}_{\text{corrected}}
\end{equation}

where the corrected time derivative enforces $\mathrm{d}F_{\text{total}}/\mathrm{d}t = 0$.

This ScalarKernel is automatically created by `AtmosphericChemistryAction` when `family_names` parameters are provided in `mode = box`.

### Usage

!syntax parameters /ScalarKernels/MCMFamilyScalarKernel

!syntax inputs /ScalarKernels/MCMFamilyScalarKernel

!syntax children /ScalarKernels/MCMFamilyScalarKernel
