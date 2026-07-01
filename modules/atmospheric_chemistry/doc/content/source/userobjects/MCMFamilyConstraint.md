# MCMFamilyConstraint

!syntax description /UserObjects/MCMFamilyConstraint

`MCMFamilyConstraint` defines chemical family conservation groups using the Differential-Algebraic Equation (DAE) method. Each family (e.g., `NOx = NO₂ + NO`, `Ox = O₃ + NO₂ + 2·NO₃`) has its total concentration conserved by designating the first member as a "DAE slack variable."

This system is the MOOSE equivalent of F0AM's family conservation mechanism in `Mass_eval.m` and `dydt_eval.m`.

### DAE Formulation

For a family with members $\{m_1, m_2, \dots, m_k\}$ and scaling factors $\{s_1, s_2, \dots, s_k\}$:

\begin{equation}
F_{\text{total}} = \sum_{i=1}^k s_i \cdot [m_i]
\end{equation}

The DAE slack variable $m_1$ absorbs the algebraic constraint:

\begin{equation}
\frac{\mathrm{d}[m_1]}{\mathrm{d}t} = -\frac{1}{s_1} \sum_{i=2}^k s_i \cdot \frac{\mathrm{d}[m_i]}{\mathrm{d}t}
\end{equation}

This enforces $\mathrm{d}F_{\text{total}}/\mathrm{d}t = 0$ while allowing non-slack members to evolve normally via the chemical ODE system.

### Usage

!syntax parameters /UserObjects/MCMFamilyConstraint

!syntax inputs /UserObjects/MCMFamilyConstraint

!syntax children /UserObjects/MCMFamilyConstraint
