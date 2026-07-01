# MCMReactionRatePostprocessor

!syntax description /Postprocessors/MCMReactionRatePostprocessor

Reports the rate of a single chemical reaction:
\begin{equation}
R_r = k_r \cdot \prod_{\text{reactants}} [C_i]
\end{equation}

in molec/cm³/s (or ppb/s if `units = ppb`). Analogous to F0AM `ExtractRates.m` and AtChem2 `productionRates.output`/`lossRates.output`.

!syntax parameters /Postprocessors/MCMReactionRatePostprocessor

!syntax inputs /Postprocessors/MCMReactionRatePostprocessor

!syntax children /Postprocessors/MCMReactionRatePostprocessor

# MCMSpeciesLossRatePostprocessor

!syntax description /Postprocessors/MCMSpeciesLossRatePostprocessor

Reports the total loss rate for a species — sum of rates of all reactions that consume the species:
\begin{equation}
L_s = \sum_{r: f_{s,r} < 0} |f_{s,r}| \cdot R_r
\end{equation}

where $f_{s,r}$ is the stoichiometric coefficient of species $s$ in reaction $r$.

!syntax parameters /Postprocessors/MCMSpeciesLossRatePostprocessor

!syntax inputs /Postprocessors/MCMSpeciesLossRatePostprocessor

!syntax children /Postprocessors/MCMSpeciesLossRatePostprocessor

# MCMSpeciesProductionRatePostprocessor

!syntax description /Postprocessors/MCMSpeciesProductionRatePostprocessor

Reports the total production rate for a species — sum of rates of all reactions that produce the species:
\begin{equation}
P_s = \sum_{r: f_{s,r} > 0} f_{s,r} \cdot R_r
\end{equation}

!syntax parameters /Postprocessors/MCMSpeciesProductionRatePostprocessor

!syntax inputs /Postprocessors/MCMSpeciesProductionRatePostprocessor

!syntax children /Postprocessors/MCMSpeciesProductionRatePostprocessor

# MCMLifetimePostprocessor

!syntax description /Postprocessors/MCMLifetimePostprocessor

Computes the chemical lifetime of a species:
\begin{equation}
\tau_s = \frac{[C_s]}{L_s}
\end{equation}

where $L_s$ is the total loss rate. Returns $10^{12}$ s as a sentinel if the loss rate is near zero (effectively infinite lifetime). Equivalent to F0AM `lifetime.m`.

!syntax parameters /Postprocessors/MCMLifetimePostprocessor

!syntax inputs /Postprocessors/MCMLifetimePostprocessor

!syntax children /Postprocessors/MCMLifetimePostprocessor
