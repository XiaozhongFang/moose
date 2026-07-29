# MAS1998 benchmark support

This directory contains the first MOOSE-native pieces for the Spee et al. global
transport-chemistry benchmark:

- final published Table 1 ground-level concentrations
- recovered 45-reaction reduced methane KPP mechanism
- final published static vertical diffusivity K(z)
- solid-body rotation wind field
- HNO3 and NO cylinder initial conditions
- lowest-layer NO emission source

The implementation intentionally follows the 1998 journal version for the
benchmark constants. The MAS-R9702 preprint at
`/mnt/d/ZoteroData/storage/CGPEEMGK/04782D.pdf` is useful background, but the
journal paper explicitly says the benchmark differs from earlier preprints. In
particular, MAS-R9702 uses a time-dependent K profile and a 22.6 km top height,
while the journal benchmark uses the static K(z) expression and 38.2 km top
height.

The 45-reaction methane chemistry mechanism in `chemistry/mas1998_methane.*`
is a C/KPP port of the recovered Wayback Machine copy of the benchmark
`chemistry.k`. CWI report NM-R9505 (`/mnt/d/ZoteroData/storage/DQYC9G5U/05028D.pdf`)
Appendix B gives the detailed CIRK reaction table, pressure convention, and
auxiliary formulas, but it is the ODE-solver version: 46 reactions between 19
species with `NOx` and `Ox` lumped species. The global benchmark page is used
as the authoritative source for the final 45-reaction, 17-variable KPP form.

The erratum `/mnt/d/ZoteroData/storage/PKRAM83L/1-s2.0-1352231096832079-main.pdf`
corrects the 1996 ODE-solver paper's `NOx`/`Ox` lumping notation and Figure 3;
it does not change the 17-variable global benchmark mechanism. CWI report
NM-R9501 (`/mnt/d/ZoteroData/storage/Z83E5FVS/05042D.pdf`) is relevant for the
IMEX/splitting context, not for replacing the benchmark chemistry constants.

The original `Ref_Sol_Benchmark_Global.text` is still unavailable. These files
are enough to build the chemistry mechanism and verify the scalar benchmark
constants, but generating a replacement high-accuracy 3D reference solution
still requires the paper-specific spherical transport discretization and solver
input deck.

Generic coupled transport support now lives in `AtmosphericChemistryCoupled`.
`transport_velocity` accepts scalar velocity component variables, `diffusivity`
accepts a constant or MOOSE Function for ordinary diffusion, and
`density_weighted_diffusivity` plus `density_weighted_air_density` add the
atmospheric `d_z(rho K d_z(c/rho))` form used by this benchmark. The MAS1998
input deck should use these Action parameters for common species
advection/diffusion kernel creation instead of duplicating per-species kernels
in the benchmark.

The spherical finite-volume horizontal transport weights are available through
`AtmosphericSphericalFVTimeDerivative` and `AtmosphericSphericalFVAdvection`.
These objects encode the physical longitude-latitude cell area and face flux
metrics for the paper's conservative spherical form. `AtmosphericSphericalFVAdvection`
also provides the MAS1998 third-order upwind flux with the published limiter.
The remaining transport gap is reduced-grid connectivity, including a closed FV
longitude seam because MOOSE periodic boundary actions do not currently apply
to FV variables, and the full coupled reference-solution input deck.
