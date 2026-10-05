# PILOT

**PILOT** — Phase field dislocation dynamics Input, Lattice & Optimization Toolkit.

PILOT combines several codes into one workflow for Phase Field Dislocation
Dynamics (PFDD) modeling of plasticity in bcc metals and alloys:

- **USFE / lattice generation** — builds alloy lattices and assigns each site an
  unstable stacking fault energy (USFE) from its local composition.
- **CRSS sweep** — locates the critical resolved shear stress (CRSS) as a
  function of temperature.
- **Stress-state visualization** — shows the resolved shear on each slip plane
  for a chosen applied stress.

Together these feed the inputs to, and interpret the results of, PFDD.