# TPU parameters (placeholder, not measured)

| Item | Value | Status |
|---|---|---|
| Model | neo-Hookean hyperelastic | placeholder; no rate dependence, no Mullins softening, no hysteresis |
| E bounds | 26 MPa and 67 MPa | bounds used in the study; confirm the exact source before citing |
| Poisson's ratio | 0.49 | assumed for a nearly incompressible elastomer; to be measured by DIC |
| Density | 1.22e-9 tonne/mm^3 | placeholder |

Abaqus inputs (units mm, tonne, s, N, MPa): C10 = mu/2 and D1 = 2/K, K = E/(3(1-2 nu)).
- E = 67 MPa, nu = 0.49: C10 = 11.242, D1 = 0.001791
- E = 26 MPa, nu = 0.49: C10 = 4.3624, D1 = 0.004615
- E = 67 MPa, nu = 0.45 (sensitivity case): C10 = 11.55, D1 = 0.008955

Next: coupon curves from Prof. Lomazzi at 5 and 500 mm/min (with one load-unload cycle) to calibrate
hyperelastic plus Prony terms. Add Mullins only if cyclic data show clear softening.

References (add to docs/references.bib): Ultimaker TPU 95A data sheet; Qi and Boyce on TPU behaviour; Gibson and Ashby, Cellular Solids.
