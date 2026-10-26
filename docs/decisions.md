# Decisions log (newest last)

- 2026-10-18  Scope agreed with Prof. Lomazzi (see scope.md). He runs all experiments; numerical work is mine.
- 2026-10-21  Tips from Prof. Lomazzi: round the internal edges of the cells (fillets); TPU is fine for now;
              use LS-DYNA rather than ANSYS Explicit Dynamics; state the element types used; next assignment is
              in-plane quasi-static compression instead of out-of-plane.
- 2026-10-22  In-plane study built in Abaqus/Explicit 2020 (solid C3D8R, enhanced hourglass; platens as discrete rigid shells).
              Honeycomb patch imported from SolidWorks (.STEP) with fillets, uniform 1 mm walls, 10 mm cell across flats, 10 mm depth.
- 2026-10-23  Quick test at 1000 mm/s average speed: ALLKE/ALLIE about 130% at mid-ramp, so the force curve is
              inertia-dominated and invalid. Rule: ALLKE/ALLIE must stay below 5% over most of the run.
- 2026-10-24  Surface-intersection warnings traced to platens flush against the patch. Fixed with gaps
              (0.1 mm top, 0.02 mm bottom); strain is computed after subtracting the gaps.
- 2026-10-24  T = 0.01 s run (about 212 mm/s average) completed; peak ALLKE 0.35 mJ. ALLIE comparison pending.
- 2026-10-26  Next: T = 0.02 and 0.04 s runs for rate convergence; then set the full-length run from the converged speed.

## Lessons recorded
- Amplitude end time must equal the step time (a 1.0 s ramp in a 0.02 s step moved the platen about 1e-6 mm).
- Rigid bodies in Abaqus/Explicit need all six DOF of the reference point constrained in every step.
- Restricting field output to the honeycomb set removes the platens from the animation (history output is unaffected).
