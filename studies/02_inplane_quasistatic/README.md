# Study 02: in-plane quasi-static compression (RHC patch)

Question: force-displacement response and first collapse of the RHC patch under in-plane compression, as a
benchmark for the experiments and the Gibson and Ashby in-plane estimates.

Setup (Abaqus/Explicit): C3D8R honeycomb (about 24,700 elements, 6 layers through the depth), neo-Hookean TPU placeholder,
discrete rigid shell platens, general contact (friction 0.3), semi-automatic mass scaling (target dt 2e-7 s).
Compression axis: Y. Strain = (|U2| - gaps) / 60.5 mm. Nominal stress = RF2 / (66.7 mm x 10 mm).

Acceptance: ALLKE/ALLIE below 5% over most of the run; the force curve must stop changing as the speed is lowered.

Workflow: edit the settings at the top of `model/setup_case.py`, run it in CAE (File > Run Script), save the .inp to `inputs/`,
export U2, RF2, ALLKE, ALLIE to `results/curves/`, and add a row to `results/cases.csv`.

`model/archive/fix_model.py` was a one-off patch script, kept for history.
