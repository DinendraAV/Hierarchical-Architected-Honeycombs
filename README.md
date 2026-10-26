# Honeycomb project: RHC vs HHC-2 under quasi-static and crush loading

Computational study of a regular hexagonal honeycomb core (RHC) against a two-level
hierarchical honeycomb (HHC-2) at equal relative density, bare core only (no skins),
extending the IAC 2025 paper. Supervisor: Prof. Luca Lomazzi (PoliMi DMEC), who runs the experiments.

**Status:** private repository. Do not make public or share data from `experiments/`
without Prof. Lomazzi's agreement.

## Scope (see docs/scope.md)
- Tier 1: periodic RVE elastic homogenisation (study 01)
- Tier 2: crush at two rates (studies 02-04)
- Tier 3, stretch: out-of-plane tension and a core-beam check (study 05)

## Layout
| Folder | Content |
|---|---|
| docs/ | scope, decisions log, meeting notes, IAC audit, references.bib |
| geometry/ | CAD (.step) and parameters for RHC and HHC-2, relative density table |
| materials/ | TPU parameters, calibration data and fits |
| studies/ | one folder per study: model/, inputs/, results/cases.csv, results/curves/ |
| experiments/ | only the data files agreed with Prof. Lomazzi |
| validation/ | numerical vs experimental comparisons |
| common/ | shared scripts (ODB to CSV, plotting) |
| paper/ | ICF16 abstract (30 Nov 2026), manuscript (target April 2027), figures |
| environment/ | solver and Python versions |

## Rules
- Every reported result has a row in the study's `results/cases.csv` and its `.inp` in `inputs/`.
- No PDFs of papers or datasheets in the repo; keep citations in `docs/references.bib`.
- Nothing from the author's employer belongs here.
- ODB files are not committed.

## Status log
See `docs/decisions.md`.
