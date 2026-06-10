# MATLAB Scripts From Supplementary Data 4

These scripts were created from `Banerjee and Eng et al Supplementary Data 4.txt` for the PSP / indigoidine growth-coupling workflow.

## Requirements

- MATLAB R2017b or compatible
- COBRA Toolbox v3.0 for the FBA/FVA and production-envelope scripts
- CellNetAnalyzer 2018.2 and CPLEX 12.8 for the cMCS calculation template
- The `.mat` model/result files named in each script header

## Files

- `calculate_coupling_mcs_ppu_template.m` records the P. putida CNA cMCS parameters from section 1.1.
- `calculateSizeOfMCS.m` computes cut-set size for each metabolite.
- `testSetcMCS_yield.m` validates reaction-vs-gene KO phenotypes for one cMCS product target.
- `addIndigoidineRxn.m` adds indigoidine synthesis, pantetheine phosphorylation, diffusion, and exchange reactions.
- `setupPutidaKT2440Model.m` applies the P. putida KT2440 constraints used in the supplement.
- `putidaAromaticRxns.m` returns the aromatic reaction list constrained in the P. putida workflow.
- `run_cmcs_validation.m` reproduces the section 1.2 cMCS validation workflow.
- `run_host_indigoidine_comparison.m` reproduces the section 2.2 host comparison workflow.
- `run_indigoidine_substrate_screen.m` reproduces the section 2.3 substrate screen.
- `run_production_envelope_yield_space.m` reproduces the section 2.4 production envelope and yield-space calculation.

## Typical Setup

From MATLAB, add this folder to the path after initializing COBRA Toolbox:

```matlab
cd DDBproject
initCobraToolbox()
changeCobraSolver('gurobi', 'all')
addpath('/Users/dbanerjee/Documents/PSP/matlab')
```

Then run the relevant `run_*.m` script from the directory containing the expected `.mat` inputs, or update the `load(...)` paths in the script header.
