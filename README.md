
## Quick start

```bash
python -m venv .venv && source .venv/bin/activate   # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt

# Run the reference scenario (writes CSV + optional plots to out/)
PYTHONPATH=src python -m Food_Metabolism.scripts.run_scenario --plot

# Rebuild plots from the latest simulation output
PYTHONPATH=src python -m Food_Metabolism.scripts.plot_results
```

Simulation behaviour is configured via `config/params.yml` (meals, initial conditions, hormone toggles, and optional parameter overrides). Reaction-rate and kinetic constants now live in `config/reaction_rates.yml`; edit that file to retune gastric emptying, absorption kinetics, pathway Vmax/Km values, glycogen handling gains, etc. The loader merges both YAML files so `params.yml` can still override any subset of those rates if desired. Results are written to `out/sim_output.csv` and plotted in `out/plots/` when `--plot` is used.

Each compartment entry in `out/sim_output.csv` now includes both the pool size (in grams) and a corresponding `_mg_dl` column, computed from the compartment’s instantaneous volume. CLI plots leverage these concentration traces so you can compare tissues on a common scale.


