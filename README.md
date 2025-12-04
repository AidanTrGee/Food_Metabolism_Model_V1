# Food Metabolism Reference Model

A physiologically grounded, simulation-ready implementation of the Food Metabolism Model described in the project assumptions document. The code tracks gastric emptying, intestinal absorption, hepatic buffering, circulatory mixing, peripheral utilisation, lymphatic fat handling, and renal glucose spill-over for a reference 90 kg adult.

Key design choices now visible in the codebase:

- **Compact compartment layout** – stomach lumen → intestinal lumen/tissue → portal vein → liver → pulmonary/arterial/venous blood loop with peripheral sinks (muscle, adipose), lymphatic fat delay, and kidney excretion.
- **Energy-driven gastric emptying** – stomach emptying depends on current caloric load with fat/protein and soluble fibre slow-down factors plus a post-meal lag.
- **Adaptive intestinal glucose transport** – SGLT1 (Michaelis–Menten), sigmoid apical GLUT2 recruitment, and capped paracellular flux, all modulated by fibre viscosity and optional insulin tweaks.
- **Hepatic buffering** – portal uptake via GLUT2 kinetics and basal hepatic glucose output that can be suppressed by insulin or stimulated by glucagon.
- **Peripheral sinks** – insulin-sensitive muscle and adipose glucose uptake plus lipid storage/slow release in adipose; lymphatic transit delays fat appearance in blood.
- **Renal threshold** – linear glucose excretion beyond ~10 mM plasma, with cumulative losses reported as a signal.

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate   # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt

# Run the reference scenario (writes CSV + optional plots to out/)
PYTHONPATH=src python -m Food_metabolism.scripts.run_scenario --plot

# Rebuild plots from the latest simulation output
PYTHONPATH=src python -m Food_metabolism.scripts.plot_results
```

Simulation behaviour is configured via `config/params.yml` (meals, initial conditions, hormone toggles, and optional parameter overrides). Reaction-rate and kinetic constants now live in `config/reaction_rates.yml`; edit that file to retune gastric emptying, absorption kinetics, pathway Vmax/Km values, glycogen handling gains, etc. The loader merges both YAML files so `params.yml` can still override any subset of those rates if desired. Results are written to `out/sim_output.csv` and plotted in `out/plots/` when `--plot` is used.

Each compartment entry in `out/sim_output.csv` now includes both the pool size (in grams) and a corresponding `_mg_dl` column, computed from the compartment’s instantaneous volume. All CLI plots and the web UI leverage these concentration traces so you can compare tissues on a common scale.

## Run the web interface

```bash
# 1. Backend API (runs on http://localhost:8000)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
PYTHONPATH=src uvicorn backend.main:app --host 0.0.0.0 --port 8000

# 2. Frontend (runs on http://localhost:5173 by default)
cd frontend
npm install
VITE_API_URL=http://localhost:8000 npm run dev
```

Open the browser at the address printed by Vite (usually <http://localhost:5173>). The React app proxies API requests to the FastAPI backend at `VITE_API_URL`; change the environment variable above if you run the backend on a different host or port.

