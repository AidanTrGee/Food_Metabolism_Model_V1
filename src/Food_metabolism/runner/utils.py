from __future__ import annotations

import csv
from typing import Dict, Iterable, List, Tuple

from ..core.base import Compartment
from ..core.units import grams_to_mg_dl


def capture_compartment_state(
    compartments: Iterable[Compartment],
) -> Tuple[
    Dict[str, float],
    Dict[str, Dict[str, float]],
    Dict[str, Dict[str, float]],
    Dict[str, float],
]:

    flat: Dict[str, float] = {}
    nested: Dict[str, Dict[str, float]] = {}
    concentrations: Dict[str, Dict[str, float]] = {}
    volumes: Dict[str, float] = {}

    for comp in compartments:
        volume = float(getattr(comp, "volume_L", 0.0) or 0.0)
        volumes[comp.name] = volume

        comp_amounts: Dict[str, float] = {}
        comp_conc: Dict[str, float] = {}

        for nutrient, pool in comp.pools.items():
            grams = float(pool.amount_g)
            key = f"{comp.name}.{nutrient}"
            flat[key] = grams
            comp_amounts[nutrient] = grams

            conc = grams_to_mg_dl(grams, volume)
            conc_key = f"{key}_mg_dl"
            flat[conc_key] = conc
            comp_conc[nutrient] = conc

        flat[f"{comp.name}.volume_L"] = volume
        nested[comp.name] = comp_amounts
        concentrations[comp.name] = comp_conc

    return flat, nested, concentrations, volumes


def write_csv(times_min: List[float], records: List[Dict[str, float]], path: str) -> None:
    keys = {"time_min"}
    for record in records:
        keys |= set(record.keys())
    columns = list(sorted(keys))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for t, record in zip(times_min, records):
            row = {k: record.get(k, 0.0) for k in columns}
            row["time_min"] = t
            writer.writerow(row)
