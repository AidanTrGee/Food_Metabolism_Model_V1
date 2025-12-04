from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Sequence

from ..core.base import Compartment, TransferFunction, Controller, Signal

Number = float

@dataclass
class SimulationLog:
    times_min: List[Number] = field(default_factory=list)
    records: List[Dict[str, Number]] = field(default_factory=list)

class SimulationEngine:
    def __init__(
        self,
        compartments: List[Compartment],
        edges: List[TransferFunction],
        controllers: Sequence[Controller],
    ):
        self.compartments = compartments
        self.edges = edges
        self.controllers = list(controllers)
        self.signals: Dict[str, Signal] = {}
        self._comp_index: Dict[str, Compartment] = {c.name: c for c in self.compartments}

    def _state_lookup(self, comp: str, nutrient: str) -> Number:
        compartment = self._comp_index.get(comp)
        if compartment is not None:
            return compartment.get(nutrient)
        return 0.0

    def get_compartment(self, name: str) -> Compartment | None:
        return self._comp_index.get(name)

    def step(self, t_min: Number, dt_min: Number):
        for ctrl in self.controllers:
            ctrl.update(t_min, dt_min, self.signals, self._state_lookup)
        for c in self.compartments: c.clear_port()
        for e in self.edges:
            gpm = e.compute_flow_gpm(t_min, dt_min, self.signals)
            dt_safe = max(dt_min, 1e-9)
            src_amt = max(e.source.get(e.nutrient), 0.0)
            tgt_amt = max(e.target.get(e.nutrient), 0.0)
            max_out = src_amt / dt_safe
            min_in = -tgt_amt / dt_safe
            if gpm > max_out:
                gpm = max_out
            elif gpm < min_in:
                gpm = min_in
            e.source.propose_flow(e.nutrient, -gpm)
            e.target.propose_flow(e.nutrient, +gpm)
        for c in self.compartments:
            for nutrient, gpm in c.port.flows_gpm.items():
                c.pools[nutrient].amount_g += gpm * dt_min
            c.step_endogenous(t_min, dt_min, self.signals)
