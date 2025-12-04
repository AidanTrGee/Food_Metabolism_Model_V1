from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Iterable, Optional
from abc import ABC, abstractmethod

Number = float

@dataclass
class Signal:
    name: str
    level: Number = 0.0

@dataclass
class NutrientPool:
    name: str
    amount_g: Number = 0.0

@dataclass
class Port:
    name: str
    flows_gpm: Dict[str, Number] = field(default_factory=dict)

class Compartment(ABC):
    def __init__(
        self,
        name: str,
        nutrients: Iterable[str],
        params: Optional[Dict[str, Number]] = None,
        volume_L: Optional[Number] = None,
    ):
        self.name = name
        self.pools: Dict[str, NutrientPool] = {n: NutrientPool(n, 0.0) for n in nutrients}
        self.params: Dict[str, Number] = params or {}
        self.port = Port(name=f"{name}:port")
        self.volume_L: Optional[Number] = volume_L

    def add(self, nutrient: str, grams: Number): self.pools[nutrient].amount_g += grams
    def get(self, nutrient: str) -> Number: return self.pools[nutrient].amount_g
    def add_many(self, additions: Dict[str, Number]):
        for nutrient, grams in additions.items():
            if nutrient in self.pools:
                self.pools[nutrient].amount_g += grams
            else:
                raise KeyError(f"{self.name} has no pool for nutrient '{nutrient}'")
    def propose_flow(self, nutrient: str, gpm: Number):
        self.port.flows_gpm[nutrient] = self.port.flows_gpm.get(nutrient, 0.0) + gpm
    def clear_port(self): self.port.flows_gpm.clear()

    def set_volume(self, volume_L: Number): self.volume_L = volume_L

    def get_concentration(self, nutrient: str) -> Number:
        if not self.volume_L or self.volume_L <= 0:
            return 0.0
        return self.get(nutrient) / self.volume_L

    @abstractmethod
    def step_endogenous(self, t_min: Number, dt_min: Number, signals: Dict[str, Signal]): ...

class TransferFunction(ABC):
    def __init__(self, name: str, source: Compartment, target: Compartment, nutrient: str, params: Optional[Dict[str, Number]] = None):
        self.name, self.source, self.target, self.nutrient = name, source, target, nutrient
        self.params: Dict[str, Number] = params or {}

    @abstractmethod
    def compute_flow_gpm(self, t_min: Number, dt_min: Number, signals: Dict[str, Signal]) -> Number: ...

class Controller(ABC):
    def __init__(self, name: str): self.name = name
    @abstractmethod
    def update(self, t_min: Number, dt_min: Number, signals: Dict[str, Signal], state_lookup): ...

def sat_mm(x: Number, Vmax: Number, Km: Number, multiplier: Number = 1.0) -> Number:
    return multiplier * (Vmax * x / (Km + x + 1e-9))
