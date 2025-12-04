from __future__ import annotations

from ..core.base import Compartment, Signal


class Lymph(Compartment):
    def __init__(self, params=None):
        params = params or {}
        super().__init__("Lymph", ["lipid"], params, volume_L=params.get("volume_L", 0.2))

    def step_endogenous(self, t_min, dt_min, signals):
        return
