from __future__ import annotations

from ..core.base import TransferFunction, Signal, sat_mm


class IntestinalFatToLymph(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        amount = self.source.get(self.nutrient)
        if amount <= 0.0:
            return 0.0
        k = float(self.params.get("k_per_min", 0.12))
        return max(0.0, k * amount)


class LymphToVenousFat(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        vmax = float(self.params.get("vmax_gpm", 2.0))
        km = float(self.params.get("km_g", 10.0))
        amount = self.source.get(self.nutrient)
        return sat_mm(amount, vmax, km)
