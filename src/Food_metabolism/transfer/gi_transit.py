from __future__ import annotations

from ..core.base import TransferFunction, Signal


class GastricEnergyEmptying(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        source_amount = self.source.get(self.nutrient)
        if source_amount <= 0.0:
            return 0.0

        rate_per_min = getattr(self.source, "emptying_rate_per_min", None)
        if rate_per_min is None:
            return 0.0

        k_ge = rate_per_min(t_min)
        multiplier = self.params.get("flow_multiplier", 1.0)
        return max(0.0, k_ge * multiplier * source_amount)
