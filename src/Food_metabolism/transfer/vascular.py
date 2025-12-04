from __future__ import annotations

from ..core.base import Signal, TransferFunction


class FirstOrderTransfer(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        amount = self.source.get(self.nutrient)
        if amount <= 0.0:
            return 0.0
        k_per_min = self.params.get("k_per_min", 1.0)
        return max(0.0, k_per_min * amount)
