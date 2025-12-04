from __future__ import annotations

from ..core.base import Signal, TransferFunction

MOLAR_MASS_GLUCOSE = 180.156  # g/mol


class RenalGlucoseExcretion(TransferFunction):

    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        blood = self.source
        amount = blood.get(self.nutrient)
        if amount <= 0.0:
            return 0.0

        volume = max(getattr(blood, "volume_L", 1.0), 1e-6)
        plasma_mM = (amount * 1000.0) / (MOLAR_MASS_GLUCOSE * volume)

        threshold = float(self.params.get("threshold_mM", 10.0))
        if plasma_mM <= threshold:
            return 0.0

        slope = float(self.params.get("slope_gpm_per_mM", 0.05))
        max_gpm = float(self.params.get("max_gpm", 3.0))
        excess = plasma_mM - threshold
        flow = slope * excess
        flow = min(flow, max_gpm)

        max_removal = amount / max(dt_min, 1e-6)
        return max(0.0, min(flow, max_removal))
