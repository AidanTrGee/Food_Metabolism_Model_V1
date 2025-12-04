from __future__ import annotations

from ..core.base import TransferFunction, Signal, sat_mm


class MuscleGlucoseUptake(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        src_amt = self.source.get(self.nutrient)
        if src_amt <= 0.0:
            return 0.0

        insulin_signal = signals.get("insulin")
        insulin_level = 0.0 if insulin_signal is None else max(0.0, min(1.0, insulin_signal.level))
        vmax_base = float(self.params.get("vmax_gpm", 3.0))
        basal_fraction = float(self.params.get("basal_fraction", 0.12))
        insulin_multiplier = float(self.params.get("insulin_multiplier", 3.0))
        vmax = vmax_base * (basal_fraction + insulin_multiplier * insulin_level)
        km = float(self.params.get("km_g", 3.0))
        return sat_mm(src_amt, vmax, km)


class AdiposeGlucoseUptake(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        src_amt = self.source.get(self.nutrient)
        if src_amt <= 0.0:
            return 0.0

        insulin_signal = signals.get("insulin")
        insulin_level = 0.0 if insulin_signal is None else max(0.0, min(1.0, insulin_signal.level))
        vmax_base = float(self.params.get("vmax_gpm", 0.8))
        basal_fraction = float(self.params.get("basal_fraction", 0.05))
        insulin_multiplier = float(self.params.get("insulin_multiplier", 2.5))
        vmax = vmax_base * (basal_fraction + insulin_multiplier * insulin_level)
        km = float(self.params.get("km_g", 3.5))
        return sat_mm(src_amt, vmax, km)
