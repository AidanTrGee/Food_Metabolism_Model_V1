from __future__ import annotations

from ..core.base import TransferFunction, Signal, sat_mm


class HepaticGlucoseUptake(TransferFunction):

    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        amount = self.source.get(self.nutrient)
        if amount <= 0.0:
            return 0.0

        insulin_signal = signals.get("insulin")
        insulin_level = 0.0 if insulin_signal is None else max(0.0, min(1.0, insulin_signal.level))
        vmax = float(self.params.get("vmax_gpm", 1.2)) * (1.0 + self.params.get("insulin_multiplier", 0.8) * insulin_level)
        km = float(self.params.get("km_g", 10.0))
        return sat_mm(amount, vmax, km)


class HepaticGlucoseOutput(TransferFunction):

    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        liver_store = self.source.get(self.nutrient)
        if liver_store <= 0.0:
            return 0.0

        insulin_signal = signals.get("insulin")
        glucagon_signal = signals.get("glucagon")
        insulin_level = 0.0 if insulin_signal is None else max(0.0, min(1.0, insulin_signal.level))
        glucagon_level = 0.0 if glucagon_signal is None else max(0.0, min(1.0, glucagon_signal.level))

        basal_output = float(self.params.get("basal_output_gpm", 0.18))
        insulin_supp = float(self.params.get("insulin_suppression", 0.6))
        glucagon_stim = float(self.params.get("glucagon_stimulation", 0.8))
        multiplier = max(0.0, 1.0 - insulin_supp * insulin_level) + glucagon_stim * glucagon_level
        multiplier = min(multiplier, float(self.params.get("max_multiplier", 3.0)))
        flow = basal_output * multiplier
        return min(liver_store / max(dt_min, 1e-6), flow)
