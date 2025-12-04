from __future__ import annotations

from ..core.base import TransferFunction, Signal, sat_mm


class AdiposeLipidUptake(TransferFunction):

    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        src_amt = self.source.get(self.nutrient)
        if src_amt <= 0.0:
            return 0.0

        insulin_signal = signals.get("insulin")
        insulin_level = 0.0 if insulin_signal is None else max(0.0, min(1.0, insulin_signal.level))
        vmax = float(self.params.get("vmax_gpm", 1.4)) * (1.0 + self.params.get("insulin_multiplier", 1.2) * insulin_level)
        km = float(self.params.get("km_g", 3.0))
        return sat_mm(src_amt, vmax, km)


class AdiposeLipidRelease(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        src_amt = self.source.get(self.nutrient)
        if src_amt <= 0.0:
            return 0.0

        insulin_signal = signals.get("insulin")
        insulin_level = 0.0 if insulin_signal is None else max(0.0, min(1.0, insulin_signal.level))
        vmax = float(self.params.get("vmax_gpm", 0.6)) / (1.0 + self.params.get("insulin_suppression", 4.0) * insulin_level)
        km = float(self.params.get("km_g", 5.0))
        return sat_mm(src_amt, vmax, km)


class AdiposeFFARelease(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        src_amt = max(self.source.get(self.nutrient), 0.0)
        if src_amt <= 0.0:
            return 0.0

        insulin_signal = signals.get("insulin")
        glucagon_signal = signals.get("glucagon")
        exercise_signal = signals.get("exercise")
        lipolysis_signal = signals.get("adipose_ffa_release_gpm")

        insulin_level = 0.0 if insulin_signal is None else max(0.0, float(insulin_signal.level))
        glucagon_level = 0.0 if glucagon_signal is None else max(0.0, float(glucagon_signal.level))
        exercise_level = 0.0 if exercise_signal is None else max(0.0, float(exercise_signal.level))
        requested = 0.0 if lipolysis_signal is None else max(0.0, float(lipolysis_signal.level))

        vmax_base = float(self.params.get("vmax_gpm", 0.6))
        km = float(self.params.get("km_g", 5.0))
        insulin_supp = float(self.params.get("insulin_suppression", 4.0))
        glucagon_gain = float(self.params.get("glucagon_stimulation", 1.2))
        exercise_gain = float(self.params.get("exercise_gain", 2.0))

        stimulation = 1.0
        stimulation += glucagon_gain * max(0.0, glucagon_level - 0.05)
        stimulation += exercise_gain * exercise_level
        stimulation = max(stimulation, 0.0)

        suppression = 1.0 + insulin_supp * max(0.0, insulin_level - 0.05)
        vmax = vmax_base * stimulation / max(suppression, 1e-9)

        rate = sat_mm(src_amt, vmax, km)
        dt_safe = max(dt_min, 1e-9)

        if requested > 0.0:
            rate = min(max(rate, requested), src_amt / dt_safe)
        else:
            rate = min(rate, src_amt / dt_safe)

        return rate
