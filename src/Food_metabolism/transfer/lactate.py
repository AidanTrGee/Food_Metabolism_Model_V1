from __future__ import annotations

from ..core.base import Signal, TransferFunction, sat_mm


class MuscleLactateExport(TransferFunction):

    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        lactate = max(self.source.get(self.nutrient), 0.0)
        if lactate <= 0.0:
            return 0.0

        params = self.params or {}
        vmax_base = float(params.get("vmax_gpm", 0.5))
        km = float(params.get("km_g", 6.0))
        basal_fraction = float(params.get("basal_fraction", 0.08))
        exercise_gain = float(params.get("exercise_gain", 4.0))
        deficit_gain = float(params.get("deficit_gain", 1.5))
        deficit_scale = max(float(params.get("deficit_scale_mmol", 1.0)), 1e-9)
        min_reg = float(params.get("min_regulation", 0.05))
        max_reg = float(params.get("max_regulation", 6.0))

        exercise = 0.0
        if signals is not None:
            signal = signals.get("exercise")
            if signal is not None:
                exercise = max(0.0, float(signal.level))
        deficit = 0.0
        if signals is not None:
            deficit_sig = signals.get("muscle_atp_deficit_mmol")
            if deficit_sig is not None:
                deficit = max(0.0, float(deficit_sig.level)) / deficit_scale

        regulation = basal_fraction + exercise_gain * exercise + deficit_gain * deficit
        regulation = max(min_reg, min(max_reg, regulation))
        vmax = max(0.0, vmax_base) * regulation
        if vmax <= 0.0:
            return 0.0

        dt_safe = max(dt_min, 1e-9)
        rate = sat_mm(lactate, vmax, km)
        rate = min(rate, lactate / dt_safe)
        return rate


class HepaticLactateUptake(TransferFunction):

    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        amount = max(self.source.get(self.nutrient), 0.0)
        if amount <= 0.0:
            return 0.0

        params = self.params or {}
        vmax_base = float(params.get("vmax_gpm", 0.45))
        km = float(params.get("km_g", 10.0))
        insulin_gain = float(params.get("insulin_gain", 2.0))
        glucagon_supp = float(params.get("glucagon_suppression", 1.5))
        reference_insulin = float(params.get("reference_insulin", 0.05))
        reference_glucagon = float(params.get("reference_glucagon", 0.05))
        min_reg = float(params.get("min_regulation", 0.05))
        max_reg = float(params.get("max_regulation", 6.0))

        insulin = reference_insulin
        glucagon = reference_glucagon
        if signals is not None:
            insulin_sig = signals.get("insulin")
            glucagon_sig = signals.get("glucagon")
            if insulin_sig is not None:
                insulin = float(insulin_sig.level)
            if glucagon_sig is not None:
                glucagon = float(glucagon_sig.level)

        regulation = 1.0
        regulation += insulin_gain * max(0.0, insulin - reference_insulin)
        regulation *= max(0.0, 1.0 - glucagon_supp * max(0.0, glucagon - reference_glucagon))
        regulation = max(min_reg, min(max_reg, regulation))

        vmax = max(0.0, vmax_base) * regulation
        if vmax <= 0.0:
            return 0.0

        dt_safe = max(dt_min, 1e-9)
        rate = sat_mm(amount, vmax, km)
        rate = min(rate, amount / dt_safe)
        return rate
