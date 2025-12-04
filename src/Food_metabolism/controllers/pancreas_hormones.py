from __future__ import annotations

from ..core.base import Controller, Signal


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


class PancreasHormoneController(Controller):

    def __init__(
        self,
        blood_compartment,
        *,
        enabled: bool = True,
        target_mg_dl: float = 90.0,
        insulin_gain: float = 0.015,
        glucagon_gain: float = 0.01,
        basal_insulin: float = 0.05,
        basal_glucagon: float = 0.05,
        max_signal: float = 0.6,
        response_alpha: float = 0.2,
        insulin_threshold_mg_dl: float = 100.0,
        glucagon_threshold_mg_dl: float = 80.0,
        glucose_filter_alpha: float = 0.15,
        max_signal_slew_per_min: float | None = 0.25,
        suppression_floor: float = 0.0,
        co_release_tolerance: float = 0.01,
    ):
        super().__init__("PancreasHormoneController")
        self.blood = blood_compartment
        self.enabled = enabled
        self.target_mg_dl = target_mg_dl
        self.insulin_gain = insulin_gain
        self.glucagon_gain = glucagon_gain
        self.basal_insulin = basal_insulin
        self.basal_glucagon = basal_glucagon
        self.max_signal = max_signal
        self.response_alpha = _clamp(response_alpha, 0.0, 1.0)
        self.insulin_threshold_mg_dl = insulin_threshold_mg_dl
        self.glucagon_threshold_mg_dl = glucagon_threshold_mg_dl
        self.glucose_filter_alpha = _clamp(glucose_filter_alpha, 0.0, 1.0)
        self.filtered_glucose_mg_dl = target_mg_dl
        self.max_signal_slew_per_min = None if max_signal_slew_per_min is None else max(0.0, max_signal_slew_per_min)
        self.suppression_floor = _clamp(suppression_floor, 0.0, max_signal)
        self.co_release_tolerance = max(0.0, co_release_tolerance)

    def update(self, t_min, dt_min, signals, state_lookup):
        insulin_sig = signals.setdefault("insulin", Signal("insulin", self.basal_insulin))
        glucagon_sig = signals.setdefault("glucagon", Signal("glucagon", self.basal_glucagon))

        if not self.enabled:
            insulin_sig.level = self.basal_insulin
            glucagon_sig.level = self.basal_glucagon
            return

        volume = max(self.blood.volume_L or 0.0, 1e-6)
        glucose_g = self.blood.get("glucose")
        mg_dl = (glucose_g * 1000.0) / (volume * 10.0)
        if self.glucose_filter_alpha > 0.0:
            self.filtered_glucose_mg_dl += self.glucose_filter_alpha * (mg_dl - self.filtered_glucose_mg_dl)
            mg_dl_control = self.filtered_glucose_mg_dl
        else:
            mg_dl_control = mg_dl

        insulin_drive = max(0.0, mg_dl_control - self.insulin_threshold_mg_dl)
        glucagon_drive = max(0.0, self.glucagon_threshold_mg_dl - mg_dl_control)
        insulin_active = insulin_drive > 0.0
        glucagon_active = glucagon_drive > 0.0

        target_insulin = _clamp(
            self.basal_insulin + self.insulin_gain * insulin_drive,
            0.0,
            self.max_signal,
        )
        target_glucagon = _clamp(
            self.basal_glucagon + self.glucagon_gain * glucagon_drive,
            0.0,
            self.max_signal,
        )

        tol = self.co_release_tolerance
        if insulin_active and target_insulin > (self.basal_insulin + tol):
            target_glucagon = min(target_glucagon, self.suppression_floor)
        if glucagon_active and target_glucagon > (self.basal_glucagon + tol):
            target_insulin = min(target_insulin, self.suppression_floor)
\
        alpha = self.response_alpha
        prev_insulin = insulin_sig.level
        prev_glucagon = glucagon_sig.level
        new_insulin = prev_insulin + alpha * (target_insulin - prev_insulin)
        new_glucagon = prev_glucagon + alpha * (target_glucagon - prev_glucagon)

        max_slew = None
        if self.max_signal_slew_per_min:
            max_slew = self.max_signal_slew_per_min * max(dt_min, 1e-6)

        if max_slew is not None:
            insulin_delta = _clamp(new_insulin - prev_insulin, -max_slew, max_slew)
            glucagon_delta = _clamp(new_glucagon - prev_glucagon, -max_slew, max_slew)
            insulin_sig.level = prev_insulin + insulin_delta
            glucagon_sig.level = prev_glucagon + glucagon_delta
        else:
            insulin_sig.level = new_insulin
            glucagon_sig.level = new_glucagon

        insulin_sig.level = _clamp(insulin_sig.level, 0.0, self.max_signal)
        glucagon_sig.level = _clamp(glucagon_sig.level, 0.0, self.max_signal)

        tol = self.co_release_tolerance
        if insulin_sig.level > self.basal_insulin + tol:
            glucagon_sig.level = min(glucagon_sig.level, self.suppression_floor)
        if glucagon_sig.level > self.basal_glucagon + tol:
            insulin_sig.level = min(insulin_sig.level, self.suppression_floor)


