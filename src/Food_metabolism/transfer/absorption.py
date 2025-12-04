from __future__ import annotations

import math

from typing import Callable, cast

from ..core.base import Signal, TransferFunction, sat_mm

MOLAR_MASS_GLUCOSE = 180.156  # g/mol


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _call_float_fn(candidate: object, default: float = 0.0) -> float:
    if callable(candidate):
        return float(cast(Callable[[], float], candidate)())
    return default


class IntestinalGlucoseAbsorption(TransferFunction):

    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        lumen = self.source
        tissue = self.target

        conc_fn = getattr(lumen, "glucose_concentration_mM", None)
        conc_lumen = _call_float_fn(conc_fn)
        if conc_lumen <= 0.0:
            return 0.0

        tissue_fn = getattr(tissue, "glucose_concentration_mM", None)
        tissue_conc = _call_float_fn(tissue_fn)

        viscosity_fn = getattr(lumen, "fiber_viscosity_modifier", None)
        fiber_modifier = _call_float_fn(viscosity_fn, default=1.0)

        insulin_signal = signals.get("insulin")
        insulin_level = 0.0 if insulin_signal is None else _clamp(insulin_signal.level, 0.0, 1.0)

        params = self.params

        sglt1_vmax = float(params.get("sglt1_vmax_gpm", 1.8))
        sglt1_km = float(params.get("sglt1_km_mM", 5.0))
        sglt1_gain = float(params.get("insulin_sglt1_boost", 0.2))
        sglt1_flux = sglt1_vmax * fiber_modifier * conc_lumen / (sglt1_km + conc_lumen)
        sglt1_flux *= 1.0 + sglt1_gain * insulin_level

        c50 = float(params.get("glut2_activation_mM", 40.0))
        slope = max(float(params.get("glut2_switch_width_mM", 8.0)), 1e-3)
        availability = 1.0 / (1.0 + math.exp(-(conc_lumen - c50) / slope))
        glut2_vmax = float(params.get("glut2_vmax_gpm", 2.5))
        glut2_km = float(params.get("glut2_km_mM", 20.0))
        glut2_insulin_att = float(params.get("insulin_glut2_attenuation", 0.2))
        gradient = max(conc_lumen - tissue_conc, 0.0)
        glut2_flux = availability * glut2_vmax * gradient / (glut2_km + gradient)
        glut2_flux *= max(0.0, 1.0 - glut2_insulin_att * insulin_level)

        active_flux = max(0.0, sglt1_flux) + max(0.0, glut2_flux)

        para_vmax = float(params.get("paracellular_vmax_gpm", 0.4))
        para_km = float(params.get("paracellular_km_mM", 40.0))
        para_cap = _clamp(float(params.get("paracellular_fraction_cap", 0.15)), 0.0, 0.5)
        para_flux_candidate = para_vmax * conc_lumen / (para_km + conc_lumen)
        max_para = para_cap * active_flux / max(1e-9, 1.0 - para_cap)
        para_flux = min(para_flux_candidate, max_para)

        total_flux = active_flux + max(0.0, para_flux)
        return max(0.0, total_flux)


class IntestinalBasolateralExport(TransferFunction):
    def compute_flow_gpm(self, t_min, dt_min, signals: dict[str, Signal]) -> float:
        amount = self.source.get(self.nutrient)
        if amount <= 0.0:
            return 0.0

        tissue_fn = getattr(self.source, "glucose_concentration_mM", None)
        tissue_conc = _call_float_fn(tissue_fn)

        portal_volume = max(getattr(self.target, "volume_L", 1.0), 1e-6)
        portal_amount = self.target.get(self.nutrient)
        portal_conc = (portal_amount * 1000.0) / (MOLAR_MASS_GLUCOSE * portal_volume)
        gradient = max(tissue_conc - portal_conc, 0.0)

        vmax = float(self.params.get("vmax_gpm", 8.0))
        km_amount = float(self.params.get("km_g", 5.0))
        gradient_ref = float(self.params.get("gradient_ref_mM", 2.0))
        gradient_factor = gradient / (gradient + gradient_ref)
        return sat_mm(amount, vmax, km_amount) * gradient_factor
