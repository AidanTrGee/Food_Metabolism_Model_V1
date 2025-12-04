from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping

from ..core.base import Signal, sat_mm

Number = float


@dataclass(frozen=True)
class GlycolysisResult:

    rate_gpm: Number
    glucose_used_g: Number
    pyruvate_g: Number
    lactate_g: Number
    regulation_factor: Number


_DEFAULTS: dict[str, dict[str, Number]] = {
    "common": {
        "vmax_gpm": 0.2,
        "km_g": 5.0,
        "basal_regulation": 1.0,
        "min_regulation": 0.05,
        "max_regulation": 6.0,
        "reference_insulin": 0.05,
        "reference_glucagon": 0.05,
        "insulin_sensitivity": 2.5,
        "glucagon_sensitivity": 0.5,
        "max_fraction_per_min": 0.05,
        "pyruvate_yield": 0.98,
        "exercise_sensitivity": 0.0,
        "lactate_threshold": 0.9,
        "max_lactate_fraction": 0.0,
    },
    "muscle": {
        "vmax_gpm": 0.35,
        "km_g": 4.0,
        "insulin_sensitivity": 2.8,
        "glucagon_sensitivity": 0.25,
        "max_fraction_per_min": 0.07,
        "exercise_sensitivity": 4.0,
        "lactate_threshold": 0.35,
        "max_lactate_fraction": 0.65,
    },
    "liver": {
        "vmax_gpm": 0.22,
        "km_g": 6.0,
        "insulin_sensitivity": 1.8,
        "glucagon_sensitivity": 3.0,
        "max_fraction_per_min": 0.04,
    },
}

_EPS = 1e-9


def _get_signal_level(signals: Mapping[str, Signal] | MutableMapping[str, Signal], name: str, default: Number) -> Number:
    signal = signals.get(name) if signals is not None else None
    if signal is None:
        return float(default)
    return float(signal.level)


def _resolve_cfg(tissue: str, params: Mapping[str, Number] | None) -> dict[str, Number]:
    cfg: dict[str, Number] = dict(_DEFAULTS["common"])
    if tissue in _DEFAULTS:
        cfg.update(_DEFAULTS[tissue])
    if params:
        cfg.update({k: float(v) for k, v in params.items()})
    return cfg


def compute_glycolysis(
    *,
    tissue: str,
    glucose_g: Number,
    volume_L: Number | None,
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal],
    params: Mapping[str, Number] | None = None,
) -> GlycolysisResult:

    if dt_min <= 0.0 or glucose_g <= 0.0:
        return GlycolysisResult(0.0, 0.0, 0.0, 0.0, 0.0)

    cfg = _resolve_cfg(tissue, params)

    insulin = _get_signal_level(signals, "insulin", cfg["reference_insulin"])
    glucagon = _get_signal_level(signals, "glucagon", cfg["reference_glucagon"])
    exercise = 0.0
    if tissue.lower() == "muscle" and cfg.get("exercise_sensitivity", 0.0) > 0.0:
        exercise = _get_signal_level(signals, "exercise", 0.0)

    regulation = cfg["basal_regulation"]
    regulation += cfg.get("insulin_sensitivity", 0.0) * max(0.0, insulin - cfg["reference_insulin"])
    regulation -= cfg.get("glucagon_sensitivity", 0.0) * max(0.0, glucagon - cfg["reference_glucagon"])

    if exercise > 0.0:
        regulation *= 1.0 + cfg.get("exercise_sensitivity", 0.0) * exercise

    regulation = max(cfg["min_regulation"], min(cfg["max_regulation"], regulation))

    vmax = max(cfg["vmax_gpm"], 0.0) * regulation
    if vmax <= 0.0:
        return GlycolysisResult(0.0, 0.0, 0.0, 0.0, regulation)

    substrate = float(glucose_g)
    km = max(cfg.get("km_g", substrate), _EPS)
    rate_gpm = sat_mm(substrate, vmax, km)

    frac_limit = max(cfg.get("max_fraction_per_min", 0.0), 0.0)
    dt_safe = max(dt_min, _EPS)
    if frac_limit > 0.0:
        rate_gpm = min(rate_gpm, (substrate * frac_limit) / dt_safe)

    rate_gpm = min(rate_gpm, substrate / dt_safe)

    glucose_used = rate_gpm * dt_min
    if glucose_used <= 0.0:
        return GlycolysisResult(0.0, 0.0, 0.0, 0.0, regulation)

    yield_coeff = max(cfg.get("pyruvate_yield", 1.0), 0.0)
    pyruvate = glucose_used * yield_coeff
    lactate = 0.0

    max_lactate_fraction = max(min(cfg.get("max_lactate_fraction", 0.0), 1.0), 0.0)
    if max_lactate_fraction > 0.0 and exercise > 0.0:
        threshold = min(max(cfg.get("lactate_threshold", 0.0), 0.0), 1.0)
        if exercise > threshold:
            intensity = (exercise - threshold) / max(1.0 - threshold, _EPS)
            intensity = max(0.0, min(1.0, intensity))
            lactate = pyruvate * intensity * max_lactate_fraction
            pyruvate = max(0.0, pyruvate - lactate)

    return GlycolysisResult(rate_gpm, glucose_used, pyruvate, lactate, regulation)
