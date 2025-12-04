from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping

from ..core.base import Signal, sat_mm

Number = float


@dataclass(frozen=True)
class LacticFermentationResult:

    rate_gpm: Number
    pyruvate_used_g: Number
    lactate_g: Number
    nad_regenerated_mmol: Number
    regulation_factor: Number


_DEFAULTS: dict[str, dict[str, Number]] = {
    "common": {
        "vmax_gpm": 0.12,
        "km_g": 2.5,
        "basal_regulation": 0.15,
        "min_regulation": 0.0,
        "max_regulation": 6.0,
        "exercise_gain": 3.0,
        "deficit_gain": 1.6,
        "deficit_scale_mmol": 1.0,
        "oxygen_suppression": 1.5,
        "reference_oxygenation": 1.0,
        "max_fraction_per_min": 0.45,
        "lactate_yield": 1.02,
        "nad_regen_per_g": 0.0,
    },
    "muscle": {
        "vmax_gpm": 0.18,
        "km_g": 2.0,
        "basal_regulation": 0.25,
        "exercise_gain": 4.5,
        "deficit_gain": 2.0,
        "deficit_scale_mmol": 0.75,
        "oxygen_suppression": 2.5,
        "max_fraction_per_min": 0.6,
        "lactate_yield": 0.98,
        "nad_regen_per_g": 2.0,
    },
}

_EPS = 1e-9


def _get_signal_level(
    signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None,
    name: str,
    default: Number,
) -> Number:
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


def compute_lactic_fermentation(
    *,
    tissue: str,
    pyruvate_g: Number,
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None,
    params: Mapping[str, Number] | None = None,
) -> LacticFermentationResult:

    if dt_min <= 0.0 or pyruvate_g <= 0.0:
        return LacticFermentationResult(0.0, 0.0, 0.0, 0.0, 0.0)

    cfg = _resolve_cfg(tissue, params)

    exercise = 0.0
    if cfg.get("exercise_gain", 0.0) > 0.0:
        exercise = max(0.0, _get_signal_level(signals, "exercise", 0.0))

    deficit_scale = max(cfg.get("deficit_scale_mmol", 1.0), _EPS)
    deficit_level = 0.0
    if cfg.get("deficit_gain", 0.0) > 0.0:
        deficit_level = max(0.0, _get_signal_level(signals, "muscle_atp_deficit_mmol", 0.0))
        deficit_level = deficit_level / deficit_scale

    oxygen = _get_signal_level(signals, "oxygenation", cfg.get("reference_oxygenation", 1.0))

    regulation = cfg.get("basal_regulation", 1.0)
    regulation += cfg.get("exercise_gain", 0.0) * exercise
    regulation += cfg.get("deficit_gain", 0.0) * deficit_level

    ox_supp = cfg.get("oxygen_suppression", 0.0)
    if ox_supp > 0.0:
        regulation *= max(0.0, 1.0 - ox_supp * max(0.0, oxygen - cfg.get("reference_oxygenation", 1.0)))

    regulation = max(cfg.get("min_regulation", 0.0), min(cfg.get("max_regulation", 1.0), regulation))

    vmax = max(cfg.get("vmax_gpm", 0.0), 0.0) * regulation
    if vmax <= 0.0:
        return LacticFermentationResult(0.0, 0.0, 0.0, 0.0, regulation)

    substrate = max(float(pyruvate_g), 0.0)
    km = max(cfg.get("km_g", substrate), _EPS)
    rate_gpm = sat_mm(substrate, vmax, km)

    frac_limit = max(cfg.get("max_fraction_per_min", 0.0), 0.0)
    dt_safe = max(dt_min, _EPS)
    if frac_limit > 0.0:
        rate_gpm = min(rate_gpm, (substrate * frac_limit) / dt_safe)

    rate_gpm = min(rate_gpm, substrate / dt_safe)

    pyruvate_used = rate_gpm * dt_min
    if pyruvate_used <= 0.0:
        return LacticFermentationResult(0.0, 0.0, 0.0, 0.0, regulation)

    lactate_yield = max(cfg.get("lactate_yield", 1.0), 0.0)
    lactate = pyruvate_used * lactate_yield

    nad_regen = max(cfg.get("nad_regen_per_g", 0.0), 0.0) * lactate

    return LacticFermentationResult(rate_gpm, pyruvate_used, lactate, nad_regen, regulation)
