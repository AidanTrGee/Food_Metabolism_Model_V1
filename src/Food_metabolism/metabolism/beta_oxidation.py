from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping

from ..core.base import Signal, sat_mm

Number = float


@dataclass(frozen=True)
class BetaOxidationResult:

    rate_gpm: Number
    lipid_used_g: Number
    acetyl_coa_g: Number
    ketone_precursor_g: Number
    regulation_factor: Number


_DEFAULTS: dict[str, dict[str, Number]] = {
    "common": {
        "vmax_gpm": 0.12,
        "km_g": 4.0,
        "basal_regulation": 1.0,
        "min_regulation": 0.05,
        "max_regulation": 6.0,
        "reference_insulin": 0.05,
        "reference_glucagon": 0.05,
        "insulin_suppression": 2.2,
        "glucagon_stimulation": 1.2,
        "exercise_stimulation": 0.0,
        "max_fraction_per_min": 0.05,
        "acetyl_yield": 0.70,
        "ketone_fraction": 0.0,
    },
    "muscle": {
        "vmax_gpm": 0.18,
        "km_g": 3.0,
        "insulin_suppression": 1.8,
        "glucagon_stimulation": 0.4,
        "exercise_stimulation": 3.5,
        "max_fraction_per_min": 0.06,
        "ketone_fraction": 0.0,
    },
    "liver": {
        "vmax_gpm": 0.14,
        "km_g": 5.0,
        "insulin_suppression": 2.6,
        "glucagon_stimulation": 2.2,
        "exercise_stimulation": 0.5,
        "max_fraction_per_min": 0.045,
        "ketone_fraction": 0.45,
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


def compute_beta_oxidation(
    *,
    tissue: str,
    lipid_g: Number,
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None,
    params: Mapping[str, Number] | None = None,
) -> BetaOxidationResult:

    if dt_min <= 0.0 or lipid_g <= 0.0:
        return BetaOxidationResult(0.0, 0.0, 0.0, 0.0, 0.0)

    cfg = _resolve_cfg(tissue, params)

    insulin = _get_signal_level(signals, "insulin", cfg["reference_insulin"])
    glucagon = _get_signal_level(signals, "glucagon", cfg["reference_glucagon"])
    exercise = 0.0
    if cfg.get("exercise_stimulation", 0.0) > 0.0:
        exercise = max(0.0, _get_signal_level(signals, "exercise", 0.0))

    regulation = cfg["basal_regulation"]
    regulation -= cfg.get("insulin_suppression", 0.0) * max(0.0, insulin - cfg["reference_insulin"])
    regulation += cfg.get("glucagon_stimulation", 0.0) * max(0.0, glucagon - cfg["reference_glucagon"])
    if exercise > 0.0:
        regulation += cfg.get("exercise_stimulation", 0.0) * exercise

    regulation = max(cfg["min_regulation"], min(cfg["max_regulation"], regulation))

    vmax = max(cfg["vmax_gpm"], 0.0) * regulation
    if vmax <= 0.0:
        return BetaOxidationResult(0.0, 0.0, 0.0, 0.0, regulation)

    substrate = max(float(lipid_g), 0.0)
    km = max(cfg.get("km_g", substrate), _EPS)
    rate_gpm = sat_mm(substrate, vmax, km)

    frac_limit = max(cfg.get("max_fraction_per_min", 0.0), 0.0)
    dt_safe = max(dt_min, _EPS)
    if frac_limit > 0.0:
        rate_gpm = min(rate_gpm, (substrate * frac_limit) / dt_safe)

    rate_gpm = min(rate_gpm, substrate / dt_safe)

    lipid_used = rate_gpm * dt_min
    if lipid_used <= 0.0:
        return BetaOxidationResult(0.0, 0.0, 0.0, 0.0, regulation)

    acetyl_yield = max(cfg.get("acetyl_yield", 0.0), 0.0)
    acetyl_coa = lipid_used * acetyl_yield

    ketone_fraction = max(0.0, min(1.0, cfg.get("ketone_fraction", 0.0)))
    ketone_precursor = acetyl_coa * ketone_fraction

    return BetaOxidationResult(rate_gpm, lipid_used, acetyl_coa, ketone_precursor, regulation)
