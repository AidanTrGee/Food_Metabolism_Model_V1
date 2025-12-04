from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, MutableMapping

from ..core.base import Signal, sat_mm

Number = float


@dataclass(frozen=True)
class GluconeogenesisResult:

    rate_gpm: Number
    glucose_g: Number
    precursor_used_g: Dict[str, Number] = field(default_factory=dict)
    regulation_factor: Number = 0.0


_DEFAULTS: dict[str, dict[str, Number]] = {
    "common": {
        "vmax_gpm": 0.12,
        "km_precursor_g": 6.0,
        "basal_regulation": 1.0,
        "min_regulation": 0.1,
        "max_regulation": 4.0,
        "insulin_suppression": 1.5,
        "glucagon_stimulation": 2.8,
        "reference_insulin": 0.05,
        "reference_glucagon": 0.05,
        "max_fraction_per_min": 0.06,
        "pyruvate_yield": 0.52,
        "lactate_yield": 0.50,
        "alanine_yield": 0.48,
    },
    "liver": {
        "vmax_gpm": 0.22,
        "glucagon_stimulation": 3.2,
        "insulin_suppression": 1.8,
        "max_fraction_per_min": 0.08,
    },
    "kidney": {
        "vmax_gpm": 0.08,
        "glucagon_stimulation": 1.5,
        "max_fraction_per_min": 0.05,
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


def compute_gluconeogenesis(
    *,
    tissue: str,
    precursors_g: Mapping[str, Number],
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal],
    params: Mapping[str, Number] | None = None,
) -> GluconeogenesisResult:

    if dt_min <= 0.0:
        return GluconeogenesisResult(0.0, 0.0, {}, 0.0)

    cfg = _resolve_cfg(tissue, params)

    insulin = _get_signal_level(signals, "insulin", cfg["reference_insulin"])
    glucagon = _get_signal_level(signals, "glucagon", cfg["reference_glucagon"])

    regulation = cfg["basal_regulation"]
    regulation -= cfg.get("insulin_suppression", 0.0) * max(0.0, insulin - cfg["reference_insulin"])
    regulation += cfg.get("glucagon_stimulation", 0.0) * max(0.0, glucagon - cfg["reference_glucagon"])
    regulation = max(cfg["min_regulation"], min(cfg["max_regulation"], regulation))

    vmax = max(cfg["vmax_gpm"], 0.0) * regulation
    if vmax <= 0.0:
        return GluconeogenesisResult(0.0, 0.0, {}, regulation)

    yields = {
        "lactate": max(cfg.get("lactate_yield", 0.0), 0.0),
        "pyruvate": max(cfg.get("pyruvate_yield", 0.0), 0.0),
        "alanine": max(cfg.get("alanine_yield", 0.0), 0.0),
    }

    total_possible = 0.0
    for name, yld in yields.items():
        if yld <= 0.0:
            continue
        total_possible += max(precursors_g.get(name, 0.0), 0.0) * yld

    if total_possible <= 0.0:
        return GluconeogenesisResult(0.0, 0.0, {}, regulation)

    km = max(cfg.get("km_precursor_g", total_possible), _EPS)
    rate_gpm = sat_mm(total_possible, vmax, km)

    frac_limit = max(cfg.get("max_fraction_per_min", 0.0), 0.0)
    dt_safe = max(dt_min, _EPS)
    if frac_limit > 0.0:
        rate_gpm = min(rate_gpm, (total_possible * frac_limit) / dt_safe)

    rate_gpm = max(0.0, min(rate_gpm, total_possible / dt_safe))

    glucose_g = rate_gpm * dt_min
    if glucose_g <= 0.0:
        return GluconeogenesisResult(0.0, 0.0, {}, regulation)

    glucose_g = min(glucose_g, total_possible)

    remaining = glucose_g
    used: Dict[str, Number] = {}
    for name in ("lactate", "pyruvate", "alanine"):
        yld = yields.get(name, 0.0)
        pool = max(precursors_g.get(name, 0.0), 0.0)
        if yld <= 0.0 or pool <= 0.0 or remaining <= 0.0:
            continue
        max_from_pool = pool * yld
        take = min(remaining, max_from_pool)
        if take <= 0.0:
            continue
        used_mass = take / max(yld, _EPS)
        used[name] = used_mass
        remaining -= take

    if remaining > _EPS:
        glucose_g -= remaining
        remaining = 0.0

    return GluconeogenesisResult(rate_gpm, glucose_g, used, regulation)
