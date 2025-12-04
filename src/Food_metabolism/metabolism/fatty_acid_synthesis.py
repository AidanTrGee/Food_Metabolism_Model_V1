from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, MutableMapping

from ..core.base import Signal, sat_mm

Number = float


@dataclass(frozen=True)
class FattyAcidSynthesisResult:
    rate_gpm: Number
    lipid_g: Number
    substrate_used_g: Dict[str, Number] = field(default_factory=dict)
    regulation_factor: Number = 0.0


_DEFAULTS: dict[str, dict[str, Number]] = {
    "common": {
        "vmax_gpm": 0.16,
        "km_g": 12.0,
        "basal_regulation": 0.2,
        "min_regulation": 0.0,
        "max_regulation": 4.5,
        "insulin_stimulation": 2.5,
        "glucagon_suppression": 1.6,
        "reference_insulin": 0.05,
        "reference_glucagon": 0.05,
        "max_fraction_per_min": 0.05,
        "glucose_lipid_yield": 0.7,
        "pyruvate_lipid_yield": 0.65,
        "lactate_lipid_yield": 0.55,
        "ffa_reesterification_yield": 1.0,
    },
    "adipose": {
        "vmax_gpm": 0.22,
        "km_g": 16.0,
        "insulin_stimulation": 3.2,
        "glucagon_suppression": 1.9,
        "max_fraction_per_min": 0.06,
        "glucose_lipid_yield": 0.75,
        "ffa_reesterification_yield": 0.95,
    },
    "liver": {
        "vmax_gpm": 0.28,
        "km_g": 18.0,
        "insulin_stimulation": 2.8,
        "glucagon_suppression": 2.3,
        "max_fraction_per_min": 0.07,
        "pyruvate_lipid_yield": 0.7,
        "lactate_lipid_yield": 0.6,
    },
}

_EPS = 1e-9


def _get_signal_level(signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None, name: str, default: Number) -> Number:
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


def compute_fatty_acid_synthesis(
    *,
    tissue: str,
    substrates_g: Mapping[str, Number],
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None,
    params: Mapping[str, Number] | None = None,
) -> FattyAcidSynthesisResult:

    if dt_min <= 0.0:
        return FattyAcidSynthesisResult(0.0, 0.0, {}, 0.0)

    total_substrate = 0.0
    for grams in substrates_g.values():
        total_substrate += max(float(grams), 0.0)
    if total_substrate <= 0.0:
        return FattyAcidSynthesisResult(0.0, 0.0, {}, 0.0)

    cfg = _resolve_cfg(tissue, params)

    insulin = _get_signal_level(signals, "insulin", cfg["reference_insulin"])
    glucagon = _get_signal_level(signals, "glucagon", cfg["reference_glucagon"])

    regulation = cfg["basal_regulation"]
    regulation += cfg.get("insulin_stimulation", 0.0) * max(0.0, insulin - cfg["reference_insulin"])
    regulation -= cfg.get("glucagon_suppression", 0.0) * max(0.0, glucagon - cfg["reference_glucagon"])
    regulation = max(cfg["min_regulation"], min(cfg["max_regulation"], regulation))

    vmax = max(cfg["vmax_gpm"], 0.0) * regulation
    if vmax <= 0.0:
        return FattyAcidSynthesisResult(0.0, 0.0, {}, regulation)

    yields = {
        "glucose": max(cfg.get("glucose_lipid_yield", 0.0), 0.0),
        "pyruvate": max(cfg.get("pyruvate_lipid_yield", 0.0), 0.0),
        "lactate": max(cfg.get("lactate_lipid_yield", 0.0), 0.0),
        "ffa": max(cfg.get("ffa_reesterification_yield", 0.0), 0.0),
    }

    potential_lipid = 0.0
    for name, grams in substrates_g.items():
        pool = max(float(grams), 0.0)
        if pool <= 0.0:
            continue
        yld = yields.get(name, 0.0)
        if yld <= 0.0:
            continue
        potential_lipid += pool * yld

    if potential_lipid <= 0.0:
        return FattyAcidSynthesisResult(0.0, 0.0, {}, regulation)

    km = max(cfg.get("km_g", potential_lipid), _EPS)
    rate_gpm = sat_mm(potential_lipid, vmax, km)

    frac_limit = max(cfg.get("max_fraction_per_min", 0.0), 0.0)
    dt_safe = max(dt_min, _EPS)
    if frac_limit > 0.0:
        rate_gpm = min(rate_gpm, (potential_lipid * frac_limit) / dt_safe)

    rate_gpm = min(rate_gpm, potential_lipid / dt_safe)

    lipid_g = rate_gpm * dt_min
    if lipid_g <= 0.0:
        return FattyAcidSynthesisResult(0.0, 0.0, {}, regulation)

    lipid_g = min(lipid_g, potential_lipid)

    remaining = lipid_g
    used: Dict[str, Number] = {}
    for name in ("glucose", "pyruvate", "lactate", "ffa"):
        pool = max(substrates_g.get(name, 0.0), 0.0)
        yld = yields.get(name, 0.0)
        if pool <= 0.0 or yld <= 0.0 or remaining <= 0.0:
            continue
        max_from_pool = pool * yld
        take = min(remaining, max_from_pool)
        if take <= 0.0:
            continue
        used_mass = take / max(yld, _EPS)
        used[name] = used_mass
        remaining -= take

    if remaining > _EPS:
        lipid_g -= remaining
        remaining = 0.0

    return FattyAcidSynthesisResult(rate_gpm, lipid_g, used, regulation)
