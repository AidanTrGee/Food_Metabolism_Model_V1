from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping

from ..core.base import Signal, sat_mm

Number = float


@dataclass(frozen=True)
class ATPDemandResult:
    rate_mmol_per_min: Number
    demand_mmol: Number
    regulation_factor: Number


@dataclass(frozen=True)
class CitricAcidCycleResult:
    rate_acetyl_coa_gpm: Number
    rate_pyruvate_gpm: Number
    acetyl_coa_used_g: Number
    pyruvate_used_g: Number
    atp_generated_mmol: Number
    nadh_generated_mmol: Number
    regulation_factor: Number


@dataclass(frozen=True)
class OxidativePhosphorylationResult:
    atp_generated_mmol: Number
    nadh_used_mmol: Number
    oxygen_used_mmol: Number
    deficit_mmol: Number
    regulation_factor: Number


_DEFAULTS: dict[str, dict[str, Number]] = {
    "common": {
        "basal_cac_regulation": 1.0,
        "cac_min_regulation": 0.1,
        "cac_max_regulation": 6.0,
        "cac_vmax_acoa_gpm": 0.06,
        "cac_km_acoa_g": 1.8,
        "cac_vmax_pyruvate_gpm": 0.05,
        "cac_km_pyruvate_g": 2.5,
        "cac_atp_per_g_acoa": 20.0,
        "cac_atp_per_g_pyruvate": 14.0,
        "cac_nadh_per_g_acoa": 6.0,
        "cac_nadh_per_g_pyruvate": 4.0,
        "cac_insulin_suppression": 1.2,
        "cac_glucagon_stimulation": 0.6,
        "cac_exercise_gain": 1.5,
        "reference_insulin": 0.05,
        "reference_glucagon": 0.05,
        "basal_atp_demand_mmol_per_min": 15.0,
        "demand_min_regulation": 0.2,
        "demand_max_regulation": 6.0,
        "demand_exercise_gain": 35.0,
        "demand_insulin_suppression": 2.5,
        "demand_glucagon_stimulation": 0.8,
        "demand_reference_exercise": 0.0,
        "op_basal_regulation": 1.0,
        "op_min_regulation": 0.2,
        "op_max_regulation": 5.0,
        "op_vmax_atp_mmol_per_min": 45.0,
        "op_coupling_efficiency": 2.8,
        "op_leak_fraction": 0.08,
        "op_reference_oxygen": 0.9,
        "op_oxygen_sensitivity": 1.7,
        "op_exercise_stimulation": 0.6,
        "op_oxygen_per_nadh": 0.5,
        "max_atp_pool_mmol": 30.0,
        "initial_atp_mmol": 15.0,
    },
    "muscle": {
        "basal_atp_demand_mmol_per_min": 24.0,
        "demand_exercise_gain": 60.0,
        "cac_vmax_acoa_gpm": 0.09,
        "cac_vmax_pyruvate_gpm": 0.08,
        "cac_exercise_gain": 2.2,
        "op_vmax_atp_mmol_per_min": 65.0,
        "op_leak_fraction": 0.1,
        "op_exercise_stimulation": 1.2,
        "op_reference_oxygen": 1.0,
        "max_atp_pool_mmol": 36.0,
        "initial_atp_mmol": 22.0,
    },
    "liver": {
        "basal_atp_demand_mmol_per_min": 12.0,
        "demand_glucagon_stimulation": 1.2,
        "cac_vmax_acoa_gpm": 0.055,
        "cac_vmax_pyruvate_gpm": 0.045,
        "cac_insulin_suppression": 0.8,
        "op_vmax_atp_mmol_per_min": 38.0,
        "op_coupling_efficiency": 2.4,
        "op_leak_fraction": 0.05,
        "max_atp_pool_mmol": 22.0,
        "initial_atp_mmol": 14.0,
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


def compute_atp_demand(
    *,
    tissue: str,
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None,
    params: Mapping[str, Number] | None = None,
) -> ATPDemandResult:

    if dt_min <= 0.0:
        return ATPDemandResult(0.0, 0.0, 0.0)

    cfg = _resolve_cfg(tissue, params)

    insulin_ref = cfg.get("reference_insulin", 0.05)
    glucagon_ref = cfg.get("reference_glucagon", 0.05)
    exercise_ref = cfg.get("demand_reference_exercise", 0.0)

    insulin = _get_signal_level(signals, "insulin", insulin_ref)
    glucagon = _get_signal_level(signals, "glucagon", glucagon_ref)
    exercise = _get_signal_level(signals, "exercise", exercise_ref)

    regulation = 1.0
    regulation += cfg.get("demand_exercise_gain", 0.0) * max(0.0, exercise - exercise_ref)
    regulation -= cfg.get("demand_insulin_suppression", 0.0) * max(0.0, insulin - insulin_ref)
    regulation += cfg.get("demand_glucagon_stimulation", 0.0) * max(0.0, glucagon - glucagon_ref)

    regulation = max(cfg.get("demand_min_regulation", 0.0), min(cfg.get("demand_max_regulation", 10.0), regulation))

    base_rate = max(cfg.get("basal_atp_demand_mmol_per_min", 0.0), 0.0)
    rate = base_rate * regulation
    demand = rate * dt_min

    return ATPDemandResult(rate, demand, regulation)


def compute_citric_acid_cycle(
    *,
    tissue: str,
    pyruvate_g: Number,
    acetyl_coa_g: Number,
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None,
    params: Mapping[str, Number] | None = None,
) -> CitricAcidCycleResult:

    if dt_min <= 0.0 or (pyruvate_g <= 0.0 and acetyl_coa_g <= 0.0):
        return CitricAcidCycleResult(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    cfg = _resolve_cfg(tissue, params)

    insulin_ref = cfg.get("reference_insulin", 0.05)
    glucagon_ref = cfg.get("reference_glucagon", 0.05)

    insulin = _get_signal_level(signals, "insulin", insulin_ref)
    glucagon = _get_signal_level(signals, "glucagon", glucagon_ref)
    exercise = _get_signal_level(signals, "exercise", 0.0)

    regulation = cfg.get("basal_cac_regulation", 1.0)
    regulation -= cfg.get("cac_insulin_suppression", 0.0) * max(0.0, insulin - insulin_ref)
    regulation += cfg.get("cac_glucagon_stimulation", 0.0) * max(0.0, glucagon - glucagon_ref)
    if exercise > 0.0:
        regulation += cfg.get("cac_exercise_gain", 0.0) * exercise

    regulation = max(cfg.get("cac_min_regulation", 0.0), min(cfg.get("cac_max_regulation", 10.0), regulation))

    vmax_acoa = max(cfg.get("cac_vmax_acoa_gpm", 0.0), 0.0) * regulation
    vmax_pyruvate = max(cfg.get("cac_vmax_pyruvate_gpm", 0.0), 0.0) * regulation

    dt_safe = max(dt_min, _EPS)

    acetyl_available = max(acetyl_coa_g, 0.0)
    pyruvate_available = max(pyruvate_g, 0.0)

    acetyl_rate = 0.0
    if vmax_acoa > 0.0 and acetyl_available > 0.0:
        km_acoa = max(cfg.get("cac_km_acoa_g", acetyl_available), _EPS)
        acetyl_rate = sat_mm(acetyl_available, vmax_acoa, km_acoa)
        acetyl_rate = min(acetyl_rate, acetyl_available / dt_safe)

    pyruvate_rate = 0.0
    if vmax_pyruvate > 0.0 and pyruvate_available > 0.0:
        km_pyruvate = max(cfg.get("cac_km_pyruvate_g", pyruvate_available), _EPS)
        pyruvate_rate = sat_mm(pyruvate_available, vmax_pyruvate, km_pyruvate)
        pyruvate_rate = min(pyruvate_rate, pyruvate_available / dt_safe)

    acetyl_used = acetyl_rate * dt_min
    pyruvate_used = pyruvate_rate * dt_min

    if acetyl_used <= 0.0 and pyruvate_used <= 0.0:
        return CitricAcidCycleResult(acetyl_rate, pyruvate_rate, 0.0, 0.0, 0.0, 0.0, regulation)

    atp_from_acetyl = acetyl_used * max(cfg.get("cac_atp_per_g_acoa", 0.0), 0.0)
    atp_from_pyruvate = pyruvate_used * max(cfg.get("cac_atp_per_g_pyruvate", 0.0), 0.0)
    atp_generated = atp_from_acetyl + atp_from_pyruvate

    nadh_from_acetyl = acetyl_used * max(cfg.get("cac_nadh_per_g_acoa", 0.0), 0.0)
    nadh_from_pyruvate = pyruvate_used * max(cfg.get("cac_nadh_per_g_pyruvate", 0.0), 0.0)
    nadh_generated = nadh_from_acetyl + nadh_from_pyruvate

    return CitricAcidCycleResult(
        acetyl_rate,
        pyruvate_rate,
        acetyl_used,
        pyruvate_used,
        atp_generated,
        nadh_generated,
        regulation,
    )


def compute_oxidative_phosphorylation(
    *,
    tissue: str,
    nadh_mmol: Number,
    atp_demand_mmol: Number,
    dt_min: Number,
    signals: Mapping[str, Signal] | MutableMapping[str, Signal] | None,
    params: Mapping[str, Number] | None = None,
) -> OxidativePhosphorylationResult:

    if dt_min <= 0.0 or atp_demand_mmol <= 0.0 or nadh_mmol <= 0.0:
        return OxidativePhosphorylationResult(0.0, 0.0, 0.0, max(atp_demand_mmol, 0.0), 0.0)

    cfg = _resolve_cfg(tissue, params)

    oxygen_ref = cfg.get("op_reference_oxygen", 0.9)
    exercise = _get_signal_level(signals, "exercise", 0.0)
    oxygen = _get_signal_level(signals, "oxygen", oxygen_ref)

    regulation = cfg.get("op_basal_regulation", 1.0)
    regulation += cfg.get("op_oxygen_sensitivity", 0.0) * max(0.0, oxygen - oxygen_ref)
    if exercise > 0.0:
        regulation += cfg.get("op_exercise_stimulation", 0.0) * exercise
    regulation = max(cfg.get("op_min_regulation", 0.0), min(cfg.get("op_max_regulation", 10.0), regulation))

    vmax_atp = max(cfg.get("op_vmax_atp_mmol_per_min", 0.0), 0.0) * regulation
    dt_safe = max(dt_min, _EPS)
    capacity = vmax_atp * dt_safe

    coupling = max(cfg.get("op_coupling_efficiency", 0.0), _EPS)
    available_from_nadh = max(nadh_mmol, 0.0) * coupling

    leak = min(max(cfg.get("op_leak_fraction", 0.0), 0.0), 0.95)

    raw_limit = min(capacity, available_from_nadh)
    if raw_limit <= 0.0:
        return OxidativePhosphorylationResult(0.0, 0.0, 0.0, max(atp_demand_mmol, 0.0), regulation)

    deliverable_limit = raw_limit * (1.0 - leak)
    atp_demand = max(atp_demand_mmol, 0.0)

    if deliverable_limit >= atp_demand:
        delivered = atp_demand
        raw_used = delivered / max(1.0 - leak, _EPS)
    else:
        delivered = deliverable_limit
        raw_used = raw_limit

    nadh_used = raw_used / coupling
    oxygen_used = nadh_used * max(cfg.get("op_oxygen_per_nadh", 0.0), 0.0)
    deficit = max(atp_demand - delivered, 0.0)

    return OxidativePhosphorylationResult(delivered, nadh_used, oxygen_used, deficit, regulation)
