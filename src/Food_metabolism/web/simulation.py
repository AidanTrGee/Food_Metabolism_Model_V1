from __future__ import annotations

import copy
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Tuple

from Food_metabolism.components.adipose import Adipose
from Food_metabolism.components.blood import ArterialBlood, Lungs, PortalVein, VenousBlood
from Food_metabolism.components.gi import IntestinalLumen, IntestinalTissue, StomachLumen
from Food_metabolism.components.kidney import Kidney
from Food_metabolism.components.liver import Liver
from Food_metabolism.components.lymph import Lymph
from Food_metabolism.components.muscle import Muscle
from Food_metabolism.configuration import load_simulation_config
from Food_metabolism.controllers.pancreas_hormones import PancreasHormoneController
from Food_metabolism.food.meals import Meal, apply_meals_to_stomach
from Food_metabolism.runner.sim import SimulationEngine
from Food_metabolism.runner.utils import capture_compartment_state
from Food_metabolism.transfer.adipose_lipids import AdiposeFFARelease, AdiposeLipidUptake
from Food_metabolism.transfer.absorption import IntestinalBasolateralExport, IntestinalGlucoseAbsorption
from Food_metabolism.transfer.gi_transit import GastricEnergyEmptying
from Food_metabolism.transfer.hepatic import HepaticGlucoseOutput, HepaticGlucoseUptake
from Food_metabolism.transfer.kidney import RenalGlucoseExcretion
from Food_metabolism.transfer.lymph_fat import IntestinalFatToLymph, LymphToVenousFat
from Food_metabolism.transfer.muscle_uptake import AdiposeGlucoseUptake, MuscleGlucoseUptake
from Food_metabolism.transfer.vascular import FirstOrderTransfer

RepoConfig = Mapping[str, Any]
ParameterOverrides = Mapping[str, float]

_PARAM_KEYS: Dict[str, Dict[str, Any]] = {
    "digestion.stomach.base_energy_emptying_kcal_min": {
        "label": "Gastric Emptying (kcal/min)",
        "step": 0.1,
        "min": 0.2,
        "max": 8.0,
    },
    "digestion.intestine.glucose_absorption.sglt1_vmax_gpm": {
        "label": "SGLT1 Vmax (g/min)",
        "step": 0.05,
        "min": 0.1,
        "max": 5.0,
    },
    "digestion.liver.output.basal_output_gpm": {
        "label": "Basal HGO (g/min)",
        "step": 0.01,
        "min": 0.05,
        "max": 0.5,
    },
    "digestion.periphery.muscle.insulin_multiplier": {
        "label": "Muscle Insulin Sensitivity (×)",
        "step": 0.1,
        "min": 0.5,
        "max": 6.0,
    },
    "kidney.excretion_k": {
        "label": "Renal Clearance Gain", "step": 0.05, "min": 0.1, "max": 2.0
    },
}


def load_base_config() -> Dict[str, Any]:
    return load_simulation_config()


def flatten_config_value(cfg: Mapping[str, Any], path: str) -> float:
    cursor: Any = cfg
    for key in path.split("."):
        if not isinstance(cursor, Mapping) or key not in cursor:
            raise KeyError(f"Config path '{path}' missing segment '{key}'")
        cursor = cursor[key]
    if isinstance(cursor, (int, float)):
        return float(cursor)
    raise TypeError(f"Config value for '{path}' is not numeric: {cursor!r}")


def apply_overrides(cfg: MutableMapping[str, Any], overrides: ParameterOverrides | None) -> None:
    if not overrides:
        return
    for path, value in overrides.items():
        keys = path.split(".")
        cursor: MutableMapping[str, Any] = cfg
        for key in keys[:-1]:
            current = cursor.get(key)
            if current is None:
                current = {}
                cursor[key] = current
            if not isinstance(current, MutableMapping):
                raise TypeError(f"Cannot set nested key '{path}' through non-mapping node '{key}'")
            cursor = current
        cursor[keys[-1]] = float(value)


def build_engine(cfg: Mapping[str, Any]) -> Tuple[SimulationEngine, List[Meal], StomachLumen]:
    digestion_cfg = cfg.get("digestion", {}) or {}
    intestine_cfg = digestion_cfg.get("intestine", {}) or {}
    liver_cfg = digestion_cfg.get("liver", {}) or {}
    periphery_cfg = digestion_cfg.get("periphery", {}) or {}
    portal_cfg = digestion_cfg.get("portal_mixing", {}) or {}
    blood_flow_cfg = cfg.get("blood_flow", {}) or {}
    kidney_cfg = cfg.get("kidney", {}) or {}
    blood_cfg = cfg.get("blood", {}) or {}

    stomach = StomachLumen(digestion_cfg.get("stomach"))
    intestinal_lumen = IntestinalLumen(intestine_cfg.get("lumen"))
    intestinal_tissue = IntestinalTissue(intestine_cfg.get("tissue"))
    portal = PortalVein()
    liver = Liver(liver_cfg.get("compartment"))
    venous = VenousBlood(blood_cfg.get("venous"))
    lungs = Lungs(blood_cfg.get("lungs"))
    arterial = ArterialBlood(blood_cfg.get("arterial"))
    muscle = Muscle(periphery_cfg.get("muscle_compartment"))
    adipose = Adipose(periphery_cfg.get("adipose_compartment"))
    lymph_cfg = digestion_cfg.get("lymph", {}) or {}
    lymph = Lymph(lymph_cfg.get("compartment", lymph_cfg))
    kidney = Kidney(kidney_cfg.get("compartment"))

    init_cfg = cfg.get("initial", {}) or {}
    venous.add("glucose", float(init_cfg.get("venous_glucose_g", 3.2)))
    arterial.add("glucose", float(init_cfg.get("arterial_glucose_g", 2.0)))
    lungs.add("glucose", float(init_cfg.get("pulmonary_glucose_g", 1.6)))
    liver.add("glucose", float(init_cfg.get("liver_glucose_g", 80.0)))
    muscle.add("glucose", float(init_cfg.get("muscle_glucose_g", 300.0)))
    adipose.add("glucose", float(init_cfg.get("adipose_glucose_g", 20.0)))
    adipose.add("lipid", float(init_cfg.get("adipose_lipid_g", 12000.0)))
    venous.add("lipid", float(init_cfg.get("venous_lipid_g", 0.3)))
    arterial.add("lipid", float(init_cfg.get("arterial_lipid_g", 0.2)))
    venous.add("ffa", float(init_cfg.get("venous_ffa_g", 0.0)))
    arterial.add("ffa", float(init_cfg.get("arterial_ffa_g", 0.0)))
    lungs.add("ffa", float(init_cfg.get("pulmonary_ffa_g", 0.0)))
    adipose.add("ffa", float(init_cfg.get("adipose_ffa_g", 0.0)))

    hormone_cfg = cfg.get("hormones", {}) or {}
    controller_cfg = hormone_cfg.get("pancreas_controller", {}) if isinstance(hormone_cfg, dict) else {}
    controllers = [
        PancreasHormoneController(
            arterial,
            enabled=hormone_cfg.get("enabled", True),
            target_mg_dl=float(hormone_cfg.get("target_mg_dl", 90.0)),
            insulin_gain=float(controller_cfg.get("insulin_gain", hormone_cfg.get("insulin_gain", 0.02))),
            glucagon_gain=float(controller_cfg.get("glucagon_gain", hormone_cfg.get("glucagon_gain", 0.015))),
            basal_insulin=float(controller_cfg.get("basal_insulin", hormone_cfg.get("basal_insulin", 0.05 if hormone_cfg.get("enabled", True) else 0.0))),
            basal_glucagon=float(controller_cfg.get("basal_glucagon", hormone_cfg.get("basal_glucagon", 0.05 if hormone_cfg.get("enabled", True) else 0.0))),
            max_signal=float(controller_cfg.get("max_signal", hormone_cfg.get("max_signal", 1.0))),
            response_alpha=float(controller_cfg.get("response_alpha", hormone_cfg.get("response_alpha", 0.25))),
            insulin_threshold_mg_dl=float(controller_cfg.get("insulin_threshold_mg_dl", hormone_cfg.get("insulin_threshold_mg_dl", 100.0))),
            glucagon_threshold_mg_dl=float(controller_cfg.get("glucagon_threshold_mg_dl", hormone_cfg.get("glucagon_threshold_mg_dl", 80.0))),
        )
    ]

    glucose_absorption_cfg = intestine_cfg.get("glucose_absorption", {}) or {}
    basolateral_cfg = intestine_cfg.get("basolateral_export", {}) or {}
    protein_abs_cfg = intestine_cfg.get("protein_absorption", {"k_per_min": 0.08})
    fat_abs_cfg = intestine_cfg.get("fat_absorption", {"k_per_min": 0.12})
    lymph_drain_cfg = lymph_cfg.get("drainage", {"vmax_gpm": 1.6, "km_g": 8.0})
    chylo_clear_cfg = lymph_cfg.get("chylomicron_clearance", lymph_drain_cfg)
    portal_glu_mix = portal_cfg.get("glucose", {"k_per_min": 4.0})
    portal_pro_mix = portal_cfg.get("protein", {"k_per_min": 4.0})

    ven_to_lung_cfg = blood_flow_cfg.get("venous_to_lungs", {"k_per_min": 1.5})
    lung_to_art_cfg = blood_flow_cfg.get("lungs_to_arterial", {"k_per_min": 4.0})
    art_to_ven_cfg = blood_flow_cfg.get("arterial_to_venous", {"k_per_min": 2.7})

    edges = [
        GastricEnergyEmptying("GE_glucose", stomach, intestinal_lumen, "glucose", {}),
        GastricEnergyEmptying("GE_lipid", stomach, intestinal_lumen, "lipid", {}),
        GastricEnergyEmptying("GE_protein", stomach, intestinal_lumen, "protein", {}),
        GastricEnergyEmptying("GE_fiber", stomach, intestinal_lumen, "soluble_fiber", {}),
        IntestinalGlucoseAbsorption("Intestinal_glucose_uptake", intestinal_lumen, intestinal_tissue, "glucose", glucose_absorption_cfg),
        IntestinalBasolateralExport("Enterocyte_to_portal", intestinal_tissue, portal, "glucose", basolateral_cfg),
        FirstOrderTransfer("Protein_absorption", intestinal_lumen, portal, "protein", protein_abs_cfg),
        IntestinalFatToLymph("Fat_to_lymph", intestinal_lumen, lymph, "lipid", fat_abs_cfg),
        LymphToVenousFat("Lymph_to_venous", lymph, venous, "lipid", chylo_clear_cfg),
        HepaticGlucoseUptake("Portal_to_liver", portal, liver, "glucose", liver_cfg.get("uptake", {})),
        FirstOrderTransfer("Portal_glucose_bypass", portal, venous, "glucose", portal_glu_mix),
        FirstOrderTransfer("Portal_protein_mix", portal, venous, "protein", portal_pro_mix),
        HepaticGlucoseOutput("Liver_output", liver, venous, "glucose", liver_cfg.get("output", {})),
        FirstOrderTransfer("Venous_to_lungs_glucose", venous, lungs, "glucose", ven_to_lung_cfg),
        FirstOrderTransfer("Venous_to_lungs_lipid", venous, lungs, "lipid", ven_to_lung_cfg),
        FirstOrderTransfer("Venous_to_lungs_ffa", venous, lungs, "ffa", ven_to_lung_cfg),
        FirstOrderTransfer("Venous_to_lungs_protein", venous, lungs, "protein", ven_to_lung_cfg),
        FirstOrderTransfer("Lungs_to_arterial_glucose", lungs, arterial, "glucose", lung_to_art_cfg),
        FirstOrderTransfer("Lungs_to_arterial_lipid", lungs, arterial, "lipid", lung_to_art_cfg),
        FirstOrderTransfer("Lungs_to_arterial_ffa", lungs, arterial, "ffa", lung_to_art_cfg),
        FirstOrderTransfer("Lungs_to_arterial_protein", lungs, arterial, "protein", lung_to_art_cfg),
        FirstOrderTransfer("Arterial_to_venous_glucose", arterial, venous, "glucose", art_to_ven_cfg),
        FirstOrderTransfer("Arterial_to_venous_lipid", arterial, venous, "lipid", art_to_ven_cfg),
        FirstOrderTransfer("Arterial_to_venous_ffa", arterial, venous, "ffa", art_to_ven_cfg),
        FirstOrderTransfer("Arterial_to_venous_protein", arterial, venous, "protein", art_to_ven_cfg),
        MuscleGlucoseUptake("Muscle_glucose", arterial, muscle, "glucose", periphery_cfg.get("muscle", {})),
        AdiposeGlucoseUptake("Adipose_glucose", arterial, adipose, "glucose", periphery_cfg.get("adipose_glucose", {})),
        AdiposeLipidUptake("Adipose_lipid_uptake", venous, adipose, "lipid", periphery_cfg.get("adipose_lipid", {})),
        AdiposeLipidUptake("Adipose_ffa_reuptake", venous, adipose, "ffa", periphery_cfg.get("adipose_ffa", periphery_cfg.get("adipose_lipid", {}))),
        AdiposeFFARelease("Adipose_ffa_release", adipose, venous, "ffa", periphery_cfg.get("adipose_release", {})),
        RenalGlucoseExcretion("Renal_glucose_excretion", arterial, kidney, "glucose", kidney_cfg),
    ]

    engine = SimulationEngine(
        compartments=[
            stomach,
            intestinal_lumen,
            intestinal_tissue,
            portal,
            liver,
            venous,
            lungs,
            arterial,
            lymph,
            muscle,
            adipose,
            kidney,
        ],
        edges=edges,
        controllers=controllers,
    )

    meals = [Meal.from_dict(m) for m in cfg.get("meals", [])]
    return engine, meals, stomach


def run_simulation(
    cfg_overrides: ParameterOverrides | None = None,
    meals: Iterable[Mapping[str, Any]] | None = None,
) -> Dict[str, Any]:
    base_cfg = load_base_config()
    cfg = copy.deepcopy(base_cfg)
    apply_overrides(cfg, cfg_overrides)
    if meals is not None:
        cfg["meals"] = [dict(m) for m in meals]

    engine, meals_cfg, stomach = build_engine(cfg)
    if meals is not None:
        meals_cfg = [Meal.from_dict(m) for m in meals]

    t_end = float(cfg.get("t_end_min", 240))
    dt = float(cfg.get("dt_min", 1.0))

    t = 0.0
    snapshots: List[Dict[str, Any]] = []

    while t <= t_end + 1e-9:
        apply_meals_to_stomach(t, dt, meals_cfg, stomach)
        _, amounts, concentrations, volumes = capture_compartment_state(engine.compartments)
        snapshot = {
            "time": t,
            "compartments": amounts,
            "concentrations": concentrations,
            "volumes_L": volumes,
            "signals": {sig.name: sig.level for sig in engine.signals.values()},
        }
        snapshots.append(snapshot)
        engine.step(t, dt)
        t += dt

    return {
        "snapshots": snapshots,
        "metadata": {
            "t_end_min": t_end,
            "dt_min": dt,
            "parameter_defaults": {key: flatten_config_value(base_cfg, key) for key in _PARAM_KEYS},
            "parameter_controls": _PARAM_KEYS,
            "blood_volume_L": float(base_cfg.get("blood", {}).get("volume_L", 5.0)),
        },
    }


def exposed_parameters() -> Dict[str, Any]:
    base_cfg = load_base_config()
    defaults = {key: flatten_config_value(base_cfg, key) for key in _PARAM_KEYS}
    return {
        "defaults": defaults,
        "controls": _PARAM_KEYS,
        "meal": base_cfg.get("meals", [{}])[0] if base_cfg.get("meals") else {
            "time_min": 10.0,
            "grams_carbohydrate": 50.0,
            "grams_fat": 10.0,
            "grams_protein": 15.0,
            "grams_soluble_fiber": 5.0,
        },
    }
