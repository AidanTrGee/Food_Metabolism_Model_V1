from __future__ import annotations

from ..core.base import Compartment, Signal, sat_mm
from ..metabolism import (
    compute_atp_demand,
    compute_beta_oxidation,
    compute_citric_acid_cycle,
    compute_glycolysis,
    compute_lactic_fermentation,
    compute_oxidative_phosphorylation,
)


class Muscle(Compartment):
    def __init__(self, params=None):
        params = params or {}
        nutrients = ["glucose", "pyruvate", "lactate", "lipid", "acetyl_coa"]
        volume = float(params.get("volume_L", 12.0))
        super().__init__("Muscle", nutrients, params, volume_L=max(volume, 1e-6))
        glycogen_capacity = float(params.get("glycogen_capacity_g", 550.0))
        initial_glycogen = float(params.get("initial_glycogen_g", 350.0))
        self.glycogen_capacity_g = max(glycogen_capacity, 0.0)
        self.glycogen_g = max(0.0, min(initial_glycogen, self.glycogen_capacity_g))
        
        energy_cfg = self._energy_cfg()
        capacity = max(float(energy_cfg.get("max_atp_pool_mmol", 0.0)), 0.0)
        initial_atp = float(energy_cfg.get("initial_atp_mmol", capacity * 0.5 if capacity > 0.0 else 0.0))
        self.atp_capacity_mmol = capacity
        self.atp_mmol = max(0.0, min(initial_atp, self.atp_capacity_mmol))
        self.atp_deficit_mmol = 0.0

    def step_endogenous(self, t_min, dt_min, signals):
        self._mobilize_glycogen(dt_min, signals)

        gly_params = None
        raw_cfg = self.params.get("glycolysis") if isinstance(self.params, dict) else None
        if isinstance(raw_cfg, dict):
            gly_params = raw_cfg

        result = compute_glycolysis(
            tissue="muscle",
            glucose_g=self.get("glucose"),
            volume_L=self.volume_L,
            dt_min=dt_min,
            signals=signals,
            params=gly_params,
        )

        lactate_request_g = 0.0
        if result.glucose_used_g > 0.0:
            self.add("glucose", -result.glucose_used_g)
            total_pyruvate = result.pyruvate_g + result.lactate_g
            if total_pyruvate > 0.0:
                self.add("pyruvate", total_pyruvate)
            lactate_request_g = max(result.lactate_g, 0.0)

        self._run_lactic_fermentation(dt_min, signals, lactate_request_g)
        self._run_beta_oxidation(dt_min, signals)
        self._update_energy_state(dt_min, signals)

        self._store_glycogen(dt_min, signals)
        self.glycogen_g = max(0.0, min(self.glycogen_g, self.glycogen_capacity_g))

        glyc_sig = signals.setdefault("muscle_glycogen_g", Signal("muscle_glycogen_g", 0.0))
        glyc_sig.level = self.glycogen_g

    def _signal_level(self, signals, name: str, default: float = 0.0) -> float:
        signal = signals.get(name) if signals is not None else None
        return float(default) if signal is None else float(signal.level)

    def _mobilize_glycogen(self, dt_min: float, signals):
        if dt_min <= 0.0 or self.glycogen_g <= 0.0:
            return 0.0

        cfg: dict[str, float] = {}
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("glycogen")
            if isinstance(raw_cfg, dict):
                cfg = {k: float(v) for k, v in raw_cfg.items()}

        exercise = max(self._signal_level(signals, "exercise", 0.0), 0.0)
        glucagon = max(self._signal_level(signals, "glucagon", 0.05), 0.0)
        insulin = max(self._signal_level(signals, "insulin", 0.05), 0.0)

        vmax = float(cfg.get("mobilize_vmax_gpm", 0.45))
        km = float(cfg.get("mobilize_km_g", 60.0))
        exercise_gain = float(cfg.get("mobilize_exercise_gain", 4.5))
        glucagon_gain = float(cfg.get("mobilize_glucagon_gain", 1.5))
        insulin_supp = float(cfg.get("mobilize_insulin_suppression", 2.2))

        multiplier = 1.0
        multiplier += exercise_gain * max(0.0, exercise)
        multiplier += glucagon_gain * max(0.0, glucagon - 0.05)
        multiplier *= max(0.0, 1.0 - insulin_supp * max(0.0, insulin - 0.05))
        vmax *= max(multiplier, 0.0)

        if vmax <= 0.0:
            return 0.0

        available = self.glycogen_g
        rate = sat_mm(available, vmax, km)
        dt_safe = max(dt_min, 1e-9)
        rate = min(rate, available / dt_safe)
        mobilized = rate * dt_min
        if mobilized <= 0.0:
            return 0.0
        self.glycogen_g -= mobilized
        self.add("glucose", mobilized)
        return mobilized

    def _store_glycogen(self, dt_min: float, signals):
        if dt_min <= 0.0:
            return 0.0

        cfg: dict[str, float] = {}
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("glycogen")
            if isinstance(raw_cfg, dict):
                cfg = {k: float(v) for k, v in raw_cfg.items()}

        insulin = max(self._signal_level(signals, "insulin", 0.05), 0.0)
        exercise = max(self._signal_level(signals, "exercise", 0.0), 0.0)

        reserve = float(cfg.get("store_glucose_reserve_g", 1.0))
        available = max(self.get("glucose") - reserve, 0.0)
        capacity = max(self.glycogen_capacity_g - self.glycogen_g, 0.0)
        if available <= 0.0 or capacity <= 0.0:
            return 0.0

        vmax = float(cfg.get("store_vmax_gpm", 0.55))
        km = float(cfg.get("store_km_g", 25.0))
        insulin_gain = float(cfg.get("store_insulin_gain", 6.0))
        exercise_supp = float(cfg.get("store_exercise_suppression", 3.0))

        multiplier = 1.0 + insulin_gain * max(0.0, insulin - 0.05)
        multiplier *= max(0.0, 1.0 - exercise_supp * exercise)
        vmax *= max(multiplier, 0.0)

        if vmax <= 0.0:
            return 0.0

        rate = sat_mm(available, vmax, km)
        dt_safe = max(dt_min, 1e-9)
        rate = min(rate, available / dt_safe, capacity / dt_safe)
        stored = rate * dt_min
        if stored <= 0.0:
            return 0.0
        self.add("glucose", -stored)
        self.glycogen_g += stored
        return stored

    def _run_beta_oxidation(self, dt_min: float, signals):
        beta_sig = signals.setdefault("muscle_beta_oxidation_gpm", Signal("muscle_beta_oxidation_gpm", 0.0))
        ketone_sig = signals.setdefault("muscle_ketone_precursor_gpm", Signal("muscle_ketone_precursor_gpm", 0.0))
        acetyl_sig = signals.setdefault("muscle_acetyl_coa_g", Signal("muscle_acetyl_coa_g", 0.0))

        if dt_min <= 0.0 or "lipid" not in self.pools:
            beta_sig.level = 0.0
            ketone_sig.level = 0.0
            acetyl_sig.level = self.get("acetyl_coa") if "acetyl_coa" in self.pools else 0.0
            return None

        cfg: dict[str, float] | None = None
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("beta_oxidation")
            if isinstance(raw_cfg, dict):
                cfg = {k: float(v) for k, v in raw_cfg.items()}

        lipid_available = self.get("lipid")
        if lipid_available <= 0.0:
            beta_sig.level = 0.0
            ketone_sig.level = 0.0
            acetyl_sig.level = self.get("acetyl_coa")
            return None

        result = compute_beta_oxidation(
            tissue="muscle",
            lipid_g=lipid_available,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )

        beta_sig.level = result.rate_gpm
        ketone_sig.level = result.ketone_precursor_g / max(dt_min, 1e-9)

        if result.lipid_used_g > 0.0:
            self.add("lipid", -result.lipid_used_g)
            self.pools["lipid"].amount_g = max(self.pools["lipid"].amount_g, 0.0)
        if "acetyl_coa" in self.pools and result.acetyl_coa_g > 0.0:
            self.add("acetyl_coa", result.acetyl_coa_g)
            self.pools["acetyl_coa"].amount_g = max(self.pools["acetyl_coa"].amount_g, 0.0)

        acetyl_sig.level = self.get("acetyl_coa")

        return result

    def _fermentation_cfg(self) -> dict[str, float]:
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("lactic_fermentation")
            if isinstance(raw_cfg, dict):
                return {k: float(v) for k, v in raw_cfg.items()}
        return {}
    
    def _energy_cfg(self) -> dict[str, float]:
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("energy")
            if isinstance(raw_cfg, dict):
                return {k: float(v) for k, v in raw_cfg.items()}
        return {}
    
    def _run_lactic_fermentation(self, dt_min: float, signals, requested_pyruvate_g: float):
        ferm_sig = signals.setdefault("muscle_lactic_fermentation_gpm", Signal("muscle_lactic_fermentation_gpm", 0.0))
        lactate_sig = signals.setdefault("muscle_lactate_g", Signal("muscle_lactate_g", 0.0))
        nad_sig = signals.setdefault("muscle_fermentation_nad_regen_mmol_per_min", Signal("muscle_fermentation_nad_regen_mmol_per_min", 0.0))

        if dt_min <= 0.0:
            ferm_sig.level = 0.0
            nad_sig.level = 0.0
            lactate_sig.level = self.get("lactate")
            return None

        available = max(self.get("pyruvate"), 0.0)
        if requested_pyruvate_g > 0.0:
            available = min(available, max(requested_pyruvate_g, 0.0))

        if available <= 0.0:
            ferm_sig.level = 0.0
            nad_sig.level = 0.0
            lactate_sig.level = self.get("lactate")
            return None

        cfg = self._fermentation_cfg()

        result = compute_lactic_fermentation(
            tissue="muscle",
            pyruvate_g=available,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )

        ferm_sig.level = result.rate_gpm
        nad_sig.level = result.nad_regenerated_mmol / max(dt_min, 1e-9)

        if result.pyruvate_used_g > 0.0:
            use_pyruvate = min(result.pyruvate_used_g, self.get("pyruvate"))
            if use_pyruvate > 0.0:
                self.add("pyruvate", -use_pyruvate)
                self.pools["pyruvate"].amount_g = max(self.pools["pyruvate"].amount_g, 0.0)

        if result.lactate_g > 0.0:
            self.add("lactate", result.lactate_g)

        lactate_sig.level = self.get("lactate")

        return result

    def _update_energy_state(self, dt_min: float, signals):
        if dt_min <= 0.0:
            return None
        
        cfg = self._energy_cfg()
        
        demand = compute_atp_demand(
            tissue="muscle",
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )
        
        pyruvate = self.get("pyruvate") if "pyruvate" in self.pools else 0.0
        acetyl = self.get("acetyl_coa") if "acetyl_coa" in self.pools else 0.0
        
        cac_result = compute_citric_acid_cycle(
            tissue="muscle",
            pyruvate_g=pyruvate,
            acetyl_coa_g=acetyl,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )
        
        if cac_result.acetyl_coa_used_g > 0.0 and "acetyl_coa" in self.pools:
            use_acetyl = min(cac_result.acetyl_coa_used_g, acetyl)
            if use_acetyl > 0.0:
                self.add("acetyl_coa", -use_acetyl)
                self.pools["acetyl_coa"].amount_g = max(self.pools["acetyl_coa"].amount_g, 0.0)
        
        if cac_result.pyruvate_used_g > 0.0 and "pyruvate" in self.pools:
            use_pyruvate = min(cac_result.pyruvate_used_g, pyruvate)
            if use_pyruvate > 0.0:
                self.add("pyruvate", -use_pyruvate)
                self.pools["pyruvate"].amount_g = max(self.pools["pyruvate"].amount_g, 0.0)
        
        remaining_demand = max(demand.demand_mmol - cac_result.atp_generated_mmol, 0.0)
        oxphos_result = compute_oxidative_phosphorylation(
            tissue="muscle",
            nadh_mmol=cac_result.nadh_generated_mmol,
            atp_demand_mmol=remaining_demand,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )
        
        produced = cac_result.atp_generated_mmol + oxphos_result.atp_generated_mmol
        required = demand.demand_mmol
        pool_before = self.atp_mmol
        total_available = produced + pool_before
        if total_available >= required:
            pool_after = total_available - required
            deficit = 0.0
        else:
            pool_after = 0.0
            deficit = required - total_available
        
        self.atp_mmol = max(0.0, min(self.atp_capacity_mmol, pool_after))
        self.atp_deficit_mmol = max(deficit, 0.0)
        
        atp_sig = signals.setdefault("muscle_atp_mmol", Signal("muscle_atp_mmol", 0.0))
        atp_sig.level = self.atp_mmol
        deficit_sig = signals.setdefault("muscle_atp_deficit_mmol", Signal("muscle_atp_deficit_mmol", 0.0))
        deficit_sig.level = self.atp_deficit_mmol
        
        production_sig = signals.setdefault("muscle_atp_production_mmol_per_min", Signal("muscle_atp_production_mmol_per_min", 0.0))
        production_sig.level = produced / max(dt_min, 1e-9)
        cac_sig = signals.setdefault("muscle_cac_rate_gpm", Signal("muscle_cac_rate_gpm", 0.0))
        cac_sig.level = cac_result.rate_acetyl_coa_gpm + cac_result.rate_pyruvate_gpm
        oxphos_sig = signals.setdefault("muscle_oxphos_atp_deficit_mmol", Signal("muscle_oxphos_atp_deficit_mmol", 0.0))
        oxphos_sig.level = oxphos_result.deficit_mmol
        
        return {
            "demand": demand,
            "citric_acid_cycle": cac_result,
            "oxidative_phosphorylation": oxphos_result,
        }
