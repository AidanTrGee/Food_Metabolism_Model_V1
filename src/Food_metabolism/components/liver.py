from __future__ import annotations

from ..core.base import Compartment, Signal, sat_mm
from ..metabolism import (
    compute_atp_demand,
    compute_beta_oxidation,
    compute_citric_acid_cycle,
    compute_fatty_acid_synthesis,
    compute_glycolysis,
    compute_gluconeogenesis,
    compute_oxidative_phosphorylation,
)


class Liver(Compartment):
    def __init__(self, params=None):
        params = params or {}
        nutrients = ["glucose", "pyruvate", "lactate", "lipid", "acetyl_coa"]
        volume = float(params.get("volume_L", 1.5))
        super().__init__("Liver", nutrients, params, volume_L=max(volume, 1e-6))
        glycogen_capacity = float(params.get("glycogen_capacity_g", 450.0))
        initial_glycogen = float(params.get("initial_glycogen_g", 300.0))
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
        self._run_gluconeogenesis(dt_min, signals)

        gly_params = None
        raw_cfg = self.params.get("glycolysis") if isinstance(self.params, dict) else None
        if isinstance(raw_cfg, dict):
            gly_params = raw_cfg

        result = compute_glycolysis(
            tissue="liver",
            glucose_g=self.get("glucose"),
            volume_L=self.volume_L,
            dt_min=dt_min,
            signals=signals,
            params=gly_params,
        )

        if result.glucose_used_g > 0.0:
            self.add("glucose", -result.glucose_used_g)
            if result.pyruvate_g > 0.0:
                self.add("pyruvate", result.pyruvate_g)
            if result.lactate_g > 0.0:
                self.add("lactate", result.lactate_g)

        self._store_glycogen(dt_min, signals)
        self.glycogen_g = max(0.0, min(self.glycogen_g, self.glycogen_capacity_g))

        glyc_sig = signals.setdefault("liver_glycogen_g", Signal("liver_glycogen_g", 0.0))
        glyc_sig.level = self.glycogen_g

        self._run_beta_oxidation(dt_min, signals)
        self._update_energy_state(dt_min, signals)

        self._run_lipogenesis(dt_min, signals)

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

        glucagon = max(self._signal_level(signals, "glucagon", 0.05), 0.0)
        insulin = max(self._signal_level(signals, "insulin", 0.05), 0.0)

        vmax = float(cfg.get("mobilize_vmax_gpm", 0.25))
        km = float(cfg.get("mobilize_km_g", 40.0))
        glucagon_gain = float(cfg.get("mobilize_glucagon_gain", 2.5))
        insulin_supp = float(cfg.get("mobilize_insulin_suppression", 1.5))

        multiplier = 1.0 + glucagon_gain * max(0.0, glucagon - 0.05)
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
        glucagon = max(self._signal_level(signals, "glucagon", 0.05), 0.0)

        reserve = float(cfg.get("store_glucose_reserve_g", 5.0))
        available = max(self.get("glucose") - reserve, 0.0)
        capacity = max(self.glycogen_capacity_g - self.glycogen_g, 0.0)
        if available <= 0.0 or capacity <= 0.0:
            return 0.0

        vmax = float(cfg.get("store_vmax_gpm", 0.4))
        km = float(cfg.get("store_km_g", 20.0))
        insulin_gain = float(cfg.get("store_insulin_gain", 4.5))
        glucagon_supp = float(cfg.get("store_glucagon_suppression", 2.5))

        multiplier = 1.0 + insulin_gain * max(0.0, insulin - 0.05)
        multiplier *= max(0.0, 1.0 - glucagon_supp * max(0.0, glucagon - 0.05))
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

    def _run_gluconeogenesis(self, dt_min: float, signals):
        cfg: dict[str, float] | None = None
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("gluconeogenesis")
            if isinstance(raw_cfg, dict):
                cfg = {k: float(v) for k, v in raw_cfg.items()}

        precursors = {
            "pyruvate": self.get("pyruvate"),
            "lactate": self.get("lactate"),
        }

        result = compute_gluconeogenesis(
            tissue="liver",
            precursors_g=precursors,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )

        if result.glucose_g <= 0.0:
            return result

        self.add("glucose", result.glucose_g)
        for name, used in result.precursor_used_g.items():
            if used > 0.0 and name in self.pools:
                self.add(name, -used)

        return result

    def _lipogenesis_cfg(self) -> dict[str, float]:
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("lipogenesis")
            if isinstance(raw_cfg, dict):
                return {k: float(v) for k, v in raw_cfg.items()}
        return {}

    def _run_lipogenesis(self, dt_min: float, signals):
        if dt_min <= 0.0:
            return None

        cfg = self._lipogenesis_cfg()

        glucose_reserve = float(cfg.get("glucose_reserve_g", 40.0))
        pyruvate_reserve = float(cfg.get("pyruvate_reserve_g", 5.0))
        lactate_reserve = float(cfg.get("lactate_reserve_g", 5.0))

        available_glucose = max(self.get("glucose") - glucose_reserve, 0.0)
        available_pyruvate = max(self.get("pyruvate") - pyruvate_reserve, 0.0)
        available_lactate = max(self.get("lactate") - lactate_reserve, 0.0)

        max_glucose = float(cfg.get("max_glucose_for_lipogenesis_g", available_glucose))
        max_pyruvate = float(cfg.get("max_pyruvate_for_lipogenesis_g", available_pyruvate))
        max_lactate = float(cfg.get("max_lactate_for_lipogenesis_g", available_lactate))

        substrates: dict[str, float] = {}
        if available_glucose > 0.0 and max_glucose > 0.0:
            substrates["glucose"] = min(available_glucose, max_glucose)
        if available_pyruvate > 0.0 and max_pyruvate > 0.0:
            substrates["pyruvate"] = min(available_pyruvate, max_pyruvate)
        if available_lactate > 0.0 and max_lactate > 0.0:
            substrates["lactate"] = min(available_lactate, max_lactate)

        lipogenesis_sig = signals.setdefault("liver_lipogenesis_gpm", Signal("liver_lipogenesis_gpm", 0.0))
        liver_lipid_sig = signals.setdefault("liver_lipid_g", Signal("liver_lipid_g", 0.0))

        if not substrates:
            lipogenesis_sig.level = 0.0
            liver_lipid_sig.level = self.get("lipid")
            return None

        result = compute_fatty_acid_synthesis(
            tissue="liver",
            substrates_g=substrates,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )

        lipogenesis_sig.level = result.rate_gpm

        if result.lipid_g <= 0.0:
            liver_lipid_sig.level = self.get("lipid")
            return result

        for name, used in result.substrate_used_g.items():
            if used <= 0.0 or name not in self.pools:
                continue
            current = self.get(name)
            self.add(name, -min(used, current))

        self.add("lipid", result.lipid_g)
        liver_lipid_sig.level = self.get("lipid")

        return result

    def _run_beta_oxidation(self, dt_min: float, signals):
        beta_sig = signals.setdefault("liver_beta_oxidation_gpm", Signal("liver_beta_oxidation_gpm", 0.0))
        ketone_sig = signals.setdefault("liver_ketone_precursor_gpm", Signal("liver_ketone_precursor_gpm", 0.0))
        acetyl_sig = signals.setdefault("liver_acetyl_coa_g", Signal("liver_acetyl_coa_g", 0.0))

        if dt_min <= 0.0:
            beta_sig.level = 0.0
            ketone_sig.level = 0.0
            acetyl_sig.level = self.get("acetyl_coa")
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
            tissue="liver",
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
        if result.acetyl_coa_g > 0.0:
            self.add("acetyl_coa", result.acetyl_coa_g)
            self.pools["acetyl_coa"].amount_g = max(self.pools["acetyl_coa"].amount_g, 0.0)

        acetyl_sig.level = self.get("acetyl_coa")

        return result

    def _energy_cfg(self) -> dict[str, float]:
        if isinstance(self.params, dict):
            raw_cfg = self.params.get("energy")
            if isinstance(raw_cfg, dict):
                return {k: float(v) for k, v in raw_cfg.items()}
        return {}

    def _update_energy_state(self, dt_min: float, signals):
        if dt_min <= 0.0:
            return None

        cfg = self._energy_cfg()

        demand = compute_atp_demand(
            tissue="liver",
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )

        pyruvate = self.get("pyruvate") if "pyruvate" in self.pools else 0.0
        acetyl = self.get("acetyl_coa") if "acetyl_coa" in self.pools else 0.0

        cac_result = compute_citric_acid_cycle(
            tissue="liver",
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
            tissue="liver",
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

        atp_sig = signals.setdefault("liver_atp_mmol", Signal("liver_atp_mmol", 0.0))
        atp_sig.level = self.atp_mmol
        deficit_sig = signals.setdefault("liver_atp_deficit_mmol", Signal("liver_atp_deficit_mmol", 0.0))
        deficit_sig.level = self.atp_deficit_mmol

        production_sig = signals.setdefault("liver_atp_production_mmol_per_min", Signal("liver_atp_production_mmol_per_min", 0.0))
        production_sig.level = produced / max(dt_min, 1e-9)
        cac_sig = signals.setdefault("liver_cac_rate_gpm", Signal("liver_cac_rate_gpm", 0.0))
        cac_sig.level = cac_result.rate_acetyl_coa_gpm + cac_result.rate_pyruvate_gpm
        oxphos_sig = signals.setdefault("liver_oxphos_atp_deficit_mmol", Signal("liver_oxphos_atp_deficit_mmol", 0.0))
        oxphos_sig.level = oxphos_result.deficit_mmol

        return {
            "demand": demand,
            "citric_acid_cycle": cac_result,
            "oxidative_phosphorylation": oxphos_result,
        }
