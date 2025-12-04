from __future__ import annotations

from ..core.base import Compartment, Signal
from ..metabolism import compute_fatty_acid_synthesis, compute_lipolysis


class Adipose(Compartment):
    def __init__(self, params=None):
        params = params or {}
        nutrients = ["glucose", "lipid", "ffa", "glycerol"]
        volume = float(params.get("volume_L", 18.0))
        super().__init__("Adipose", nutrients, params, volume_L=max(volume, 1e-6))

    def step_endogenous(self, t_min, dt_min, signals):
        self._run_lipolysis(dt_min, signals)
        self._run_lipogenesis(dt_min, signals)

    def _section_cfg(self, name: str) -> dict[str, float]:
        if not isinstance(self.params, dict):
            return {}
        section = self.params.get(name)
        if isinstance(section, dict):
            return {k: float(v) for k, v in section.items()}
        return {}

    def _run_lipolysis(self, dt_min: float, signals):
        if dt_min <= 0.0:
            return None

        cfg = self._section_cfg("lipolysis")
        stored = max(self.get("lipid"), 0.0)
        result = compute_lipolysis(
            tissue="adipose",
            stored_lipid_g=stored,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )

        lipolysis_sig = signals.setdefault("adipose_lipolysis_gpm", Signal("adipose_lipolysis_gpm", 0.0))
        ffa_sig = signals.setdefault("adipose_ffa_release_gpm", Signal("adipose_ffa_release_gpm", 0.0))
        glycerol_sig = signals.setdefault("adipose_glycerol_release_gpm", Signal("adipose_glycerol_release_gpm", 0.0))

        lipolysis_sig.level = result.rate_gpm

        if result.lipid_hydrolyzed_g <= 0.0:
            ffa_sig.level = 0.0
            glycerol_sig.level = 0.0
            return result

        dt_safe = max(dt_min, 1e-9)
        ffa_sig.level = result.ffa_g / dt_safe if result.ffa_g > 0.0 else 0.0
        glycerol_sig.level = result.glycerol_g / dt_safe if result.glycerol_g > 0.0 else 0.0

        if result.lipid_hydrolyzed_g > 0.0:
            self.add("lipid", -min(result.lipid_hydrolyzed_g, stored))
        if result.ffa_g > 0.0:
            self.add("ffa", result.ffa_g)
        if result.glycerol_g > 0.0:
            self.add("glycerol", result.glycerol_g)

        adipose_ffa_sig = signals.setdefault("adipose_ffa_pool_g", Signal("adipose_ffa_pool_g", 0.0))
        adipose_ffa_sig.level = self.get("ffa")

        return result

    def _run_lipogenesis(self, dt_min: float, signals):
        if dt_min <= 0.0:
            return None

        cfg = self._section_cfg("lipogenesis")

        glucose_reserve = float(cfg.get("glucose_reserve_g", 4.0))
        available_glucose = max(self.get("glucose") - glucose_reserve, 0.0)

        ffa_reesterify_fraction = float(cfg.get("ffa_reesterify_fraction", 0.4))
        ffa_reesterify_fraction = max(0.0, min(1.0, ffa_reesterify_fraction))
        available_ffa = max(self.get("ffa"), 0.0) * ffa_reesterify_fraction

        substrates: dict[str, float] = {}
        if available_glucose > 0.0:
            substrates["glucose"] = available_glucose
        if available_ffa > 0.0:
            substrates["ffa"] = available_ffa

        if not substrates:
            signals.setdefault("adipose_lipogenesis_gpm", Signal("adipose_lipogenesis_gpm", 0.0)).level = 0.0
            return None

        result = compute_fatty_acid_synthesis(
            tissue="adipose",
            substrates_g=substrates,
            dt_min=dt_min,
            signals=signals,
            params=cfg,
        )

        lipogenesis_sig = signals.setdefault("adipose_lipogenesis_gpm", Signal("adipose_lipogenesis_gpm", 0.0))
        lipogenesis_sig.level = result.rate_gpm

        if result.lipid_g <= 0.0:
            return result

        for name, used in result.substrate_used_g.items():
            if used <= 0.0 or name not in self.pools:
                continue
            current = self.get(name)
            self.add(name, -min(used, current))

        self.add("lipid", result.lipid_g)

        adipose_lipid_sig = signals.setdefault("adipose_lipid_g", Signal("adipose_lipid_g", 0.0))
        adipose_lipid_sig.level = self.get("lipid")

        adipose_ffa_sig = signals.setdefault("adipose_ffa_pool_g", Signal("adipose_ffa_pool_g", 0.0))
        adipose_ffa_sig.level = self.get("ffa")

        return result
