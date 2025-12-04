from __future__ import annotations

from ..core.base import Compartment, Signal
from ..metabolism import compute_gluconeogenesis


class Kidney(Compartment):

    def __init__(self, params=None):
        params = params or {}
        volume = float(params.get("volume_L", 0.3))
        super().__init__(
            "Kidney",
            ["glucose", "pyruvate", "lactate"],
            params,
            volume_L=max(volume, 1e-6),
        )
        self.excreted_total_g = 0.0
        self.cortical_glucose_g = max(float(params.get("initial_cortical_glucose_g", 0.0)), 0.0)
        if self.cortical_glucose_g > 0.0:
            self.add("glucose", self.cortical_glucose_g)

    def step_endogenous(self, t_min, dt_min, signals):
        gng = self._run_gluconeogenesis(dt_min, signals)
        if gng.glucose_g > 0.0:
            self.cortical_glucose_g += gng.glucose_g
            signal = signals.setdefault(
                "renal_gluconeogenesis_gpm",
                Signal("renal_gluconeogenesis_gpm", 0.0),
            )
            signal.level = gng.rate_gpm

        amount = self.get("glucose")
        if amount <= 0.0:
            return

        reserve = min(self.cortical_glucose_g, amount)
        excreted = max(amount - reserve, 0.0)
        if excreted > 0.0:
            self.excreted_total_g += excreted
        self.pools["glucose"].amount_g = reserve
        self.cortical_glucose_g = min(self.cortical_glucose_g, reserve)
        signal = signals.setdefault(
            "renal_glucose_excreted_total_g",
            Signal("renal_glucose_excreted_total_g", 0.0),
        )
        signal.level = self.excreted_total_g
        cortex_sig = signals.setdefault(
            "renal_cortical_glucose_g",
            Signal("renal_cortical_glucose_g", 0.0),
        )
        cortex_sig.level = self.cortical_glucose_g

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
            tissue="kidney",
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
