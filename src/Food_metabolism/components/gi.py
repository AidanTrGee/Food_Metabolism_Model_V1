from __future__ import annotations

import math
from typing import Any, Dict, Mapping, Optional

from ..core.base import Compartment, Signal

MOLAR_MASS_GLUCOSE = 180.156  # g/mol

ENERGY_KCAL_PER_G: Dict[str, float] = {
	"glucose": 4.0,
	"lipid": 9.0,
	"protein": 4.0,
}

SPECIFIC_VOLUME_L_PER_G: Dict[str, float] = {
	"glucose": 1.0e-3,
	"lipid": 1.1e-3,
	"protein": 1.0e-3,
	"soluble_fiber": 1.2e-3,
	"water": 1.0e-3,
}


def _clamp(value: float, lower: float, upper: float) -> float:
	return max(lower, min(upper, value))


def _numeric_subset(config: Mapping[str, Any]) -> Dict[str, float]:
	return {
		str(name): float(amount)
		for name, amount in config.items()
		if isinstance(name, str) and isinstance(amount, (int, float))
	}


def _get_numeric_mapping(config: Mapping[str, Any], key: str) -> Dict[str, float]:
	value = config.get(key, {})
	if isinstance(value, Mapping):
		return _numeric_subset(value)
	return {}


def _get_float(config: Mapping[str, Any], key: str, default: float) -> float:
	value = config.get(key, default)
	if isinstance(value, (int, float)):
		return float(value)
	return default


class StomachLumen(Compartment):

	def __init__(self, params: Optional[Dict[str, float]] = None):
		raw_params: Dict[str, Any] = dict(params or {})
		nutrients = ["glucose", "lipid", "protein", "soluble_fiber", "water"]
		#should we add fructose for this first model or wait until we have results first. 

		numeric_params = _numeric_subset(raw_params)
		residual_volume = _get_float(
			raw_params,
			"residual_volume_L",
			_get_float(raw_params, "volume_L", 0.075),
		)
		residual_volume = max(residual_volume, 1e-6)
		super().__init__("StomachLumen", nutrients, numeric_params, volume_L=residual_volume)

		macro_cfg = _get_numeric_mapping(raw_params, "macro_modifiers")
		fiber_cfg = _get_numeric_mapping(raw_params, "fiber")
		lag_cfg = _get_numeric_mapping(raw_params, "lag")

		self.base_energy_emptying_kcal_min = _get_float(raw_params, "base_energy_emptying_kcal_min", 2.3)
		self.max_emptying_rate_per_min = _get_float(raw_params, "max_emptying_rate_per_min", 2.5)
		self.min_energy_kcal = _get_float(raw_params, "min_energy_kcal", 0.05)
		self.fat_coeff = float(macro_cfg.get("fat_coeff", 0.7))
		self.protein_coeff = float(macro_cfg.get("protein_coeff", 0.2))

		self.fiber_reference_g = float(fiber_cfg.get("reference_g", 4.0))
		self.fiber_beta = float(fiber_cfg.get("beta", fiber_cfg.get("alpha", 0.25)))
		self.fiber_min_modifier = float(fiber_cfg.get("min_modifier", 0.35))

		self.lag_base_min = float(lag_cfg.get("base_min", 10.0))
		self.lag_fat_coeff_min = float(lag_cfg.get("fat_coeff_min", 8.0))
		self.lag_fiber_coeff_min = float(lag_cfg.get("fiber_coeff_min", 4.0))

		self.residual_volume_L = residual_volume
		self.min_volume_L = _get_float(raw_params, "min_volume_L", 5e-3)
		self.specific_volume_L_per_g = dict(SPECIFIC_VOLUME_L_PER_G)
		spec_cfg = _get_numeric_mapping(raw_params, "specific_volume_L_per_g")
		for nutrient, value in spec_cfg.items():
			self.specific_volume_L_per_g[nutrient] = max(value, 0.0)

		self.last_meal_time = -1e9
		self.lag_ready_time = 0.0

		self._update_volume()

	def step_endogenous(self, t_min: float, dt_min: float, signals: Dict[str, Signal]):
		self._update_volume()

	def ingest_meal(self, t_event_min: float, additions_g: Dict[str, float]):
		self.add_many(additions_g)
		self.last_meal_time = t_event_min
		self.lag_ready_time = max(self.lag_ready_time, t_event_min + self._compute_lag_duration())

	def energy_content_kcal(self) -> float:
		return sum(ENERGY_KCAL_PER_G.get(n, 0.0) * pool.amount_g for n, pool in self.pools.items())

	def _energy_fraction(self, nutrient: str) -> float:
		total = self.energy_content_kcal()
		if total <= 1e-9:
			return 0.0
		return ENERGY_KCAL_PER_G.get(nutrient, 0.0) * self.get(nutrient) / total

	def fiber_equivalent(self) -> float:
		reference = max(self.fiber_reference_g, 1e-6)
		return self.get("soluble_fiber") / reference

	def viscosity_modifier(self) -> float:
		ratio = self.fiber_equivalent()
		modifier = 1.0 / (1.0 + self.fiber_beta * ratio)
		return _clamp(modifier, self.fiber_min_modifier, 1.0)

	def emptying_rate_per_min(self, t_min: float) -> float:
		energy = self.energy_content_kcal()
		if energy <= self.min_energy_kcal:
			return 0.0
		if t_min < self.lag_ready_time:
			return 0.0

		fat_frac = self._energy_fraction("lipid")
		prot_frac = self._energy_fraction("protein")
		macro_modifier = math.exp(-self.fat_coeff * fat_frac - self.protein_coeff * prot_frac)
		macro_modifier = _clamp(macro_modifier, 0.0, 1.0)

		energy_flow_kcal_min = self.base_energy_emptying_kcal_min * macro_modifier * self.viscosity_modifier()
		rate = energy_flow_kcal_min / max(energy, self.min_energy_kcal)
		return _clamp(rate, 0.0, self.max_emptying_rate_per_min)

	def _compute_lag_duration(self) -> float:
		fat_frac = self._energy_fraction("lipid")
		fiber_ratio = self.fiber_equivalent()
		lag = self.lag_base_min + self.lag_fat_coeff_min * fat_frac + self.lag_fiber_coeff_min * fiber_ratio
		return max(0.0, lag)

	def _food_volume_L(self) -> float:
		total = self.residual_volume_L
		for nutrient, pool in self.pools.items():
			spec = self.specific_volume_L_per_g.get(nutrient, SPECIFIC_VOLUME_L_PER_G.get(nutrient, 1.0e-3))
			total += max(pool.amount_g, 0.0) * spec
		return max(self.min_volume_L, total)

	def _update_volume(self) -> None:
		self.volume_L = self._food_volume_L()

	def add(self, nutrient: str, grams: float):
		super().add(nutrient, grams)
		self._update_volume()

	def add_many(self, additions: Dict[str, float]):
		super().add_many(additions)
		self._update_volume()


class StomachTissue(Compartment):

	def __init__(self, params: Optional[Dict[str, float]] = None):
		params = params or {}
		super().__init__("StomachTissue", ["glucose"], params, volume_L=float(params.get("volume_L", 0.015)))

	def step_endogenous(self, t_min: float, dt_min: float, signals: Dict[str, Signal]):
		return


class IntestinalLumen(Compartment):

	def __init__(self, params: Optional[Dict[str, float]] = None):
		raw_params: Dict[str, Any] = dict(params or {})
		nutrients = ["glucose", "lipid", "protein", "soluble_fiber"]
		numeric_params = _numeric_subset(raw_params)
		baseline_volume = _get_float(
			raw_params,
			"baseline_volume_L",
			_get_float(raw_params, "volume_L", _get_float(raw_params, "lumen_volume_L", 0.35)),
		)
		baseline_volume = max(baseline_volume, 1e-6)
		super().__init__("IntestinalLumen", nutrients, numeric_params, volume_L=baseline_volume)

		fiber_cfg = _get_numeric_mapping(raw_params, "fiber")
		self.fiber_reference_g = float(fiber_cfg.get("reference_g", 4.0))
		self.fiber_alpha = float(fiber_cfg.get("alpha", 0.25))
		self.fiber_min_modifier = float(fiber_cfg.get("min_modifier", 0.4))
		#need to fix this

		self.residual_volume_L = baseline_volume
		self.min_volume_L = _get_float(raw_params, "min_volume_L", 5e-3)
		self.specific_volume_L_per_g = dict(SPECIFIC_VOLUME_L_PER_G)
		spec_cfg = _get_numeric_mapping(raw_params, "specific_volume_L_per_g")
		for nutrient, value in spec_cfg.items():
			self.specific_volume_L_per_g[nutrient] = max(value, 0.0)

		self._update_volume()

	def step_endogenous(self, t_min: float, dt_min: float, signals: Dict[str, Signal]):
		self._update_volume() 

	def _lumen_volume_L(self) -> float:
		total = self.residual_volume_L
		for nutrient, pool in self.pools.items():
			spec = self.specific_volume_L_per_g.get(nutrient, SPECIFIC_VOLUME_L_PER_G.get(nutrient, 1.0e-3))
			total += max(pool.amount_g, 0.0) * spec
		return max(self.min_volume_L, total)

	def _update_volume(self) -> None:
		self.volume_L = self._lumen_volume_L()

	def add(self, nutrient: str, grams: float):
		super().add(nutrient, grams)
		self._update_volume()

	def add_many(self, additions: Dict[str, float]):
		super().add_many(additions)
		self._update_volume()

	def fiber_viscosity_modifier(self) -> float:
		fiber = self.get("soluble_fiber")
		if fiber <= 0.0:
			return 1.0
		ratio = fiber / max(self.fiber_reference_g, 1e-6)
		modifier = 1.0 / (1.0 + self.fiber_alpha * ratio)
		return _clamp(modifier, self.fiber_min_modifier, 1.0)

	def glucose_concentration_mM(self) -> float:
		volume = self.volume_L if self.volume_L is not None else 0.0
		volume = max(volume, 1e-6)
		glucose_g = self.get("glucose")
		conc_mM = (glucose_g * 1000.0) / (MOLAR_MASS_GLUCOSE * volume) 
		return max(0.0, conc_mM)


class IntestinalTissue(Compartment):

	def __init__(self, params: Optional[Dict[str, float]] = None):
		params = params or {}
		nutrients = ["glucose"]
		volume = float(params.get("volume_L", 0.025))
		super().__init__("IntestinalTissue", nutrients, params, volume_L=max(volume, 1e-6))

	def step_endogenous(self, t_min: float, dt_min: float, signals: Dict[str, Signal]):
		return

	def glucose_concentration_mM(self) -> float:
		volume = self.volume_L if self.volume_L is not None else 0.0
		volume = max(volume, 1e-6)
		glucose_g = self.get("glucose")
		return (glucose_g * 1000.0) / (MOLAR_MASS_GLUCOSE * volume)

