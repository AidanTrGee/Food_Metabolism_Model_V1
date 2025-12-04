from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Any, Dict


@dataclass
class Meal:
    time_min: float
    grams_carbohydrate: float
    grams_fat: float
    grams_protein: float
    grams_soluble_fiber: float = 0.0
    grams_water: float = 0.0
    label: str | None = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]):
        return cls(
            time_min=data["time_min"],
            grams_carbohydrate=data.get("grams_carbohydrate", data.get("grams_glucose", 0.0)),
            grams_fat=data.get("grams_fat", 0.0),
            grams_protein=data.get("grams_protein", 0.0),
            grams_soluble_fiber=data.get("grams_soluble_fiber", 0.0),
            grams_water=data.get("grams_water", 0.0),
            label=data.get("label") if data.get("label") is None or isinstance(data.get("label"), str) else str(data.get("label")),
        )

    def as_stomach_additions(self) -> Dict[str, float]:
        return {
            "glucose": self.grams_carbohydrate,
            "lipid": self.grams_fat,
            "protein": self.grams_protein,
            "soluble_fiber": self.grams_soluble_fiber,
            "water": self.grams_water,
        }


def apply_meals_to_stomach(t_min, dt_min, meals, stomach):
    for meal in meals:
        if t_min <= meal.time_min < t_min + dt_min + 1e-9:
            stomach.ingest_meal(meal.time_min, meal.as_stomach_additions())
