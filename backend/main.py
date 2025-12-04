from __future__ import annotations

from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from Food_metabolism.web.simulation import exposed_parameters, run_simulation


class MealPayload(BaseModel):
    time_min: float = Field(..., ge=0)
    grams_carbohydrate: float = Field(..., ge=0)
    grams_fat: float = Field(0.0, ge=0)
    grams_protein: float = Field(0.0, ge=0)
    grams_soluble_fiber: float = Field(0.0, ge=0)
    grams_water: float = Field(0.0, ge=0)
    label: Optional[str] = None


class SimulationRequest(BaseModel):
    parameter_overrides: Optional[Dict[str, float]] = None
    meals: Optional[List[MealPayload]] = None


app = FastAPI(title="Food Metabolism Simulator API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health() -> Dict[str, str]:
    return {"status": "ok"}


@app.get("/config")
def get_config() -> Dict[str, object]:
    return exposed_parameters()


@app.post("/simulate")
def simulate(request: SimulationRequest) -> Dict[str, object]:
    try:
        meals = None
        if request.meals is not None:
            meals = [meal.model_dump() for meal in request.meals]
        result = run_simulation(request.parameter_overrides, meals)
        return result
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(status_code=400, detail=str(exc)) from exc
