from __future__ import annotations
from dataclasses import dataclass
from ..core.base import Controller, Signal

@dataclass
class ExerciseBlock:
    start_min: float; end_min: float; intensity: float  # 0..1

class ExerciseController(Controller):
    def __init__(self, blocks): super().__init__("ExerciseController"); self.blocks = blocks
    def update(self, t_min, dt_min, signals, state_lookup):
        level = 0.0
        for b in self.blocks:
            if b.start_min <= t_min < b.end_min: level = max(level, b.intensity)
        signals.setdefault("exercise", Signal("exercise", 0.0)).level = level
