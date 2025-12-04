from __future__ import annotations

from typing import Optional

Number = float


def grams_to_mg_dl(amount_g: Number, volume_L: Optional[Number]) -> Number:
    if volume_L is None:
        return 0.0
    volume = float(volume_L)
    if volume <= 0.0:
        return 0.0
    mg = float(amount_g) * 1000.0
    deciliters = volume * 10.0
    if deciliters <= 0.0:
        return 0.0
    return mg / deciliters
