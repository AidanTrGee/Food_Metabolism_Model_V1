from __future__ import annotations

from ..core.base import Compartment, Signal


class VenousBlood(Compartment):
    def __init__(self, params=None):
        params = params or {}
        nutrients = ["glucose", "lipid", "ffa", "protein", "lactate"]
        volume = float(params.get("volume_L", 3.5))
        super().__init__("VenousBlood", nutrients, params, volume_L=max(volume, 1e-6))

    def step_endogenous(self, t_min, dt_min, signals):
        return


class ArterialBlood(Compartment):
    def __init__(self, params=None):
        params = params or {}
        nutrients = ["glucose", "lipid", "ffa", "protein", "lactate"]
        volume = float(params.get("volume_L", 2.0))
        super().__init__("ArterialBlood", nutrients, params, volume_L=max(volume, 1e-6))

    def step_endogenous(self, t_min, dt_min, signals):
        return


class Lungs(Compartment):

    def __init__(self, params=None):
        params = params or {}
        nutrients = ["glucose", "lipid", "ffa", "protein", "lactate"]
        volume = float(params.get("volume_L", 1.3))
        super().__init__("Lungs", nutrients, params, volume_L=max(volume, 1e-6))

    def step_endogenous(self, t_min, dt_min, signals):
        return


class PortalVein(Compartment):
    def __init__(self, params=None):
        params = params or {}
        volume = float(params.get("volume_L", 0.4))
        super().__init__("PortalVein", ["glucose", "protein"], params, volume_L=max(volume, 1e-6))

    def step_endogenous(self, t_min, dt_min, signals):
        return

Blood = VenousBlood
