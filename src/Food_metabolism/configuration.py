from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
_CONFIG_DIR = _REPO_ROOT / "config"
_PARAMS_BASENAME = "params.yml"
_REACTION_RATES_BASENAME = "reaction_rates.yml"


def _load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Configuration file {path} must contain a mapping at the top level.")
    return data


def _deep_merge_dicts(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {k: v for k, v in base.items()}
    for key, override_value in override.items():
        base_value = result.get(key)
        if isinstance(base_value, Mapping) and isinstance(override_value, Mapping):
            result[key] = _deep_merge_dicts(base_value, override_value)
        else:
            result[key] = override_value
    return result


def load_simulation_config(
    *,
    params_path: str | Path | None = None,
    reaction_rates_path: str | Path | None = None,
) -> dict[str, Any]:
    params_file = Path(params_path) if params_path else _CONFIG_DIR / _PARAMS_BASENAME
    rates_file = Path(reaction_rates_path) if reaction_rates_path else _CONFIG_DIR / _REACTION_RATES_BASENAME

    params_cfg = _load_yaml(params_file)
    reaction_cfg = _load_yaml(rates_file)
    if not reaction_cfg:
        return params_cfg
    if not params_cfg:
        return reaction_cfg
    return _deep_merge_dicts(reaction_cfg, params_cfg)
