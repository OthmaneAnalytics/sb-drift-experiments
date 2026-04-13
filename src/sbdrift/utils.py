from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import yaml


def as_array(x: Sequence[float] | float | np.ndarray) -> np.ndarray:
    """Convert input to a float numpy array."""
    return np.asarray(x, dtype=float)


def ensure_dir(path: str | Path) -> Path:
    """Create a directory if it does not exist and return it as a Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML file into a dictionary."""
    p = Path(path)
    with p.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data if data is not None else {}


def save_json(data: Any, path: str | Path) -> None:
    """Save data as pretty JSON."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True, default=_json_default)


def _json_default(obj: Any) -> Any:
    """JSON serializer for numpy / pathlib objects."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")
