from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "01_rate.py"

SPEC = importlib.util.spec_from_file_location(
    "rate_script_for_cache_tests",
    SCRIPT_PATH,
)
assert SPEC is not None
assert SPEC.loader is not None

RATE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RATE
SPEC.loader.exec_module(RATE)


def base_inputs() -> dict[str, object]:
    return {
        "cfg": {
            "family": "gg",
            "model_id": "test_model",
            "dim": 2,
            "box": [[-1.0, 1.0], [-1.0, 1.0]],
            "eval_box": [[-0.5, 0.5], [-0.5, 0.5]],
            "interval": [0.2, 1.0],
            "A": [[1.0, 0.0], [0.0, 1.0]],
        },
        "model_id": "test_model",
        "t0": 0.6,
        "xi0": np.array([0.0, 0.0]),
        "x_grid": np.array(
            [
                [-0.5, -0.5],
                [-0.5, 0.5],
                [0.5, -0.5],
                [0.5, 0.5],
            ],
            dtype=float,
        ),
        "truth_grid_2d": 121,
    }


def digest(**overrides: object) -> str:
    values = base_inputs()
    values.update(overrides)
    metadata = RATE.truth_cache_metadata(**values)
    return str(metadata["digest"])


def test_truth_cache_digest_tracks_truth_defining_inputs() -> None:
    baseline = base_inputs()
    baseline_digest = digest()

    assert digest(t0=0.61) != baseline_digest

    assert digest(
        xi0=np.array([0.1, 0.0]),
    ) != baseline_digest

    changed_grid = np.asarray(
        baseline["x_grid"],
        dtype=float,
    ).copy()
    changed_grid[0, 0] += 1e-8

    assert digest(x_grid=changed_grid) != baseline_digest
    assert digest(truth_grid_2d=241) != baseline_digest
    assert digest(model_id="another_model") != baseline_digest

    changed_cfg = copy.deepcopy(baseline["cfg"])
    assert isinstance(changed_cfg, dict)
    changed_cfg["A"][0][0] = 0.9

    assert digest(cfg=changed_cfg) != baseline_digest


class DummyModel:
    dim = 1


class DummyEngine:
    def __init__(self) -> None:
        self.model = DummyModel()
        self.calls = 0

    def a_star(
        self,
        t: float,
        x: np.ndarray,
        xi: np.ndarray,
        grid_points_2d: int,
    ) -> np.ndarray:
        self.calls += 1
        x_value = float(np.asarray(x).reshape(-1)[0])
        xi_value = float(np.asarray(xi).reshape(-1)[0])

        return np.array(
            [
                x_value
                + xi_value
                + float(t)
                + float(grid_points_2d)
            ]
        )


def one_dimensional_inputs() -> dict[str, object]:
    return {
        "cfg": {
            "family": "dummy",
            "model_id": "dummy_1d",
            "dim": 1,
            "box": [[-1.0, 1.0]],
            "eval_box": [[-0.5, 0.5]],
            "interval": [0.2, 1.0],
        },
        "model_id": "dummy_1d",
        "t0": 0.6,
        "xi0": np.array([0.0]),
        "x_grid": np.array(
            [[-0.5], [0.0], [0.5]],
            dtype=float,
        ),
        "truth_grid_2d": 121,
    }


def test_truth_cache_hit_and_input_change(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(RATE, "ROOT", tmp_path)

    engine = DummyEngine()
    inputs = one_dimensional_inputs()

    first = RATE.get_truth(
        engine=engine,
        cache_mode="use",
        **inputs,
    )
    calls_after_first = engine.calls

    second = RATE.get_truth(
        engine=engine,
        cache_mode="use",
        **inputs,
    )

    assert calls_after_first == len(inputs["x_grid"])
    assert engine.calls == calls_after_first
    np.testing.assert_array_equal(first, second)

    changed = dict(inputs)
    changed["truth_grid_2d"] = 241

    RATE.get_truth(
        engine=engine,
        cache_mode="use",
        **changed,
    )

    assert engine.calls == (
        calls_after_first + len(inputs["x_grid"])
    )

    cache_files = list(
        (
            tmp_path
            / "results"
            / "processed"
            / "rate"
            / "truth_cache"
        ).glob("v2_*.npz")
    )

    assert len(cache_files) == 2


def test_invalid_cache_shape_is_recomputed(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(RATE, "ROOT", tmp_path)

    engine = DummyEngine()
    inputs = one_dimensional_inputs()

    metadata = RATE.truth_cache_metadata(**inputs)
    cache_path = RATE.make_truth_cache_path(metadata)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        cache_path,
        truth=np.zeros((1, 1)),
        metadata_json=np.array(
            json.dumps(
                metadata,
                sort_keys=True,
                separators=(",", ":"),
            )
        ),
    )

    truth = RATE.get_truth(
        engine=engine,
        cache_mode="use",
        **inputs,
    )

    assert truth.shape == (3, 1)
    assert engine.calls == 3


def test_cache_mode_off_neither_reads_nor_writes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(RATE, "ROOT", tmp_path)

    engine = DummyEngine()
    inputs = one_dimensional_inputs()

    RATE.get_truth(
        engine=engine,
        cache_mode="off",
        **inputs,
    )
    RATE.get_truth(
        engine=engine,
        cache_mode="off",
        **inputs,
    )

    assert engine.calls == 2 * len(inputs["x_grid"])

    cache_root = (
        tmp_path
        / "results"
        / "processed"
        / "rate"
        / "truth_cache"
    )

    assert not cache_root.exists()
