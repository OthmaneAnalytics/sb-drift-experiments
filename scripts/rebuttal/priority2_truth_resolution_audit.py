#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sbdrift.models import GGModel, MMModel, load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.utils import load_yaml


RESOLUTIONS_2D = [81, 121, 161, 201, 241]
REFERENCE_RESOLUTION = 241
CHUNK_SIZE = 16

RUN_ID = datetime.now().strftime(
    "priority2_truth_%Y%m%d_%H%M%S"
)
OUT = ROOT / "results" / "rebuttal" / "priority2" / RUN_ID
OUT.mkdir(parents=True, exist_ok=False)


def evaluation_grid(
    cfg: dict,
    points_per_dimension: int,
) -> tuple[list[np.ndarray], np.ndarray]:
    box = np.asarray(cfg["eval_box"], dtype=float)

    axes = [
        np.linspace(
            box[k, 0],
            box[k, 1],
            points_per_dimension,
        )
        for k in range(len(box))
    ]

    mesh = np.meshgrid(*axes, indexing="ij")
    points = np.stack(
        [component.reshape(-1) for component in mesh],
        axis=1,
    )

    return axes, points


def trapezoid_weights(
    low: float,
    high: float,
    n: int,
) -> tuple[np.ndarray, np.ndarray]:
    axis = np.linspace(low, high, n)
    spacing = float(axis[1] - axis[0])

    weights = np.full(n, spacing, dtype=float)
    weights[0] *= 0.5
    weights[-1] *= 0.5

    return axis, weights


def vectorized_conditional_pdf(
    model,
    xi: np.ndarray,
    y_points: np.ndarray,
) -> np.ndarray:
    if isinstance(model, GGModel):
        distribution = model._conditional_dist(xi)
        values = (
            distribution._mvn.pdf(y_points)
            / distribution._norm_const
        )

    elif isinstance(model, MMModel):
        mixing_weight = model.gate(xi)
        distribution_1 = model._conditional_component_1(xi)
        distribution_2 = model._conditional_component_2(xi)

        values = (
            mixing_weight
            * distribution_1._mvn.pdf(y_points)
            / distribution_1._norm_const
            + (1.0 - mixing_weight)
            * distribution_2._mvn.pdf(y_points)
            / distribution_2._norm_const
        )

    else:
        raise TypeError(
            f"Unsupported model type: {type(model).__name__}"
        )

    values = np.asarray(values, dtype=float).reshape(-1)

    if not np.all(np.isfinite(values)):
        raise RuntimeError(
            "Conditional density contains nonfinite values."
        )

    return values


def truth_field_2d_vectorized(
    model,
    t: float,
    xi: np.ndarray,
    x_grid: np.ndarray,
    resolution: int,
    chunk_size: int = CHUNK_SIZE,
) -> np.ndarray:
    if model.dim != 2:
        raise ValueError("This function requires a 2D model.")

    axis_0, weights_0 = trapezoid_weights(
        float(model.low[0]),
        float(model.high[0]),
        resolution,
    )
    axis_1, weights_1 = trapezoid_weights(
        float(model.low[1]),
        float(model.high[1]),
        resolution,
    )

    mesh_0, mesh_1 = np.meshgrid(
        axis_0,
        axis_1,
        indexing="ij",
    )
    y_points = np.stack(
        [mesh_0.reshape(-1), mesh_1.reshape(-1)],
        axis=1,
    )

    tensor_weights = (
        weights_0[:, None] * weights_1[None, :]
    ).reshape(-1)

    density = vectorized_conditional_pdf(
        model,
        xi,
        y_points,
    )

    delta_t = float(model.u - t)
    total_interval = float(model.u - model.s)

    xi_array = np.asarray(xi, dtype=float).reshape(2)

    initial_tilt = np.sum(
        (y_points - xi_array[None, :]) ** 2,
        axis=1,
    ) / (2.0 * total_interval)

    base_weights = tensor_weights * density

    output = np.empty_like(x_grid, dtype=float)

    for start in range(0, len(x_grid), chunk_size):
        stop = min(start + chunk_size, len(x_grid))
        x_chunk = np.asarray(
            x_grid[start:stop],
            dtype=float,
        )

        terminal_distance = np.sum(
            (
                y_points[None, :, :]
                - x_chunk[:, None, :]
            )
            ** 2,
            axis=2,
        )

        exponent = (
            -terminal_distance / (2.0 * delta_t)
            + initial_tilt[None, :]
        )

        tilted_weights = (
            np.exp(exponent)
            * base_weights[None, :]
        )

        denominator = np.sum(
            tilted_weights,
            axis=1,
        )
        numerator = tilted_weights @ y_points

        if np.any(denominator <= 0.0):
            raise RuntimeError(
                "Nonpositive population denominator."
            )

        output[start:stop] = (
            numerator / denominator[:, None]
            - x_chunk
        ) / delta_t

    if not np.all(np.isfinite(output)):
        raise RuntimeError(
            "Computed 2D truth contains nonfinite values."
        )

    return output


def truth_field_1d_reference(
    model,
    t: float,
    xi: np.ndarray,
    x_grid: np.ndarray,
) -> np.ndarray:
    low = float(model.low[0])
    high = float(model.high[0])
    xi_array = np.asarray(xi, dtype=float).reshape(1)

    conditional_pdf = model.conditional_pdf_fn(xi_array)

    delta_t = float(model.u - t)
    total_interval = float(model.u - model.s)

    output = np.empty((len(x_grid), 1), dtype=float)

    for index, x in enumerate(x_grid):
        x_value = float(np.asarray(x).reshape(-1)[0])

        def tilted_density(y: float) -> float:
            y_array = np.array([y], dtype=float)

            exponent = (
                -((y - x_value) ** 2)
                / (2.0 * delta_t)
                + ((y - xi_array[0]) ** 2)
                / (2.0 * total_interval)
            )

            return (
                float(np.exp(exponent))
                * float(conditional_pdf(y_array))
            )

        denominator, _ = quad(
            tilted_density,
            low,
            high,
            epsabs=1e-12,
            epsrel=1e-12,
            limit=500,
        )

        numerator, _ = quad(
            lambda y: y * tilted_density(y),
            low,
            high,
            epsabs=1e-12,
            epsrel=1e-12,
            limit=500,
        )

        output[index, 0] = (
            numerator / denominator - x_value
        ) / delta_t

    return output


def vector_field_ise(
    estimate: np.ndarray,
    reference: np.ndarray,
    axes: list[np.ndarray],
) -> float:
    squared_norm = np.sum(
        (estimate - reference) ** 2,
        axis=1,
    )

    if len(axes) == 1:
        return float(
            np.trapezoid(
                squared_norm,
                axes[0],
            )
        )

    shaped = squared_norm.reshape(
        len(axes[0]),
        len(axes[1]),
    )

    integrated_axis_1 = np.trapezoid(
        shaped,
        axes[1],
        axis=1,
    )

    return float(
        np.trapezoid(
            integrated_axis_1,
            axes[0],
            axis=0,
        )
    )


def comparison_metrics(
    estimate: np.ndarray,
    reference: np.ndarray,
    axes: list[np.ndarray],
) -> dict[str, float]:
    difference = estimate - reference
    vector_errors = np.linalg.norm(
        difference,
        axis=1,
    )

    reference_norms = np.linalg.norm(
        reference,
        axis=1,
    )

    sup_vector = float(np.max(vector_errors))
    sup_component = float(np.max(np.abs(difference)))
    rms_vector = float(
        np.sqrt(np.mean(vector_errors**2))
    )
    reference_scale = max(
        float(np.max(reference_norms)),
        1e-15,
    )

    return {
        "sup_vector_error": sup_vector,
        "sup_component_error": sup_component,
        "relative_sup_vector_error":
            sup_vector / reference_scale,
        "rms_vector_error": rms_vector,
        "ise": vector_field_ise(
            estimate,
            reference,
            axes,
        ),
    }


def legacy_cache_path(
    model_id: str,
    t: float,
    xi: np.ndarray,
    n_grid: int,
) -> Path:
    xi_tag = "_".join(
        f"{value:+.3f}"
        for value in np.asarray(xi).reshape(-1)
    )

    return (
        ROOT
        / "results"
        / "processed"
        / "rate"
        / "truth_cache"
        / (
            f"{model_id}_t{t:.3f}_"
            f"xi{xi_tag}_n{n_grid}.npz"
        )
    )


def probe_indices(
    points_per_dimension: int,
) -> list[int]:
    n = points_per_dimension
    coordinates = [
        (0, 0),
        (0, n // 2),
        (0, n - 1),
        (n // 2, 0),
        (n // 2, n // 2),
        (n // 2, n - 1),
        (n - 1, 0),
        (n - 1, n // 2),
        (n - 1, n - 1),
    ]

    return [i * n + j for i, j in coordinates]


summary_rows: list[dict[str, object]] = []
runtime_rows: list[dict[str, object]] = []
field_archive: dict[str, np.ndarray] = {}
implementation_errors: list[str] = []


# ------------------------------------------------------------------
# Full 1D paper grids
# ------------------------------------------------------------------
for config_name in ["gg_1d.yaml", "mm_1d.yaml"]:
    config_path = ROOT / "configs" / config_name
    cfg = load_yaml(config_path)
    model = load_model_from_config(cfg)
    engine = TruthEngine(model)

    model_id = str(cfg["model_id"])
    t = float(cfg["rate_time"])
    xi = np.asarray(cfg["xi0"], dtype=float)
    axes, x_grid = evaluation_grid(cfg, 200)

    start = time.perf_counter()
    engine_field = np.asarray(
        [
            engine.a_star(
                t,
                x=x,
                xi=xi,
            )
            for x in x_grid
        ],
        dtype=float,
    )
    engine_seconds = time.perf_counter() - start

    start = time.perf_counter()
    reference_field = truth_field_1d_reference(
        model,
        t,
        xi,
        x_grid,
    )
    reference_seconds = time.perf_counter() - start

    metrics = comparison_metrics(
        engine_field,
        reference_field,
        axes,
    )

    summary_rows.append(
        {
            "model_id": model_id,
            "dimension": 1,
            "comparison": "TruthEngine_vs_tighter_quad",
            "estimate_resolution": "adaptive_quad_1e-10",
            "reference_resolution": "adaptive_quad_1e-12",
            **metrics,
        }
    )

    runtime_rows.extend(
        [
            {
                "model_id": model_id,
                "calculation": "TruthEngine",
                "resolution": "adaptive_quad_1e-10",
                "seconds": engine_seconds,
            },
            {
                "model_id": model_id,
                "calculation": "independent_reference",
                "resolution": "adaptive_quad_1e-12",
                "seconds": reference_seconds,
            },
        ]
    )

    field_archive[
        f"{model_id}_engine"
    ] = engine_field
    field_archive[
        f"{model_id}_reference"
    ] = reference_field

    legacy_path = legacy_cache_path(
        model_id,
        t,
        xi,
        200,
    )

    if legacy_path.exists():
        with np.load(
            legacy_path,
            allow_pickle=False,
        ) as data:
            legacy = np.asarray(
                data["truth"],
                dtype=float,
            )

        if legacy.shape != reference_field.shape:
            implementation_errors.append(
                f"{model_id}: legacy cache shape "
                f"{legacy.shape}, expected "
                f"{reference_field.shape}"
            )
        else:
            legacy_metrics = comparison_metrics(
                legacy,
                reference_field,
                axes,
            )

            summary_rows.append(
                {
                    "model_id": model_id,
                    "dimension": 1,
                    "comparison":
                        "legacy_cache_vs_reference",
                    "estimate_resolution":
                        "legacy_tracked_cache",
                    "reference_resolution":
                        "adaptive_quad_1e-12",
                    **legacy_metrics,
                }
            )


# ------------------------------------------------------------------
# Complete 21 x 21 two-dimensional paper grids
# ------------------------------------------------------------------
for config_name in ["gg_2d.yaml", "mm_2d.yaml"]:
    config_path = ROOT / "configs" / config_name
    cfg = load_yaml(config_path)
    model = load_model_from_config(cfg)
    engine = TruthEngine(model)

    model_id = str(cfg["model_id"])
    t = float(cfg["rate_time"])
    xi = np.asarray(cfg["xi0"], dtype=float)
    axes, x_grid = evaluation_grid(cfg, 21)

    fields_by_resolution: dict[int, np.ndarray] = {}

    for resolution in RESOLUTIONS_2D:
        print(
            f"[{model_id}] computing full "
            f"21x21 field at resolution "
            f"{resolution}",
            flush=True,
        )

        start = time.perf_counter()

        field = truth_field_2d_vectorized(
            model,
            t,
            xi,
            x_grid,
            resolution,
        )

        elapsed = time.perf_counter() - start

        fields_by_resolution[resolution] = field
        field_archive[
            f"{model_id}_grid{resolution}"
        ] = field

        runtime_rows.append(
            {
                "model_id": model_id,
                "calculation":
                    "vectorized_tensor_trapezoid",
                "resolution": resolution,
                "seconds": elapsed,
            }
        )

        print(
            f"[{model_id}] resolution "
            f"{resolution} completed in "
            f"{elapsed:.2f} seconds",
            flush=True,
        )

    reference = fields_by_resolution[
        REFERENCE_RESOLUTION
    ]

    for resolution in RESOLUTIONS_2D:
        if resolution == REFERENCE_RESOLUTION:
            continue

        metrics = comparison_metrics(
            fields_by_resolution[resolution],
            reference,
            axes,
        )

        summary_rows.append(
            {
                "model_id": model_id,
                "dimension": 2,
                "comparison":
                    f"grid{resolution}_vs_grid"
                    f"{REFERENCE_RESOLUTION}",
                "estimate_resolution": resolution,
                "reference_resolution":
                    REFERENCE_RESOLUTION,
                **metrics,
            }
        )

    # Confirm the vectorized implementation reproduces TruthEngine.
    indices = probe_indices(21)
    engine_probe = np.asarray(
        [
            engine.a_star(
                t,
                x=x_grid[index],
                xi=xi,
                grid_points_2d=121,
            )
            for index in indices
        ],
        dtype=float,
    )
    vectorized_probe = fields_by_resolution[121][
        indices
    ]

    probe_difference = np.max(
        np.abs(
            engine_probe - vectorized_probe
        )
    )

    summary_rows.append(
        {
            "model_id": model_id,
            "dimension": 2,
            "comparison":
                "vectorized121_vs_TruthEngine121_probes",
            "estimate_resolution": 121,
            "reference_resolution": 121,
            "sup_vector_error": float(
                np.max(
                    np.linalg.norm(
                        vectorized_probe - engine_probe,
                        axis=1,
                    )
                )
            ),
            "sup_component_error": float(
                probe_difference
            ),
            "relative_sup_vector_error": float(
                np.max(
                    np.linalg.norm(
                        vectorized_probe - engine_probe,
                        axis=1,
                    )
                )
                / max(
                    float(
                        np.max(
                            np.linalg.norm(
                                engine_probe,
                                axis=1,
                            )
                        )
                    ),
                    1e-15,
                )
            ),
            "rms_vector_error": float(
                np.sqrt(
                    np.mean(
                        np.linalg.norm(
                            vectorized_probe
                            - engine_probe,
                            axis=1,
                        )
                        ** 2
                    )
                )
            ),
            "ise": np.nan,
        }
    )

    if probe_difference > 1e-8:
        implementation_errors.append(
            f"{model_id}: vectorized audit and "
            f"TruthEngine differ by {probe_difference:.3e}"
        )

    legacy_path = legacy_cache_path(
        model_id,
        t,
        xi,
        21,
    )

    if legacy_path.exists():
        with np.load(
            legacy_path,
            allow_pickle=False,
        ) as data:
            legacy = np.asarray(
                data["truth"],
                dtype=float,
            )

        if legacy.shape != reference.shape:
            implementation_errors.append(
                f"{model_id}: legacy cache shape "
                f"{legacy.shape}, expected "
                f"{reference.shape}"
            )
        else:
            for reference_resolution in [121, 241]:
                metrics = comparison_metrics(
                    legacy,
                    fields_by_resolution[
                        reference_resolution
                    ],
                    axes,
                )

                summary_rows.append(
                    {
                        "model_id": model_id,
                        "dimension": 2,
                        "comparison":
                            "legacy_cache_vs_grid"
                            f"{reference_resolution}",
                        "estimate_resolution":
                            "legacy_tracked_cache",
                        "reference_resolution":
                            reference_resolution,
                        **metrics,
                    }
                )


summary = pd.DataFrame(summary_rows)
runtimes = pd.DataFrame(runtime_rows)

summary.to_csv(
    OUT / "truth_resolution_summary.csv",
    index=False,
)
runtimes.to_csv(
    OUT / "truth_resolution_runtimes.csv",
    index=False,
)

np.savez_compressed(
    OUT / "truth_fields.npz",
    **field_archive,
)

metadata = {
    "run_id": RUN_ID,
    "resolutions_2d": RESOLUTIONS_2D,
    "reference_resolution": REFERENCE_RESOLUTION,
    "chunk_size": CHUNK_SIZE,
    "paper_grid_1d": 200,
    "paper_grid_2d": 21,
}

(
    OUT / "audit_metadata.json"
).write_text(
    json.dumps(
        metadata,
        indent=2,
        sort_keys=True,
    ),
    encoding="utf-8",
)


print()
print("=" * 120)
print("TRUTH RESOLUTION SUMMARY")
print("=" * 120)

with pd.option_context(
    "display.max_rows",
    100,
    "display.max_columns",
    30,
    "display.width",
    300,
):
    print(summary.to_string(index=False))


print()
print("=" * 120)
print("RUNTIMES")
print("=" * 120)
print(runtimes.to_string(index=False))


print()
print("=" * 120)
print("COMPUTATION STATUS")
print("=" * 120)

if implementation_errors:
    print("FAIL")

    for error in implementation_errors:
        print(f"- {error}")

    raise SystemExit(1)

print("PASS")
print(f"Outputs: {OUT}")
