from pathlib import Path
import runpy

import numpy as np
import pytest
from scipy.integrate import quad

from sbdrift.estimator import DriftEstimator
from sbdrift.kernels import ProductKernel
from sbdrift.models import load_model_from_config
from sbdrift.truth_engine import TruthEngine
from sbdrift.utils import load_yaml


ROOT = Path(__file__).resolve().parents[1]
RATE_SCRIPT = runpy.run_path(str(ROOT / "scripts" / "01_rate.py"))
make_bandwidth_grid = RATE_SCRIPT["make_bandwidth_grid"]


def load_model(config_name: str):
    cfg = load_yaml(ROOT / "configs" / config_name)
    return load_model_from_config(cfg)


def test_epanechnikov_kernel_normalizes_in_1d():
    kernel = ProductKernel()
    grid = np.linspace(-1.0, 1.0, 20_001)
    integral = np.trapezoid(kernel.base(grid), grid)

    assert integral == pytest.approx(1.0, abs=1e-8)
    assert kernel.base(np.array([-1.1, 1.1])).tolist() == [0.0, 0.0]


@pytest.mark.parametrize(
    ("dim", "factor"),
    [
        (1, 81.0),
        (2, 9.0),
    ],
)
def test_bandwidth_grid_respects_effective_sample_floor(dim, factor):
    M = 1000
    hs = make_bandwidth_grid(
        M=M,
        dim=dim,
        h0=1.2,
        q=2 ** (-0.5),
        min_h_factor=factor,
    )

    expected_min = factor * M ** (-1.0 / dim)

    assert np.all(np.isfinite(hs))
    assert np.all(hs > 0)
    assert np.all(np.diff(hs) > 0)
    assert hs[0] == pytest.approx(expected_min, abs=1e-12)
    assert np.all(M * hs**dim >= 81.0 - 1e-9)


@pytest.mark.parametrize(
    "config_name",
    [
        "gg_1d.yaml",
        "gg_2d.yaml",
        "mm_1d.yaml",
        "mm_2d.yaml",
    ],
)
def test_sampling_is_reproducible_and_inside_support(config_name):
    model = load_model(config_name)

    xs1, xu1 = model.sample(24, np.random.default_rng(991))
    xs2, xu2 = model.sample(24, np.random.default_rng(991))

    assert xs1.shape == (24, model.dim)
    assert xu1.shape == (24, model.dim)
    np.testing.assert_array_equal(xs1, xs2)
    np.testing.assert_array_equal(xu1, xu2)

    assert np.all(xs1 >= model.low)
    assert np.all(xs1 <= model.high)
    assert np.all(xu1 >= model.low)
    assert np.all(xu1 <= model.high)


@pytest.mark.parametrize(
    ("config_name", "xi"),
    [
        ("gg_1d.yaml", np.array([0.0])),
        ("mm_1d.yaml", np.array([0.8])),
    ],
)
def test_conditional_density_integrates_to_one(config_name, xi):
    model = load_model(config_name)
    pdf = model.conditional_pdf_fn(xi)

    value, _ = quad(
        lambda y: pdf(np.array([y])),
        float(model.low[0]),
        float(model.high[0]),
        epsabs=1e-9,
        epsrel=1e-9,
    )

    assert value == pytest.approx(1.0, abs=1e-7)


def test_estimator_matches_direct_weighted_ratio():
    model = load_model("gg_1d.yaml")
    xs, xu = model.sample(64, np.random.default_rng(1234))
    estimator = DriftEstimator(model=model, xs=xs, xu=xu)

    t = 0.6
    x = np.array([0.2])
    xi = np.array([0.0])
    h = 0.7

    weights = estimator.kernel_weights(xi, h)

    dt = model.u - t
    delta = model.u - model.s

    F = np.exp(
        -((xu[:, 0] - x[0]) ** 2) / (2.0 * dt)
        + ((xu[:, 0] - xi[0]) ** 2) / (2.0 * delta)
    )

    g1 = np.mean(F * weights)
    g2 = np.mean(xu[:, 0] * F * weights)
    expected = (g2 / g1 - x[0]) / dt

    observed = estimator.a_hat_point(
        t=t,
        x=x,
        xi=xi,
        h=h,
    )[0]

    assert g1 > 0
    assert observed == pytest.approx(expected, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize(
    ("config_name", "x", "xi"),
    [
        ("gg_1d.yaml", np.array([0.2]), np.array([0.0])),
        ("mm_1d.yaml", np.array([0.3]), np.array([0.8])),
    ],
)
def test_population_truth_is_finite(config_name, x, xi):
    model = load_model(config_name)
    engine = TruthEngine(model)

    denominator = engine.D_star(t=0.6, x=x, xi=xi)
    drift = engine.a_star(t=0.6, x=x, xi=xi)

    assert np.isfinite(denominator)
    assert denominator > 0
    assert drift.shape == (model.dim,)
    assert np.all(np.isfinite(drift))
