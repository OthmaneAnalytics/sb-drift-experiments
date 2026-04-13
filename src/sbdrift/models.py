from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Sequence

import numpy as np
from scipy.special import expit
from scipy.stats import multivariate_normal, norm

from .utils import as_array


Array = np.ndarray


def _box_arrays(box: Sequence[Sequence[float]], dim: int) -> tuple[Array, Array]:
    box_arr = np.asarray(box, dtype=float)
    if box_arr.shape != (dim, 2):
        raise ValueError(f"Box must have shape ({dim}, 2), got {box_arr.shape}")
    low = box_arr[:, 0]
    high = box_arr[:, 1]
    if np.any(low >= high):
        raise ValueError("Each box interval must satisfy low < high.")
    return low, high


def _rectangle_prob(mean: Array, cov: Array, low: Array, high: Array) -> float:
    dim = len(mean)
    if dim == 1:
        sigma = float(np.sqrt(cov[0, 0]))
        return float(norm.cdf(high[0], loc=mean[0], scale=sigma) - norm.cdf(low[0], loc=mean[0], scale=sigma))

    mvn = multivariate_normal(mean=mean, cov=cov, allow_singular=False)
    total = 0.0
    for bits in product([0, 1], repeat=dim):
        corner = np.where(np.array(bits, dtype=int) == 1, high, low)
        sign = (-1) ** (dim - int(np.sum(bits)))
        total += sign * float(mvn.cdf(corner))
    return max(total, 1e-15)


@dataclass(frozen=True)
class TruncatedGaussian:
    mean: Array
    cov: Array
    low: Array
    high: Array

    def __post_init__(self) -> None:
        object.__setattr__(self, "mean", as_array(self.mean).reshape(-1))
        cov = np.asarray(self.cov, dtype=float)
        if cov.ndim == 0:
            cov = np.array([[float(cov)]], dtype=float)
        if cov.ndim == 1:
            cov = np.diag(cov)
        object.__setattr__(self, "cov", cov)
        dim = len(self.mean)
        if cov.shape != (dim, dim):
            raise ValueError("Covariance shape does not match mean dimension.")
        low = as_array(self.low).reshape(-1)
        high = as_array(self.high).reshape(-1)
        if len(low) != dim or len(high) != dim:
            raise ValueError("Box dimension does not match mean dimension.")
        object.__setattr__(self, "low", low)
        object.__setattr__(self, "high", high)
        norm_const = _rectangle_prob(self.mean, self.cov, self.low, self.high)
        object.__setattr__(self, "_norm_const", norm_const)
        object.__setattr__(self, "_mvn", multivariate_normal(mean=self.mean, cov=self.cov, allow_singular=False))

    @property
    def dim(self) -> int:
        return len(self.mean)

    def in_box(self, x: Array) -> bool:
        x = as_array(x).reshape(-1)
        return bool(np.all(x >= self.low) and np.all(x <= self.high))

    def pdf(self, x: Sequence[float] | Array) -> float:
        x = as_array(x).reshape(-1)
        if not self.in_box(x):
            return 0.0
        return float(self._mvn.pdf(x) / self._norm_const)

    def sample(self, n: int, rng: np.random.Generator) -> Array:
        samples: list[Array] = []
        needed = n
        batch = max(512, 4 * n)
        while needed > 0:
            draw = rng.multivariate_normal(self.mean, self.cov, size=batch)
            keep = np.all((draw >= self.low) & (draw <= self.high), axis=1)
            kept = draw[keep]
            if kept.size:
                take = kept[:needed]
                samples.append(take)
                needed -= len(take)
        return np.vstack(samples)


class BaseModel:
    name: str
    dim: int
    low: Array
    high: Array
    s: float
    u: float

    def xs_density(self, xi: Sequence[float] | Array) -> float:
        raise NotImplementedError

    def xu_conditional_density(self, y: Sequence[float] | Array, xi: Sequence[float] | Array) -> float:
        raise NotImplementedError

    def conditional_pdf_fn(self, xi: Sequence[float] | Array):
        raise NotImplementedError

    def sample(self, n: int, rng: np.random.Generator) -> tuple[Array, Array]:
        raise NotImplementedError


@dataclass
class GGModel(BaseModel):
    name: str
    dim: int
    s: float
    u: float
    low: Array
    high: Array
    xs_dist: TruncatedGaussian
    A: Array
    b: Array
    eps_cov: Array

    @classmethod
    def from_config(cls, cfg: dict) -> "GGModel":
        dim = int(cfg["dim"])
        low, high = _box_arrays(cfg["box"], dim)
        xs_dist = TruncatedGaussian(cfg["m_s"], cfg["Sigma_s"], low, high)
        A = np.asarray(cfg["A"], dtype=float)
        if dim == 1:
            A = np.array([[float(A)]]) if np.ndim(A) == 0 else np.asarray(A, dtype=float).reshape(1, 1)
        b = as_array(cfg["b"]).reshape(dim)
        eps_cov = np.asarray(cfg["Sigma_eps"], dtype=float)
        if eps_cov.ndim == 0:
            eps_cov = np.array([[float(eps_cov)]])
        if eps_cov.ndim == 1:
            eps_cov = np.diag(eps_cov)
        return cls(
            name=cfg.get("name", f"GG{dim}D"),
            dim=dim,
            s=float(cfg["interval"][0]),
            u=float(cfg["interval"][1]),
            low=low,
            high=high,
            xs_dist=xs_dist,
            A=A,
            b=b,
            eps_cov=eps_cov,
        )

    def xs_density(self, xi: Sequence[float] | Array) -> float:
        return self.xs_dist.pdf(xi)

    def _conditional_dist(self, xi: Sequence[float] | Array) -> TruncatedGaussian:
        xi = as_array(xi).reshape(self.dim)
        mean = self.A @ xi + self.b
        return TruncatedGaussian(mean, self.eps_cov, self.low, self.high)

    def xu_conditional_density(self, y: Sequence[float] | Array, xi: Sequence[float] | Array) -> float:
        return self._conditional_dist(xi).pdf(y)

    def conditional_pdf_fn(self, xi: Sequence[float] | Array):
        dist = self._conditional_dist(xi)
        return dist.pdf

    def sample(self, n: int, rng: np.random.Generator) -> tuple[Array, Array]:
        xs = self.xs_dist.sample(n, rng)
        ys = np.zeros_like(xs)
        for i, xi in enumerate(xs):
            ys[i] = self._conditional_dist(xi).sample(1, rng)[0]
        return xs, ys


@dataclass
class MMModel(BaseModel):
    name: str
    dim: int
    s: float
    u: float
    low: Array
    high: Array
    xs_components: tuple[TruncatedGaussian, TruncatedGaussian]
    A1: Array
    b1: Array
    Sigma1: Array
    A2: Array
    b2: Array
    Sigma2: Array
    alpha0: float
    alpha: Array

    @classmethod
    def from_config(cls, cfg: dict) -> "MMModel":
        dim = int(cfg["dim"])
        low, high = _box_arrays(cfg["box"], dim)
        xs_components = (
            TruncatedGaussian(cfg["m1"], cfg["S"], low, high),
            TruncatedGaussian(cfg["m2"], cfg["S"], low, high),
        )

        def _mat(x: Sequence[Sequence[float]] | float) -> Array:
            arr = np.asarray(x, dtype=float)
            if arr.ndim == 0:
                return np.array([[float(arr)]])
            if arr.ndim == 1 and dim == 1:
                return arr.reshape(1, 1)
            return arr

        def _cov(x: Sequence[Sequence[float]] | Sequence[float] | float) -> Array:
            arr = np.asarray(x, dtype=float)
            if arr.ndim == 0:
                return np.array([[float(arr)]])
            if arr.ndim == 1:
                return np.diag(arr)
            return arr

        gate = cfg["gate"]
        return cls(
            name=cfg.get("name", f"MM{dim}D"),
            dim=dim,
            s=float(cfg["interval"][0]),
            u=float(cfg["interval"][1]),
            low=low,
            high=high,
            xs_components=xs_components,
            A1=_mat(cfg["A1"]),
            b1=as_array(cfg["b1"]).reshape(dim),
            Sigma1=_cov(cfg["Sigma1"]),
            A2=_mat(cfg["A2"]),
            b2=as_array(cfg["b2"]).reshape(dim),
            Sigma2=_cov(cfg["Sigma2"]),
            alpha0=float(gate.get("alpha0", 0.0)),
            alpha=as_array(gate["alpha"]).reshape(dim),
        )

    def gate(self, xi: Sequence[float] | Array) -> float:
        xi = as_array(xi).reshape(self.dim)
        return float(expit(self.alpha0 + float(self.alpha @ xi)))

    def xs_density(self, xi: Sequence[float] | Array) -> float:
        xi = as_array(xi)
        return 0.5 * self.xs_components[0].pdf(xi) + 0.5 * self.xs_components[1].pdf(xi)

    def _conditional_component_1(self, xi: Sequence[float] | Array) -> TruncatedGaussian:
        xi = as_array(xi).reshape(self.dim)
        return TruncatedGaussian(self.A1 @ xi + self.b1, self.Sigma1, self.low, self.high)

    def _conditional_component_2(self, xi: Sequence[float] | Array) -> TruncatedGaussian:
        xi = as_array(xi).reshape(self.dim)
        return TruncatedGaussian(self.A2 @ xi + self.b2, self.Sigma2, self.low, self.high)

    def xu_conditional_density(self, y: Sequence[float] | Array, xi: Sequence[float] | Array) -> float:
        pi = self.gate(xi)
        return pi * self._conditional_component_1(xi).pdf(y) + (1.0 - pi) * self._conditional_component_2(xi).pdf(y)

    def conditional_pdf_fn(self, xi: Sequence[float] | Array):
        pi = self.gate(xi)
        d1 = self._conditional_component_1(xi)
        d2 = self._conditional_component_2(xi)
        return lambda y: pi * d1.pdf(y) + (1.0 - pi) * d2.pdf(y)

    def sample(self, n: int, rng: np.random.Generator) -> tuple[Array, Array]:
        comp = rng.binomial(1, 0.5, size=n)
        xs = np.zeros((n, self.dim), dtype=float)
        for k in [0, 1]:
            idx = np.where(comp == k)[0]
            if len(idx):
                xs[idx] = self.xs_components[k].sample(len(idx), rng)
        ys = np.zeros_like(xs)
        for i, xi in enumerate(xs):
            if rng.uniform() < self.gate(xi):
                ys[i] = self._conditional_component_1(xi).sample(1, rng)[0]
            else:
                ys[i] = self._conditional_component_2(xi).sample(1, rng)[0]
        return xs, ys


def load_model_from_config(cfg: dict) -> BaseModel:
    family = cfg["family"].lower()
    if family == "gg":
        return GGModel.from_config(cfg)
    if family == "mm":
        return MMModel.from_config(cfg)
    raise ValueError(f"Unknown family: {cfg['family']}")
