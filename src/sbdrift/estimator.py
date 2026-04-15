from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .kernels import ProductKernel, make_kernel
from .models import BaseModel
from .utils import as_array

Array = np.ndarray


@dataclass
class DriftEstimator:
    model: BaseModel
    xs: Array
    xu: Array
    kernel: ProductKernel | None = None

    def __post_init__(self) -> None:
        self.xs = np.asarray(self.xs, dtype=float)
        self.xu = np.asarray(self.xu, dtype=float)
        if self.xs.ndim == 1:
            self.xs = self.xs.reshape(-1, 1)
        if self.xu.ndim == 1:
            self.xu = self.xu.reshape(-1, 1)
        if self.xs.shape != self.xu.shape:
            raise ValueError("xs and xu must have the same shape")
        if self.xs.shape[1] != self.model.dim:
            raise ValueError("sample dimension does not match model dimension")
        if self.kernel is None:
            self.kernel = make_kernel("epanechnikov")

    @property
    def n(self) -> int:
        return self.xs.shape[0]

    @property
    def dim(self) -> int:
        return self.model.dim

    def delta_t(self, t: float) -> float:
        return float(self.model.u - t)

    def kernel_weights(self, xi: Sequence[float] | Array, h: float) -> Array:
        xi_arr = as_array(xi).reshape(1, self.dim)
        return self.kernel.Kh(self.xs - xi_arr, h)

    def f_hat(self, xi: Sequence[float] | Array, h: float) -> float:
        return float(np.mean(self.kernel_weights(xi, h)))

    def _F_matrix(self, t: float, x_grid: Array, xi: Sequence[float] | Array) -> Array:
        x_grid = np.asarray(x_grid, dtype=float)
        if x_grid.ndim == 1:
            x_grid = x_grid.reshape(-1, self.dim)
        xi_arr = as_array(xi).reshape(1, 1, self.dim)
        xu = self.xu.reshape(1, self.n, self.dim)
        xg = x_grid.reshape(-1, 1, self.dim)
        dt = self.delta_t(t)
        delta = float(self.model.u - self.model.s)
        term1 = -np.sum((xu - xg) ** 2, axis=2) / (2.0 * dt)
        term2 = np.sum((xu - xi_arr) ** 2, axis=2) / (2.0 * delta)
        return np.exp(term1 + term2)

    def a_hat_grid(
        self,
        t: float,
        x_grid: Sequence[Sequence[float]] | Array,
        xi: Sequence[float] | Array,
        h: float,
        f_floor: float = 1e-12,
        d_floor: float = 1e-12,
        return_details: bool = False,
    ) -> Array | tuple[Array, dict[str, Array | float]]:
        x_grid = np.asarray(x_grid, dtype=float)
        if x_grid.ndim == 1:
            x_grid = x_grid.reshape(-1, self.dim)
        xi_arr = as_array(xi).reshape(self.dim)

        w = self.kernel_weights(xi_arr, h)
        fhat = max(float(np.mean(w)), f_floor)

        Fmat = self._F_matrix(t, x_grid, xi_arr)
        weighted = Fmat * w[None, :]
        g1 = np.mean(weighted, axis=1)
        g2 = np.mean(weighted[:, :, None] * self.xu[None, :, :], axis=1)

        Dhat = np.maximum(g1 / fhat, d_floor)
        Nhat = g2 / fhat
        ahat = (Nhat / Dhat[:, None] - x_grid) / self.delta_t(t)

        if return_details:
            return ahat, {
                "f_hat": fhat,
                "g1_hat": g1,
                "g2_hat": g2,
                "D_hat": Dhat,
                "N_hat": Nhat,
                "kernel_weights": w,
                "F_matrix": Fmat,
            }
        return ahat

    def point_details(
        self,
        t: float,
        x: Sequence[float] | Array,
        xi: Sequence[float] | Array,
        h: float,
        f_floor: float = 1e-12,
        d_floor: float = 1e-12,
    ) -> dict[str, Array | float]:
        x_arr = as_array(x).reshape(1, self.dim)
        ahat, details = self.a_hat_grid(t=t, x_grid=x_arr, xi=xi, h=h, f_floor=f_floor, d_floor=d_floor, return_details=True)
        details = dict(details)
        details["a_hat"] = np.asarray(ahat, dtype=float).reshape(self.dim)
        details["x"] = x_arr.reshape(self.dim)
        details["xi"] = as_array(xi).reshape(self.dim)
        details["t"] = float(t)
        return details

    def conditional_moment(self, xi: Sequence[float] | Array, h: float, values: Array, f_floor: float = 1e-12) -> Array:
        values = np.asarray(values, dtype=float)
        w = self.kernel_weights(xi, h)
        fhat = max(float(np.mean(w)), f_floor)
        if values.ndim == 1:
            return np.array([float(np.mean(values * w) / fhat)])
        return np.mean(values * w[:, None], axis=0) / fhat

    def a_hat_point(
        self,
        t: float,
        x: Sequence[float] | Array,
        xi: Sequence[float] | Array,
        h: float,
        f_floor: float = 1e-12,
        d_floor: float = 1e-12,
    ) -> Array:
        return np.asarray(self.point_details(t=t, x=x, xi=xi, h=h, f_floor=f_floor, d_floor=d_floor)["a_hat"], dtype=float)
