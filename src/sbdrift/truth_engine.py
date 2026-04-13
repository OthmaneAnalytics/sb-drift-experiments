from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from scipy.integrate import quad

from .models import BaseModel
from .utils import as_array


Array = np.ndarray


@dataclass
class TruthEngine:
    model: BaseModel

    def delta_t(self, t: float) -> float:
        return float(self.model.u - t)

    def F(self, t: float, xi: Sequence[float] | Array, x: Sequence[float] | Array, y: Sequence[float] | Array) -> float:
        xi = as_array(xi).reshape(self.model.dim)
        x = as_array(x).reshape(self.model.dim)
        y = as_array(y).reshape(self.model.dim)
        dt = self.delta_t(t)
        delta = self.model.u - self.model.s
        term1 = -np.sum((y - x) ** 2) / (2.0 * dt)
        term2 = np.sum((y - xi) ** 2) / (2.0 * delta)
        return float(np.exp(term1 + term2))

    def D_star(self, t: float, x: Sequence[float] | Array, xi: Sequence[float] | Array, grid_points_2d: int = 201) -> float:
        dim = self.model.dim
        if dim == 1:
            low = float(self.model.low[0])
            high = float(self.model.high[0])
            xi_arr = as_array(xi).reshape(1)
            x_arr = as_array(x).reshape(1)

            qpdf = self.model.conditional_pdf_fn(xi_arr)

            def integrand(y: float) -> float:
                return self.F(t, xi_arr, x_arr, np.array([y])) * qpdf(np.array([y]))

            val, _ = quad(integrand, low, high, epsabs=1e-10, epsrel=1e-10, limit=200)
            return float(val)

        if dim == 2:
            return float(self._integrate_2d_scalar(lambda y: self.F(t, xi, x, y), xi, grid_points_2d))

        raise NotImplementedError("TruthEngine currently supports d=1 or d=2 for deterministic truth integration.")

    def N_star(self, t: float, x: Sequence[float] | Array, xi: Sequence[float] | Array, grid_points_2d: int = 201) -> Array:
        dim = self.model.dim
        if dim == 1:
            low = float(self.model.low[0])
            high = float(self.model.high[0])
            xi_arr = as_array(xi).reshape(1)
            x_arr = as_array(x).reshape(1)

            qpdf = self.model.conditional_pdf_fn(xi_arr)

            def integrand(y: float) -> float:
                yy = np.array([y])
                return y * self.F(t, xi_arr, x_arr, yy) * qpdf(yy)

            val, _ = quad(integrand, low, high, epsabs=1e-10, epsrel=1e-10, limit=200)
            return np.array([float(val)])

        if dim == 2:
            res = self._integrate_2d_vector(lambda y: y * self.F(t, xi, x, y), xi, grid_points_2d)
            return res

        raise NotImplementedError("TruthEngine currently supports d=1 or d=2 for deterministic truth integration.")

    def a_star(self, t: float, x: Sequence[float] | Array, xi: Sequence[float] | Array, grid_points_2d: int = 201) -> Array:
        x_arr = as_array(x).reshape(self.model.dim)
        D = self.D_star(t, x_arr, xi, grid_points_2d=grid_points_2d)
        N = self.N_star(t, x_arr, xi, grid_points_2d=grid_points_2d)
        return (N / D - x_arr) / self.delta_t(t)

    def _integrate_2d_scalar(self, func: Callable[[Array], float], xi: Sequence[float] | Array, grid_points: int) -> float:
        low = self.model.low
        high = self.model.high
        ys0 = np.linspace(low[0], high[0], grid_points)
        ys1 = np.linspace(low[1], high[1], grid_points)
        d0 = ys0[1] - ys0[0]
        d1 = ys1[1] - ys1[0]
        out = np.zeros((grid_points, grid_points), dtype=float)
        xi_arr = as_array(xi).reshape(2)
        qpdf = self.model.conditional_pdf_fn(xi_arr)
        for i, y0 in enumerate(ys0):
            for j, y1 in enumerate(ys1):
                y = np.array([y0, y1], dtype=float)
                out[i, j] = func(y) * qpdf(y)
        return float(np.trapezoid(np.trapezoid(out, ys1, axis=1), ys0, axis=0))

    def _integrate_2d_vector(self, func: Callable[[Array], Array], xi: Sequence[float] | Array, grid_points: int) -> Array:
        low = self.model.low
        high = self.model.high
        ys0 = np.linspace(low[0], high[0], grid_points)
        ys1 = np.linspace(low[1], high[1], grid_points)
        xi_arr = as_array(xi).reshape(2)
        qpdf = self.model.conditional_pdf_fn(xi_arr)
        out0 = np.zeros((grid_points, grid_points), dtype=float)
        out1 = np.zeros((grid_points, grid_points), dtype=float)
        for i, y0 in enumerate(ys0):
            for j, y1 in enumerate(ys1):
                y = np.array([y0, y1], dtype=float)
                val = func(y) * qpdf(y)
                out0[i, j] = val[0]
                out1[i, j] = val[1]
        v0 = np.trapezoid(np.trapezoid(out0, ys1, axis=1), ys0, axis=0)
        v1 = np.trapezoid(np.trapezoid(out1, ys1, axis=1), ys0, axis=0)
        return np.array([float(v0), float(v1)])
