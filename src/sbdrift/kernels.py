from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

Array = np.ndarray
KernelName = Literal["epanechnikov"]


@dataclass(frozen=True)
class ProductKernel:
    name: KernelName = "epanechnikov"

    def base(self, u: Array) -> Array:
        if self.name != "epanechnikov":
            raise ValueError(f"Unsupported kernel: {self.name}")
        out = 0.75 * (1.0 - u**2)
        out = np.where(np.abs(u) <= 1.0, out, 0.0)
        return out

    def K(self, z: Array) -> Array:
        """Product kernel K(z) for z of shape (..., d)."""
        z = np.asarray(z, dtype=float)
        if z.ndim == 0:
            z = z.reshape(1, 1)
        if z.ndim == 1:
            z = z.reshape(-1, 1)
        vals = self.base(z)
        return np.prod(vals, axis=-1)

    def Kh(self, z: Array, h: float) -> Array:
        """Scaled kernel h^{-d} K(z/h) for z of shape (..., d)."""
        z = np.asarray(z, dtype=float)
        if z.ndim == 0:
            z = z.reshape(1, 1)
        if z.ndim == 1:
            z = z.reshape(-1, 1)
        d = z.shape[-1]
        return (h ** (-d)) * self.K(z / h)


def make_kernel(name: KernelName = "epanechnikov") -> ProductKernel:
    return ProductKernel(name=name)
