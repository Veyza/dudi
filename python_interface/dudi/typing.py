"""
Typing aliases for the DUDI Python interface.

- FloatArray: np.ndarray of float64 (any shape)
- Vec3:       np.ndarray of float64 with shape (3,) — used for positions/vectors in meters
              (moon-centered). Shape is enforced at runtime where relevant (e.g. via as_vec3).
"""

from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
Vec3 = NDArray[np.float64]
