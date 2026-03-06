"""
Data models for the DUDI Python interface.

These mirror the Fortran derived types in src/define_types.f90:
- Point:    (r, alpha, beta, rvector, r_scaled, compute)
- Source:   (r, alphaM, betaM, rrM, zeta, eta, symmetry_axis, ejection_angle_distr, ud,
             sd, production_fun, production_rate, is_jet)

Units (DUDI convention):
- Distances in meters (moon-centered).
- Angles in radians.
- Ejection speeds in m/s.
- 3-vectors are numpy float64 arrays with shape (3,).
"""

from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from .typing import Vec3


def as_vec3(x, name: str = "vector") -> Vec3:
    v = np.asarray(x, dtype=np.float64)
    if v.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {v.shape}.")
    return v


def normalize(v: Vec3, atol: float = 1e-15) -> Vec3:
    v = as_vec3(v, "normalize(v)")
    n = float(np.linalg.norm(v))
    if not math.isfinite(n) or n <= atol:
        raise ValueError("Cannot normalize near-zero vector.")
    return (v / n).astype(np.float64)


def spherical_to_cartesian(r: float, alpha: float, beta: float) -> Vec3:
    sa, ca = math.sin(alpha), math.cos(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    return np.array([r * sa * cb, r * sa * sb, r * ca], dtype=np.float64)


@dataclass
class EjectionSpeedProperties:
    """Matches Fortran type(ejection_speed_properties).

    ud_shape : int
        Selector for ejection speed PDF.
    umin : float
        Minimum ejection speed [m/s], umin >= 0.
    umax : float
        Maximum ejection speed [m/s], umax >= umin.
    """

    ud_shape: int
    umin: float  # m/s
    umax: float  # m/s

    def __post_init__(self) -> None:
        if not isinstance(self.ud_shape, int):
            raise ValueError("ud_shape must be int.")
        if not (math.isfinite(self.umin) and self.umin >= 0.0):
            raise ValueError("umin must be finite and >= 0.")
        if not (math.isfinite(self.umax) and self.umax >= self.umin):
            raise ValueError("umax must be finite and >= umin.")


@dataclass
class Point:
    """
    Location where the dust density is evaluated (moon-centered).

    r : float
        Radial distance from the center of the moon [m], r >= 0.
    alpha : float
        Polar angle [rad].
    beta : float
        Eastern longitude [rad].
    rvector : Vec3
        Cartesian coordinates [m] in moon-centered frame, shape (3,).

    Matches Fortran position_in_space (r_scaled and compute are derived in the bridge).
    """

    r: float  # m
    alpha: float  # rad
    beta: float  # rad
    rvector: Vec3  # m, moon-centered

    def __post_init__(self) -> None:
        if not (math.isfinite(self.r) and self.r >= 0.0):
            raise ValueError("Point.r must be finite and >= 0.")
        if not (math.isfinite(self.alpha) and math.isfinite(self.beta)):
            raise ValueError("Point angles must be finite.")
        _ = as_vec3(self.rvector, "Point.rvector")


@dataclass
class Source:
    """
    Dust source definition (moon-centered).

    r, alphaM, betaM : float
        Source spherical coordinates [m, rad, rad].
    rrM : Vec3
        Source Cartesian position [m] in moon-centered frame, shape (3,).
    zeta, eta : float
        Zenith and azimuth [rad] of the ejection symmetry axis (see main README).
    symmetry_axis : Vec3
        Unit 3-vector along the ejection symmetry axis (moon-centered).
    ejection_angle_distr : int
        Selector for the ejection direction distribution (e.g. 1=Gaussian, 2=uniform cone, 3=parabola).
    ud : EjectionSpeedProperties
        Ejection speed PDF parameters (ud_shape, umin, umax) in m/s.
    sd : int
        Selector for dust size distribution.
    production_fun : int
        Selector for production rate function (<= 0 for constant).
    production_rate : float
        Production rate parameter (normalization), used by production_fun.
    is_jet : bool
        True if ejection is concentrated (e.g. omega < 0.1).

    Matches Fortran source_properties (ui, Gu_precalc are computed in the bridge).
    """

    r: float
    alphaM: float
    betaM: float
    rrM: Vec3
    zeta: float
    eta: float
    symmetry_axis: Vec3
    ejection_angle_distr: int
    ud: EjectionSpeedProperties
    sd: int
    production_fun: int
    production_rate: float
    is_jet: bool

    def __post_init__(self) -> None:
        if not (math.isfinite(self.r) and self.r >= 0.0):
            raise ValueError("Source.r must be finite and >= 0.")
        for name, val in (
            ("alphaM", self.alphaM),
            ("betaM", self.betaM),
            ("zeta", self.zeta),
            ("eta", self.eta),
        ):
            if not math.isfinite(val):
                raise ValueError(f"Source.{name} must be finite.")
        _ = as_vec3(self.rrM, "Source.rrM")
        ax = as_vec3(self.symmetry_axis, "Source.symmetry_axis")
        n = float(np.linalg.norm(ax))
        if not (math.isfinite(n) and abs(n - 1.0) <= 1e-9):
            raise ValueError("Source.symmetry_axis must be unit length.")
        if not isinstance(self.ejection_angle_distr, int):
            raise ValueError("ejection_angle_distr must be int.")
        if not isinstance(self.sd, int):
            raise ValueError("sd must be int.")
        if not isinstance(self.production_fun, int):
            raise ValueError("production_fun must be int.")
        if not math.isfinite(self.production_rate):
            raise ValueError("production_rate must be finite.")
        if not isinstance(self.is_jet, bool):
            raise ValueError("is_jet must be bool.")
