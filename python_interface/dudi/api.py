"""
Thin Python API for DUDI (moon-centered dust density and mean velocity).

Maps typed dataclasses (Point, Source) to the ctypes bridge. All coordinates
and distances in meters, angles in rad, time in seconds.
"""

from __future__ import annotations

from typing import Sequence
import numpy as np

from .models import Point, Source
from ._bridge_ctypes import (
    call_density as _call_density,
    call_mean_velocity as _call_mean_velocity,
    call_batch_points as _call_batch_points,
    call_batch_sources as _call_batch_sources,
    call_batch_sources_points as _call_batch_sources_points,
    call_mean_velocity_batch_points as _call_mean_velocity_batch_points,
    call_mean_velocity_batch_sources as _call_mean_velocity_batch_sources,
    call_mean_velocity_batch_sources_points as _call_mean_velocity_batch_sources_points,
)


def _source_kwargs(source: Source) -> dict:
    """Bridge keyword args for a single source (scalars)."""
    return {
        "src_r": float(source.r),
        "src_alphaM": float(source.alphaM),
        "src_betaM": float(source.betaM),
        "src_zeta": float(source.zeta),
        "src_eta": float(source.eta),
        "src_eject_distr": int(source.ejection_angle_distr),
        "src_ud_shape": int(source.ud.ud_shape),
        "src_umin": float(source.ud.umin),
        "src_umax": float(source.ud.umax),
        "src_sd": int(source.sd),
        "src_production_fun": int(source.production_fun),
        "src_production_rate": float(source.production_rate),
        "src_is_jet": bool(source.is_jet),
    }


def _point_kwargs(point: Point) -> dict:
    """Bridge keyword args for a single point (scalars)."""
    return {
        "point_r": float(point.r),
        "point_alpha": float(point.alpha),
        "point_beta": float(point.beta),
    }


# ----------------------------------------------------------------------
# Scalar API
# ----------------------------------------------------------------------


def density(point: Point, source: Source, tnow: float) -> tuple[float, float]:
    """
    Dust number density at a single point from a single source.

    Returns
    -------
    (density_bound, density_unbound) : tuple of float
    """
    return _call_density(
        tnow=float(tnow),
        **_point_kwargs(point),
        **_source_kwargs(source),
    )


def mean_velocity(point: Point, source: Source, tnow: float) -> np.ndarray:
    """
    Mean velocity and number density at a single point from a single source.

    Returns
    -------
    integral : ndarray, shape (8,)
        [v_up_x, v_up_y, v_up_z, v_down_x, v_down_y, v_down_z, n_up, n_down]
    """
    return _call_mean_velocity(
        tnow=float(tnow),
        **_point_kwargs(point),
        **_source_kwargs(source),
    )


# ----------------------------------------------------------------------
# Batched API: density
# ----------------------------------------------------------------------


def batch_over_points(
    points: Sequence[Point],
    source: Source,
    tnow: float,
) -> np.ndarray:
    """
    Density at many points for a single source.

    Returns
    -------
    densities : ndarray, shape (N, 2)
        [bound, unbound] per point.
    """
    pts = list(points)
    if not pts:
        return np.empty((0, 2), dtype=np.float64)

    point_r = np.array([float(p.r) for p in pts], dtype=np.float64)
    point_alpha = np.array([float(p.alpha) for p in pts], dtype=np.float64)
    point_beta = np.array([float(p.beta) for p in pts], dtype=np.float64)

    return _call_batch_points(
        point_r=point_r,
        point_alpha=point_alpha,
        point_beta=point_beta,
        tnow=float(tnow),
        **_source_kwargs(source),
    )


def batch_over_sources(
    point: Point,
    sources: Sequence[Source],
    tnow: float,
) -> np.ndarray:
    """
    Density at a single point from many sources (one result per source).

    Returns
    -------
    densities : ndarray, shape (N, 2)
        [bound, unbound] per source.
    """
    srcs = list(sources)
    if not srcs:
        return np.empty((0, 2), dtype=np.float64)

    src_r = np.array([float(s.r) for s in srcs], dtype=np.float64)
    src_alphaM = np.array([float(s.alphaM) for s in srcs], dtype=np.float64)
    src_betaM = np.array([float(s.betaM) for s in srcs], dtype=np.float64)
    src_zeta = np.array([float(s.zeta) for s in srcs], dtype=np.float64)
    src_eta = np.array([float(s.eta) for s in srcs], dtype=np.float64)
    src_eject_distr = np.array(
        [int(s.ejection_angle_distr) for s in srcs], dtype=np.int32
    )
    src_ud_shape = np.array([int(s.ud.ud_shape) for s in srcs], dtype=np.int32)
    src_umin = np.array([float(s.ud.umin) for s in srcs], dtype=np.float64)
    src_umax = np.array([float(s.ud.umax) for s in srcs], dtype=np.float64)
    src_sd = np.array([int(s.sd) for s in srcs], dtype=np.int32)
    src_production_fun = np.array([int(s.production_fun) for s in srcs], dtype=np.int32)
    src_production_rate = np.array(
        [float(s.production_rate) for s in srcs], dtype=np.float64
    )
    src_is_jet = np.array([1 if s.is_jet else 0 for s in srcs], dtype=np.int32)

    return _call_batch_sources(
        point_r=float(point.r),
        point_alpha=float(point.alpha),
        point_beta=float(point.beta),
        src_r=src_r,
        src_alphaM=src_alphaM,
        src_betaM=src_betaM,
        src_zeta=src_zeta,
        src_eta=src_eta,
        src_eject_distr=src_eject_distr,
        src_ud_shape=src_ud_shape,
        src_umin=src_umin,
        src_umax=src_umax,
        src_sd=src_sd,
        src_production_fun=src_production_fun,
        src_production_rate=src_production_rate,
        src_is_jet=src_is_jet,
        tnow=float(tnow),
    )


def batch_over_points_sources(
    points: Sequence[Point],
    sources: Sequence[Source],
    tnow: float,
) -> np.ndarray:
    """
    Total density at many points from many sources (contributions summed per point).

    Returns
    -------
    densities : ndarray, shape (n_points, 2)
        [bound, unbound] per point (sum over all sources).
    """
    pts = list(points)
    srcs = list(sources)
    if not pts:
        return np.empty((0, 2), dtype=np.float64)
    if not srcs:
        return np.zeros((len(pts), 2), dtype=np.float64)

    n_points = len(pts)
    n_sources = len(srcs)

    point_r = np.array([float(p.r) for p in pts], dtype=np.float64)
    point_alpha = np.array([float(p.alpha) for p in pts], dtype=np.float64)
    point_beta = np.array([float(p.beta) for p in pts], dtype=np.float64)

    src_r = np.array([float(s.r) for s in srcs], dtype=np.float64)
    src_alphaM = np.array([float(s.alphaM) for s in srcs], dtype=np.float64)
    src_betaM = np.array([float(s.betaM) for s in srcs], dtype=np.float64)
    src_zeta = np.array([float(s.zeta) for s in srcs], dtype=np.float64)
    src_eta = np.array([float(s.eta) for s in srcs], dtype=np.float64)
    src_eject_distr = np.array(
        [int(s.ejection_angle_distr) for s in srcs], dtype=np.int32
    )
    src_ud_shape = np.array([int(s.ud.ud_shape) for s in srcs], dtype=np.int32)
    src_umin = np.array([float(s.ud.umin) for s in srcs], dtype=np.float64)
    src_umax = np.array([float(s.ud.umax) for s in srcs], dtype=np.float64)
    src_sd = np.array([int(s.sd) for s in srcs], dtype=np.int32)
    src_production_fun = np.array([int(s.production_fun) for s in srcs], dtype=np.int32)
    src_production_rate = np.array(
        [float(s.production_rate) for s in srcs], dtype=np.float64
    )
    src_is_jet = np.array([1 if s.is_jet else 0 for s in srcs], dtype=np.int32)

    return _call_batch_sources_points(
        point_r=point_r,
        point_alpha=point_alpha,
        point_beta=point_beta,
        src_r=src_r,
        src_alphaM=src_alphaM,
        src_betaM=src_betaM,
        src_zeta=src_zeta,
        src_eta=src_eta,
        src_eject_distr=src_eject_distr,
        src_ud_shape=src_ud_shape,
        src_umin=src_umin,
        src_umax=src_umax,
        src_sd=src_sd,
        src_production_fun=src_production_fun,
        src_production_rate=src_production_rate,
        src_is_jet=src_is_jet,
        tnow=float(tnow),
    )


# ----------------------------------------------------------------------
# Batched API: mean velocity
# ----------------------------------------------------------------------


def mean_velocity_batch_points(
    points: Sequence[Point],
    source: Source,
    tnow: float,
) -> np.ndarray:
    """
    Mean velocity and number density at many points for a single source.

    Returns
    -------
    integral : ndarray, shape (N, 8)
        [v_up(3), v_down(3), n_up, n_down] per point.
    """
    pts = list(points)
    if not pts:
        return np.empty((0, 8), dtype=np.float64)

    point_r = np.array([float(p.r) for p in pts], dtype=np.float64)
    point_alpha = np.array([float(p.alpha) for p in pts], dtype=np.float64)
    point_beta = np.array([float(p.beta) for p in pts], dtype=np.float64)

    return _call_mean_velocity_batch_points(
        point_r=point_r,
        point_alpha=point_alpha,
        point_beta=point_beta,
        tnow=float(tnow),
        **_source_kwargs(source),
    )


def mean_velocity_batch_sources(
    point: Point,
    sources: Sequence[Source],
    tnow: float,
) -> np.ndarray:
    """
    Mean velocity and number density at a single point from many sources.

    Returns
    -------
    integral : ndarray, shape (N, 8)
        [v_up(3), v_down(3), n_up, n_down] per source.
    """
    srcs = list(sources)
    if not srcs:
        return np.empty((0, 8), dtype=np.float64)

    n = len(srcs)
    src_r = np.array([float(s.r) for s in srcs], dtype=np.float64)
    src_alphaM = np.array([float(s.alphaM) for s in srcs], dtype=np.float64)
    src_betaM = np.array([float(s.betaM) for s in srcs], dtype=np.float64)
    src_zeta = np.array([float(s.zeta) for s in srcs], dtype=np.float64)
    src_eta = np.array([float(s.eta) for s in srcs], dtype=np.float64)
    src_eject_distr = np.array(
        [int(s.ejection_angle_distr) for s in srcs], dtype=np.int32
    )
    src_ud_shape = np.array([int(s.ud.ud_shape) for s in srcs], dtype=np.int32)
    src_umin = np.array([float(s.ud.umin) for s in srcs], dtype=np.float64)
    src_umax = np.array([float(s.ud.umax) for s in srcs], dtype=np.float64)
    src_sd = np.array([int(s.sd) for s in srcs], dtype=np.int32)
    src_production_fun = np.array([int(s.production_fun) for s in srcs], dtype=np.int32)
    src_production_rate = np.array(
        [float(s.production_rate) for s in srcs], dtype=np.float64
    )
    src_is_jet = np.array([1 if s.is_jet else 0 for s in srcs], dtype=np.int32)

    return _call_mean_velocity_batch_sources(
        point_r=float(point.r),
        point_alpha=float(point.alpha),
        point_beta=float(point.beta),
        src_r=src_r,
        src_alphaM=src_alphaM,
        src_betaM=src_betaM,
        src_zeta=src_zeta,
        src_eta=src_eta,
        src_eject_distr=src_eject_distr,
        src_ud_shape=src_ud_shape,
        src_umin=src_umin,
        src_umax=src_umax,
        src_sd=src_sd,
        src_production_fun=src_production_fun,
        src_production_rate=src_production_rate,
        src_is_jet=src_is_jet,
        tnow=float(tnow),
    )


def mean_velocity_batch_points_sources(
    points: Sequence[Point],
    sources: Sequence[Source],
    tnow: float,
) -> np.ndarray:
    """
    Mean velocity and number density at many points from many sources (summed per point).

    Returns
    -------
    integral : ndarray, shape (n_points, 8)
        [v_up(3), v_down(3), n_up, n_down] per point (sum over all sources).
    """
    pts = list(points)
    srcs = list(sources)
    if not pts:
        return np.empty((0, 8), dtype=np.float64)
    if not srcs:
        return np.zeros((len(pts), 8), dtype=np.float64)

    point_r = np.array([float(p.r) for p in pts], dtype=np.float64)
    point_alpha = np.array([float(p.alpha) for p in pts], dtype=np.float64)
    point_beta = np.array([float(p.beta) for p in pts], dtype=np.float64)

    src_r = np.array([float(s.r) for s in srcs], dtype=np.float64)
    src_alphaM = np.array([float(s.alphaM) for s in srcs], dtype=np.float64)
    src_betaM = np.array([float(s.betaM) for s in srcs], dtype=np.float64)
    src_zeta = np.array([float(s.zeta) for s in srcs], dtype=np.float64)
    src_eta = np.array([float(s.eta) for s in srcs], dtype=np.float64)
    src_eject_distr = np.array(
        [int(s.ejection_angle_distr) for s in srcs], dtype=np.int32
    )
    src_ud_shape = np.array([int(s.ud.ud_shape) for s in srcs], dtype=np.int32)
    src_umin = np.array([float(s.ud.umin) for s in srcs], dtype=np.float64)
    src_umax = np.array([float(s.ud.umax) for s in srcs], dtype=np.float64)
    src_sd = np.array([int(s.sd) for s in srcs], dtype=np.int32)
    src_production_fun = np.array([int(s.production_fun) for s in srcs], dtype=np.int32)
    src_production_rate = np.array(
        [float(s.production_rate) for s in srcs], dtype=np.float64
    )
    src_is_jet = np.array([1 if s.is_jet else 0 for s in srcs], dtype=np.int32)

    return _call_mean_velocity_batch_sources_points(
        point_r=point_r,
        point_alpha=point_alpha,
        point_beta=point_beta,
        src_r=src_r,
        src_alphaM=src_alphaM,
        src_betaM=src_betaM,
        src_zeta=src_zeta,
        src_eta=src_eta,
        src_eject_distr=src_eject_distr,
        src_ud_shape=src_ud_shape,
        src_umin=src_umin,
        src_umax=src_umax,
        src_sd=src_sd,
        src_production_fun=src_production_fun,
        src_production_rate=src_production_rate,
        src_is_jet=src_is_jet,
        tnow=float(tnow),
    )
