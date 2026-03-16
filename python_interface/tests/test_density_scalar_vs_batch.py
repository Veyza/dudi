"""
Test 1: scalar vs batched density for identical inputs.

Checks that calling `dudi.api.density` point-by-point and calling
`dudi.api.batch_over_points` with the same `Point`/`Source` data
produce numerically consistent results.
"""

from __future__ import annotations

import numpy as np

from dudi.api import (
    batch_over_points,
    batch_over_points_sources,
    batch_over_sources,
    density,
)
from dudi.models import EjectionSpeedProperties, Point, Source, spherical_to_cartesian


RM = 1.8216e6  # must match const.f90 rm


def _make_source(src_alphaM: float, src_betaM: float) -> Source:
    src_r = RM
    src_rrM = spherical_to_cartesian(src_r, src_alphaM, src_betaM)
    symmetry_axis = src_rrM / np.linalg.norm(src_rrM)
    return Source(
        r=src_r,
        alphaM=src_alphaM,
        betaM=src_betaM,
        rrM=src_rrM,
        zeta=0.1,
        eta=0.2,
        symmetry_axis=symmetry_axis,
        ejection_angle_distr=1,
        ud=EjectionSpeedProperties(ud_shape=1, umin=50.0, umax=200.0),
        sd=2,
        production_fun=0,
        production_rate=1e10,
        is_jet=True,
    )


def test_density_scalar_and_batch_match() -> None:
    """Scalar `density` vs `batch_over_points` should give the same numbers."""
    src = _make_source(0.5, 0.7)
    # three arbitrary points above the surface
    rs = RM * np.array([1.5, 1.7, 2.0])
    alphas = np.array([0.3, 1.0, 1.2])
    betas = np.array([0.1, 0.5, 1.0])
    pts = [
        Point(
            r=float(r),
            alpha=float(a),
            beta=float(b),
            rvector=spherical_to_cartesian(float(r), float(a), float(b)),
        )
        for r, a, b in zip(rs, alphas, betas)
    ]
    tnow = 600.0

    # scalar per-point
    scalar = np.array([density(p, src, tnow) for p in pts], dtype=float)
    # batched in Fortran
    batched = batch_over_points(pts, src, tnow)

    np.testing.assert_allclose(batched, scalar, rtol=1e-6, atol=0.0)
    print("Density scalar and batch match")


def test_density_sources_and_points_match() -> None:
    """
    `batch_over_sources` (per-source contributions) and `batch_over_points_sources`
    (sum over sources per point) should be consistent with scalar `density` calls.
    """
    # Single test point above the surface
    r = 1.9 * RM
    alpha = 0.4
    beta = 0.9
    point = Point(
        r=float(r),
        alpha=float(alpha),
        beta=float(beta),
        rvector=spherical_to_cartesian(float(r), float(alpha), float(beta)),
    )

    # Several sources with different locations
    alphasM = np.array([0.3, 0.8, 1.1])
    betasM = np.array([0.2, 0.9, 1.3])
    sources = [_make_source(float(aM), float(bM)) for aM, bM in zip(alphasM, betasM)]

    tnow = 700.0

    # Scalar per-source densities at this point, then sum over sources
    scalar_per_source = np.vstack([density(point, src, tnow) for src in sources])
    summed_from_scalars = scalar_per_source.sum(axis=0)

    # Batched over sources: per-source contributions at this point
    batched_sources = batch_over_sources(point, sources, tnow)

    # Batched over points and sources, summing contributions per point
    batched_point_total = batch_over_points_sources([point], sources, tnow)[0]

    # Per-source consistency
    np.testing.assert_allclose(batched_sources, scalar_per_source, rtol=1e-6, atol=0.0)

    # Sum over sources: all three ways of computing the total must agree
    summed_from_batched = batched_sources.sum(axis=0)
    np.testing.assert_allclose(
        summed_from_batched, summed_from_scalars, rtol=1e-6, atol=0.0
    )
    np.testing.assert_allclose(
        batched_point_total, summed_from_scalars, rtol=1e-6, atol=0.0
    )
