"""
Test 2: mean-velocity batched API shape and finiteness.

Checks that `mean_velocity` and its batched counterpart
`mean_velocity_batch_points` return arrays of the documented shape
and with finite (non-NaN, non-inf) entries.
"""

from __future__ import annotations

import numpy as np

from dudi.api import mean_velocity, mean_velocity_batch_points
from dudi.models import EjectionSpeedProperties, Point, Source, spherical_to_cartesian


RM = 1.8216e6  # must match const.f90 rm


def _make_source() -> Source:
    src_r = RM
    src_alphaM = 0.8
    src_betaM = 1.2
    src_rrM = spherical_to_cartesian(src_r, src_alphaM, src_betaM)
    symmetry_axis = src_rrM / np.linalg.norm(src_rrM)
    return Source(
        r=src_r,
        alphaM=src_alphaM,
        betaM=src_betaM,
        rrM=src_rrM,
        zeta=0.2,
        eta=0.3,
        symmetry_axis=symmetry_axis,
        ejection_angle_distr=1,
        ud=EjectionSpeedProperties(ud_shape=1, umin=100.0, umax=300.0),
        sd=2,
        production_fun=0,
        production_rate=5e9,
        is_jet=True,
    )


def test_mean_velocity_scalar_and_batch_shapes_finite() -> None:
    """Batched mean_velocity returns (N, 8) and agrees with scalar call."""
    src = _make_source()
    rs = RM * np.array([1.3, 1.6, 1.9, 2.2])
    alphas = np.array([0.4, 0.6, 1.0, 1.2])
    betas = np.array([0.0, 0.4, 0.8, 1.2])
    pts = [
        Point(
            r=float(r),
            alpha=float(a),
            beta=float(b),
            rvector=spherical_to_cartesian(float(r), float(a), float(b)),
        )
        for r, a, b in zip(rs, alphas, betas)
    ]
    tnow = 400.0

    scalar = np.vstack([mean_velocity(p, src, tnow) for p in pts])
    batched = mean_velocity_batch_points(pts, src, tnow)

    assert batched.shape == (len(pts), 8)
    np.testing.assert_allclose(batched, scalar, rtol=1e-6, atol=0.0)
    assert np.isfinite(batched).all()
