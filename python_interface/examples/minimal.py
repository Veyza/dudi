"""
Minimal example: initialize a point and source, call DUDI density and mean_velocity.

Units: distances [m], angles [rad], time [s], speeds [m/s].

Prereqs:
  1) Build the Fortran ctypes bridge:
       bash python_interface/fortran_bridge/build_ctypes_bridge.sh
  2) (optional) Install in editable mode:
       pip install -e .

Run from repo root:
  python python_interface/examples/minimal.py
"""

import numpy as np

from dudi.api import density, mean_velocity, batch_over_points, batch_over_sources
from dudi.models import Point, Source, EjectionSpeedProperties, spherical_to_cartesian

# Moon radius (e.g. Enceladus), meters
RM = 1.8216e6

# Observation point: above the pole, 2 moon radii
point_r = 1.5 * RM
point_alpha = 0.1
point_beta = 0.1
point_rvector = spherical_to_cartesian(point_r, point_alpha, point_beta)
pt = Point(r=point_r, alpha=point_alpha, beta=point_beta, rvector=point_rvector)

# Source on the surface at the pole; symmetry axis = outward radial (0, 0, 1)
src_r = RM
src_alphaM = 0.12
src_betaM = 0.12
src_rrM = spherical_to_cartesian(src_r, src_alphaM, src_betaM)
src_zeta = 0.0
src_eta = 0.0
symmetry_axis = src_rrM / src_r  # unit vector
src = Source(
    r=src_r,
    alphaM=src_alphaM,
    betaM=src_betaM,
    rrM=src_rrM,
    zeta=src_zeta,
    eta=src_eta,
    symmetry_axis=symmetry_axis,
    ejection_angle_distr=1,  # 1=Gaussian, 2=uniform cone, 3=parabola
    ud=EjectionSpeedProperties(ud_shape=1, umin=50.0, umax=500.0),  # m/s
    sd=2,
    production_fun=0,  # constant production
    production_rate=1e10,
    is_jet=True,
)

tnow = 600.0  # seconds

# Scalar calls
d_bound, d_unbound = density(pt, src, tnow)
print("density(point, source, tnow):")
print(f"  bound   = {d_bound}")
print(f"  unbound = {d_unbound}")

integral = mean_velocity(pt, src, tnow)
v_up = integral[0:3]
v_down = integral[3:6]
n_up, n_down = integral[6], integral[7]
print("mean_velocity(point, source, tnow):")
print(f"  v_up   = {v_up}")
print(f"  v_down = {v_down}")
print(f"  n_up   = {n_up}, n_down = {n_down}")
