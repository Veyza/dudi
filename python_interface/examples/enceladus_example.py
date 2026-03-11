"""
Reproduce the Fortran example `examples/enceladus_example.f90` and the plot from
`scripts/e2plot.py` using the DUDI Python interface.

Prereqs:
  1) Build the Fortran ctypes bridge:
       bash python_interface/fortran_bridge/build_ctypes_bridge.sh
  2) (optional) Install in editable mode:
       pip install -e .

Run from the repo root:
  python python_interface/examples/enceladus.py
"""

from __future__ import annotations

import os
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from dudi.api import batch_over_points
from dudi.models import (
    EjectionSpeedProperties,
    Point,
    Source,
    spherical_to_cartesian,
)


# Constants matching those used in the Fortran code (const.f90)
PI = float(np.pi)
HALFPI = 0.5 * PI
DEG2RAD = PI / 180.0
RM = 252e3  # moon radius [m]
BG = 0.01  # background density added in cassini_flyby_out


def load_sources(path: str = "./input_data_files/Enceladus_jet.dat") -> List[Source]:
    """
    Read Enceladus jet source parameters from file, mimicking read_sources_params.

    File columns (per source):
        alphaM_deg, betaM_deg, zeta_deg, eta_deg,
        production_fun, production_rate,
        ud_shape, umin, umax,
        ejection_angle_distr, sd
    """
    arr = np.loadtxt(path)
    arr = np.atleast_2d(arr)

    sources: List[Source] = []
    for row in arr:
        (
            alpha_deg,
            beta_deg,
            zeta_deg,
            eta_deg,
            production_fun,
            production_rate,
            ud_shape,
            umin,
            umax,
            ead,
            sd,
        ) = row

        alphaM = HALFPI - alpha_deg * DEG2RAD
        betaM = beta_deg * DEG2RAD
        zeta = zeta_deg * DEG2RAD
        eta = eta_deg * DEG2RAD

        r = RM
        rrM = spherical_to_cartesian(r, alphaM, betaM)

        # The Fortran bridge recomputes symmetry_axis and Gu_precalc internally,
        # so we only need a unit vector placeholder here.
        symmetry_axis = rrM / np.linalg.norm(rrM)

        src = Source(
            r=float(r),
            alphaM=float(alphaM),
            betaM=float(betaM),
            rrM=rrM,
            zeta=float(zeta),
            eta=float(eta),
            symmetry_axis=symmetry_axis,
            ejection_angle_distr=int(ead),
            ud=EjectionSpeedProperties(
                ud_shape=int(ud_shape),
                umin=float(umin),
                umax=float(umax),
            ),
            sd=int(sd),
            production_fun=int(production_fun),
            production_rate=float(production_rate),
            is_jet=True,
        )
        sources.append(src)

    return sources


def load_cassini_E2(
    path: str = "./input_data_files/Cassini_E2_flyby.dat",
) -> Tuple[np.ndarray, List[Point]]:
    """
    Read Cassini E2 trajectory, mimicking read_Cassini_E2.

    File columns per line:
        ttab, r, alpha_deg, beta_deg
    """
    ttab, r, alpha_deg, beta_deg = np.loadtxt(path, unpack=True)

    alpha = HALFPI - alpha_deg * DEG2RAD
    beta = beta_deg * DEG2RAD

    points: List[Point] = []
    for ri, ai, bi in zip(r, alpha, beta):
        rv = spherical_to_cartesian(float(ri), float(ai), float(bi))
        points.append(
            Point(
                r=float(ri),
                alpha=float(ai),
                beta=float(bi),
                rvector=rv,
            )
        )

    return ttab, points


def main() -> None:
    # Load source(s) and trajectory
    sources = load_sources()
    ttab, points = load_cassini_E2()

    if not sources:
        raise RuntimeError("No sources loaded from Enceladus_jet.dat")
    if not points:
        raise RuntimeError("No points loaded from Cassini_E2_flyby.dat")

    tnow = 0.0  # as in enceladus_example.f90

    # Sum contributions from all sources (Ns is small; loop in Python is fine)
    nt = len(points)
    density_total = np.zeros((nt, 2), dtype=np.float64)
    for src in sources:
        dens_src = batch_over_points(points, src, tnow)
        density_total += dens_src

    # Match cassini_flyby_out: write ttab and sum(density) + BG
    total_number_density = density_total.sum(axis=1) + BG

    results_dir = "./results"
    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, "E2_profile.dat")
    np.savetxt(
        out_path,
        np.column_stack([ttab, total_number_density]),
        fmt=["%9.3f", "%11.4e"],
    )
    print(f"Wrote model profile to {out_path}")

    # Reproduce scripts/e2plot.py: compare model with HRD measurements
    model_t = ttab
    model_dens = total_number_density

    hrd_data = np.loadtxt("./input_data_files/E2_1.6.txt", usecols=range(2))
    hrd_t = hrd_data[:, 0]
    hrd_dens = hrd_data[:, 1]

    plt.figure(1)
    plt.ylabel("number density of grains $> 1.6\\ \\mu m$")
    plt.xlabel("seconds from the moment of the closest approach")
    plt.ylim(0, 0.12)
    plt.suptitle("HRD number density profile of E2 flyby")
    plt.plot(model_t, model_dens, "k-", label="model number density")
    plt.plot(hrd_t, hrd_dens, "bo", label="HRD measurements")
    plt.legend(loc="upper left")
    plt.show()


if __name__ == "__main__":
    main()
