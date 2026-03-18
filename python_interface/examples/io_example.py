"""
Reproduce the Fortran example `examples/io_example.f90` (Io volcano images)
and the plots from `scripts/volcano_image.py` using the DUDI Python interface.

Prereqs:
  1) Build the Fortran ctypes bridge:
       bash python_interface/fortran_bridge/build_ctypes_bridge.sh
  2) (optional) Install in editable mode:
       pip install -e .

Run from the repo root:
  python python_interface/examples/io_volcano.py
"""

from __future__ import annotations

import os
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from dudi.api import batch_over_points
from dudi.models import EjectionSpeedProperties, Point, Source


# ----------------------------------------------------------------------
# Constants (match const.f90 and image_construction.f90)
# ----------------------------------------------------------------------

PI = float(np.pi)
TWOPI = 2.0 * PI
HALFPI = 0.5 * PI
SQRT2 = float(np.sqrt(2.0))
DEG2RAD = PI / 180.0

RM = 1.8216e6  # moon radius [m]

# image_construction constants
Vrad = 4.5e-2
Hrad = 4.5e-2
Vpix = 64
Hpix = 64
ni = 20
sampdist_large = 3.0e2  # meters

BG_IMAGE = 1e-15  # background added in io_example


def myatan1(re0: float, im: float) -> float:
    """Modified atan returning angle in [0, 2pi], see help.f90/myatan1."""
    re = float(re0)
    # avoid problems with zero re
    if 0.0 < re < 1e-12:
        re = 1e-12
    if -1e-12 < re < 0.0:
        re = -1e-12
    ang = float(np.arctan2(im, re))
    if re < 0.0:
        ang += PI
    if re > 0.0 and im < 0.0:
        ang += TWOPI
    return ang


def dist_between_2lines(
    M1: np.ndarray, s1: np.ndarray, M2: np.ndarray, s2: np.ndarray
) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Python port of help.f90/dist_between_2lines.

    Returns (d, K1, K2) where K1 and K2 are closest points on the 1st and 2nd lines,
    and d is the shortest distance between them.
    """
    s = np.cross(s1, s2)
    if float(np.linalg.norm(s)) < 1e-3:
        # nearly parallel lines
        K2 = M2.copy()
        vtmp = M2 - M1
        K1 = M1 + s1 * float(np.linalg.norm(vtmp))
        vtmp0 = np.cross(s1, vtmp)
        d = float(np.linalg.norm(vtmp0))
        return d, K1, K2

    t1 = np.cross(s1, s)
    t2 = np.cross(s2, s)

    # First intersection (K2)
    A = np.zeros((3, 3), dtype=float)
    B = np.zeros(3, dtype=float)

    A[0, :] = t1
    B[0] = float(np.dot(t1, M1))

    A[1, 0] = s2[1]
    A[1, 1] = -s2[0]
    A[1, 2] = 0.0
    B[1] = -(s2[0] * M2[1] - s2[1] * M2[0])

    A[2, 0] = s2[2]
    A[2, 1] = 0.0
    A[2, 2] = -s2[0]
    B[2] = -(s2[0] * M2[2] - s2[2] * M2[0])

    Ainv = np.linalg.inv(A)
    K2 = Ainv @ B

    # Second intersection (K1)
    A[0, :] = t2
    B[0] = float(np.dot(t2, M2))

    A[1, 0] = -s1[1]
    A[1, 1] = s1[0]
    A[1, 2] = 0.0
    B[1] = s1[0] * M1[1] - s1[1] * M1[0]

    A[2, 0] = -s1[2]
    A[2, 1] = 0.0
    A[2, 2] = s1[0]
    B[2] = s1[0] * M1[2] - s1[2] * M1[0]

    Ainv = np.linalg.inv(A)
    K1 = Ainv @ B

    vtmp = K1 - K2
    d = float(np.linalg.norm(vtmp))
    return d, K1, K2


# ----------------------------------------------------------------------
# Volcano source parameters (get_volcano_params)
# ----------------------------------------------------------------------


def get_volcano_source() -> Source:
    """
    Python version of inputdata.get_volcano_params for Ns = 1.
    """
    alphaM = 5.0 * DEG2RAD
    betaM = 315.0 * DEG2RAD
    zeta = 3.0 * DEG2RAD
    eta = PI

    production_fun = 2
    production_rate = 1.0e14
    ud_shape = 2
    umin = 720.0
    umax = 750.0
    ejection_angle_distr = 1
    sd = 2
    r = RM

    # rrM from spherical (r, alphaM, betaM)
    sa, ca = np.sin(alphaM), np.cos(alphaM)
    cb, sb = np.cos(betaM), np.sin(betaM)
    rrM = np.array([r * sa * cb, r * sa * sb, r * ca], dtype=float)
    symmetry_axis = rrM / float(np.linalg.norm(rrM))
    return Source(
        r=float(r),
        alphaM=float(alphaM),
        betaM=float(betaM),
        rrM=rrM,
        zeta=float(zeta),
        eta=float(eta),
        symmetry_axis=symmetry_axis,
        ejection_angle_distr=int(ejection_angle_distr),
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


# ----------------------------------------------------------------------
# Line of sight construction (line_of_sight)
# ----------------------------------------------------------------------


def line_of_sight(source: Source) -> np.ndarray[Any, np.dtype[object]]:
    """
    Python reimplementation of Fortran `line_of_sight` from image_construction.f90.

    Returns
    -------
    points : ndarray of Point, shape (2*ni+1, 2*Hpix+1, 2*Vpix+1)
        Python analogue of Fortran points(-ni:ni,-Hpix:Hpix,-Vpix:Vpix).
    """
    Hcamscale = Hrad / float(Hpix)
    Vcamscale = Vrad / float(Vpix)

    scpos = np.array(
        [RM * SQRT2 * 3.0, RM * SQRT2 * 3.0, RM * 1.12],
        dtype=float,
    )

    xcam = np.array([0.0, 0.0, 1.0], dtype=float)
    ycam = np.array([-1.0 / SQRT2, 1.0 / SQRT2, 0.0], dtype=float)
    zcam = np.array([-1.0 / SQRT2, -1.0 / SQRT2, 0.0], dtype=float)

    vtmp = scpos - source.rrM
    r0 = np.linalg.norm(vtmp)  # unused but kept for clarity
    _ = r0

    arcrm = RM / np.linalg.norm(scpos)
    mooncenterdir = np.array(
        [-scpos[2] / RM, 0.0, scpos[0] * SQRT2 / RM],
        dtype=float,
    )
    mooncenterdir /= np.linalg.norm(mooncenterdir)

    # Allocate array of Point objects
    nk = 2 * ni + 1
    nh = 2 * Hpix + 1
    nv = 2 * Vpix + 1
    points = np.empty((nk, nh, nv), dtype=object)

    for ii in range(-Vpix, Vpix + 1):
        # vertical pixel coordinate
        pixdircam_y = 0.0
        if ii < 0:
            pixdircam_y = (ii + 0.5) * Vcamscale
        elif ii > 0:
            pixdircam_y = (ii - 0.5) * Vcamscale

        for i in range(-Hpix, Hpix + 1):
            pixdircam_x = 0.0
            if i < 0:
                pixdircam_x = (i + 0.5) * Hcamscale
            elif i > 0:
                pixdircam_x = (i - 0.5) * Hcamscale

            pixdircam_z = np.sqrt(1.0 - pixdircam_x**2 - pixdircam_y**2)
            pixdircam = np.array([pixdircam_x, pixdircam_y, pixdircam_z], dtype=float)

            # angle between pixdircam and mooncenterdir
            arctmp = np.arccos(np.clip(np.dot(pixdircam, mooncenterdir), -1.0, 1.0))

            if arctmp >= arcrm and i != 0 and ii != 0:
                # Transform pixdircam to moon-centered CS
                vtmp = np.array([xcam[0], ycam[0], zcam[0]], dtype=float)
                pixdir_x = np.dot(pixdircam, vtmp)
                vtmp = np.array([xcam[1], ycam[1], zcam[1]], dtype=float)
                pixdir_y = np.dot(pixdircam, vtmp)
                vtmp = np.array([xcam[2], ycam[2], zcam[2]], dtype=float)
                pixdir_z = np.dot(pixdircam, vtmp)
                pixdir = np.array([pixdir_x, pixdir_y, pixdir_z], dtype=float)

                _, K1, _ = dist_between_2lines(
                    scpos, pixdir, source.rrM, source.symmetry_axis
                )

                for iii in range(-ni, ni + 1):
                    k = iii + ni
                    h = i + Hpix
                    v = ii + Vpix
                    rvec = K1 + pixdir * iii * sampdist_large
                    r = np.linalg.norm(rvec)
                    alpha = np.arccos(rvec[2] / r)
                    beta = np.arctan2(rvec[1], rvec[0])
                    compute = (r < 1.98e6) and (r > RM)
                    points[k, h, v] = Point(
                        r=r,
                        alpha=alpha,
                        beta=beta,
                        rvector=rvec,
                        compute=compute,
                    )
            else:
                for iii in range(-ni, ni + 1):
                    k = iii + ni
                    h = i + Hpix
                    v = ii + Vpix
                    points[k, h, v] = Point(
                        r=0.0,
                        alpha=0.0,
                        beta=0.0,
                        rvector=np.zeros(3, dtype=float),
                        compute=False,
                    )

    return points


# ----------------------------------------------------------------------
# Integral over line of sight (Integral_over_LoS)
# ----------------------------------------------------------------------


def integral_over_LoS(dens0: np.ndarray) -> np.ndarray:
    """
    Python version of image_construction.Integral_over_LoS.

    dens0 : shape (2*ni+1, 2*Hpix+1, 2*Vpix+1)
        Total number density along each line of sight.
    Returns
    -------
    image : shape (2*Hpix+1, 2*Vpix+1)
    """
    nk, nh, nv = dens0.shape
    assert nk == 2 * ni + 1
    assert nh == 2 * Hpix + 1
    assert nv == 2 * Vpix + 1

    image = np.zeros((nh, nv), dtype=float)

    for ii in range(-Vpix, Vpix + 1):
        for i in range(-Hpix, Hpix + 1):
            if i == 0 or ii == 0:
                continue
            h = i + Hpix
            v = ii + Vpix
            tmpTr = 0.0
            for k in range(-ni, ni):
                kk = k + ni
                kk1 = kk + 1
                tmpTr += dens0[kk, h, v] + dens0[kk1, h, v]
            image[h, v] = float(sampdist_large * 0.5 * tmpTr)

    return image


# ----------------------------------------------------------------------
# Main driver
# ----------------------------------------------------------------------


def main() -> None:
    source = get_volcano_source()

    # Precompute LOS geometry (does not depend on time)
    points = line_of_sight(source)

    nk, nh, nv = points.shape

    results_dir = "./results"
    os.makedirs(results_dir, exist_ok=True)

    # Moments of time (same as volcano_image.py)
    moments = [200 * (i + 1) for i in range(9)]

    for m_idx, tnow in enumerate(moments, start=1):
        print(f"Computing image {m_idx} at t = {tnow} s")

        dens0 = np.zeros((nk, nh, nv), dtype=float)

        # Collect all points where we need density, then call the Fortran
        # batch_over_points routine once to keep the innermost loop in Fortran.
        eval_points: List[Point] = []
        eval_indices: List[Tuple[int, int, int]] = []

        for ii in range(-Vpix, Vpix + 1):
            for i in range(-Hpix, Hpix + 1):
                h = i + Hpix
                v = ii + Vpix
                for k in range(-ni, ni + 1):
                    kk = k + ni
                    pt = points[kk, h, v]
                    if not pt.compute:
                        continue
                    eval_points.append(pt)
                    eval_indices.append((kk, h, v))

        if eval_points:
            dens_batch = batch_over_points(eval_points, source, float(tnow))
            # dens_batch shape: (N, 2) [bound, unbound]
            for (kk, h, v), dv in zip(eval_indices, dens_batch):
                dens0[kk, h, v] = float(dv[0] + dv[1])

        image = integral_over_LoS(dens0)

        # Add background where line of sight does not cross the moon's disc
        for ii in range(-Vpix, Vpix + 1):
            for i in range(-Hpix, Hpix + 1):
                h = i + Hpix
                v = ii + Vpix
                # Use k=0 slice (iii=0 in Fortran)
                pt_center = points[ni, h, v]
                if pt_center.r >= RM:
                    image[h, v] += BG_IMAGE

        # Write image to ./results/{m}.dat in the same layout as result_image_out
        fname = os.path.join(results_dir, f"{m_idx}.dat")
        # Build rows in the order ii=-Vpix..-1, then 1..Vpix, skipping ii=0,
        # and concatenating image(-Hpix:-1,ii), image(1:Hpix,ii).
        rows = []
        for ii in list(range(-Vpix, 0)) + list(range(1, Vpix + 1)):
            v = ii + Vpix
            left = image[(-Hpix) + Hpix : 0 + Hpix, v]  # indices -Hpix..-1
            right = image[1 + Hpix : Hpix + 1 + Hpix, v]  # indices 1..Hpix
            row = np.concatenate([left, right])
            rows.append(row)
        data = np.vstack(rows)
        np.savetxt(fname, data, fmt="%.6e")
        print(f"Wrote image data to {fname}")

    # ------------------------------------------------------------------
    # Plotting, similar to scripts/volcano_image.py
    # ------------------------------------------------------------------
    plt.gray()
    for i, tnow in enumerate(moments):
        dat_path = os.path.join(results_dir, f"{i + 1}.dat")
        d = np.loadtxt(dat_path, usecols=range(2 * Hpix))
        d = np.rot90(d)
        d = np.where(d > 1e-20, d, 1e-19)
        d = np.log10(d)
        plt.figure(i + 1)
        imgplot = plt.imshow(d, vmin=-17, vmax=-4)
        plt.colorbar()
        plt.text(75, 117, f"time = {tnow} s", color="w")
        imname = os.path.join(results_dir, f"volcano_{i + 1}.png")
        plt.savefig(imname)
        print(f"Saved plot to {imname}")


if __name__ == "__main__":
    main()
