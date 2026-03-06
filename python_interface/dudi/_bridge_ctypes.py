"""
ctypes bindings to the DUDI Fortran bridge (libpy_dudi_bridge.so).

Exposes: call_density, call_mean_velocity, call_batch_points, call_batch_sources,
call_batch_sources_points, and the mean_velocity batch variants.
All coordinates and distances in meters (moon-centered), angles in rad, time in s.
"""

from __future__ import annotations
import ctypes as C
from pathlib import Path
import numpy as np

_lib_path = Path(__file__).with_name("libpy_dudi_bridge.so")
if not _lib_path.exists():
    raise ImportError(
        f"Shared library not found: {_lib_path}. "
        "Run: bash python_interface/fortran_bridge/build_ctypes_bridge.sh"
    )
_lib = C.CDLL(str(_lib_path))

# ctypes array types
Vec3d = np.ctypeslib.ndpointer(dtype=np.float64, shape=(3,), flags="C_CONTIGUOUS")
Vec1d = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
Int1d = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")


def _as_vec3(a) -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"expected shape (3,), got {arr.shape}")
    return np.ascontiguousarray(arr)


def _as_1d_f64(a, name: str = "array") -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D, got shape {arr.shape}")
    return np.ascontiguousarray(arr)


def _as_2d_f64(a, name: str = "array") -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be 2D, got shape {arr.shape}")
    return np.ascontiguousarray(arr)


def _as_1d_i32(a, name: str = "array") -> np.ndarray:
    arr = np.asarray(a, dtype=np.int32)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D, got shape {arr.shape}")
    return np.ascontiguousarray(arr)


def _c_int(b: bool | int) -> int:
    return 1 if bool(b) else 0


# ==============  py_dudi_density  ==============
# void py_dudi_density(double point_r, point_alpha, point_beta,
#   double src_r, src_alphaM, src_betaM, src_zeta, src_eta,
#   int src_eject_distr, src_ud_shape, double src_umin, src_umax,
#   int src_sd, src_production_fun, double src_production_rate,
#   int src_is_jet, double tnow, double density_out[2])
_lib.py_dudi_density.argtypes = [
    C.c_double,
    C.c_double,
    C.c_double,  # point
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,  # source geom
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_double,  # ud
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_int,  # sd, production_fun, rate, is_jet
    C.c_double,  # tnow
    np.ctypeslib.ndpointer(dtype=np.float64, shape=(2,), flags="C_CONTIGUOUS"),
]
_lib.py_dudi_density.restype = None


def call_density(
    *,
    point_r: float,
    point_alpha: float,
    point_beta: float,
    src_r: float,
    src_alphaM: float,
    src_betaM: float,
    src_zeta: float,
    src_eta: float,
    src_eject_distr: int,
    src_ud_shape: int,
    src_umin: float,
    src_umax: float,
    src_sd: int,
    src_production_fun: int,
    src_production_rate: float,
    src_is_jet: bool | int,
    tnow: float,
) -> tuple[float, float]:
    """Single-point density. Returns (density_bound, density_unbound)."""
    out = np.empty(2, dtype=np.float64)
    _lib.py_dudi_density(
        C.c_double(point_r),
        C.c_double(point_alpha),
        C.c_double(point_beta),
        C.c_double(src_r),
        C.c_double(src_alphaM),
        C.c_double(src_betaM),
        C.c_double(src_zeta),
        C.c_double(src_eta),
        C.c_int(src_eject_distr),
        C.c_int(src_ud_shape),
        C.c_double(src_umin),
        C.c_double(src_umax),
        C.c_int(src_sd),
        C.c_int(src_production_fun),
        C.c_double(src_production_rate),
        C.c_int(_c_int(src_is_jet)),
        C.c_double(tnow),
        out,
    )
    return float(out[0]), float(out[1])


# ==============  py_dudi_mean_velocity  ==============
# same inputs + integral_out[8]
_lib.py_dudi_mean_velocity.argtypes = [
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_double,
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_int,
    C.c_double,
    np.ctypeslib.ndpointer(dtype=np.float64, shape=(8,), flags="C_CONTIGUOUS"),
]
_lib.py_dudi_mean_velocity.restype = None


def call_mean_velocity(
    *,
    point_r: float,
    point_alpha: float,
    point_beta: float,
    src_r: float,
    src_alphaM: float,
    src_betaM: float,
    src_zeta: float,
    src_eta: float,
    src_eject_distr: int,
    src_ud_shape: int,
    src_umin: float,
    src_umax: float,
    src_sd: int,
    src_production_fun: int,
    src_production_rate: float,
    src_is_jet: bool | int,
    tnow: float,
) -> np.ndarray:
    """Single-point mean velocity. Returns array of 8: v_up[3], v_down[3], n_up, n_down."""
    out = np.empty(8, dtype=np.float64)
    _lib.py_dudi_mean_velocity(
        C.c_double(point_r),
        C.c_double(point_alpha),
        C.c_double(point_beta),
        C.c_double(src_r),
        C.c_double(src_alphaM),
        C.c_double(src_betaM),
        C.c_double(src_zeta),
        C.c_double(src_eta),
        C.c_int(src_eject_distr),
        C.c_int(src_ud_shape),
        C.c_double(src_umin),
        C.c_double(src_umax),
        C.c_int(src_sd),
        C.c_int(src_production_fun),
        C.c_double(src_production_rate),
        C.c_int(_c_int(src_is_jet)),
        C.c_double(tnow),
        out,
    )
    return out


# ==============  py_dudi_batch_points  ==============
# (n_points, density_out[2*n_points], point_r[n], point_alpha[n], point_beta[n],
#  src scalars..., tnow)
_lib.py_dudi_batch_points.argtypes = [
    C.c_int,
    Vec1d,  # n_points, density_out
    Vec1d,
    Vec1d,
    Vec1d,  # point_r, point_alpha, point_beta
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_double,
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_int,
    C.c_double,  # tnow
]
_lib.py_dudi_batch_points.restype = None


def call_batch_points(
    *,
    point_r: np.ndarray,
    point_alpha: np.ndarray,
    point_beta: np.ndarray,
    src_r: float,
    src_alphaM: float,
    src_betaM: float,
    src_zeta: float,
    src_eta: float,
    src_eject_distr: int,
    src_ud_shape: int,
    src_umin: float,
    src_umax: float,
    src_sd: int,
    src_production_fun: int,
    src_production_rate: float,
    src_is_jet: bool | int,
    tnow: float,
) -> np.ndarray:
    """Batch over points, one source. Returns (n_points, 2): [bound, unbound] per point."""
    r = _as_1d_f64(point_r, "point_r")
    alpha = _as_1d_f64(point_alpha, "point_alpha")
    beta = _as_1d_f64(point_beta, "point_beta")
    if not (r.size == alpha.size == beta.size):
        raise ValueError("point_r, point_alpha, point_beta must have same length")
    n = int(r.size)
    density_out = np.empty(2 * n, dtype=np.float64)
    _lib.py_dudi_batch_points(
        n,
        density_out,
        r,
        alpha,
        beta,
        C.c_double(src_r),
        C.c_double(src_alphaM),
        C.c_double(src_betaM),
        C.c_double(src_zeta),
        C.c_double(src_eta),
        C.c_int(src_eject_distr),
        C.c_int(src_ud_shape),
        C.c_double(src_umin),
        C.c_double(src_umax),
        C.c_int(src_sd),
        C.c_int(src_production_fun),
        C.c_double(src_production_rate),
        C.c_int(_c_int(src_is_jet)),
        C.c_double(tnow),
    )
    return density_out.reshape(n, 2)


# ==============  py_dudi_batch_sources  ==============
_lib.py_dudi_batch_sources.argtypes = [
    C.c_int,
    Vec1d,  # n_sources, density_out
    C.c_double,
    C.c_double,
    C.c_double,  # point
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,  # src_r, alphaM, betaM, zeta, eta
    Int1d,
    Int1d,
    Vec1d,
    Vec1d,  # eject_distr, ud_shape, umin, umax
    Int1d,
    Int1d,
    Vec1d,
    Int1d,  # sd, production_fun, production_rate, is_jet
    C.c_double,
]
_lib.py_dudi_batch_sources.restype = None


def call_batch_sources(
    *,
    point_r: float,
    point_alpha: float,
    point_beta: float,
    src_r: np.ndarray,
    src_alphaM: np.ndarray,
    src_betaM: np.ndarray,
    src_zeta: np.ndarray,
    src_eta: np.ndarray,
    src_eject_distr: np.ndarray,
    src_ud_shape: np.ndarray,
    src_umin: np.ndarray,
    src_umax: np.ndarray,
    src_sd: np.ndarray,
    src_production_fun: np.ndarray,
    src_production_rate: np.ndarray,
    src_is_jet: np.ndarray,
    tnow: float,
) -> np.ndarray:
    """Batch over sources, one point. Returns (n_sources, 2): [bound, unbound] per source."""
    r_arr = _as_1d_f64(src_r, "src_r")
    n = int(r_arr.size)
    for name, arr in [
        ("src_alphaM", src_alphaM),
        ("src_betaM", src_betaM),
        ("src_zeta", src_zeta),
        ("src_eta", src_eta),
        ("src_umin", src_umin),
        ("src_umax", src_umax),
        ("src_production_rate", src_production_rate),
    ]:
        a = (
            _as_1d_f64(arr, name)
            if np.issubdtype(np.asarray(arr).dtype, np.floating)
            else arr
        )
        if a.size != n:
            raise ValueError(f"{name} length must be {n}")
    for name, arr in [
        ("src_eject_distr", src_eject_distr),
        ("src_ud_shape", src_ud_shape),
        ("src_sd", src_sd),
        ("src_production_fun", src_production_fun),
        ("src_is_jet", src_is_jet),
    ]:
        a = (
            _as_1d_i32(arr, name)
            if hasattr(arr, "size")
            else np.asarray(arr, dtype=np.int32)
        )
        if a.size != n:
            raise ValueError(f"{name} length must be {n}")

    alphaM = _as_1d_f64(src_alphaM, "src_alphaM")
    betaM = _as_1d_f64(src_betaM, "src_betaM")
    zeta = _as_1d_f64(src_zeta, "src_zeta")
    eta = _as_1d_f64(src_eta, "src_eta")
    umin = _as_1d_f64(src_umin, "src_umin")
    umax = _as_1d_f64(src_umax, "src_umax")
    rate = _as_1d_f64(src_production_rate, "src_production_rate")
    eject = _as_1d_i32(src_eject_distr, "src_eject_distr")
    udsh = _as_1d_i32(src_ud_shape, "src_ud_shape")
    sd = _as_1d_i32(src_sd, "src_sd")
    pfun = _as_1d_i32(src_production_fun, "src_production_fun")
    isjet = _as_1d_i32(np.asarray(src_is_jet, dtype=np.int32), "src_is_jet")

    density_out = np.empty(2 * n, dtype=np.float64)
    _lib.py_dudi_batch_sources(
        n,
        density_out,
        C.c_double(point_r),
        C.c_double(point_alpha),
        C.c_double(point_beta),
        r_arr,
        alphaM,
        betaM,
        zeta,
        eta,
        eject,
        udsh,
        umin,
        umax,
        sd,
        pfun,
        rate,
        isjet,
        C.c_double(tnow),
    )
    return density_out.reshape(n, 2)


# ==============  py_dudi_batch_sources_points  ==============
_lib.py_dudi_batch_sources_points.argtypes = [
    C.c_int,
    C.c_int,
    Vec1d,  # n_points, n_sources, density_out
    Vec1d,
    Vec1d,
    Vec1d,  # point_r, point_alpha, point_beta
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Int1d,
    Int1d,
    Vec1d,
    Vec1d,
    Int1d,
    Int1d,
    Vec1d,
    Int1d,
    C.c_double,
]
_lib.py_dudi_batch_sources_points.restype = None


def call_batch_sources_points(
    *,
    point_r: np.ndarray,
    point_alpha: np.ndarray,
    point_beta: np.ndarray,
    src_r: np.ndarray,
    src_alphaM: np.ndarray,
    src_betaM: np.ndarray,
    src_zeta: np.ndarray,
    src_eta: np.ndarray,
    src_eject_distr: np.ndarray,
    src_ud_shape: np.ndarray,
    src_umin: np.ndarray,
    src_umax: np.ndarray,
    src_sd: np.ndarray,
    src_production_fun: np.ndarray,
    src_production_rate: np.ndarray,
    src_is_jet: np.ndarray,
    tnow: float,
) -> np.ndarray:
    """Batch over points and sources; density summed over sources. Returns (n_points, 2)."""
    r_pt = _as_1d_f64(point_r, "point_r")
    n_points = int(r_pt.size)
    alpha_pt = _as_1d_f64(point_alpha, "point_alpha")
    beta_pt = _as_1d_f64(point_beta, "point_beta")
    if alpha_pt.size != n_points or beta_pt.size != n_points:
        raise ValueError("point arrays must have same length")

    r_src = _as_1d_f64(src_r, "src_r")
    n_sources = int(r_src.size)
    alphaM = _as_1d_f64(src_alphaM, "src_alphaM")
    betaM = _as_1d_f64(src_betaM, "src_betaM")
    zeta = _as_1d_f64(src_zeta, "src_zeta")
    eta = _as_1d_f64(src_eta, "src_eta")
    umin = _as_1d_f64(src_umin, "src_umin")
    umax = _as_1d_f64(src_umax, "src_umax")
    rate = _as_1d_f64(src_production_rate, "src_production_rate")
    eject = _as_1d_i32(src_eject_distr, "src_eject_distr")
    udsh = _as_1d_i32(src_ud_shape, "src_ud_shape")
    sd = _as_1d_i32(src_sd, "src_sd")
    pfun = _as_1d_i32(src_production_fun, "src_production_fun")
    isjet = _as_1d_i32(np.asarray(src_is_jet, dtype=np.int32), "src_is_jet")
    for arr in [
        alphaM,
        betaM,
        zeta,
        eta,
        umin,
        umax,
        rate,
        eject,
        udsh,
        sd,
        pfun,
        isjet,
    ]:
        if arr.size != n_sources:
            raise ValueError(f"source array length must be n_sources={n_sources}")

    density_out = np.empty(2 * n_points, dtype=np.float64)
    _lib.py_dudi_batch_sources_points(
        n_points,
        n_sources,
        density_out,
        r_pt,
        alpha_pt,
        beta_pt,
        r_src,
        alphaM,
        betaM,
        zeta,
        eta,
        eject,
        udsh,
        umin,
        umax,
        sd,
        pfun,
        rate,
        isjet,
        C.c_double(tnow),
    )
    return density_out.reshape(n_points, 2)


# ==============  py_dudi_mean_velocity_batch_points  ==============
_lib.py_dudi_mean_velocity_batch_points.argtypes = [
    C.c_int,
    Vec1d,  # n_points, integral_out[8*n_points]
    Vec1d,
    Vec1d,
    Vec1d,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_double,
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_double,
    C.c_int,
    C.c_int,
    C.c_double,
    C.c_int,
    C.c_double,
]
_lib.py_dudi_mean_velocity_batch_points.restype = None


def call_mean_velocity_batch_points(
    *,
    point_r: np.ndarray,
    point_alpha: np.ndarray,
    point_beta: np.ndarray,
    src_r: float,
    src_alphaM: float,
    src_betaM: float,
    src_zeta: float,
    src_eta: float,
    src_eject_distr: int,
    src_ud_shape: int,
    src_umin: float,
    src_umax: float,
    src_sd: int,
    src_production_fun: int,
    src_production_rate: float,
    src_is_jet: bool | int,
    tnow: float,
) -> np.ndarray:
    """Mean-velocity batch over points. Returns (n_points, 8): [v_up(3), v_down(3), n_up, n_down] per point."""
    r = _as_1d_f64(point_r, "point_r")
    n = int(r.size)
    alpha = _as_1d_f64(point_alpha, "point_alpha")
    beta = _as_1d_f64(point_beta, "point_beta")
    if alpha.size != n or beta.size != n:
        raise ValueError("point arrays must have same length")
    integral_out = np.empty(8 * n, dtype=np.float64)
    _lib.py_dudi_mean_velocity_batch_points(
        n,
        integral_out,
        r,
        alpha,
        beta,
        C.c_double(src_r),
        C.c_double(src_alphaM),
        C.c_double(src_betaM),
        C.c_double(src_zeta),
        C.c_double(src_eta),
        C.c_int(src_eject_distr),
        C.c_int(src_ud_shape),
        C.c_double(src_umin),
        C.c_double(src_umax),
        C.c_int(src_sd),
        C.c_int(src_production_fun),
        C.c_double(src_production_rate),
        C.c_int(_c_int(src_is_jet)),
        C.c_double(tnow),
    )
    return integral_out.reshape(n, 8)


# ==============  py_dudi_mean_velocity_batch_sources  ==============
_lib.py_dudi_mean_velocity_batch_sources.argtypes = [
    C.c_int,
    Vec1d,  # n_sources, integral_out[8*n_sources]
    C.c_double,
    C.c_double,
    C.c_double,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Int1d,
    Int1d,
    Vec1d,
    Vec1d,
    Int1d,
    Int1d,
    Vec1d,
    Int1d,
    C.c_double,
]
_lib.py_dudi_mean_velocity_batch_sources.restype = None


def call_mean_velocity_batch_sources(
    *,
    point_r: float,
    point_alpha: float,
    point_beta: float,
    src_r: np.ndarray,
    src_alphaM: np.ndarray,
    src_betaM: np.ndarray,
    src_zeta: np.ndarray,
    src_eta: np.ndarray,
    src_eject_distr: np.ndarray,
    src_ud_shape: np.ndarray,
    src_umin: np.ndarray,
    src_umax: np.ndarray,
    src_sd: np.ndarray,
    src_production_fun: np.ndarray,
    src_production_rate: np.ndarray,
    src_is_jet: np.ndarray,
    tnow: float,
) -> np.ndarray:
    """Mean-velocity batch over sources. Returns (n_sources, 8)."""
    r_src = _as_1d_f64(src_r, "src_r")
    n = int(r_src.size)
    alphaM = _as_1d_f64(src_alphaM, "src_alphaM")
    betaM = _as_1d_f64(src_betaM, "src_betaM")
    zeta = _as_1d_f64(src_zeta, "src_zeta")
    eta = _as_1d_f64(src_eta, "src_eta")
    umin = _as_1d_f64(src_umin, "src_umin")
    umax = _as_1d_f64(src_umax, "src_umax")
    rate = _as_1d_f64(src_production_rate, "src_production_rate")
    eject = _as_1d_i32(src_eject_distr, "src_eject_distr")
    udsh = _as_1d_i32(src_ud_shape, "src_ud_shape")
    sd = _as_1d_i32(src_sd, "src_sd")
    pfun = _as_1d_i32(src_production_fun, "src_production_fun")
    isjet = _as_1d_i32(np.asarray(src_is_jet, dtype=np.int32), "src_is_jet")
    for arr in [
        alphaM,
        betaM,
        zeta,
        eta,
        umin,
        umax,
        rate,
        eject,
        udsh,
        sd,
        pfun,
        isjet,
    ]:
        if arr.size != n:
            raise ValueError(f"source array length must be {n}")

    integral_out = np.empty(8 * n, dtype=np.float64)
    _lib.py_dudi_mean_velocity_batch_sources(
        n,
        integral_out,
        C.c_double(point_r),
        C.c_double(point_alpha),
        C.c_double(point_beta),
        r_src,
        alphaM,
        betaM,
        zeta,
        eta,
        eject,
        udsh,
        umin,
        umax,
        sd,
        pfun,
        rate,
        isjet,
        C.c_double(tnow),
    )
    return integral_out.reshape(n, 8)


# ==============  py_dudi_mean_velocity_batch_sources_points  ==============
_lib.py_dudi_mean_velocity_batch_sources_points.argtypes = [
    C.c_int,
    C.c_int,
    Vec1d,  # n_points, n_sources, integral_out[8*n_points]
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Vec1d,
    Int1d,
    Int1d,
    Vec1d,
    Vec1d,
    Int1d,
    Int1d,
    Vec1d,
    Int1d,
    C.c_double,
]
_lib.py_dudi_mean_velocity_batch_sources_points.restype = None


def call_mean_velocity_batch_sources_points(
    *,
    point_r: np.ndarray,
    point_alpha: np.ndarray,
    point_beta: np.ndarray,
    src_r: np.ndarray,
    src_alphaM: np.ndarray,
    src_betaM: np.ndarray,
    src_zeta: np.ndarray,
    src_eta: np.ndarray,
    src_eject_distr: np.ndarray,
    src_ud_shape: np.ndarray,
    src_umin: np.ndarray,
    src_umax: np.ndarray,
    src_sd: np.ndarray,
    src_production_fun: np.ndarray,
    src_production_rate: np.ndarray,
    src_is_jet: np.ndarray,
    tnow: float,
) -> np.ndarray:
    """Mean-velocity batch over points and sources (summed). Returns (n_points, 8)."""
    r_pt = _as_1d_f64(point_r, "point_r")
    n_points = int(r_pt.size)
    alpha_pt = _as_1d_f64(point_alpha, "point_alpha")
    beta_pt = _as_1d_f64(point_beta, "point_beta")
    if alpha_pt.size != n_points or beta_pt.size != n_points:
        raise ValueError("point arrays must have same length")

    r_src = _as_1d_f64(src_r, "src_r")
    n_sources = int(r_src.size)
    alphaM = _as_1d_f64(src_alphaM, "src_alphaM")
    betaM = _as_1d_f64(src_betaM, "src_betaM")
    zeta = _as_1d_f64(src_zeta, "src_zeta")
    eta = _as_1d_f64(src_eta, "src_eta")
    umin = _as_1d_f64(src_umin, "src_umin")
    umax = _as_1d_f64(src_umax, "src_umax")
    rate = _as_1d_f64(src_production_rate, "src_production_rate")
    eject = _as_1d_i32(src_eject_distr, "src_eject_distr")
    udsh = _as_1d_i32(src_ud_shape, "src_ud_shape")
    sd = _as_1d_i32(src_sd, "src_sd")
    pfun = _as_1d_i32(src_production_fun, "src_production_fun")
    isjet = _as_1d_i32(np.asarray(src_is_jet, dtype=np.int32), "src_is_jet")
    for arr in [
        alphaM,
        betaM,
        zeta,
        eta,
        umin,
        umax,
        rate,
        eject,
        udsh,
        sd,
        pfun,
        isjet,
    ]:
        if arr.size != n_sources:
            raise ValueError(f"source array length must be n_sources={n_sources}")

    integral_out = np.empty(8 * n_points, dtype=np.float64)
    _lib.py_dudi_mean_velocity_batch_sources_points(
        n_points,
        n_sources,
        integral_out,
        r_pt,
        alpha_pt,
        beta_pt,
        r_src,
        alphaM,
        betaM,
        zeta,
        eta,
        eject,
        udsh,
        umin,
        umax,
        sd,
        pfun,
        rate,
        isjet,
        C.c_double(tnow),
    )
    return integral_out.reshape(n_points, 8)
