Derivations of the formulae implemented in the DUDI (DUst DIstribution) model
are detailed in:
Ershova, A. & Schmidt, J. (2021). Two-body model for the spatial
distribution of dust ejected from an atmosphereless body. Astronomy and
Astrophysics, 650, A186.

Please, cite this paper if using this software in your work.

# DUDI Python Interface

## 1. What is this?

This directory provides the **Python interface**
to the **DUDI** model for dust number density and related quantities
around an atmosphereless spherical moon.
The original DUDI implementation is entirely in Fortran-95, and no
changes were made to the physical or numerical core for this Python release.

To support calling DUDI directly from Python, a small Fortran-2003
layer is added using:

- `bind(C)` interfaces,
- C-compatible derived types,
- a compiled shared library exposed to Python through `ctypes`.

This Python interface provides:

- A documented and object-oriented high-level API
  (`Point`, `Source`, convenience helpers for calling `DUDI` and
  `DUDI_mean_velocity`, and batching routines).
- Full access to all scientifically relevant features of DUDI without writing
  in Fortran.
- 1-to-1 Python equivalents of the core Fortran example programs
  (`enceladus_example`, `europa_example`, `io_example`, plus small
  diagnostic scripts).


## 2. Architecture Overview

Python (your scripts, notebooks)
        |
        v      high-level API + data models + batching utilities
  python_interface/dudi/api.py
  python_interface/dudi/models.py
        |
        v      thin ctypes wrapper (NumPy <-> raw C arrays)
  python_interface/dudi/ctypes_bridge.py
        |
        v      shared library exposing Fortran routines as C functions
  python_interface/dudi/libpy_dudi_bridge.so
        |
        v      Fortran-2003 wrapper module (bind(C))
  python_interface/fortran_bridge/py_bridge.f90 and helpers
        |
        v      Fortran-95 DUDI scientific core (unchanged)


Python never touches Fortran derived types directly. The C-bindable Fortran
wrappers (`py_bridge.f90`) receive only C-friendly scalars/arrays, reconstruct
the Fortran types (`position_in_space`, `source_properties`), call the
real routines in the `integrator` module (e.g. `DUDI`, `DUDI_mean_velocity`),
and return real(8) values to Python.

**Types & precision**
The DUDI routines compute densities as real(4). Our wrappers convert them
to real(8) (double) just before returning so Python gets a normal float.
Python vectors are passed as NumPy float64 contiguous arrays of shape (3,).
Fortran LOGICAL inputs at the boundary are passed as C integers (0/1).

All computational kernels and all **OpenMP parallelism** remain entirely in
**Fortran**, where they are most efficient. The Python layer is thin and
designed only for data preparation, convenience, and analysis.


## 3. Installation and Build (Linux)

The model has been tested and validated on **Linux** systems.
Source installation requires a working Fortran toolchain.

### 3.1. System prerequisites

Install:

- `gfortran` (Fortran 95/2003 compiler)
- OpenMP (`libgomp` – normally included with GCC)
- Python ≥ 3.9
- `pip`
- GitHub CLI or `git`

Example (Ubuntu-like):

sudo apt-get install gfortran libgomp1 python3 python3-pip git

### 3.2. Clone the repository

git clone https://github.com/Veyza/dudi.git
cd dudi

### 3.3 Set up the moon for modeling

Set the moon's mass (moon_mass [kg]) and radius (rm [meters])
in the module const.f90.

### 3.4. Build the Fortran ctypes bridge

bash python_interface/fortran_bridge/build_ctypes_bridge.sh

When finished, you should see:

python_interface/dudi/libpy_dudi_bridge.so

### 3.5. Install the Python package
Create an environment to avoid conflicts and system complaints:

python3 -m venv .venv
source .venv/bin/activate

Editable mode:

pip install -e .

### Tests

To run the Python interface tests locally (from the repo root):

1. Build the Fortran ctypes bridge (only needed once per build):

```bash
bash python_interface/fortran_bridge/build_ctypes_bridge.sh
```

2. (Recommended) Create and activate a virtual environment and install the package with test dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
```

3. Run the tests:

```bash
python3 -m pytest -q
```

The tests live under `tests/` and include:

- `test_density_scalar_vs_batch.py` – compares scalar `density` calls to `batch_over_points`.
- `test_mean_velocity_shapes.py` – checks `mean_velocity` and batched shapes/values.
- `test_examples_enceladus_smoke.py` – runs the `enceladus_example` script end-to-end.

## 4. Model input and numerical method

### 4.1 General model input
For the physical meaning of all input quantities, users should refer to
the main DUDI package `README`, in particular to the sections
“Specifying Parameters”, “Specifying Functions Describing the Dust Ejection”,
and “Utilizing the DUDI Subroutine”. These describe the required model
parameters, units, and conventions used by DUDI. The Python classes defined
in `python_interface/dudi/models.py` (`Point`, `Source`) are direct mirrors
of the Fortran derived types defined in the DUDI core (`position_in_space`,
`source_properties`) and follow the same data structure and validation logic.
Likewise, the high-level API functions exposed in `dudi.api` use the same
notation and argument structure as the corresponding Fortran routines in the
`integrator` module, enabling users to rely on the main README and the
original Fortran documentation when preparing input for the Python interface.

### 4.2 Distributions
The DUDI Python interface uses the parametric size, ejection-speed, and
ejection-direction distributions implemented in `src/distributions_fun.f90`.
Distribution families are selected via the integer flags `sd`,
`ud_shape`, and `ejection_angle_distr` in the `Source` object, matching
the Fortran code.

## 5. Python examples

The Python examples can be found in the folder `python_interface/examples`
and may serve as templates for your own applications. Apart from
`minimal.py`, which is a small Python-only demonstration, each script is a
direct translation of one of the Fortran examples in `examples/` and
reproduces the same physical setup and output structure.

### `minimal.py`

`minimal.py` is a compact, end-to-end sanity check and usage example. It
builds a single `Point` and `Source`, then calls the high-level wrappers
around `DUDI` and `DUDI_mean_velocity` once and prints the resulting
densities and average velocities to stdout. It is intended as a minimal
test that your Fortran bridge, Python installation, and basic calling
convention are working correctly.

### `enceladus_example.py`

`enceladus_example.py` is the Python analogue of
`examples/enceladus_example.f90`. It reads the jet and trajectory input
files (`input_data_files/Enceladus_jet.dat`,
`input_data_files/Cassini_E2_flyby.dat`), constructs the corresponding
array of `Source` and `Point` objects, and evaluates dust number density
along the Cassini E2 flyby using `DUDI`. The script writes the resulting
profile to the `results/` directory in the same format as the Fortran
example so that it can be plotted with `scripts/e2plot.py`.

### `europa_example.py`

`europa_example.py` is the Python analogue of `examples/europa_example.f90`.
It configures several sources and deposition points on the surface of
Europa (mirroring `get_europa_input` in `src/inputdata.f90`), calls
`DUDI` in batches over these locations, and writes the resulting mass
flux profiles to the `results/` directory. The output layout matches the
Fortran example and can be plotted with `scripts/deposition.py`.

### `io_example.py`

`io_example.py` is the Python analogue of `examples/io_example.f90`. It
constructs the grid of lines-of-sight for the synthetic volcano images
(as in `src/image_construction.f90`), evaluates the dust cross-section
at all sample points using `DUDI`, integrates along each line of sight,
and writes a sequence of 2D images to `results/`. These can be visualized
with `scripts/volcano_image.py`.

## 6. Maintenance notes

### 6.1 Rebuilding the Fortran bridge (release and debug modes)

The shared library `libpy_dudi_bridge.so` is generated by the build script:

`python_interface/fortran_bridge/build_ctypes_bridge.sh`

By default, the script compiles the Fortran core with optimization flags (-O2)
and OpenMP enabled.

To build the bridge in debug mode, which enables runtime checks, disables
optimizations, and may activate additional debugging output inside Fortran,
set the environment variable `DEBUG=1`:

`DEBUG=1 bash python_interface/fortran_bridge/build_ctypes_bridge.sh`

This produces a slower but more traceable version of the shared library, useful
when diagnosing issues such as NaNs, unexpected densities, mismatched array
shapes, or memory-related bugs.

After rebuilding the library, reinstall the Python package:
`pip install -e .`

### 6.2 Rebuilding or modifying Fortran sources

If **any** Fortran source file is changed, the shared library must be rebuilt.

### 6.3 Continuous Integration

The GitHub Actions workflows located in .github/workflows/ run tests on every
push or pull request.
A full CI run:

  1. Build the Fortran ctypes bridge,

  2. Install the Python package,

  3. Run pytest.

If CI is extended to include macOS or Windows, update the classifiers in
pyproject.toml accordingly and ensure the build script is adapted.
