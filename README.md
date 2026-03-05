## DUDI – Specialized Workflow Branch

This branch contains an experimental set of workflows for the DUDI model.
It is **not** intended to be merged into `main`. This `README.md` documents the branch‑specific executables and how to run them.

The three primary workflows are:

- **`flyby_profiles`** – dust number–density profiles along Cassini flybys.
- **`plume_vertical_structure`** – vertical slices through the plume (dust mass, gas, composition).
- **`horizontal_structure`** – horizontal slices of dust mass density above the south pole.

All paths below are relative to the repository root.

---

## Build, dependencies, and global parameters

### Build

- **Compiler**: `gfortran` (or compatible) with OpenMP.
- **Build tool**: GNU `make`.

```bash
make flyby_profile          # builds ./bin/flyby_profile
make vertical_structure     # builds ./bin/vertical_structure
make horizontal_structure   # builds ./bin/horizontal_structure
```

### R dependencies (for plotting)

R is optional but recommended for plotting results.

Non‑standard R packages used in this branch:

- **`akima`** – interpolation for horizontal slices.
- **`fields`** – `image.plot`, color bars, etc.
- **`RColorBrewer`** – discrete colour palettes.
- **`latex2exp`** – LaTeX‑style axis labels for some plots.

You can install them (from the repo root) via:

```r
source("scripts/R-packages_list.R")
```

or manually:

```r
install.packages(c("akima", "fields", "RColorBrewer", "latex2exp"))
```

### Global parameter `p` in `const.f90`

The module `const` in `src/const.f90` defines the integer parameter `p`:

```fortran
! p = 0 -- number density is computed, 1 -- mean radius,
! 2 -- cross section, 3 -- mass density
integer, parameter :: p = 0
```

Different workflows require **specific values of `p`** before you compile and run:

- **Number‑density workflows**:
  - `flyby_profiles`
  - `plume_vertical_structure` for **gas** and **composition**
  - **Required**: `p = 0`

- **Mass‑density workflows**:
  - `plume_vertical_structure` for **dust mass**
  - `horizontal_structure`
  - **Required**: `p = 3`

If `p` is wrong, the programs will stop with a clear error (e.g. “wrong p parameter, in const.f90 must be p = 3”).

Typical workflow:

1. Edit `src/const.f90` and set `p` to `0` or `3` as required.
2. Rebuild the relevant executable(s) with `make`.
3. Run the workflow.

---

## Workflow 1: `flyby_profiles`

### Overview

`flyby_profiles` computes **dust number density** along Cassini flyby trajectories for a variety of Enceladus encounters.

- **Source**: `examples/flyby_profile.f90`
- **Binary**: `./bin/flyby_profile`
- **Outputs**: `./results/E<code>_profile.dat`
- **Plotting**: `scripts/flyby_plots.R` (uses `latex2exp`)

### Required settings

- **Set** `p = 0` in `src/const.f90` (number density).
- Build:

```bash
make flyby_profile
```

### Arguments (`fnum`)

The program takes one numeric argument `fnum` that selects a particular flyby/profile definition. It is interpreted by `define_flyby_params` in `src/inputdata.f90`, which:

- selects the Cassini SPICE‑derived trajectory file,
- sets dust size‑distribution integration limits (`rmin`, `rmax`, etc.),
- sets a scaling factor `varfact`.

Examples of supported `fnum` values:

- **E5** cases: `5`, `5.2`, `5.01`, `5.02`, `5.21`, `5.22`
- **E17** cases: `17`, `17.0`, `17.17`, `17.27`, `17.37`, `17.2`, `17.3`
- **E7** cases: `7.1`, `7.2`, `7.3`
- **E21** cases: `21.1`, `21.2`, `21.3`
- Other profiles: `3.1`, `3.2`, `3.3`, `4.1`, `4.2`, etc.

### Running

```bash
# Example: E5 profile
./bin/flyby_profile 5

# Example: E17 profile
./bin/flyby_profile 17
```

A convenience target runs several standard profiles:

```bash
make run-flyby
```

This creates profile files in `./results/` such as `E5_profile.dat`.
Use `scripts/flyby_plots.R` to visualize the profiles.

---

## Workflow 2: `plume_vertical_structure`

### Overview

`plume_vertical_structure` computes **vertical slices through the plume** in a fixed plane, for:

- **dust mass density** (dust plume),
- **gas number density** (gas plume),
- **composition** (salt‑rich fraction in dust).

- **Source**: `examples/plume_vert_slice.f90`
- **Binary**: `./bin/vertical_structure`
- **Outputs**: several matrix files in `./results/`
- **Plotting**: `scripts/plot_vertical_structure.R`

### Command‑line argument (mode)

The program takes a single argument selecting the quantity:

```bash
./bin/vertical_structure mass   # dust mass distribution
./bin/vertical_structure gas    # gas distribution
./bin/vertical_structure comp   # composition (salt-rich fraction)
```

Internally this maps to a numeric `qid`:

- `"mass"` → `qid = 0.1`
- `"gas"`  → `qid = 0.2`
- `"comp"` → `qid = 0.4`

Numeric arguments `0.1`, `0.2`, `0.4` are also accepted for backward compatibility.

### Required `p` value

- **Dust mass (`mass` / `0.1`)**:
  - **Set `p = 3`**.
- **Gas (`gas` / `0.2`)**:
  - **Set `p = 0`**.
- **Composition (`comp` / `0.4`)**:
  - **Set `p = 0`**.

The program checks `p` for each mode and stops if it is inconsistent.

### Spatial resolution: `cellsize` and `nt`

At the top of `examples/plume_vert_slice.f90`:

```fortran
integer, parameter :: nt      = 300        ! grid size (nt x nt)
real(8), parameter :: cellsize = 3d3       ! cell size in meters
```

- **`nt`** controls the number of grid points in each direction within the vertical plane.
- **`cellsize`** sets the physical grid spacing in meters; it is passed to:
  - `get_flyby_plane` in `src/inputdata.f90` (plane geometry),
  - `vertical_slicematrix_out` and `composition_matrix_out` in `src/dataoutmod.f90` (file naming).

You can change `cellsize` in `plume_vert_slice.f90` and recompile to get coarser/finer vertical structure.

### Output filenames and `_a-panel` / `_b-panel`

In `src/dataoutmod.f90`:

- **Vertical density slices** (`vertical_slicematrix_out`):

  - Dust (`qid = 0.1`):

    - if `cellsize < 2.9d3`:
      `./results/dens_in_plane_dust_a-panel.dat`
    - else:
      `./results/dens_in_plane_dust_b-panel.dat`

  - Gas (`qid = 0.2`):

    - if `cellsize < 2.9d3`:
      `./results/dens_in_plane_gas_a-panel.dat`
    - else:
      `./results/dens_in_plane_gas_b-panel.dat`

- **Composition matrices** (`composition_matrix_out`, `qid = 0.4`):

  - if `cellsize < 2.9d3`:

    - `./results/salt_poor_plane_a-panel.dat`
    - `./results/salt_rich_plane_a-panel.dat`

  - else:

    - `./results/salt_poor_plane_b-panel.dat`
    - `./results/salt_rich_plane_b-panel.dat`

The R plotting script (`scripts/plot_vertical_structure.R`) expects these filenames for “a‑panel” (close) and “b‑panel” (far) plots.

### Running examples

```bash
# Dust mass vertical slice (near plane, requires p = 3)
# 1) Set p = 3 in src/const.f90
make vertical_structure
# 2) Run:
./bin/vertical_structure mass

# Gas vertical slice (requires p = 0)
# 1) Set p = 0 in src/const.f90
make vertical_structure
# 2) Run:
./bin/vertical_structure gas

# Composition slice (requires p = 0)
./bin/vertical_structure comp
```

Then, from R:

```r
source("scripts/plot_vertical_structure.R")
# call the plotting functions defined there to generate PNGs
```

---

## Workflow 3: `horizontal_structure`

### Overview

`horizontal_structure` computes **horizontal slices** of **dust mass density** above Enceladus’ south pole at specified altitudes.

- **Source**: `examples/plume_horizontal_structure.f90`
- **Binary**: `./bin/horizontal_structure`
- **Outputs**: `./results/surface_dust_mass_density_<AltKm>km.dat`
- **Plotting**: `scripts/plot_horizontal_structure.R` → `./results/horizontal_slices.png`

### Required `p` value

This workflow uses the dust mass normalization:

- **Set `p = 3` in `src/const.f90`**.

If `p` is not `3`, the program stops with an error (“const.f90 must have p = 3 for mass density.”).

### Geometry and surface grid

Key parameters in `plume_horizontal_structure.f90`:

```fortran
integer, parameter :: Njets      = 100
integer, parameter :: Ndsources  = 160

real(8), parameter :: griddist  = 3d3      ! spacing along surface [m]
real(8), parameter :: maxalphaM = pi       ! 180 deg -> latitude -90 deg
real(8), parameter :: minalphaM = 155d0 * deg2rad   ! latitude -65 deg
```

- The code constructs a quasi‑uniform set of surface points, covering **latitudes from −65° down to −90°** around the south pole.
- `griddist = 3 km` controls the approximate spacing between neighbouring surface points.

The grid is generated by `get_surface_points` in `src/inputdata.f90`, using the specified altitude.

### Command‑line argument (altitude)

`horizontal_structure` expects **one argument: altitude in meters**:

```bash
./bin/horizontal_structure 25000     # 25 km
./bin/horizontal_structure 50000     # 50 km
./bin/horizontal_structure 75000     # 75 km
./bin/horizontal_structure 100000    # 100 km
```

For an altitude `A` meters, the output file is:

```text
./results/surface_dust_mass_density_<A_km>km.dat
```

where `<A_km> = nint(A / 1000)`, e.g.:

- `surface_dust_mass_density_25km.dat`
- `surface_dust_mass_density_50km.dat`, etc.

Each line in the file has:

```text
latitude_deg   longitude_deg   dust_mass_density_kg_per_m3
```

### Source setup

`plume_horizontal_structure.f90` reuses the same **dust mass** source model as the `"mass"` mode of `plume_vertical_structure`:

- Salt‑poor jets (`sd = 1`) from `vertical_jets.dat`.
- Salt‑rich dust from jets (`sd = 3`) and diffuse sources from `diffuse_sources.dat`.

It accumulates mass density from:

- all jets,
- all diffuse sources,

and then applies the same global scaling factor `varfact`.

### Plotting horizontal slices

From R:

```r
source("scripts/plot_horizontal_structure.R")

# default altitudes: 25, 50, 75, 100 km
plot_horizontal_structure()
```

This reads:

- `surface_dust_mass_density_25km.dat`
- `surface_dust_mass_density_50km.dat`
- `surface_dust_mass_density_75km.dat`
- `surface_dust_mass_density_100km.dat`

converts density from **kg/m³** to **g/m³**, and generates `./results/horizontal_slices.png`, a 2×2 panel figure with:

- smooth filled‑contour backgrounds,
- bright markers for jets and diffuse sources,
- annotated “Saturn” and “Orbit” arrows,
- a color bar (g/m³) for each panel.

---

## Quick reference: required `p` per workflow

- **`flyby_profiles`** (`bin/flyby_profile`):
  - **Set** `p = 0`
  - Argument: numeric `fnum` (flyby/profile code).
- **`plume_vertical_structure`** (`bin/vertical_structure`):
  - `"mass"` / `0.1` (dust mass): **`p = 3`**
  - `"gas"` / `0.2` (gas number density): **`p = 0`**
  - `"comp"` / `0.4` (composition): **`p = 0`**
- **`horizontal_structure`** (`bin/horizontal_structure`):
  - **Set** `p = 3`
  - Argument: altitude in meters.

This README is specific to this branch and documents only the workflows and conventions that differ from the main branch.
