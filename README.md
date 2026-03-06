## DUDI – Specialized Workflow Branch

This branch contains the setup of the DUDI package used for modeling the dust plume of Enceladus. The model implementation follows the approach described in
Ershova, A., Schmidt, J., Postberg, F., Khawaja, N., Nölle, L., Srama, R., Kempf, S., & Southworth, B. (2024). Modeling the Enceladus dust plume based on in situ measurements performed with the Cassini Cosmic Dust Analyzer. Astronomy & Astrophysics, 689, A114. https://doi.org/10.1051/0004-6361/202450429

The input data include lists of dust sources with their corresponding physical parameters and the ephemerides of the Cassini spacecraft, which are used to simulate Cassini flybys of Enceladus and compute the resulting dust environment along the spacecraft trajectory, as well as the precomputed file "background_size_distribution.dat" used to take into account the E ring dust.

The file "vertical_jets.dat" contains latitudes and Eastern longitudes of the jets from Porco et al, 2014. However, the jets' zenith angles and azimuts are set to zeros. The file "diffuse_sources.dat" has the same structure. The files with Cassini ephemeridae have 4 columns: seconds from the moment of the closest approach to Enceladus, radial distance to the moon center, latitude and Eastern longitude.

The three primary workflows and their output:

- **`flyby_profiles`** – dust number–density profiles along Cassini flybys used to obtain the plots in Figs. 12 - 15 of Ershova et al., 2024.
The output of this workflow is a 5-column table: seconds from the closest approach, number density of dust particles within the detectable size range of the Cassini instrument (HRD or CA), proportion of type1, type2, and type3 dust. In the case of E7 and E21 flybys type1 and type2 proportions do not matter as separate quantities but should be summed to obtain the proportion of salt-poor dust.

- **`plume_vertical_structure`** – vertical slice of the plume (dust mass, gas, composition) shown in Figs. 19-20 of Ershova et al., 2024.
This workflow outputs a text file with a matrix (or 2 matrices) 300x300 numbers corresponding to the computed quantity in the equidistant grid nodes. Depending on the keywords, the output can be gas density, dust mass density (salt-rich and salt-poor dust together), or dust number density - separately salt-rich and salt-poor particles. The dust mass of number density are computed for the particles with sizes between 0.1 um and 15 um.

- **`horizontal_structure`** – dust mass density distribution at the fixed altitude above the south pole.
This workflow output has 3 columns: latitude, Eastern longitude, and the dust mass density at the given point (altitude above Enceladus surface is fixed and given in the name of the output file).

**Units.** Number density in the output files is always in m^{-3}, mass density is in kg/m^3, angles are written in degrees.

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

### The branch contains also R scripts to plot the results. The plots color scheme differs from the ones in Ershova et al., 2024. You can install the packages needed for plotting using the script.

```r
source("scripts/R-packages_list.R")
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

If `p` is wrong, the programs will stop with a warning.

Typical workflow:

1. Edit `src/const.f90` and set `p` to `0` or `3` as required.
2. Rebuild the relevant executable(s) with `make`.
3. Run the workflow.

---

## Workflow 1: `flyby_profiles`

`flyby_profiles` computes **dust number density** along Cassini flyby trajectories for a variety of Enceladus encounters.

- **Source**: `examples/flyby_profile.f90`
- **Binary**: `./bin/flyby_profile`
- **Outputs**: `./results/EN_profile.dat`, N - code related to the flyby number
- **Plotting in R**:
HRDplot3in1(7) or HRDplot3in1(21) to plot model profiles from Figs. 12 and 13 respectively.
flybyplot(5) or flybyplot(17) to plot model profiles as in Figs. 14 and 15 respectively.
The repository currently does not contain the CDA measurer profiles shown in the plots in Ershova et al., 2024 along with the model fits.

### Required settings

- **Set** `p = 0` in `src/const.f90` (number density).
### Compile:

```bash
make flyby_profile
```

### Arguments (`fnum`)

The program takes one numeric argument `fnum` that selects a particular flyby/profile definition.

Examples of supported `fnum` values:
- **E7** cases: `7.1`, `7.2`, `7.3`        - profiles of the E7 flyby with the three size thresholds of the HRD measurements.
- **E21** cases: `21.1`, `21.2`, `21.3`    - profiles of the E21 flyby with the three size thresholds of the HRD measurements.

**A typo was found in the code originally used to produce the plots in Figs. 14 and 15 of Ershova et al., 2024.** The mistake is related to the definition of the parameters R₁ and R₂, and the corresponding coefficients a and b (see Sect. 4.4 of Ershova et al., 2024). The error has only a minor effect on the fits, and these parameters do not affect any of the conclusions of the paper. Therefore, we did not publish an erratum. Instead, we provide the corrected code for the compositional profiles in this branch to avoid confusion. **The correct values of the parameters R₁ and R₂ must be 0.9 um and 1.1 um for the E17 flyby and 0.3 um and 0.5 um for the E5.**
- **E5** cases: `5`, `5.2`
- **E17** cases: `17`, `17.17`

There is a target in `makefile` to run all the flyby profiles (runs for several seconds).

```bash
make run-flyby
```
---

## Workflow 2: `plume_vertical_structure`
**This workflow includes a possibility to calculate the gas density. DUDI is a dust code and the approach to calculations of gas density differs from DUDI's standard algorithm. The author of the code currently suspects a mistake in the gas density calculations and does not recommend using the gas density distribution obtained with this model until the matter is properly investigated.**

### Compile:
```bash
make vertical_structure
```

- **dust mass density** (dust plume),
- **gas number density** (gas plume),
- **composition** (salt‑rich fraction in dust).

The program takes a single argument selecting the quantity:

```bash
./bin/vertical_structure "mass"   # dust mass distribution
./bin/vertical_structure "gas"    # gas distribution
./bin/vertical_structure "comp"   # composition (salt-rich fraction)
```

### Required `p` value

- **Dust mass (`mass`)**:
  - **Set `p = 3`**.
- **Gas (`gas`)**:
  - **Set `p = 0`**.
- **Composition (`comp`)**:
  - **Set `p = 0`**.

### Spatial resolution: `cellsize` and `nt`

At the top of `examples/plume_vert_slice.f90`:

```fortran
integer, parameter :: nt      = 300        ! grid size (nt x nt)
real(8), parameter :: cellsize = 3d3       ! cell size in meters, for the panel (a) plot `cellsize = 1d3`
```

- **`nt`** controls the number of grid points in each direction within the vertical plane.
- **`cellsize`** sets the physical grid spacing in meters; it is passed to:
  - `get_flyby_plane` in `src/inputdata.f90` (plane geometry),
  - `vertical_slicematrix_out` and `composition_matrix_out` in `src/dataoutmod.f90` (file naming).

You can change `cellsize` in `plume_vert_slice.f90` and recompile to get coarser/finer vertical structure.
The program always computes a square matrix of values though in the panel (a) plot only a sub-region of it is shown.

### Output filenames and `_a-panel` / `_b-panel`

In `src/dataoutmod.f90`:

- **Density distributions (mass)** (`vertical_slicematrix_out`):

  - Dust:

    - if `cellsize < 2.9d3`:
      `./results/dens_in_plane_dust_a-panel.dat`
    - else:
      `./results/dens_in_plane_dust_b-panel.dat`

  - Gas:

    - if `cellsize < 2.9d3`:
      `./results/dens_in_plane_gas_a-panel.dat`
    - else:
      `./results/dens_in_plane_gas_b-panel.dat`

- **Number density distributions** (`composition_matrix_out`, `qid = 0.4`):

  - if `cellsize < 2.9d3`:

    - `./results/salt_poor_plane_a-panel.dat`
    - `./results/salt_rich_plane_a-panel.dat`

  - else:

    - `./results/salt_poor_plane_b-panel.dat`
    - `./results/salt_rich_plane_b-panel.dat`

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
Then plot in R:
```r
source("scripts/plot_vertical_structure.R")
composition_in_plane(close=T)     # Fig. 20, panel (a)
composition_in_plane()            # Fig. 20, panel (b)
dust2gas_in_plane(close = T)      # Fig. 19, panel (a)
dust2gas_in_plane()               # Fig. 19, panel (b)
density_in_plane("gas", T)        # plots gass density like in a-panel
density_in_plane("gas")           # plots gass density like in b-panel
density_in_plane("mass", T)       # plots dust mass density like in a-panel
density_in_plane("mass")          # plots dust mass density like in b-panel
density_in_plane("number", T)     # plots dust number density like in a-panel
density_in_plane("number")        # plots dust number density like in b-panel
```
---

## Workflow 3: `horizontal_structure`

### Compile:
```bash
make horizontal_structure
```
```
- The code constructs a quasi‑uniform set of surface points, covering **latitudes from −65° down to −90°** around the south pole.
- `griddist = 3 ` [km] controls the approximate spacing between neighbouring surface points.

### Command‑line argument (altitude)

`horizontal_structure` expects **one argument: altitude in meters**. To obtain the equivalent of Fig. 16 from Ershova et al., 2024 run:

```bash
./bin/horizontal_structure 25d3     # 25 km
./bin/horizontal_structure 50d3     # 50 km
./bin/horizontal_structure 75d3     # 75 km
./bin/horizontal_structure 100d3    # 100 km
```

For an altitude `A` meters, the output file is:

```text
./results/surface_dust_mass_density_<A_km>km.dat
```

where `<A_km> = nint(A / 1000)`, e.g.:

- `surface_dust_mass_density_25km.dat`
- `surface_dust_mass_density_50km.dat`, etc.

Then plot in R:

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

converts density from **kg/m³** to **g/m³**, and generates `./results/horizontal_slices.png` and creates 2×2 panel figure.

---
This README is specific to this branch and documents only the workflows and conventions that differ from the main branch.
