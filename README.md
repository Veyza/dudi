= README for DUDI (DUst DIstribution) Package

DUDI is a Fortran-95 software package developed to simulate the two-body
model for the distribution of dust ejected from the surface of an
atmosphereless celestial body. This tool is based on the research by
Anastasiia Ershova & Juergen Schmidt, as detailed in the paper:
Ershova, A. & Schmidt, J. (2021). Two-body model for the spatial
distribution of dust ejected from an atmosphereless body. Astronomy and
Astrophysics, 650.

== License

DUDI is distributed under GNU GENERAL PUBLIC LICENSE Version 3.

= Table of Contents

1. Prerequisites
2. Basic Usage
3. Specifying Parameters
4. Specifying Functions Describing the Dust Ejection
5. Utilizing the DUDI Subroutine
6. Applying the DUDI Subroutine
7. Output of the Results
8. Controlling the Accuracy of Calculations
9. Compilation
10. Running Examples
11. Possible Issues
12. Changes in the Latest Version


== 1. Prerequisites

To compile and run DUDI, the following is required:
* gfortran - Fortran compiler

For generating plots from the example applications described in this
README, the following additional software is needed:
* Python3, NumPy, and Matplotlib

Note: If these are not installed, DUDI will still produce results, but
plots will not be generated.



== 2. Basic Usage

DUDI requires setting up specific parameters before running simulations.
These parameters can be configured in various `..f90` files as described
in the later sections of this README. For details on modifying these
parameters and using the software to simulate dust distributions, refer
to the sections 'Specify the key parameters' and 'Running Examples'.



== 3. Specifying Parameters

Configure model parameters in `src/const.f90` for physical properties and
simulation settings:

* `moon_mass`   (real(8)): Mass of the moon (kg).
* `rm`          (real(8)): Radius of the moon (meters).
* `rho`         (real(8)): Density of dust particles (kg/m^3).
* `flux`        (logical): Set `.TRUE.` for dust flux, `.FALSE.` for density.
* `p`           (integer): Output type (0 = number density, 1 = mean radius, 2 =
                            cross section, 3 = mass density).
* `rmin`        (real(8)): Lower size limit of dust particles (microns).
* `rmax`        (real(8)): Upper size limit of dust particles (microns).
* `GRN`         (integer): Number of Gu(Rmin,Rmax) values precomputed.
* `order_R`     (integer): Integration order for Gu(Rmin,Rmax) values.
* `order_v_el`  (integer): Integration order for elliptic trajectories.
* `order_v_hy`  (integer): Integration order for hyperbolic trajectories.



== 4. Specifying Functions Describing the Dust Ejection

Customize dust ejection characteristics in `src/distributions_fun.f90`. This file
includes functions that define the size, speed, and direction distribution of
dust ejection, and the dust production rate. Each function utilizes Fortran's
`select case` operator, where a selector specifies a different
distribution to compute. Several distributions are already implemented and are
used in the example applications.

* `size_distribution(R, sd, fR)`:
    - `R` (real(8)): Particle radius (microns)
    - `sd` (integer): Selector for the size distribution type
    - Returns `fR` (real(8)): Probability density function (PDF) of size
      distribution evaluated at R.
 -- Write your own distribution at`case(4)` or further.

* `ejection_speed_distribution(ud, u, R, fu)`:
    - `ud` (type(ejection_speed_properties)): Structure containing parameters
      of the ejection speed distribution, defined as follows:
        * `ud%ud_shape` (integer): Selector for the ejection speed PDF.
        * `ud%umin` (real(8)): Minimum ejection speed (m/s).
        * `ud%umax` (real(8)): Maximum ejection speed (m/s).
    - `u` (real(8)): Ejection speed (m/s)
    - `R` (real(8)): Particle radius (microns)
    - Returns `fu` (real(8)): PDF of ejection speed distribution.
 -- Write your own distribution at`case(3)` or further.

* `ejection_direction_distribution(distribution_shape, wpsi, psi, lambdaM,
  zeta, eta, fpsi)`:
    - istribution parameters defining the shape and orientation
    - Returns `fpsi` (real(8)): PDF of ejection direction distribution.
 -- Write your own distribution at`case(4)` or further.

* `production_rate(t, gamma0, ratefun, gammarate)`:
    - `t` (real(8)): Time of ejection (seconds)
    - `gamma0` (real(8)): Parameter used in defining the production rate
    - `ratefun` (integer): Selector for the rate function
    - Writes result to `gammarate` (real(8)).
 -- Write your own distribution at`case(3)` or further.



== 5. Utilizing the DUDI Subroutine

The `DUDI(density, point, source, tnow)` subroutine performs the calculations.
The subroutine is located in the `src/integrator.f90` file. It writes the result to
the `density` variable and takes as input:

- `density`: An array of two real numbers, which the subroutine updates:
  - First number: Density of particles on bound orbits.
  - Second number: Density of particles on unbound orbits.
- `tnow`: real(8) variable specifying the moment in time at which the density is
  calculated.
- `source`: A structure containing detailed information about the dust source.
- `point`: A structure specifying the spacecraft's position in space.

Both `source` and `point` structures include multiple precomputed quantities used
at various stages of the calculations. These structures do not have all independent
values, requiring careful setup and validation.

The `point` is a structure of a derived type `position_in_space` which defines the
point in space where the density is calculated.

*Components of the `position_in_space` structure:*
- `point%r` (real(8)): Radial distance of the point from the center of the moon
   (meters).
- `point%r_scaled` (real(8)): Radial distance scaled to the moon's radius.
   This is calculated as:
   point%r_scaled = point%r / rm
   where `rm` is the radius of the moon set in module const..f90.
- `point%alpha` (real(8)): Polar angle in the moon's centered coordinate system
   (radians).
- `point%beta` (real(8)): Eastern longitude in the moon's centered coordinate
   system (radians).
- `point%rvector` (real(8) 3D-vector): Cartesian coordinates of the point
   in the moon-centered coordinate system, can be calculated using:
   point%rvector(1) = point%r * sin(point%alpha) * cos(point%beta)
   point%rvector(2) = point%r * sin(point%alpha) * sin(point%beta)
   point%rvector(3) = point%r * cos(point%alpha)
- `point%compute` (logical): Indicates whether to calculate density at this
   point. Useful for excluding unnecessary areas from computations, enhancing
   efficiency.


The `source` is a structure of a derived type `source_properties` which contains
parameters describing the dust ejection.

*Components of the `source_properties` structure:*
- `source%r` (real(8)): Radial distance of the point source from the moon center
   (meters).
- `source%alphaM` (real(8)): Polar angle of the point source (radians).
- `source%betaM` (real(8)): Eastern longitude of the point source (radians).
- `source%rrM` (real(8) 3D-vector): Cartesian coordinates of the point source
   in the moon-centered coordinate system. Calculated as:
   source%rrM(1) = source%r * sin(source%alphaM) * cos(source%betaM)
   source%rrM(2) = source%r * sin(source%alphaM) * sin(source%betaM)
   source%rrM(3) = source%r * cos(source%alphaM)
- `source%zeta` (real(8)): Zenith angle of the axis around which ejection
   is symmetrical (radians).
- `source%eta` (real(8)): Azimuth of the axis
   (counted from the local North, clockwise) (radians).
- `source%symmetry_axis` (real(8) 3D-vector): Unit vector in moon-centered
   coordinate system pointing to the direction of the ejection symmetry axis.
   module input_data contains the subroutine jet_direction so that
   the following line of the code:

     call jet_direction(source%alphaM, source%betaM, source%zeta, source%eta, &
                         source%rrM, source%symmetry_axis)

   writes the correct values into the source%symmetry_axis

- `source%ejection_angle_distr` (integer): Selector parameter for different
  types of ejection angle distributions.
- `source%ud` (type(ejection_speed_properties)): Structure containing parameters
   of the ejection speed distribution, including:
   - `ud_shape` (integer): Selector for the ejection speed distribution PDF.
   - `umin` (real(8)): Minimum possible ejection speed (m/s).
   - `umax` (real(8)): Maximum possible ejection speed (m/s).
   - `source%sd` (integer): Selector parameter for different size distributions.
- `source%ui` (real(8) array of GRN elements): Interpolation grid
   for `Gu(Rmin,Rmax)`.
- `source%Gu_precalc` (real(8) array of GRN elements): Precalculated values
   of `Gu(Rmin,Rmax)`
   module gu contains the subroutine Gu_integral
   so that the following line of code

call Gu_integral(source%ui, source%Gu_precalc, source%sd, source%ud, rmin, rmax)

   writes the correct values into source%ui and source%Gu_precalc.
- `source%production_fun` (integer): Parameter defining the function used as the
   time-dependent dust production rate.
- `source%production_rate` (real(8)): The dust production rate in case
   of stationary ejection, used as a parameter for the function in
   non-stationary cases.
- `source%is_jet` (logical): True if the ejection is concentrated. Recommended
   to be set when `source%omega < 0.1`.



== 6. Applying the DUDI Subroutine

The `DUDI` subroutine, located in `src/integrator.f90`, performs the numerical
integration required to compute dust density. This subroutine must be called
from a main program, which manages the input and output of the data.

A template for such a main program, `examples/main_program.f90`, is provided within the
package to assist users in setting up their simulations quickly and efficiently.

The `main_program..f90` file is provided as a template for the main program. It
uses subroutines from the `input_data` module to create:
- An array of dust sources.
- An array of points in space for density calculations.

The program employs nested loops over sources and points, invoking the `DUDI`
subroutine from the `integrator` module at each point. The resulting density at
each point is the aggregate of densities from all sources.

OpenMP is utilized to accelerate the calculations, leveraging the independence
of calculations for each source-point pair. Ultimately, a function from
the `dataout` module is called to write the results to a file.


*Customization:*
- Users are encouraged to modify the `examples/main_program.f90` template and the related
  input/output subroutines according to their specific research needs.



== 7. Output of the Results

The `result_out` function in the `src/dataoutmod.f90` module handles the output of
simulation results. It writes data to the file `twobody_model_result.dat` in the
`./results/` directory with the following column structure:
- 1st column: Total density at the point.
- 2nd column: Density of particles on elliptic orbits.
- 3rd column: Density of particles on hyperbolic orbits (1st column is the sum
  of the 2nd and 3rd columns).
- 4th column: Spacecraft's radial distance in meters.
- 5th column: Spacecraft's latitude in degrees.
- 6th column: Spacecraft's eastern longitude in degrees.

Users have the flexibility to customize the output by writing their own function
within the `src/dataoutmod.f90` module.



== 8. Controlling the Accuracy of Calculations

The program monitors numerical accuracy. If a poorly conditioned case occurs
where accuracy falls below standard levels, it generates a file named `fort.666`.
This file logs details about the issues, aiding in pinpointing the problems.

Even with several indications of poor accuracy in `fort.666`, this does not
completely discredit the overall solution. Numerical difficulties may appear
at some integration steps for each "source and spacecraft position" pair, but
not necessarily at all steps.

Furthermore, if the number of warnings in `fort.666` exceeds the limit set by
`maxNofWarnings` in the `src/const.f90` module, the program will stop and print a
warning message to the command line.



== 9. Compilation

A makefile is provided for easy compilation using the gfortran compiler.

- To compile the program, enter the following command in the command line:

                   make
This command produces an executable dudi

- The command:

                   make clean

removes all files *.o, *.mod, and *dudi (including the executable produced by
make-command)

The makefile also allows one to obtain model results and generate plots for
the examples described in the paper by Ershova & Schmidt, 2021. Below are
detailed instructions on how to run the example code.



== 10. Running Examples

DUDI includes several examples to demonstrate its capabilities. Each example
uses specific parameters and input files.

*Example 1: The Number Density Profile
of the E2 Flyby of the Cassini Spacecraft at Enceladus*
  - Main Program: `examples/enceladus_example.f90`
    - Performs a loop over 100 points along the Cassini spacecraft trajectory.
    - Computes the number density of dust from a tilted jet representing the
      Enceladus dust plume.
  - Parameter Settings in `const..f90`:
    * `moon_mass = 1.08022d+20` // Enceladus's mass in kg
    * `rm = 252d+3`             // Enceladus's mean radius in meters
    * `flux = .FALSE.`          // Interested in number density, not flux
    * `p = 0`                   // Outputs a number density profile
    * `rmin = 1.6d0`            // Lower threshold of HRD sensitivity (microns)
    * `rmax = 6.0d0`            // Upper boundary (microns), assuming a small
                                // probability to detect larger particles
    * `GRN = 2000`              // Number of precalculated values of Gu integral
    * `order_R = 30`            // Order of integration for Gu(Rmin,Rmax) over R
    * `order_v_el = 50`         // Gaussian quadrature order for elliptic
                                // trajectories
    * `order_v_hy = 20`         // Gaussian quadrature order for hyperbolic
                                // trajectories
  - Input Files:
    * `Enceladus_jet.dat`        // Contains source properties
    * `Cassini_E2_flyby.dat`     // Contains spacecraft coordinates
      - These files are located in the directory `input_data_files`.
      - Coordinates, zenith angles, and azimuth of the jet adopted from Porco
        et al., 2014.
      - Coordinates obtained using the NAIF SPICE package.
  - Command to Run:
   - Execute: `make enceladus`
   - This compiles the program, producing `enceladus_model`, and runs it.
  - Output:
   - Results are saved in `E2_profile.dat` within the `results` directory.
   - Output format: 2 columns, with time in seconds from 2005-07-14 19:49:21
     (closest approach) and the dust number-density.
   - If Python3, NumPy, and Matplotlib are installed, a plot comparing the model
     profile with observational data will be displayed.


*Example 2: Dust Deposition on the Surface of Europa*
  - Main Program: `examples/europa_example.f90`
    - Calculates the dust flux onto Europa's surface for four distinct cases,
      each representing a combination of one of two different size distributions
      and one of two ejection direction distributions.
    - The total mass production rate is computed separately for each of these
      four scenarios.
  - Parameter Settings in `const..f90`:
    * `moon_mass = 4.8d+22`  // Europa's mass in kg

 The plot shown in the paper (Ershova & Schmidt, 2021) was generated with an
 erroneous value of rm = 3.12e+6, Europa's mean diameter (instead of radius)
 in meters. We retain  this value here to allow reproduction of Fig. 13 from
 (Ershova & Schmidt, 2021). The case of surface deposition on Europa was purely
 illustrative and did not lead to scientifically significant conclusions.
 Therefore, we do not consider it necessary to publish an erratum to address
 this issue.

    * `rm = 3.12d+6`         // Europa's mean diameter in meters
    * `rho = 920d0`          // Dust grains density in kg/m3 (assuming water
                             // ice composition)
    * `flux = .TRUE.`        // Calculating flux
    * `p = 3`                // Outputs a mass flux
    * `rmin = 0.2d0`         // Lower limit of size distributions (microns)
    * `rmax = 20d0`          // Upper limit of size distributions (microns)
    * `GRN = 2000`           // Number of precalculated values of Gu integral
    * `order_R = 10`         // Order of integration for Gu(Rmin,Rmax) over R
    * `order_v_el = 30`      // Gaussian quadrature order for elliptic
                             // trajectories
    * `order_v_hy = 5`       // Order does not matter as all particles are on
                             // elliptic trajectories
  - Source and Location Settings:
    - No input files used.
    - Parameters for the sources and deposition points are set in the subroutine
      `get_europa_input` in the module `inputdata..f90`.
  - Command to Run:
    - Execute: `make europa`
    - Compiles the program, producing `europa_model`, and runs it.
  - Output:
    - Results are saved in files `narrow_jet_shallow_sd.dat`,
      `diffuse_source_steep_sd.dat`, `diffuse_source_shallow_sd.dat`, and
      `narrow_jet_steep_sd.dat` in the `./results` directory.
    - Each output file includes two columns: the first column shows the distance
      from the source in km, and the second column shows the mass flux in
      kg/m^2/s.
    - The total mass production rate for each size distribution is printed in
      the terminal.
    - If Python3, NumPy, and Matplotlib are installed, a plot will be produced
      and saved as `mass_deposition.png` in the `./results` directory.


*Example 3: The Images of a Fictive Volcano Erupted on Io*
  - Main Program: `examples/io_example.f90`
    - Utilizes additional modules including `image_construction..f90`.
    - Constructs images simulating a CCD camera with 128x128 pixels,
      calculating the line of sight for each pixel across a grid of 41 points.
    - Computes the 2nd moment of the number density at each point and
      integrates along the lines of sight to construct an image.
    - The moon's disc is visible in the image; pixels covering it have zero
      brightness. Points along lines of sight crossing the moon's disk are
      excluded from calculations (point%compute = .FALSE.).
    - Creates 9 images representing non-stationary ejection at 9 specific
      moments.
  - Parameter Settings in `const..f90`:
    * `moon_mass = 8.94d+22`  // Io's mass in kg
    * `rm = 1.8216d+6`        // Io's mean radius in meters
    * `flux = .FALSE.`        // Not calculating flux
    * `p = 2`                 // Outputs a cross section covered by the dust
    * `rmin = 0.2d0`          // Lower size limit (microns), within optical
                              // wavelength interval
    * `rmax = 0.4d0`       // Upper size limit (microns), twice less than rmin
    * `GRN = 3`            // Minimal value for GRN due to uniform distribution
    * `order_R = 10`       // Order of integration Gu(Rmin,Rmax) over R
    * `order_v_el = 5`     // Small quadrature order for short velocity interval
    * `order_v_hy = 5`     // All particles are on elliptic trajectories, this
                              parameter is not important
  - Source and Location Settings:
    - No input files used.
    - Source parameters are configured in `get_volcano_params` subroutine
      within `inputdata..f90`.
    - Grid for calculating particle cross section is formed in the `line_of_sight`
      subroutine in `image_construction..f90`.
  - Command to Run:
    - Execute: `make io`
    - Compiles the program, producing `io_model`, and runs it.
  - Output:
    - Results are saved in files named `1.dat` to `9.dat` in the `results`
      directory, corresponding to each moment of image capture.
    - If Python3, NumPy, and Matplotlib are installed, 9 images of the volcano
      are produced with the time after the start of ejection noted in the
      bottom-right corner. These images are named `volcano_1.png` to
      `volcano_9.png` and stored in the `results` directory.



== 11. Possible Issues

If you encounter a "Segmentation fault" when running a program compiled with
the `-openmp` key (as specified in our Makefile), it may be due to the stack
size limit.
To address this:

* In `bash`, increase the stack size by running:
  `ulimit -s unlimited`
* In `tcshell`, use:
  `limit stacksize unlimited`



== 12. Changes in the Latest Version

=== Version 1.0.1

- **Warning Management:**
  DUDI now tracks warnings logged to `fort.666`. If the count exceeds the
  `maxNofWarnings` set in `const..f90`, the program will terminate and alert
  the user via the command line. This prevents the creation of large warning
  files due to parameter errors.

- **Update to `size_distribution`:**
  The `size_distribution` function in `distributions_fun..f90` no longer
  uses the "p" parameter from `const..f90`. It now returns only the PDF for
  the specified grain radius. Required calculations based on "p" are handled
  by functions in `gu..f90`.

- **Ejection Velocity Correction:**
  Fixed an issue in `integrator..f90` where the integration was
  incorrectly handled when minimum ejection velocity exceeded the moon's
  escape velocity.

- **Io Example Optimization:**
  Multidimensional arrays are reshaped and nested loops reordered in the Io
  example, enhancing performance.



=== Version 1.1.0

- **Subroutine `Gu_integral` enhancements:**
  This subroutine evaluates the function G_u^p(R_min, R_max) (see Eq. 14 of
  Ershova & Schmidt, 2021) for a grid of u-values. In previous versions,
  `Gu_integral` used integration limits predefined in the module const and
  always returned an array of evaluated integrals.

  Now, `Gu_integral` accepts the integration limits, `rlim1` and `rlim2`, as
  input variables. Additionally, if `rlim1 = rlim2`, `Gu_integral` returns the
  product of the size distribution and ejection direction distribution,
  evaluated at `rlim1` and the grid of u-values. This enhancement allows for
  the calculation of the differential number density (i.e., without integration
  over R in Eq. 20). The evaluated value corresponds to
  n(r, \alpha, \beta, R) instead of n(r, \alpha, \beta, R_min < R < R_max).

- **Source file extension changed to .f90:**
  This change allows the Intel compiler to be used for DUDI compilation
  without any additional consequences.



=== Version 1.2.0

  A feature has been added to compute the average velocity vector of dust grains
  at a given point and time.

- **New Subroutine `DUDI_mean_velocity`:**
  The subroutine `DUDI_mean_velocity` has been added to the module
  `integrator.f90`. It computes the average velocity vector of dust grains at
  a given point and time.

  **Input:**
  - The input parameters are identical to those of the main subroutine `DUDI`.

  **Output:**
  - Returns an array of 8 real numbers:
    1. A 3D vector representing the average velocity of dust grains moving upward.
    2. A 3D vector representing the average velocity of dust grains moving downward.
    3. The number density of dust grains moving upward.
    4. The number density of dust grains moving downward.

  **New Subroutine `Integrand_mean_flux`:**
  - The subroutine `DUDI_mean_velocity` calls a newly introduced subroutine
    `Integrand_mean_flux` from the module `twobody_fun.f90`.

  - A detailed description of the algorithm can be found in the file
    **`DUDI_average_velocity_vectors.pdf´** available in this repository.


=== Version 1.2.1
  - Restructured repository with clear src/, examples/, scripts/, bin/, build/, and results/ folders.
  - Updated Makefile and .gitignore accordingly
