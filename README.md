This repository contains DUDI (DUst DIstribution) the code that implements 
the two-body model for the configuration of dust ejected from the surface 
of an atmosphereless body by Anastasiia Ershova & Juergen Schmidt. 
The derivation of the model, along with the formulae implemented in the code,
as well as detailed examples of the model capabilities 
can be found in the paper:
*reference to the paper*

DUDI is written in Fortran-95. A makefile is provided for compilation.
The user must have the compiler gfortran installed. There are no other 
prerequisites for running the program itself.

There are three example applications of the model described in the paper. 
For the convenience of the user implementations of these examples are also 
provided with this package. To obtain the plots corresponding to each example
the user must have installed python3, numpy, and matplotlib. The functions 
that make the plots are also contained in this repository.

The code is distributed under the terms of GNU GENERAL PUBLIC LICENSE Version 3

Bellow find the instructions how to work with DUDI.
The outline of the instructions is following:

Specify the key parameters 
--------------------------
Specify the functions describing the dust ejection
--------------------------------------------------
Supply input data
-----------------
Call the subroutine DUDI(num_dens, point, source, tnow)
-------------------------------------------------------
Output the result
-----------------
Controlling the accuracy of calculations
----------------------------------------
Compilation
-----------
Example 1. The number density profile of the E2 flyby of the Cassini spacecraft
at Enceladus.
-------------------------------------------------------------------------------
Example 2. Dust deposition on the surface of Europa
---------------------------------------------------
Example 3. The images of a fictive volcano erupted on Io
--------------------------------------------------------
Possible issues
---------------




Specify the key parameters 
-------------------------------------------------------------------------------
    open the file const.f95 and set values to the following parameters:
	

    moon_mass                  mass of the moon, kg

    rm                         radius of the moon, meters


    rho                        density of the dust particles' material 
                               in kg/m^3,(needed if one wants to compute
                               mass density or mass fluxes)

    flux                       if set .TRUE. then dust flux through the
                               surface parallel to the moon surface is computed
                               instead of density

    p                          0 -- number density is computed, 1 -- mean
                               radius, 2 -- cross section, 3 -- mass density

    rmin                       lower boundary for function Gu(rmin, rmax),
                               microns

    rmax                       upper boundary for function Gu(rmin, rmax),
                               microns

    GRN                        number of Gu(Rmin,Rmax) values that are
                               precalculated other values are obtained
                               by interpolation from precalculated values

    order_R                    order of Gauss-Legendre quadrature formula
                               used for integration over particle radius R
                               (possible values are 5, 10, 20, 30)

    order_v_el                 order of Gauss-Legendre quadrature formula
                               used for integration over velocity to obtain
                               the particles density separately for particles
                               on bound (elliptic) trajectories
                               (possible values are 5, 10, 20, 30, 40, 50)
                               
    order_v_hy                 order of Gauss-Legendre quadrature formula
                               used for integration over velocity to obtain
                               the particles density separately for particles
                               on unbound (hyperbolic) trajectories
                               (possible values are 5, 10, 20, 30, 40, 50)                               
                               


Specify the functions describing the dust ejection
------------------------------------------------------------------------------
Size, ejection speed, and ejection direction distributions, as well as
the time-dependent dust production rate are represented by functions listed
in the file "distributions_fun.f95". There are several examples for such 
distributions already implemented in the code, but the user can provide 
his/her own expressions in the body of the corresponding functions. 
Choices are implemented in the provided code in terms of the "select case" 
operator. The functions in the file "distributions_fun.f95" that likely 
need to be modified are the following:

   ____________________________________________________________________________
    size_distribution(R, sd, mom)

      R             particle radius in microns
      sd            integer number used to select the PDF
      mom           parameter defining the obtained quantity: 
                    0 -- number density, 1 -- mean radius, 
                    2 -- cross section, 3 -- volume 
                    (same as the parameter p in Formula 15)

      the value returned by this function is stored in the variable fR

      to use your own distribution you should go to the lines:

      case(4)
      ! HERE IS THE PLACE FOR WRITING YOUR OWN PDF
        fR = 0d0

      and replace "fR = 0d0" with fR = the desired PDF


   ____________________________________________________________________________
    ejection_speed_distribution(ud, u, R)

      ud            a structure which contains parameters of the ejection speed
                    distribution designated as ud%ud_shape, ud%umin, ud%umax,
                    where
      ud%ud_shape   integer used to choose the PDF
      umin          minimum ejection speed in m/s
      umax          maximum ejection speed in m/s
      u             ejection speed in m/s
      R             particle radius in microns

      the value of the function is stored in the variable fu

      to use your own distribution you should got to the lines:

      case(3)
      ! HERE IS THE PLACE FOR WRITING YOUR OWN PDF
        fu = 0d0

      and replace "fu = 0d0" with fu = the desired PDF


   ____________________________________________________________________________
     ejection_direction_distribution(distribution_shape, wpsi, psi, lambdaM,
                                                                     zeta, eta)

      distribution_shape 	an integer parameter defining the distribution
                                shape
      wpsi                      polar angle in the coordinate system where 
                                the distribution is axisymmetric (radians)
      psi                       polar angle in the horizontal coordinate system
      lambdaM                   azimuth in the horizontal CS (radians)
      zeta                      zenith angle of the distribution symmetry axis 
                                in the horizontal CS (radians)
      eta                       azimuth of the distribution symmetry axis in
                                the horizontal CS (radians)

      (see Figure 1 of the paper)

      the value of the function is stored in the variable fpsi

      The PDF is written for the variable wpsi in the coordinate system where
      the ejection is axisymmetric. In this system azimuth of ejection is 
      assumed to be distributed as 1/2pi. 

      to use your own distribution you should go to the lines:

      case(4)
      ! HERE IS THE PLACE FOR WRITING YOUR OWN PDF
        fpsi = 0d0

      and replace "fpsi = 0d0" with fpsi = the desired PDF


   ____________________________________________________________________________
    function production_rate(t, gamma0, ratefun)

      t                        moment of ejection (seconds)

      gamma0                   parameter which can be used in the definition
                               of the function. In the implemented examples
                               gamma0 is the production rate at maximum
                               (particles/second)

      ratefun                  parameter used to choose the expression for
                               the production rate


      to use your own expression for the production rate you should go 
      to the lines:

      case(4)
      ! HERE IS THE PLACE TO WRITE YOUR OWN FUNCTION FOR THE PRODUCTION RATE
        gammarate = 0d0

      and replace "gammarate = 0d0" with gammarate = the desired function





Supply input data
-------------------------------------------------------------------------------
All the calculations are performed by the subroutine 
IntegrateNumberDensity(num_dens, point, source, tnow) 
which will return the result in the variable num_res

    num_res	        array of 2 real numbers. The first number is
                        the density of particles on bound orbits,
                        the second number is the density of particles
                        on escaping trajectories 

    tnow                moment of time at which the density is calculated

    source              structure containing information about the source
                        ejecting dust 

    point               structure containing the information about
                        the spacecraft position in space. 


    Both structures contain multiple quantities that must be precomputed 
    as they are used at different stages of the calculations. 
    Not all the values are independent. 
    Bellow you find the description of the structures "source" and "point"
    and the instructions how to set the correct values for their components.

_______________________________________________________________________________
point is a structure of a derived type position_in_space which defines 
the point in space where the density is calculated. 
position_in_space contains the following elements:

    real(8)                  point%r              spacecraft's radial distance 
                                                  from the center of the moon 
                                                  in meters

    real(8)                  point%r_scaled       spacecraft's radial distance 
                                                  from the center of the moon 
                                                  in units of the moon radius : 
                                                  *****************************
                                                  point%r_scaled = point%r / rm
                                                  *****************************

    real(8)                  point%alpha          spacecraft's polar angle in 
                                                  the moon's centered coordinate 
                                                  system (radians)

    real(8)                  point%beta           spacecraft's eastern longitude
                                                  in the moon's centered 
                                                  coordinate system (radians)

    3d-vector of real(8)     point%rvector        Cartesian coordinates of the 
                                                  spacecraft in the moon-centered
                                                  coordinate system :

    ***************************************************************  
    point%rvector(1) = point%r * sin(point%alpha) * cos(point%beta)
    point%rvector(2) = point%r * sin(point%alpha) * sin(point%beta)
    point%rvector(3) = point%r * cos(point%alpha)
    ***************************************************************
    
    logical                  point%compute        TRUE if the density in this 
                                                  point must be calculated (allows 
                                                  to exclude regions where the 
                                                  calculations are unnecessary, 
                                                  see Example 3)


_______________________________________________________________________________
source is a structure of a derived type source_properties which contains 
parameters describing the dust ejection

    real(8)     source%alphaM                   polar angle of the point 
                                                source (radians)

    real(8)     source%betaM                    eastern longitude of the point 
                                                source (radians)

    real(8)     source%rrM                      (3d-vector) Cartesian 
                                                coordinates of the point 
                                                source in the moon-centered CS:
                                                        
    *********************************************************                                                    
    source%rrM(1) = rm * sin(point%alphaM) * cos(point%betaM)
    source%rrM(2) = rm * sin(point%alphaM) * sin(point%betaM)
    source%rrM(3) = rm * cos(point%alphaM)
    *********************************************************
    
    real(8)     source%zeta                     zenith angle of the axis around
                                                which ejection is symmetrical

    real(8)     source%eta                      azimuth of this axis (counted 
                                               from the local North, clockwise)

    real(8)     source%symmetry_axis            (3d-vector) unit vector in moon
                                                -centered coordinate system 
                                                pointing to the direction of 
                                                the axis around which ejection 
                                                is symmetrical
    ***************************************************************************	
    module input_data contains the subroutine jet_direction so that 
    the following line of the code:
    ___________________________________________________________________________           

    call jet_direction(source%alphaM, source%betaM, source%zeta, source%eta, 
                                              source%rrM, source%symmetry_axis)
    ___________________________________________________________________________
																
    writes the correct values into the source%symmetry_axis
    ***************************************************************************

    integer          source%ejection_angle_distr     parameter to select from 
                                                     the given types 
                                                     of distributions

    type(ejection_speed_properties)     source%ud    parameters of the ejection
                                                     speed distribution 
                                                     which are 

                                 ud_shape      integer parameter used to select 
                                               the expression for the 
                                               distribution PDF
                                 umin          minimal possible ejection speed
                                               (m/s)
                                 umax          maximal possible ejection speed
                                               (m/s)


    integer         source%sd               parameter to select from the given 
                                            size-distributions

    real(8)         source%ui               (array of GRN elements) 
                                            interpolation grid for Gu(Rmin,Rmax)

    real(8)         source%Gu_precalc       (array of GRN elements) 
                                            Gu(Rmin,Rmax) precalculated

    ***************************************************************************
    module gu contains the subroutine Gu_integral 
    so that the following line of code
    ___________________________________________________________________________
   
    call Gu_integral(source%ui, source%Gu_precalc, source%sd, source%ud)
    ___________________________________________________________________________

   writes the correct values into source%ui and source%Gu_precalc. 
   Here p is the parameter defining the dimension of the function Gu. 
   (see Formula 15)
   rmin, rmax and p are defined in module const
   ****************************************************************************
   
    integer        source%production_fun      parameter defining the function 
                                              used as time-dependent dust 
                                              production rate. Any value <= 0 
                                              corresponds to the stationary 
                                              case

    real(8)        source%production_rate     the dust production rate in case 
                                              of the stationary ejection, in 
                                              a non-stationary case this 
                                              variable can be used as 
                                              a parameter for the function 
                                              production_rate

    logical        source%is_jet              .TRUE. if the ejection is 
                                              concentrated: the recommended 
                                              way to define source%is_jet:
                                           source%is_jet = (source%omega < 0.1)





The user can utilize the subroutines constructing the structures "source" and
"point" that are provided with the examples, or write her own ones in the 
module "inputdata.f95"



Call the subroutine DUDI(num_dens, point, source, tnow)
-------------------------------------------------------------------------------
Subroutine DUDI in the module integrand.f95 is the piece of code doing the main
work. It performes numerical integration to compute density at the given point
in space for the source with given properties at the given moment of time tnow.


The file "main_program.f95" is provided as a template of the main program. 
It calls the subroutines from the module input_data to obtain an array of 
sources ejecting dust and an array of the points in space where one wants to 
compute density. Then the program makes a double loop over the sources and 
over the points calling the subroutine DUDI from the module integrator. 
In each point the result density is the sum of densities of the 
dust from all the sources. OpenMP is used to speed up the calculations as for
 each pair of source and point the calculations are independent. Finally, 
the main program calls the function from module dataout to write the result 
into the file.

Modify the file "main_program.f95" and the subroutines for input and output 
for your own needs.



Output of the result
------------------------------------------------------------------------------
The function result_out in the module "dataoutmod.f95" is provided for output
of the result. It writes the result into the file "twobody_model_result.dat" 
in the directory "./results/" as follows:

    1st column is the density at the point

    2nd column is the density of particles on elliptic orbits

    3rd column is the density of particles on hyperbolic orbits

    (1st column is sum of the 2nd and the 3rd)
		
    4th column is the spacecraft's radial distance in meters

    5th column is the spacecraft's latitude in degrees

    6th column is the spacecraft's eastern longitude in degrees



If needed, the user can write his own function for output in the module 
"dataoutmod.f95".



Controlling the accuracy of calculations
-------------------------------------------------------------------------------
If the program encounters a numerically poorly conditioned case, in which 
accuracy is bellow the usual level the file "fort.666" is created into which
information about the difficulty is printed out. This allows one to locate 
where the difficulty arises and how bad it is. Normally, there will be no such
file. Having several indications of a bad accuracy also does not completely 
discredit the whole solution, because for each pair "source and spacecraft 
position" the numerically difficult cases appear at some integration steps, 
not at all of them.



Compilation
-------------------------------------------------------------------------------
The makefile is provided for compilation. 
It utilizes the common compiler gfortran

_______________________________________________________________________________
Compile the program typing in command line:

    make

the command produces an executable file named "twobody_model"

_______________________________________________________________________________
Command
		
    make clean

will remove all the executable files and files with extension .o and .mod
  

The makefile can be also used to obtain the model results and the corresponding 
plots that are described in the paper as the examples. Bellow you find the 
detailed instructions how to run the example code.
								




Example 1. The number density profile of the E2 flyby of the Cassini spacecraft
at Enceladus
-------------------------------------------------------------------------------
The main program enceladus_example.f95 performs a loop over 100 points along 
the spacecraft trajectory and computes in these points the number density of 
dust from one tilted jet representing the Enceladus dust plume.

Set the following values to the parameters in the module const.f95:

  moon_mass = 1.08022d+20 ! Enceladus's mass in kg
  rm = 252d+3             ! Enceladus's mean radius in meters

  flux = .FALSE.          ! we are currently not interested in flux
	
  p = 0                   ! this will give us a number density profile
  rmin = 1.6d0            ! micron, the lower threshold of the HRD sensitivity
  rmax = 6.0d0            ! micron, a fairly large upper boundary, it is
                          ! assumed that the probability to detect a particle
                          ! larger than rmax is negligibly small
		
  GRN = 2000              ! number of precalculated values of Gu integral
  order_R = 30            ! order of integration Gu(Rmin,Rmax) over R
                          ! to obtain Gu 
  order_v_el = 50         ! order of Gaussian quadrature for integration of 
                          ! n(r, alpha, beta, v, theta, lambda) over velocity 
                          ! interval corresponding to elliptic trajectories
  order_v_hy = 20         ! order of Gaussian quadrature for integration of 
                          ! n(r, alpha, beta, v, theta, lambda) over velocity 
                          ! interval corresponding to hyperbolic trajectories

The directory input_data_files contains files "Enceladus_jet.dat" and 
"Cassini_E2_flyby.dat" with the input data. 
"Enceladus_jet.dat" contains a line with the source properties.
The coordinates, zenith angles and azimuth of the jet were adopted 
from the paper of Porco et al., 2014.

"Cassini_E2_flyby.dat" contains the coordinates of the Cassini spacecraft at 
the required moments. The coordinates were obtained with the SPICE package 
and stored into this file for simplicity.

The following command

make enceladus

will compile the program, producing the executable "enceladus_model", 
and run it. 
The model result will be output into the file "E2_profile.dat"
in the directory "results".
The output consists of 2 columns. The first column is the time in seconds 
counted from 2005-07-14 19:49:21 (the moment of the closest approach),
the second column is the dust number-density.
If your computer has python3, numpy and matplotlib installed then a plot 
with comparison of the model profile with the observational data will 
be shown on the screen.


Example 2. Dust deposition on the surface of Europa
-------------------------------------------------------------------------------
The main program europa_example.f95 performs calculation of the dust flux 
onto the moon surface. The total mass production rate is computed for 2 
different size distributions. The surface deposition is calculated 
separately for 4 sources different in size and ejection direction distributions

Set the following values to the parameters in the module const.f95:

  moon_mass = 4.8d+22     ! Europa's mass in kg
  rm = 1.56d+6          ! Europa's mean radius in meters
	
  rho = 920d0             ! dust grains density in kg/m^3 
                          ! (assuming pure water ice composition)
  flux = .TRUE.           ! this time we are calculating flux
	
  p = 3                   ! this will give us a mass flux
  rmin = 0.2d0            ! the lower limit of the applied size distributions,
                          ! microns
  rmax = 20d0             ! the upper limit of the applied size distributions,
                          ! microns
		
  GRN = 2000              ! number of precalculated values of Gu integral
  order_R = 10            ! order of integration Gu(Rmin,Rmax) over R
                          ! to obtain Gu 
  order_v_el = 30         ! order of Gaussian quadrature for integration of 
                          ! n(r, alpha, beta, v, theta, lambda) over velocity 
                          ! interval corresponding to elliptic trajectories

  order_v_hy = 5          ! does not matter in this example because all
                          ! the particles are on elliptic trajectories

No input files are used. The parameters of the sources and the points on 
the surface, where the flux is calculated, are set in the subroutine 
get_europa_input in the module "inputdata.f95"

The following command

	make europa

will compile the program, produce the executable "europa_model", and run it. 
The result will be put in the files "narrow_jet_shallow_sd.dat", 
"diffuse_source_steep_sd.dat", "diffuse_source_shallow_sd.dat", and 
"narrow_jet_steep_sd.dat" in the directory "./results". 
They correspond to the four sources with different combinations of the size
distribution and the ejection direction distribution. Each output file 
consists of two columns where the first column will be the distance from 
the source in km, and the second column will be the mass flux in kg/m^2/s. 
The total mass production rate for each size distribution will be printed 
in the terminal. If you have python3, numpy and matplotlib installed 
on your computer a plot will be produced and saved to the file
"mass_deposition.png" in the directory "./results".



Example 3. The images of a fictive volcano erupted on Io
-------------------------------------------------------------------------------
The main program io_example.f95 constructs images using in addition to 
the usual list of necessary modules the module "image_construction.f95". 
For a hypothetical position and pointing of the CCD camera with the 
dimension of 128x128 pixels, the program defines the direction of the 
line of sight of each pixel and a grid of 41 points on each of these lines 
of sight. The program computes the 2nd moment of the number density 
at each point and then integrates along the lines of sight constructing 
an image. The moon's disc is visible in the image, the pixels covered by 
it have zero value of brightness. Therefore, the points along the lines 
of sights crossing the moon disk are excluded from the calculations by setting
for them parameter point%compute = .FALSE. at the stage of the formation of 
the points grid. The ejection is considered non-stationary and 
the program creates 9 images as if they were taken at 9 specified moments.

set the following values to the parameters in the module const.f95:

  moon_mass = 8.94d+22    ! Io's mass in kg
  rm = 1.8216d+6          ! Io's mean radius in meters

  flux = .FALSE.          ! we are currently not interested in flux

  p = 2                   ! this will give us a cross section covered by 
                          ! the dust particles
  rmin = 0.2d0            ! micron, we assume that we observe the particles
                          ! with sizes inside 
                          ! the optical wavelength interval
  rmax = 0.4d0            ! therefore the radii are twice less

  GRN = 3                 ! in this example we use a uniform distribution 
                          ! of ejection
                          ! velocity, therefore Gu = const so we set the 
                          ! minimal possible value for GRN.
  order_R = 10            ! size interval of integration Gu(Rmin,Rmax) over
                          ! R to obtain Gu 
                          ! is relatively short, so the order of quadrature
                          ! can be smaller too.
  order_v_el = 5          ! the interval for the integration over velocity
                          ! is short, so we use 
                          ! a small order of the quadrature
                          
  order_v_hy = 5          ! does not matter in this example because all the
                          ! particles are on elliptic trajectories

No input files are used. The parameters of the source are set up directly 
inside the subroutine get_volcano_params in the module "inputdata.f95". 
The grid of the points, where the particles cross section is calculated, 
is formed in the subroutine line_of_sight in the module image_construction.f95. 
The same module contains the parameters defining the image dimension in radians
and in pixels. 

The command

make io

will compile the program, producing the executable "io_model", and run it. 
The result will be outputted into the files "*.dat", (* stands for numbers from
1 to 9)  in the directory "results". If the user has python3, numpy and 
matplotlib installed there will be also produced 9 images of the volcano 
with the corresponding time in seconds after the start of ejection written 
in the bottom-right corner. The plots will be in the directory "results", 
they will be named "volcano_*.png" (* stands for numbers from 1 to 9).

	

Possible issues
-------------------------------------------------------------------------------
If you obtain a "Segmentation fault" message when running the program compiled
with the -openmp key (also used in our Makefile), try setting the core dump 
file size unlimited. In bash it is done by the command:

                  ulimit -s unlimited
In tcshell use comand:
                  limit stacksize unlimited
