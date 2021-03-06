! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! **REFERENCE TO THE PAPER IN THE RIGHT FORMAT SO THAT IT COULD BE JUST COPIED AND PASTED TO THE REFERENCE LIST**

! File: const.f95
! Description: The fundamental constants and the numerical parameters
!              that are used by many subroutines

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module const
	implicit none
		
		! fundamental constants and auxilary numbers
		real(8), parameter :: pi = 3.141592653589793d0
		real(8), parameter :: halfpi = pi / 2d0
		real(8), parameter :: sqrtpi = sqrt(pi)
		real(8), parameter :: twopi = 2d0 * pi
		real(8), parameter :: sqrt2d0 = sqrt(2d0)
		real(8), parameter :: deg2rad = pi / 1.8d+2
		real(8), parameter :: rad2deg = 1.8d+2 / pi
		real(8), parameter :: gravity_constant = 6.674d-11				! m^3 / kg / s^2
		
		! parameters defining the moon
		real(8), parameter :: moon_mass = 8.94d+22		! kg
		real(8), parameter :: gm = gravity_constant * moon_mass
		real(8), parameter :: rm = 1.8216d+6		! meters
		real(8), parameter :: vesc = sqrt(gm * 2d0 / rm)
				
		! other parameters defining the quantity of interest
		! density of the dust particles (if one wants to compute mass)
		! if the quantity which we want to compute is dust flux through
		! the surface parallel to the moon surface
		! parameter flux should be .TRUE., if it is .FALSE. then density is computed
		real(8), parameter :: rho = 920d0		! kg/m^3 
		logical, parameter :: flux = .FALSE.
		
		! parameters of Gu function
		! p = 0 -- number density is computed, 1 -- mean radius,
		! 2 -- cross section, 3 -- mass density
		integer, parameter :: p = 2
		! lower boundary for function Gu(rmin, rmax), microns
		real(8), parameter :: rmin = 0.2d0
		! upper boundary for function Gu(rmin, rmax), microns
		real(8), parameter :: rmax = 0.4d0
		
		! parameters controling accuracy
		! number of precalculated values of GR(u) integral
		! <=> u/u_gas goes from 0 to 1 with step 1/GRN
		integer, parameter :: GRN = 3
		! order of integration G(R,u) over R to obtain GR(u)
		integer, parameter :: order_R = 10
		! order of Gaussian quadrature for integration
		! of n(r, alpha, beta, v, theta, lambda) over velocity (v)
		! separately for bound and unbound particles
		integer, parameter :: order_v_el = 5
		integer, parameter :: order_v_hy = 5
		
	

			
end module const
