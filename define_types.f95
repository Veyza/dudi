! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! File: define_types.f95
! Description: Definition of the structures used by the routines

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module define_types
	use const
	implicit none
	
	type ejection_speed_properties
		integer ud_shape                ! parameter defining which distribution is used for ejection speed
		real(8) umax                    ! gas velocity
		real(8) umin                    ! a parameter in velocity distribution
	end type ejection_speed_properties

	type source_properties              ! parameters of dust ejection
		real(8) rrM(3)                    ! Cartesian coordinates of a point
		                                  ! source in the moon-centered coordinate system
		real(8) r
		real(8) alphaM                        ! polar angle of the point source
		real(8) betaM                         ! eastern longitude of the point source
		real(8) zeta                          ! zenith angle of the axis around which ejection is symmetrical
		real(8) eta                           ! azimuth of this axis (counted from the local North, clockwise)
		real(8) symmetry_axis(3)              ! unit vector in moon-centered coordinate system pointing to the direction of the axis around which ejection is symmetrical
		type(ejection_speed_properties) ud    ! parameters of ejection speed distribution
		integer ejection_angle_distr          ! parameter defining which ejection angle distribution is used; 1 -- Gaussian, 2 -- uniform cone, 3 -- parabola inside a cone
		integer sd                            ! parameter to select ejected dust size distribution
		real(8) ui(GRN)                       ! interpolation grid for GRm(u,Rmin,Rmax) precalculation
		real(8) Gu_precalc(GRN)               ! Gu(Rmin,Rmax)
		integer production_fun                ! parameter used to select a function for production rate (if <= 0, the production rate is constant)
		real(8) production_rate               ! parameter used in definition of the function for production rate (usually the normalization factor)
		logical is_jet                        ! .TRUE. if the ejection is concentrated (omega is small)
	end type source_properties
	
	type position_in_space             ! where the dust density is to be calculated
		real(8) r                        ! distance from the center of the moon
		real(8) r_scaled                 ! distance from the center of the moon divided by the moon radius
		real(8) alpha                    ! polar angle
		real(8) beta                     ! eastern longitude
		real(8) rvector(3)               ! Cartesian coordinates of the point 
		logical compute	                 ! marker that density should or should not be conputed in this point
	end type position_in_space


end module define_types
