! This file is a part of DUDI, the Fortran-90 implementation 
! of the two-body model for dust dynamics
! Version 1.2.1
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! File: europa_example.f90
! Description: The program computes the surface deposition of dust for the points
!              at different distances from the source.
!              The calculations areperformed for 4 different sources.
!              The mass of dust produced during the specified time period
!              is computed for 2 different size distributions.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

program europa_example
	use const
	use define_types
	use distributions_fun
	use integrator
	use inputdata
	use gu
	use dataoutmod
	USE OMP_LIB
	
	implicit none
	real(8), parameter :: r1 = 0.2d0
	real(8), parameter :: r2 = 20d0
	real(8), parameter :: europa_orbital_period = 3.06822d+5		! seconds
	! number of sources
	integer, parameter :: Ns = 4
	! number of points on the SC trajectory for which the number density is to be calculated
	integer, parameter :: nt = 80
	integer i_s, i
	real massflux(nt,2)
	real(8), parameter :: tnow = 0d0
	real(8) dphi(nt), ddphi, tstep
	real(8) mass_shallow, mass_steep, m1, m2
	type(source_properties) source(Ns)
	type(position_in_space) point(nt)
	
	call get_europa_input(Ns, source, nt, point, dphi)

	call mass_production(m1, source(1)%sd, r1, r2)
	call mass_production(m2, source(2)%sd, r1, r2)
	! integrate the production rate over time

	mass_steep = source(1)%production_rate * m2
	mass_shallow = source(1)%production_rate * m1
	
	write(*,*) '   '
	write(*,'(A84,e10.3,x,A2)') 'with the shallow size distribution the total &
								mass produced in a second', mass_shallow, 'kg'
	write(*,'(A84,e10.3,x,A2)') 'with the steep size distribution the total &
								mass produced in a second  ', mass_steep, 'kg'
	write(*,*) '   '
								
	do i_s = 1, Ns
		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, source, massflux)
		!$OMP DO
		do i = 1, nt
			call DUDI(massflux(i,:), point(i), source(i_s), tnow)
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		call surface_deposition_out(i_s, massflux(:,1), nt, dphi)
	enddo

	
end
