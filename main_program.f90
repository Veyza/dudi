! This file is a part of DUDI, the Fortran-90 implementation 
! of the two-body model for dust dynamics
! Version 1.1.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! File: main_program.f90
! Description: A suggested template for application of the model

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

program main_program
	use const
	use define_types
	use integrator
	use inputdata
	use dataoutmod
	USE OMP_LIB
	
	implicit none
	! number of sources
	integer, parameter :: Ns = 1
	! number of points on the SC trajectory
	! for which the number density is to be calculated
	integer, parameter :: nt = 1
	real(8), parameter :: tnow = 600d0
	integer i_s, i, ii
	real density(nt,2), tmp_res(nt,2)
	type(source_properties) source(Ns)
	type(position_in_space) point(nt)
	
	call read_sources_params(source, Ns, './input_data_files/sources.dat')
	call read_spacecraft_coordinates(point, nt, './input_data_files/points.dat')
	
	density = 0d0

	do i_s = 1, Ns
		
		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, source, tmp_res)
		!$OMP DO
		do i = 1, nt
			call DUDI(tmp_res(i,:), point(i), source(i_s), tnow)
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		density = density + tmp_res
	enddo

	call result_out(density, nt, point)
	
	
	
end
