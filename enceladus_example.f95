! This file is a part of the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.1
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! File: enceladus_example.f95
! Description: The program calculates the number density of the dust from the Enceladus plume
!              at the points along the trajectory of the Cassini spacecraft
!              during its E2 flyby at Enceladus

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com
program enceladus_example
	use const
	use define_types
	use integrator
	use inputdata
	use dataoutmod
	USE OMP_LIB
	
	implicit none
	! number of sources
	integer, parameter :: Ns = 1
	! number of points on the SC trajectory for which the number density is to be calculated
	integer, parameter :: nt = 100
	real(8), parameter :: tnow = 0d0
	real, parameter :: bg = 0.01
	integer i_s, i, ii
	real density(nt,2), tmp_res(nt,2), ttab(nt)
	type(source_properties) source(Ns)
	type(position_in_space) point(nt)
	character(len = 44) :: fname = './input_data_files/Enceladus_jet.dat'
	
	call read_sources_params(source, Ns, fname)
	call read_Cassini_E2(point, ttab, nt)
	
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
		! add to the result the density of the dust from each source
		density = density + tmp_res
	enddo
	
	call cassini_flyby_out(density, ttab, bg, nt)
	
	
	
end
