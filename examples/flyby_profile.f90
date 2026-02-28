! This file is a part of DUDI, the Fortran-90 implementation
! of the two-body model for dust dynamics
! Version 1.2.3
! This is free software. You can use and redistribute it
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt,
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186

! File: flyby_profile.f90
! Description: The fundamental constants and the numerical parameters
!              that are used by many subroutines

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

program main_program
	use const
	use define_types
	use integrator
	use inputdata
	use dataoutmod
	use bgmod
	USE OMP_LIB

	implicit none
	! number of sources
	integer, parameter :: Njets = 100
	integer, parameter :: Nds = 160
	! number of points on the SC trajectory
	! for which the number density is to be calculated
	integer, parameter :: nt = 100
	real(8), parameter :: tnow = 600d0
	real(8) varfact, rlim1, rlim2, rlim1salt, rlim2salt, bg
	real ttab(nt)
	integer i_s, i, ii, ty
	real tmp_res(nt,2), type1(nt), type2(nt), type3(nt)
	type(source_properties) difsources(Nds), jets(Njets)
	type(position_in_space) point(nt)
	character(len = 5) chfnum
	character(len = 40) resname, fname
	real fnum

	call getarg(1, chfnum)
	read(chfnum,*) fnum

	call define_flyby_params(fnum, fname, rlim1, rlim2, rlim1salt, rlim2salt, varfact)
	call read_Cassini_flyby(point, ttab, nt, fname)

	! calculating type 1 number density
	if(int(fnum) == 5 .or. int(fnum) == 17) then
		ty = int(fnum) * 10 + 1
		write(*,*) 'sd', ty
	else
		ty = 1
	endif
	call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', ty, varfact, rlim1, rlim2)

	type1 = 0.0

	do i_s = 1, Njets

		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, jets, tmp_res)
		!$OMP DO
		do i = 1, nt
			call DUDI(tmp_res(i,:), point(i), jets(i_s), tnow)
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		type1 = type1 + tmp_res(:,1) + tmp_res(:,2)
	enddo

	! calculating type 2 number density
	type2 = 0.0
	if(int(fnum) == 5 .or. int(fnum) == 17) then
		ty = int(fnum) * 10 + 2
		write(*,*) 'sd', ty
		call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', ty, varfact, rlim1, rlim2)

		do i_s = 1, Njets

			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, jets, tmp_res)
			!$OMP DO
			do i = 1, nt
				call DUDI(tmp_res(i,:), point(i), jets(i_s), tnow)
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type2 = type2 + tmp_res(:,1) + tmp_res(:,2)
		enddo
	endif
	! calculating type 3 number density
	call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', 3, varfact, rlim1salt, rlim2salt)

	type3 = 0.0

	do i_s = 1, Njets

		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, jets, tmp_res)
		!$OMP DO
		do i = 1, nt
			call DUDI(tmp_res(i,:), point(i), jets(i_s), tnow)
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		type3 = type3 + tmp_res(:,1) + tmp_res(:,2)
	enddo

	call get_diffuse_sources(difsources, Nds, './input_data_files/diffuse_sources.dat', varfact, rlim1salt, rlim2salt)

	do i_s = 1, Nds

		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, difsources, tmp_res)
		!$OMP DO
		do i = 1, nt
			call DUDI(tmp_res(i,:), point(i), difsources(i_s), tnow)
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		type3 = type3 + tmp_res(:,1) + tmp_res(:,2)
	enddo

	! taking into account the E ring background
   	bg = integrateBG(rlim1, rlim2, 1.6d0, Rg_upperlim, 0.03d0)
   	type1 = type1 + real(bg * 0.6d0, kind=kind(type1))
   	type2 = type2 + real(bg * 0.33d0, kind=kind(type2))
   	type3 = type3 + real(bg * 0.07d0, kind=kind(type3))

	resname = './results/E' // trim(chfnum) // '_profile.dat'

	call cassini_flyby_out(type1, type2, type3, ttab, nt, resname)



end
