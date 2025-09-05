! This file is a part of DUDI, the Fortran-90 implementation
! of the two-body model for dust dynamics
! Version 1.2.2
! This is free software. You can use and redistribute it
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt,
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186

! File: io_example.f90
! Description: The program constructs 4 images of a fictive volcano eruption on Io,
!              the images correspond to different stages of the eruption

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

program io_example
	use const
	use define_types
	use image_construction
	use integrator
	use inputdata
	use distributions_fun
	use Gu
	use help
	USE OMP_LIB
	implicit none
	integer, parameter :: Ns = 1                  ! number of sources
	real(8) tnow
	integer i, ii, iii, i_s, nn, itmp, m
	real density(-ni:ni,-Hpix:Hpix, -Vpix:Vpix)
	real image(-Hpix:Hpix, -Vpix:Vpix), source_image(-Hpix:Hpix, -Vpix:Vpix)
	real tmp_res(2, -ni:ni, -Hpix:Hpix, -Vpix:Vpix)
	type(position_in_space) points(-ni:ni, -Hpix:Hpix, -Vpix:Vpix)
	type(source_properties) :: sources(Ns)
	real, parameter :: bg = 1d-15
	real(8) moments(4)
	! moments of time at that the images are constructed
	call get_volcano_params(sources, Ns)

	tnow = 0d0
	do m = 1, 9
		tnow = tnow + 200d0
		image = 0.0
		do i_s = 1, Ns
			density = 0.0
			tmp_res = 0.0

			! form the integration grid for this specific plume
			call line_of_sight(points, sources(i_s))

			do iii = -Vpix, Vpix
			do i = -Hpix, Hpix

			!$OMP PARALLEL PRIVATE(ii) &
			!$OMP SHARED(nn, i, iii, sources, points)
			!$OMP DO
				! computing the 2nd moment of the dust density
				! in the nodes along the line of sight
				do ii = -ni, ni
					! if the point isn't on the LOS crossing the moon and
					! not on the LOS having 0-index like (i,0) or (0,iii)
					if(points(ii,i,iii)%compute) then
						! calculate the density
						call DUDI(tmp_res(:,ii,i,iii), points(ii,i,iii), sources(i_s), tnow)
					else
						tmp_res(:,ii,i,iii) = 0d0
					endif

				enddo		! ii = -nn, nn (over the line of sight)
			!$OMP END DO
			!$OMP END PARALLEL

			enddo	! i = -Hpix, Hpix
			enddo	! iii = -Vpix, Vpix

			! bound and unbound from the i_s source together
			density = tmp_res(1,:,:,:) + tmp_res(2,:,:,:)
			tmp_res = 0d0
			! create an image of this source
			call Integral_over_LoS(density, source_image)
			! add it to the resulting image
			image = image + source_image

		enddo		! i_s = Ns, 1, -1 (sources)
		do ii = -Vpix, Vpix
			do i = -Hpix, Hpix
				if(points(0,i,ii)%r >= rm) then
					image(i,ii) = image(i,ii) + bg
				endif
			enddo
		enddo

		call result_image_out(image, m)
	enddo


end
