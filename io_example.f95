! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! **REFERENCE TO THE PAPER IN THE RIGHT FORMAT SO THAT IT COULD BE JUST COPIED AND PASTED TO THE REFERENCE LIST**

! File: io_example.f95
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
	real density(-Hpix:Hpix, -Vpix:Vpix, -ni:ni)
	real image(-Hpix:Hpix, -Vpix:Vpix), source_image(-Hpix:Hpix, -Vpix:Vpix)
	real tmp_res(-Hpix:Hpix, -Vpix:Vpix, -ni:ni, 2)
	type(position_in_space) points(-Hpix:Hpix, -Vpix:Vpix, -ni:ni)
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
				
			do i = -Hpix, Hpix
			do iii = -Vpix, Vpix
			
			!$OMP PARALLEL PRIVATE(ii) &
			!$OMP SHARED(nn, i, iii, sources, points)
			!$OMP DO
				! computing the 2nd moment of the dust density
				! in the nodes along the line of sight
				do ii = -ni, ni
					! if the point isn't on the LOS crossing the moon and
					! not on the LOS having 0-index like (i,0) or (0,iii)
					if(points(i,iii,ii)%compute) then
						! calculate the density
						call DUDI(tmp_res(i,iii,ii,:), points(i,iii,ii), sources(i_s), tnow)
					else
						tmp_res(i,iii,ii,:) = 0d0
					endif

				enddo		! ii = -nn, nn (over the line of sight)
			!$OMP END DO
			!$OMP END PARALLEL

			enddo	! iii = -Vpix, Vpix
			enddo	! i = -Hpix, Hpix
			
			! bound and unbound from the i_s source together
			density = tmp_res(:,:,:,1) + tmp_res(:,:,:,2)
			tmp_res = 0d0
			! create an image of this source
			call Integral_over_LoS(density, source_image)
			! add it to the resulting image
			image = image + source_image

		enddo		! i_s = Ns, 1, -1 (sources)
		do i = -Hpix, Hpix
			do ii = -Vpix, Vpix
				if(points(i,ii,0)%r >= rm) then
					image(i,ii) = image(i,ii) + bg
				endif
			enddo
		enddo

		call result_image_out(image, m)
	enddo
	
	 
end
