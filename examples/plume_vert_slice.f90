program plume_vert_slice
	use const
	use help
	use define_types
	use integrator
	use inputdata
	use dataoutmod
	use bgmod
	use gu
	USE OMP_LIB

	implicit none
	! number of sources
	integer, parameter :: Ns = 460
	integer, parameter :: Ndsources = 160
	integer, parameter :: Njets = (Ns - Ndsources) / 3
	real(8), parameter :: griddist = 2d3
	real(8), parameter :: latdist = asin(griddist/rm)
	real(8), parameter :: maxalphaM = pi
	real(8), parameter :: minalphaM = 155d0 * deg2rad
	real, parameter :: H2Omass = 2.998e-26
	! number of points on the SC trajectory for which the number density is to be calculated
	integer, parameter :: nt = 300
	real(8), parameter :: tnow = 0d0
	integer i_s, i, ii, nn
	real, allocatable, dimension(:) :: ttab
	real, allocatable, dimension(:,:) :: type1, type2, type3
	real, allocatable, dimension(:,:,:) :: tmp_res
	type(source_properties) source(Ns)
	type(position_in_space), allocatable, dimension(:) :: fltr
	type(position_in_space), allocatable, dimension(:,:) :: point
	character(len = 60) :: fname = './input_data_files/Enceladus_plume_Porco_jets2.dat'
	character(len = 21) :: resname = "./results/porco_vert/"
	character(len = 40) outfile
	character(len = 42) :: flybytr
	real(8) rmin, rmax, varfact, rminsalt, rmaxsalt
	real fnum
	character(len = 5) chfnum
	real(8) bg, production_salt_poor, production_salt_rich, production_salt_rich_jets
	real(8) mass_salt_poor, mass_salt_rich, mass, mass_salt_rich_jets, k, tmpmin, tmpmax
	real(8) sdnorm(3), alphatmp, alphatmp2, beta1, beta2

	call getarg(1, chfnum)
	read(chfnum,*) fnum


	alphatmp = pi * 0.8d0 !(90d0 + 56.506d0) * deg2rad !pi * 0.8d0
	beta1 = 45d0 * deg2rad !16.843 * deg2rad !
	alphatmp2 = alphatmp !(90d0 + 84.909d0) * deg2rad
	beta2 = 225d0 * deg2rad !12.593 * deg2rad !225d0 * deg2rad


	allocate(type1(nt,nt), type3(nt,nt), tmp_res(nt,nt,2), point(nt,nt))

	production_salt_poor = 0d0
	production_salt_rich = 0d0
	production_salt_rich_jets = 0d0
	sdnorm = 1d0

	call get_flyby_data(fnum, flybytr, outfile, rmin, rmax, rminsalt, rmaxsalt, varfact)

	if(fnum == 0.1) then
		write(*,*) "calculating dust distribution"
		if(p /= 3) then
			write(*,*) "wrong p parameter"
			stop
		endif
		if((rmin < Rg_upperlim)) then
			write(*,*) "check rmin, rmax params"
			stop
		endif
		call read_sources_params(source, Ns, fname, rmin, rmax, rminsalt, rmaxsalt, varfact)
        call get_jets(source(1:Njets), Njets, './input_data_files/vertical_jets.dat', &
                       1, varfact, rmin, rmax)
        call get_diffuse_sources(source(Njets+1:Ns), Ns-Njets, './input_data_files/diffuse_sources.dat', &
                            varfact, rmin, rmax)

		call get_flyby_plane(point, nt, &
		(/sin(alphatmp) * cos(beta1), sin(alphatmp) * sin(beta1), cos(alphatmp)/), &
		(/sin(alphatmp2) * cos(beta2), sin(alphatmp2) * sin(beta2), cos(alphatmp2)/), fnum)
	endif
	if(fnum == 0.4) then
		write(*,*) "calculating composition distribution"
		if(p /= 0) then
			write(*,*) "wrong p parameter"
			stop
		endif
		if((rmin < Rg_upperlim)) then
			write(*,*) "check rleft, rright params"
			stop
		endif
		allocate(type2(nt,nt))
		type2 = 0d0
		call read_sources_params(source, Ns, fname, rmin, rmax, rminsalt, rmaxsalt, varfact)
!~ 		call vert_sect(point, nt)
!~ 		call get_flyby_plane(point, nt, source(36)%rrM, source(60)%rrM, fnum)
		call get_flyby_plane(point, nt, &
		(/sin(alphatmp) * cos(beta1), sin(alphatmp) * sin(beta1), cos(alphatmp)/), &
		(/sin(alphatmp2) * cos(beta2), sin(alphatmp2) * sin(beta2), cos(alphatmp2)/), fnum)
	endif
	if(fnum == 0.2) then
		write(*,*) "calculating gas distribution"
		if((rmin < Rg_upperlim)) then
			write(*,*) "check rleft, rright params"
			stop
		endif
		if(p /= 0) then
			write(*,*) "wrong p parameter"
			stop
		endif
		call get_gas_plume(source, Ns, fname, rmin, rmax, sdnorm)
!~ 		call vert_sect(point, nt)
!~ 		call get_flyby_plane(point, nt, source(36)%rrM, source(60)%rrM, fnum)
		call get_flyby_plane(point, nt, &
		(/sin(alphatmp) * cos(beta1), sin(alphatmp) * sin(beta1), cos(alphatmp)/), &
		(/sin(alphatmp2) * cos(beta2), sin(alphatmp2) * sin(beta2), cos(alphatmp2)/), fnum)
	endif

	tmpmin = 0d0
	ii = 0
	do i = 1, Njets
		production_salt_poor = production_salt_poor + source(i)%production_rate
	enddo
	do i = Njets*2+1, Ns
		production_salt_rich = production_salt_rich + source(i)%production_rate
	enddo
	do i = Njets*2+1, Ns-Ndsources
		production_salt_rich_jets = production_salt_rich_jets + source(i)%production_rate
	enddo
	write(*,'(4(e10.4, x) f8.4)') production_salt_poor, production_salt_rich, &
	production_salt_rich_jets, production_salt_rich - production_salt_rich_jets

	if(fnum == 0.1 .or. fnum == 0.4) then
		call mass_production(mass, source(1)%sd, 0.1d0, Rg_upperlim)
		mass_salt_poor = mass
		call mass_production(mass, source(Njets+1)%sd, 0.1d0, Rg_upperlim)
		mass_salt_poor = mass_salt_poor + mass
		write(*,*) mass
		mass_salt_poor = mass_salt_poor * production_salt_poor * varfact

		call mass_production(mass, source(Njets*3+1)%sd, rmin, Rg_upperlim)
		mass_salt_rich = mass * production_salt_rich * varfact

		write(*,*) 'mass salt poor', mass_salt_poor
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-poor', mass_salt_poor, &
												'[kg/s]'
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich', mass_salt_rich, &
												'[kg/s]'
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich in jets', &
										mass * production_salt_rich_jets * varfact,  '[kg/s]'
	endif
	if(fnum == 0.2) then
		write(*,*) 'production rate of gas'
		write(*,*) 'jets', production_salt_poor * H2Omass / sdnorm(1), 'kg/s'
		write(*,*) 'diffuse sources', production_salt_rich * H2Omass , 'kg/s'
	endif

	type1 = 0d0 ; type3 = 0d0

	! salt-poor dust from jets
	do i_s = 1, Njets
		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, source, tmp_res)
		!$OMP DO
		do i = 1, nt
		do ii = 1, nt
			tmp_res(i,ii,:) = 0d0
			if(point(i,ii)%compute) then
				call DUDI(tmp_res(i,ii,:), point(i,ii), source(i_s), tnow)
			endif
			if(tmp_res(i,ii,1) /= tmp_res(i,ii,1)) then
				write(*,*) 'NaN obtained', tmp_res(i,ii,1)
				stop
			endif
		enddo
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		! add to the result the density of the dust from each source
		type1 = type1 + tmp_res(:,:,1) + tmp_res(:,:,2)
	enddo
	if(fnum == 0.4) then
		! same sources but different size-distribution; type2
		do i_s = Njets+1, 2*Njets
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, source, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), source(i_s), tnow)
				endif
				if(tmp_res(i,ii,1) /= tmp_res(i,ii,1)) then
					write(*,*) 'NaN obtained', tmp_res(i,ii,1)
					stop
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			! add to the result the density of the dust from each source
			type2 = type2 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
	endif
	if(fnum == 0.1 .or. fnum == 0.4) then
		! same sources but different size-distribution: type 3
		do i_s = 2*Njets+1, 3*Njets
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, source, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), source(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			! add to the result the density of the dust from each source
			type3 = type3 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
	endif
!~ 	! salt-rich dust from diffuse sources
	do i_s = Ns-Ndsources+1, Ns
		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, source, tmp_res)
		!$OMP DO
		do i = 1, nt
		do ii = 1, nt
			tmp_res(i,ii,:) = 0d0
			if(point(i,ii)%compute) then
				call DUDI(tmp_res(i,ii,:), point(i,ii), source(i_s), tnow)
			endif
		enddo
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		! add to the result the density of the dust from each source
		type3 = type3 + tmp_res(:,:,1) + tmp_res(:,:,2)
	enddo

	type1 = type1 * real(varfact, kind=kind(type1))
	type3 = type3 * real(varfact, kind=kind(type1))

	if(fnum == 0.2) then
		call vertical_slicematrix_out((type1+type3) * H2Omass, nt, fnum)
	endif
	if(fnum == 0.1) then
		call vertical_slicematrix_out(type1+type3, nt, fnum)
	endif
	if(fnum == 0.4) then
		call composition_matrix_out(type1, type2, type3, nt, fnum)
		deallocate(type2)
	endif

	deallocate(type1, type3, tmp_res, point)


end
