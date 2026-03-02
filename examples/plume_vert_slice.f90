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
	! number of sources (same as flyby_profile)
	integer, parameter :: Njets = 100
	integer, parameter :: Ndsources = 160
	real(8), parameter :: griddist = 2d3
	real(8), parameter :: latdist = asin(griddist/rm)
	real(8), parameter :: maxalphaM = pi
	real(8), parameter :: minalphaM = 155d0 * deg2rad
	real, parameter :: H2Omass = 2.998e-26
	! number of points in the plane
	integer, parameter :: nt = 300
	real(8), parameter :: tnow = 0d0
	integer i_s, i, ii
	real, allocatable, dimension(:,:) :: type1, type3
	real, allocatable, dimension(:,:,:) :: tmp_res
	type(source_properties) jets(Njets), difsources(Ndsources)
	type(position_in_space), allocatable, dimension(:,:) :: point
	character(len = 42) :: flybytr
	real(8) rmin, rmax, varfact, rminsalt, rmaxsalt
	real fnum
	character(len = 5) chfnum
	real(8) production_salt_poor, production_salt_rich, production_salt_rich_jets
	real(8) mass_salt_poor, mass_salt_rich, mass
	real(8) alphatmp, alphatmp2, beta1, beta2

	! what are we calculating? 0.1: dust mass distribution, 0.2: gas distribution, 0.4: composition distribution
	call getarg(1, chfnum)
	read(chfnum,*) fnum

	allocate(type1(nt,nt), type3(nt,nt), tmp_res(nt,nt,2), point(nt,nt))

	! setup flyby/plane parameters (rmin, rmax, etc.)
	call define_flyby_params(fnum, flybytr, rmin, rmax, rminsalt, rmaxsalt, varfact)

	! polar angle and longitude of the plane where we calculate the density
	alphatmp = pi * 0.8d0
	beta1 = 45d0 * deg2rad
	alphatmp2 = alphatmp
	beta2 = 225d0 * deg2rad
	call get_flyby_plane(point, nt, &
	(/sin(alphatmp) * cos(beta1), sin(alphatmp) * sin(beta1), cos(alphatmp)/), &
	(/sin(alphatmp2) * cos(beta2), sin(alphatmp2) * sin(beta2), cos(alphatmp2)/), fnum)

	type1 = 0d0
	type3 = 0d0
	production_salt_poor = 0d0
	production_salt_rich = 0d0
	production_salt_rich_jets = 0d0

	if(fnum == 0.1) then
		write(*,*) "calculating dust distribution"
		if(p /= 3) then
			write(*,*) "wrong p parameter"
			stop
		endif
		! type1: salt-poor dust from jets (sd=1)
		call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', 1, varfact, rmin, rmax)
		production_salt_poor = 0d0
		do i = 1, Njets
			production_salt_poor = production_salt_poor + jets(i)%production_rate
		enddo
		call mass_production(mass, 1, 0.1d0, Rg_upperlim)
		mass_salt_poor = mass * production_salt_poor * varfact
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-poor', mass_salt_poor, '[kg/s]'
		do i_s = 1, Njets
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, jets, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), jets(i_s), tnow)
				endif
				if(tmp_res(i,ii,1) /= tmp_res(i,ii,1)) then
					write(*,*) 'NaN obtained', tmp_res(i,ii,1)
					stop
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type1 = type1 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
		! type3: salt-rich dust from jets (sd=3) + diffuse sources
		call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', 3, varfact, rminsalt, rmaxsalt)
		call get_diffuse_sources(difsources, Ndsources, './input_data_files/diffuse_sources.dat', varfact, rminsalt, rmaxsalt)

		do i = 1, Njets
			production_salt_rich_jets = production_salt_rich_jets + jets(i)%production_rate
		enddo
		production_salt_rich = production_salt_rich_jets
		do i = 1, Ndsources
			production_salt_rich = production_salt_rich + difsources(i)%production_rate
		enddo
		call mass_production(mass, 3, rmin, Rg_upperlim)
		mass_salt_rich = mass * production_salt_rich * varfact
		write(*,*) mass
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich', mass_salt_rich, '[kg/s]'
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich in jets', &
			mass * production_salt_rich_jets * varfact, '[kg/s]'
		do i_s = 1, Njets
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, jets, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), jets(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type3 = type3 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
		do i_s = 1, Ndsources
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, difsources, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), difsources(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type3 = type3 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo

	else if(fnum == 0.4) then
		write(*,*) "calculating composition distribution"
		if(p /= 0) then
			write(*,*) "wrong p parameter"
			stop
		endif
		! type1: salt-poor dust from jets (sd=1)
		call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', 1, varfact, rmin, rmax)
		production_salt_poor = 0d0
		do i = 1, Njets
			production_salt_poor = production_salt_poor + jets(i)%production_rate
		enddo
		call mass_production(mass, 1, 0.1d0, Rg_upperlim)
		mass_salt_poor = mass * production_salt_poor * varfact
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-poor', mass_salt_poor, '[kg/s]'
		do i_s = 1, Njets
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, jets, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), jets(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type1 = type1 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
		! type3: salt-rich dust from jets (sd=3) + diffuse (no type2)
		call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', 3, varfact, rminsalt, rmaxsalt)
		call get_diffuse_sources(difsources, Ndsources, './input_data_files/diffuse_sources.dat', varfact, rminsalt, rmaxsalt)
		do i = 1, Njets
			production_salt_rich_jets = production_salt_rich_jets + jets(i)%production_rate
		enddo
		production_salt_rich = production_salt_rich_jets
		do i = 1, Ndsources
			production_salt_rich = production_salt_rich + difsources(i)%production_rate
		enddo
		call mass_production(mass, 3, rmin, Rg_upperlim)
		mass_salt_rich = mass * production_salt_rich * varfact
		write(*,*) mass
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich', mass_salt_rich, '[kg/s]'
		write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich in jets', &
			mass * production_salt_rich_jets * varfact, '[kg/s]'
		do i_s = 1, Njets
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, jets, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), jets(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type3 = type3 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
		do i_s = 1, Ndsources
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, difsources, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), difsources(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type3 = type3 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo

	else if(fnum == 0.2) then
		write(*,*) "calculating gas distribution"
		if(p /= 0) then
			write(*,*) "wrong p parameter"
			stop
		endif
		call get_gas_jets(jets, Njets, './input_data_files/vertical_jets.dat', rmin, rmax)
		call get_gas_diffuse_sources(difsources, Ndsources, './input_data_files/diffuse_sources.dat', rmin, rmax)
		do i = 1, Njets
			production_salt_poor = production_salt_poor + jets(i)%production_rate
		enddo
		do i = 1, Ndsources
			production_salt_rich = production_salt_rich + difsources(i)%production_rate
		enddo
		write(*,*) 'production rate of gas'
		write(*,*) 'jets', production_salt_poor * H2Omass, 'kg/s'
		write(*,*) 'diffuse sources', production_salt_rich * H2Omass, 'kg/s'
		do i_s = 1, Njets
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, jets, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), jets(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type1 = type1 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
		do i_s = 1, Ndsources
			!$OMP PARALLEL PRIVATE(i) &
			!$OMP SHARED(i_s, point, difsources, tmp_res)
			!$OMP DO
			do i = 1, nt
			do ii = 1, nt
				tmp_res(i,ii,:) = 0d0
				if(point(i,ii)%compute) then
					call DUDI(tmp_res(i,ii,:), point(i,ii), difsources(i_s), tnow)
				endif
			enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			type3 = type3 + tmp_res(:,:,1) + tmp_res(:,:,2)
		enddo
	endif

	type1 = type1 * real(varfact, kind=kind(type1))
	type3 = type3 * real(varfact, kind=kind(type1))

	if(fnum == 0.2) then
		call vertical_slicematrix_out((type1+type3) * H2Omass, nt, fnum)
	else if(fnum == 0.1) then
		call vertical_slicematrix_out(type1+type3, nt, fnum)
	else if(fnum == 0.4) then
		call composition_matrix_out(type1, type3, nt, fnum)
	endif

	deallocate(type1, type3, tmp_res, point)

end program plume_vert_slice
