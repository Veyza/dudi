program plume_horizontal_structure
	use const
	use help
	use define_types
	use integrator
	use inputdata
	use gu
	use bgmod
	use omp_lib

	implicit none

	integer, parameter :: Njets = 100
	integer, parameter :: Ndsources = 160

	real(8), parameter :: griddist = 3d3                ! spacing along surface [m]
	real(8), parameter :: maxalphaM = pi                ! 180 deg -> latitude -90 deg
	real(8), parameter :: minalphaM = 155d0 * deg2rad   ! latitude -65 deg

	integer :: nt, nlat, i, itmp, i_s
	real(8) :: latdist, altitude
	type(position_in_space), allocatable :: points(:)
	type(source_properties) :: jets(Njets), difsources(Ndsources)
	real, allocatable :: dens(:)
	real :: tmp_res(2)

	real(8) :: rmin, rmax, rminsalt, rmaxsalt, varfact
	real(8) :: production_salt_poor, production_salt_rich, production_salt_rich_jets
	real(8) :: mass_salt_poor, mass_salt_rich, mass
	real(8), parameter :: tnow = 0d0

	character(len = 64) :: arg
	integer :: alt_km
	character(len = 256) :: outname
	real(8) :: latdeg, londeg

	!----------------------------------------------------------------------
	! Read altitude (in meters) from command line
	!----------------------------------------------------------------------
	call getarg(1, arg)
	if(len_trim(arg) == 0) then
		write(*,*) 'Usage: ./bin/horizontal_structure <altitude_in_meters>'
		stop
	endif
	read(arg,*) altitude

	!----------------------------------------------------------------------
	! Check that we are in dust mass mode
	!----------------------------------------------------------------------
	if(p /= 3) then
		write(*,*) 'Error: const.f90 must have p = 3 for mass density.'
		stop
	endif

	!----------------------------------------------------------------------
	! Geometry of the surface grid: latitudes from -65 deg to -90 deg
	!----------------------------------------------------------------------
	latdist = griddist / rm

	nlat = int((maxalphaM - minalphaM) * rm / griddist) + 1
	nt = 1
	do i = 1, nlat
		itmp = int(twopi * sin(latdist * dble(i)) * rm / griddist) + 1
		nt = nt + itmp
	enddo

	allocate(points(nt))
	allocate(dens(nt))

	call get_surface_points(points, nt, griddist, latdist, maxalphaM, minalphaM, altitude)

	dens = 0.0

	!----------------------------------------------------------------------
	! Source setup: same as vertical_structure "mass" workflow (qid = 0.1)
	! Here we use the same parameter values that define_flyby_params would
	! set for qid = 0.1 (vertical plane case).
	!----------------------------------------------------------------------
	rmin     = 0.1d0
	rmax     = Rg_upperlim
	rminsalt = 0.1d0
	rmaxsalt = Rg_upperlim
	varfact  = 1d0

	production_salt_poor = 0d0
	production_salt_rich = 0d0
	production_salt_rich_jets = 0d0

	! type1: salt-poor dust from jets (sd=1)
	call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', 1, varfact, rmin, rmax)
	do i_s = 1, Njets
		production_salt_poor = production_salt_poor + jets(i_s)%production_rate
	enddo
	call mass_production(mass, 1, 0.1d0, Rg_upperlim)
	mass_salt_poor = mass * production_salt_poor * varfact
	write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-poor', mass_salt_poor, '[kg/s]'

	! Accumulate contribution from salt-poor jets
	do i_s = 1, Njets
		write(*,*) 'jet', i_s, '/', Njets
		!$OMP PARALLEL DO PRIVATE(i, tmp_res) SHARED(points, jets, dens, i_s)
		do i = 1, nt
			tmp_res = 0.0
			if(points(i)%compute) then
				call DUDI(tmp_res, points(i), jets(i_s), tnow)
			endif
			if(tmp_res(1) /= tmp_res(1)) then
				write(*,*) 'NaN obtained', tmp_res(1)
				stop
			endif
			dens(i) = dens(i) + tmp_res(1) + tmp_res(2)
		enddo
		!$OMP END PARALLEL DO
	enddo

	! type3: salt-rich dust from jets (sd=3) + diffuse sources
	call get_jets(jets, Njets, './input_data_files/vertical_jets.dat', 3, varfact, rminsalt, rmaxsalt)
	call get_diffuse_sources(difsources, Ndsources, './input_data_files/diffuse_sources.dat', varfact, rminsalt, rmaxsalt)

	do i_s = 1, Njets
		production_salt_rich_jets = production_salt_rich_jets + jets(i_s)%production_rate
	enddo
	production_salt_rich = production_salt_rich_jets
	do i_s = 1, Ndsources
		production_salt_rich = production_salt_rich + difsources(i_s)%production_rate
	enddo

	call mass_production(mass, 3, rmin, Rg_upperlim)
	mass_salt_rich = mass * production_salt_rich * varfact
	write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich', mass_salt_rich, '[kg/s]'
	write(*,'(A35, 2x, f10.3, 2x, A6)') 'overall mass of salt-rich in jets', &
		mass * production_salt_rich_jets * varfact, '[kg/s]'

	! Accumulate contribution from salt-rich jets
	do i_s = 1, Njets
		write(*,*) 'jet', i_s, '/', Njets
		!$OMP PARALLEL DO PRIVATE(i, tmp_res) SHARED(points, jets, dens, i_s)
		do i = 1, nt
			tmp_res = 0.0
			if(points(i)%compute) then
				call DUDI(tmp_res, points(i), jets(i_s), tnow)
			endif
			dens(i) = dens(i) + tmp_res(1) + tmp_res(2)
		enddo
		!$OMP END PARALLEL DO
	enddo

	! Accumulate contribution from diffuse sources
	do i_s = 1, Ndsources
		write(*,*) 'diffuse source', i_s, '/', Ndsources
		!$OMP PARALLEL DO PRIVATE(i, tmp_res) SHARED(points, difsources, dens, i_s)
		do i = 1, nt
			tmp_res = 0.0
			if(points(i)%compute) then
				call DUDI(tmp_res, points(i), difsources(i_s), tnow)
			endif
			dens(i) = dens(i) + tmp_res(1) + tmp_res(2)
		enddo
		!$OMP END PARALLEL DO
	enddo

	! Apply overall scaling factor, as in vertical_structure "mass"
	dens = dens * real(varfact, kind=kind(dens))

	!----------------------------------------------------------------------
	! Output: latitude [deg], longitude [deg], dust mass density
	!----------------------------------------------------------------------
	alt_km = nint(altitude / 1d3)
	write(outname,'(A,I0,A)') './results/surface_dust_mass_density_', alt_km, 'km.dat'

	open(111, file = trim(outname), status = 'replace')
	do i = 1, nt
		latdeg = 90d0 - points(i)%alpha / deg2rad
		londeg = points(i)%beta / deg2rad
		write(111,*) latdeg, londeg, dens(i)
	enddo
	close(111)

	deallocate(points, dens)

end program plume_horizontal_structure
