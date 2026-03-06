! This file is a part of DUDI, the Fortran-90 implementation
! of the two-body model for dust dynamics
! Version 1.2.3
! This is free software. You can use and redistribute it
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt,
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186

! File: inputdata.f90
! Description: The subroutines suplying the parameters of sources, and
!              spacecraft positions, and the auxiliary functions used by them

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module inputdata
	implicit none
		contains

			! input the spacecraft positions at the defined moments
			! from the file fname
			! the file was prepared in advance using SPICE software
			subroutine read_Cassini_flyby(points, ttab, nt, fname)
				use const
				use define_types
				implicit none
				integer, intent(in) :: nt
				type(position_in_space), intent(out) :: points(nt)
				real, intent(out) :: ttab(nt)
				integer i
				character(*), intent(in) :: fname

				open(200, file = fname, status = 'old')
					do i = 1, nt
						read(200,*) ttab(i), points(i)%r, points(i)%alpha, points(i)%beta
						points(i)%alpha = points(i)%alpha * deg2rad
						points(i)%alpha = halfpi - points(i)%alpha
						points(i)%beta = points(i)%beta * deg2rad
						points(i)%rvector(1) = points(i)%r * sin(points(i)%alpha) &
						                                   * cos(points(i)%beta)
						points(i)%rvector(2) = points(i)%r * sin(points(i)%alpha) &
						                                   * sin(points(i)%beta)
						points(i)%rvector(3) = points(i)%r * cos(points(i)%alpha)
						points(i)%r_scaled = points(i)%r / rm
						points(i)%compute = .TRUE.
					enddo
				close(200)

			end subroutine read_Cassini_flyby




			! input the parameters of the sources from the file with the name
			! stored in the variable fname
			subroutine get_diffuse_sources(sources, Ns, fname, k, rlim1, rlim2)
				use const
				use define_types
				use gu
				implicit none
				integer, intent(in) :: Ns
				type(source_properties), intent(out) :: sources(Ns)
				real(8), intent(in) :: k
				real(8), intent(in) :: rlim1, rlim2
				integer i
				character(*), intent(in) :: fname
				character(len = 5) header

				open(100, file = fname, status = 'old')
					read(100,*) header
					do i = 1, Ns
						read(100,*) sources(i)%alphaM, sources(i)%betaM, &
									sources(i)%zeta, sources(i)%eta
						sources(i)%r = rm
						sources(i)%alphaM = sources(i)%alphaM * deg2rad
						sources(i)%betaM = sources(i)%betaM * deg2rad
						sources(i)%alphaM = halfpi - sources(i)%alphaM

						sources(i)%zeta = sources(i)%zeta * deg2rad
						sources(i)%eta = sources(i)%eta * deg2rad

						sources(i)%rrM(1) = rm * sin(sources(i)%alphaM) &
						                       * cos(sources(i)%betaM)
						sources(i)%rrM(2) = rm * sin(sources(i)%alphaM) &
						                       * sin(sources(i)%betaM)
						sources(i)%rrM(3) = rm * cos(sources(i)%alphaM)

						sources(i)%is_jet = .FALSE.

						sources(i)%sd = 3

						sources(i)%ud%ud_shape = 3
						sources(i)%ud%umin = 0d0		! m/s
						sources(i)%ud%umax = 400d0		! m/s

						sources(i)%ejection_angle_distr = 3

						sources(i)%production_fun = 0
						sources(i)%production_rate = 4d11 * k

						call jet_direction(sources(i)%betaM, sources(i)%zeta, &
						     sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)

						call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, &
						               sources(i)%sd, sources(i)%ud, rlim1, rlim2)

					enddo
				close(100)

			end subroutine get_diffuse_sources


			! gas jets: 100 jets with gas-specific sd, ud, and production rate
			subroutine get_gas_jets(jets, Njets, fname)
				use const
				use define_types
				use gu
				implicit none
				integer, intent(in) :: Njets
				type(source_properties), intent(out) :: jets(Njets)
				integer i
				character(*), intent(in) :: fname
				character(5) header

				open(100, file = fname, status = 'old')
					read(100,*) header
					do i = 1, Njets
						read(100,*) jets(i)%alphaM, jets(i)%betaM, &
									jets(i)%zeta, jets(i)%eta
						jets(i)%r = rm
						jets(i)%alphaM = jets(i)%alphaM * deg2rad
						jets(i)%betaM = jets(i)%betaM * deg2rad
						jets(i)%alphaM = halfpi - jets(i)%alphaM

						jets(i)%rrM(1) = rm * sin(jets(i)%alphaM) &
						* cos(jets(i)%betaM)
						jets(i)%rrM(2) = rm * sin(jets(i)%alphaM) &
												* sin(jets(i)%betaM)
						jets(i)%rrM(3) = rm * cos(jets(i)%alphaM)

						jets(i)%zeta = jets(i)%zeta * deg2rad
						jets(i)%eta = jets(i)%eta * deg2rad

						jets(i)%is_jet = .TRUE.

						jets(i)%ud%umax = 800d0		! m/s

						jets(i)%ejection_angle_distr = 1

						jets(i)%ud%ud_shape = 0
						jets(i)%ud%umin = jets(i)%ud%umax - 20d0
						jets(i)%sd = 0    ! uniform, non-normalized fR = 1d0
						jets(i)%production_rate = 1.000667d26

						call Gu_integral(jets(i)%ui, jets(i)%Gu_precalc, &
						               jets(i)%sd, jets(i)%ud, 1d0, 1d0)

						! normalization to account for dimension when treating gas molecules
						! as dust particles with a size distribution between 0.1 and Rg_upperlim
						jets(i)%Gu_precalc = jets(i)%Gu_precalc * (Rg_upperlim - 0.1d0)

						call jet_direction(jets(i)%betaM, jets(i)%zeta, &
						     jets(i)%eta, jets(i)%rrM, jets(i)%symmetry_axis)
					enddo
					close(100)
			end subroutine get_gas_jets



			! gas diffuse sources: 160 sources with gas-specific sd, ud, and production rate
			subroutine get_gas_diffuse_sources(sources, Ndsources, fname)
				use const
				use define_types
				use gu
				implicit none
				integer, intent(in) :: Ndsources
				type(source_properties), intent(out) :: sources(Ndsources)
				integer i
				character(*), intent(in) :: fname
				character(len = 5) header

				open(100, file = fname, status = 'old')
					read(100,*) header
					do i = 1, Ndsources
						read(100,*) sources(i)%alphaM, sources(i)%betaM, &
									sources(i)%zeta, sources(i)%eta
						sources(i)%r = rm

						sources(i)%alphaM = sources(i)%alphaM * deg2rad
						sources(i)%betaM = sources(i)%betaM * deg2rad
						sources(i)%alphaM = halfpi - sources(i)%alphaM

						sources(i)%rrM(1) = rm * sin(sources(i)%alphaM) &
						* cos(sources(i)%betaM)
						sources(i)%rrM(2) = rm * sin(sources(i)%alphaM) &
												* sin(sources(i)%betaM)
						sources(i)%rrM(3) = rm * cos(sources(i)%alphaM)

						sources(i)%is_jet = .FALSE.

						sources(i)%sd = 0	! uniform, non-normalized fR = 1d0
						sources(i)%ud%umax = 400d0		! m/s

						sources(i)%ejection_angle_distr = 3


						sources(i)%ud%ud_shape = 0
						sources(i)%ud%umin = sources(i)%ud%umax - 20d0
						sources(i)%production_rate = 6.25417d24

						call jet_direction(sources(i)%betaM, sources(i)%zeta, &
						     sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)

						call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, &
						               sources(i)%sd, sources(i)%ud, 1d0, 1d0)
						! normalization to account for dimension when treating gas molecules
						! as dust particles with a size distribution between 0.1 and Rg_upperlim
						sources(i)%Gu_precalc = sources(i)%Gu_precalc * (Rg_upperlim - 0.1d0)
					enddo
				close(100)
			end subroutine get_gas_diffuse_sources



			subroutine get_flyby_plane(points, nt, r1, r2, fnum, cellsize)
				use const
				use help
				use define_types
				implicit none
				integer, intent(in) :: nt
				type(position_in_space), intent(out) :: points(nt, nt)
				real, intent(in) :: fnum
				real(8), intent(in) :: r1(3), r2(3), cellsize
				integer i, ii
				real(8) xtmp(3), ytmp(3), ntmp(3)

				if(fnum == 0.2 .or. fnum == 0.1 .or. fnum == 0.4) then
					! vector normal to the plane containing r1 and r2
					ntmp = -vector_product(r1, r2)
					ntmp = ntmp / norma3d(ntmp)
					! projection of the polar axis to the plane containing r1 and r2
					ytmp = (/0d0, 0d0, 1d0/) - dot_product((/0d0, 0d0, 1d0/), ntmp) * ntmp
					ytmp = -ytmp / norma3d(ytmp) * cellsize
					xtmp = vector_product(ytmp, ntmp)
					do i = 1, nt
					do ii = 1, nt
						if(cellsize < 2.9d3) then
							points(i,ii)%rvector = (nt/2-i-0.5) * xtmp + ytmp/cellsize*222.0d3 &
							+ (ii) * ytmp + ntmp * 0d0
						else
							points(i,ii)%rvector = (nt/2-i-0.5) * xtmp &
							+ (ii) * ytmp + ntmp * 0d0
						endif

						points(i,ii)%r = norma3d(points(i,ii)%rvector)
						points(i,ii)%alpha = acos(points(i,ii)%rvector(3) / points(i,ii)%r)
						points(i,ii)%beta = myatan1(points(i,ii)%rvector(1), points(i,ii)%rvector(2))
						points(i,ii)%r_scaled = points(i,ii)%r / rm
						points(i,ii)%compute = points(i,ii)%r_scaled >= 1d0
					enddo
					enddo
				endif


			end subroutine get_flyby_plane




			subroutine get_surface_points(points, nt, griddist, latdist, maxalphaM, minalphaM, altitude)
				use const
				use help
				use define_types
				implicit none
				real(8), intent(in) :: griddist, latdist, maxalphaM, minalphaM, altitude
				integer, intent(in) :: nt
				type(position_in_space), intent(out) :: points(nt)
				real(8) mincos, maxcos, difcos
				real(8) tmp, tmp2
				integer i, ii, iii, itmp

				points(1)%alpha = 179d0 * deg2rad
				points(1)%beta = 0d0
				iii = 2
				write(*,*) 'n lat circles =', int((maxalphaM - minalphaM) * rm / griddist) + 1
				write(*,*) 'nn =', nt, points(nt)%alpha
				do i = 1, int((maxalphaM - minalphaM) * rm / griddist) + 1
					itmp = int(twopi * sin(latdist * i) * rm / griddist) + 1
!~ 					write(*,*) i, itmp, iii
					do ii = 1, itmp
						points(iii)%alpha = pi - latdist * i
						points(iii)%beta = ii * twopi / itmp
						iii = iii + 1
					enddo
				enddo

				do i = 1, nt
					points(i)%rvector(3) = cos(points(i)%alpha)
					points(i)%rvector(1) = sin(points(i)%alpha) * cos(points(i)%beta)
					points(i)%rvector(2) = sin(points(i)%alpha) * sin(points(i)%beta)
					points(i)%rvector = points(i)%rvector * (rm + altitude)
					points(i)%r = norma3d(points(i)%rvector)
					points(i)%r_scaled = points(i)%r / rm
					points(i)%compute = .TRUE.
				enddo

			end subroutine get_surface_points



			subroutine define_flyby_params(fnum, flybytr, rmin, rmax, rminsalt, rmaxsalt, varfact)
				use const
				use distributions_fun
				real, intent(in) :: fnum
				character(*), intent(out) :: flybytr
				real(8), intent(out) :: rmin, rmax, rminsalt, rmaxsalt, varfact
				real(8) fR

				if(fnum == 0.1 .or. fnum == 0.2 .or. fnum == 0.4) then
					flybytr = './input_data_files/vertical_plane.dat'
					rmin = 0.1d0
					rmax = Rg_upperlim
					rminsalt = 0.1d0
					rmaxsalt = Rg_upperlim
					varfact = 1d0
				endif
				if(fnum == 5.2) then
					flybytr = './input_data_files/Cassini_E5_flyby.dat'
					rmin = 0.2d0
					rmax = 1.7d0
					rminsalt = 0.1d0
					rmaxsalt = 1.1d0
					varfact = 2d0
					Rleft = 0.3d0
					Rright = 0.5d0
					fR = size_distribution(Rright, 1)
					acoeflin = fR / (Rright - Rleft)
					bcoeflin = -acoeflin * Rleft
					write(*,*) 'Rleft=', Rleft
					write(*,*) 'Rright=', Rright
					write(*,*) 'a=', acoeflin
					write(*,*) 'b=', bcoeflin
					write(*,*) 'E5 profile corrected, k =', varfact
				endif
				if(fnum == 5) then
					flybytr = './input_data_files/Cassini_E5_flyby.dat'
					rmin = 0.2d0
					rmax = 1.7d0
					rminsalt = 0.1d0
					rmaxsalt = 1.1d0
					varfact = 1d0
					Rleft = 0.3d0
					Rright = 0.5d0
					fR = size_distribution(Rright, 1)
					acoeflin = fR / (Rright - Rleft)
					bcoeflin = -acoeflin * Rleft
					write(*,*) 'Rleft=', Rleft
					write(*,*) 'Rright=', Rright
					write(*,*) 'a=', acoeflin
					write(*,*) 'b=', bcoeflin
					write(*,*) 'E5 profile corrected, k =', varfact
				endif
				if(fnum == 17.17) then
					rmin = 0.8d0
					rmax = 6.2d0
					rminsalt = 0.5d0
					rmaxsalt = 4.1d0
					flybytr = './input_data_files/Cassini_E17_flyby.dat'
					varfact = 1.7d0
					Rleft = 0.9d0
					Rright = 1.1d0
					fR = size_distribution(Rright, 1)
					acoeflin = fR / (Rright - Rleft)
					bcoeflin = -acoeflin * Rleft
					write(*,*) 'Rleft=', Rleft
					write(*,*) 'Rright=', Rright
					write(*,*) 'a=', acoeflin
					write(*,*) 'b=', bcoeflin
					write(*,*) 'E17 profile corrected, k =', varfact
				endif
				if(fnum == 17) then
					rmin = 0.8d0
					rmax = 6.2d0
					rminsalt = 0.5d0
					rmaxsalt = 4.1d0
					flybytr = './input_data_files/Cassini_E17_flyby.dat'
					varfact = 1.0d0
					Rleft = 0.9d0
					Rright = 1.1d0
					fR = size_distribution(Rright, 1)
					acoeflin = fR / (Rright - Rleft)
					bcoeflin = -acoeflin * Rleft
					write(*,*) 'Rleft=', Rleft
					write(*,*) 'Rright=', Rright
					write(*,*) 'a=', acoeflin
					write(*,*) 'b=', bcoeflin
					write(*,*) 'E17 profile corrected, k =', varfact
				endif
				if(fnum == 7.1) then
					flybytr = './input_data_files/Cassini_E7_flyby.dat'
					rmin = 1.75d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 7.2) then
					flybytr = './input_data_files/Cassini_E7_flyby.dat'
					rmin = 3.1d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 7.3) then
					flybytr = './input_data_files/Cassini_E7_flyby.dat'
					rmin = 6.7d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 3.1) then
					flybytr = './input_data_files/Cassini_E3_flyby.dat'
					rmin = 0.1d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 3.2) then
					flybytr = './input_data_files/Cassini_E3_flyby.dat'
					rmin = 0.5d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 3.3) then
					flybytr = './input_data_files/Cassini_E3_flyby.dat'
					rmin = 1d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 21.1) then
					flybytr = './input_data_files/Cassini_E21_flyby.dat'
					rmin = 1.6d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 21.2) then
					flybytr = './input_data_files/Cassini_E21_flyby.dat'
					rmin = 2.85d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif
				if(fnum == 21.3) then
					flybytr = './input_data_files/Cassini_E21_flyby.dat'
					rmin = 6.2d0
					rmax = Rg_upperlim
					rminsalt = rmin
					rmaxsalt = rmax
					varfact = 1d0
				endif

			end subroutine define_flyby_params




			function cylindrical2point(rh, theta, z) result(point)
				use const
				use help
				use define_types
				type(position_in_space) point
				real rh, theta, z

				point%r = sqrt(rh*rh + z * z)
				point%alpha = acos(z / point%r)
				point%beta = theta
				point%rvector(1) = point%r * sin(point%alpha) &
						                   * cos(point%beta)
				point%rvector(2) = point%r * sin(point%alpha) &
						                   * sin(point%beta)
				point%rvector(3) = z
				point%r_scaled = point%r / rm
				point%compute = .TRUE.


			end function cylindrical2point


			function cartesian2point(x, y, z) result(point)
				use const
				use help
				use define_types
				type(position_in_space) point
				real x, y, z

				point%r = sqrt(x * x + y * y + z * z)
				point%alpha = acos(z / point%r)
				point%beta = atan(y, x)
				point%rvector(1) = x
				point%rvector(2) = y
				point%rvector(3) = z
				point%r_scaled = point%r / rm
				point%compute = .TRUE.


			end function cartesian2point




			! input the parameters of the sources from the file with the name
			! stored in the variable fname
			subroutine get_jets(sources, Ns, fname, ty, k, rlim1, rlim2)
				use const
				use define_types
				use gu
				implicit none
				integer, intent(in) :: Ns
				type(source_properties), intent(out) :: sources(Ns)
				integer, intent(in) :: ty
				real(8), intent(in) :: k
				real(8), intent(in) :: rlim1, rlim2
				integer i
				character(*), intent(in) :: fname
				character(5) header

				open(100, file = fname, status = 'old')
					read(100,*) header
					do i = 1, Ns
						read(100,*) sources(i)%alphaM, sources(i)%betaM, &
									sources(i)%zeta, sources(i)%eta
						sources(i)%r = rm
						sources(i)%alphaM = sources(i)%alphaM * deg2rad
						sources(i)%betaM = sources(i)%betaM * deg2rad
						sources(i)%alphaM = halfpi - sources(i)%alphaM

						sources(i)%zeta = sources(i)%zeta * deg2rad
						sources(i)%eta = sources(i)%eta * deg2rad

						sources(i)%rrM(1) = rm * sin(sources(i)%alphaM) &
						                       * cos(sources(i)%betaM)
						sources(i)%rrM(2) = rm * sin(sources(i)%alphaM) &
						                       * sin(sources(i)%betaM)
						sources(i)%rrM(3) = rm * cos(sources(i)%alphaM)

						sources(i)%is_jet = .TRUE.

						sources(i)%sd = ty

						sources(i)%ud%ud_shape = 1
						sources(i)%ud%umin = 0d0		! m/s
						sources(i)%ud%umax = 800d0		! m/s

						sources(i)%ejection_angle_distr = 1

						sources(i)%production_fun = 0
						! salt-poor size distributions
						if(ty == 1 .or. ty == 51 .or. ty == 52 &
						     .or. ty == 171 .or. ty == 172) then
							sources(i)%production_rate = 1.9d13 * k
							call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, &
						               sources(i)%sd, sources(i)%ud, rlim1, rlim2)
						endif
						! salt-rich size distribution
						if(ty == 3) then
							sources(i)%production_rate = 1.5d12 * k
							call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, &
						               sources(i)%sd, sources(i)%ud, rlim1, rlim2)
						endif

						call jet_direction(sources(i)%betaM, sources(i)%zeta, &
						     sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)


					enddo
				close(100)

			end subroutine get_jets




			! obtains Cartesian coordinates of a unit vector
			! in the moon-centered coordinate system
			! aligned with the ejection symmetry axis
			! the source's position (rrM - Cartesian coordinates, betaM - eastern longitude),
			! the jet's zenith angle (zeta) and azimuth (eta) are known
			subroutine jet_direction(betaM, zeta, eta, rrM, jetdir)
				use const
				use help
				real(8), intent(in) :: betaM, zeta, eta, rrM(3)
				real(8), intent(out) :: jetdir(3)
				real(8) xj(3), yj(3), rtmp(3)
				real(8) xout, yout, zout, tmpang

				rtmp = rrM / rm
				tmpang = 3d0 * halfpi-betaM
				if(zeta /= 0d0) then
					call eulrot(0d0, 0d0, tmpang, rtmp(1), rtmp(2), rtmp(3), &
					            xout, yout, zout, .FALSE.)
					rtmp(1) = xout ; rtmp(2) = yout ; rtmp(3) = zout

					xj(1) = 0d0
					xj(2) = sign(1d0, rtmp(3)) * abs(rtmp(3))
					xj(3) = -sign(1d0, rtmp(2)) * abs(rtmp(2))
					xj = xj / norma3d(xj)

					yj = vector_product(rtmp,xj)

					jetdir = sin(zeta) * cos(eta) * xj &
							- sin(zeta) * sin(eta) * yj &
							+ cos(zeta) * rtmp
					jetdir = jetdir / norma3d(jetdir)

					call eulrot(0d0, 0d0, tmpang, rtmp(1), rtmp(2), rtmp(3), &
					            xout, yout, zout, .TRUE.)
					rtmp(1) = xout ; rtmp(2) = yout ; rtmp(3) = zout

					call eulrot(0d0, 0d0, tmpang, jetdir(1), jetdir(2), jetdir(3), &
					            xout, yout, zout, .TRUE.)
					jetdir(1) = xout ; jetdir(2) = yout ; jetdir(3) = zout
				else
					jetdir = rtmp
				endif


			end subroutine jet_direction



end module inputdata
