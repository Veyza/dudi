! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! File: inputdata.f95
! Description: The subroutines suplying the parameters of sources, and 
!              spacecraft positions, and the auxiliary functions used by them

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module inputdata
	implicit none
		contains

			! input the spacecraft coordinates from the file with the name 
			! stored in the variable fname
			! in the file fname the coordinates are written as follows
			! radial distance from the moon center [m], latitude [deg], eastern longitude [deg]
			subroutine read_spacecraft_coordinates(points, nt, fname)
				use const
				use define_types
				implicit none
				integer, intent(in) :: nt
				type(position_in_space), intent(out) :: points(nt)
				integer i
				character(len = 29), intent(in) :: fname
				
				open(200, file = fname, status = 'old')
					do i = 1, nt
						read(200,*) points(i)%r, points(i)%alpha, points(i)%beta
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
			
			end subroutine read_spacecraft_coordinates
			
			
			
			! input the spacecraft positions from the file './input_data_files/Cassini_E2_flyby.dat'
			! the file contains the coordinates of the Cassini spacecraft
			! during its E2 flyby at Enceladus
			! the file was prepared in advance using SPICE software
			subroutine read_Cassini_E2(points, ttab, nt)
				use const
				use define_types
				implicit none
				integer, intent(in) :: nt
				type(position_in_space), intent(out) :: points(nt)
				real, intent(out) :: ttab(nt)
				integer i
				character(len = 41) :: fname = './input_data_files/Cassini_E2_flyby.dat'
				
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
			
			end subroutine read_Cassini_E2		




			! the subroutine provides the input data for the example
			! calculations of the Europa surface depositions
			subroutine get_europa_input(Ns, sources, nt, points, dphi)
				use define_types
				use gu
				implicit none
				integer, intent(in) :: Ns, nt
				real(8), intent(out) :: dphi(nt)
				type(source_properties), intent(out) :: sources(Ns)
				type(position_in_space), intent(out) :: points(nt)
				integer i, ead_choice(4), sd_choice(4)
				
				ead_choice = (/1, 3, 3, 1/)
				sd_choice = (/2, 3, 2, 3/)
				
				! define 4 sources with the same coordinates and verticle axis
				! of symmetry but different size- and ejection direction distributions
				do i = 1, Ns
					sources(i)%alphaM = halfpi
					sources(i)%betaM = 0d0
					sources(i)%zeta = 0d0
					sources(i)%eta = 0d0
					sources(i)%production_fun = 0
					sources(i)%production_rate = 1d14
					sources(i)%ud%ud_shape = 1
					sources(i)%ud%umin = 0d0
					sources(i)%ud%umax = 500d0
					sources(i)%ejection_angle_distr = ead_choice(i)
					sources(i)%sd = sd_choice(i)
					sources(i)%r = rm

									
					sources(i)%rrM(1) = rm * sin(sources(i)%alphaM) &
					                       * cos(sources(i)%betaM)
					sources(i)%rrM(2) = rm * sin(sources(i)%alphaM) &
					                       * sin(sources(i)%betaM)	
					sources(i)%rrM(3) = rm * cos(sources(i)%alphaM)
						
					sources(i)%is_jet = .TRUE.
					
					call jet_direction(sources(i)%betaM, sources(i)%zeta, sources(i)%eta, &
					                   sources(i)%rrM, sources(i)%symmetry_axis)
					call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, sources(i)%sd, &
					                 sources(i)%ud)
				enddo
				
				! the points are equidistantly placed on a 10deg arc
				! having the source at one of the arc's ends
				forall(i = 1:nt) dphi(i) = 2d0 * dble(i) / dble(nt) * deg2rad
				
				do i = 1, nt
					points(i)%r = rm
					points(i)%alpha = halfpi 
					points(i)%beta = dphi(i)
					points(i)%rvector(1) = points(i)%r * sin(points(i)%alpha) &
					                                   * cos(points(i)%beta)
					points(i)%rvector(2) = points(i)%r * sin(points(i)%alpha) &
					                                   * sin(points(i)%beta)
					points(i)%rvector(3) = points(i)%r * cos(points(i)%alpha)
					points(i)%r_scaled = 1d0
					points(i)%compute = .TRUE.
				enddo				
			
			
			end subroutine get_europa_input



			
			
			subroutine get_volcano_params(source, Ns)
				use const
				use define_types
				use gu
				implicit none
				integer, intent(in) :: Ns
				type(source_properties), intent(out) :: source(Ns)
				integer i
				
				do i = 1, Ns
					source(i)%alphaM = 5d0 * deg2rad
					source(i)%betaM = 315d0 * deg2rad
					source(i)%zeta = 3d0 * deg2rad
					source(i)%eta = pi
					source(i)%production_fun = 2
					source(i)%production_rate = 1d14
					source(i)%ud%ud_shape = 2
					source(i)%ud%umin = 720d0
					source(i)%ud%umax = 750d0		
					source(i)%ejection_angle_distr = 1
					source(i)%sd = 2
					source(i)%r = rm

									
					source(i)%rrM(1) = rm * sin(source(i)%alphaM) &
					                      * cos(source(i)%betaM)
					source(i)%rrM(2) = rm * sin(source(i)%alphaM) &
					                      * sin(source(i)%betaM)
					source(i)%rrM(3) = rm * cos(source(i)%alphaM)
						
					source(i)%is_jet = .TRUE.
					
					call jet_direction(source(i)%betaM, source(i)%zeta, &
					         source(i)%eta, source(i)%rrM, source(i)%symmetry_axis)
					call Gu_integral(source(i)%ui, source(i)%Gu_precalc, &
					                 source(i)%sd, source(i)%ud)
				enddo
			
			end subroutine get_volcano_params
			
			
			
			! input the parameters of the sources from the file with the name
			! stored in the variable fname
			subroutine read_sources_params(sources, Ns, fname)
				use const
				use define_types
				use gu
				implicit none
				integer, intent(in) :: Ns
				type(source_properties), intent(out) :: sources(Ns)
				integer i
				character(*), intent(in) :: fname
				
				open(100, file = fname, status = 'old')
					do i = 1, Ns
						read(100,*) sources(i)%alphaM, sources(i)%betaM, &
									sources(i)%zeta, sources(i)%eta, &
									sources(i)%production_fun, sources(i)%production_rate, &
									sources(i)%ud%ud_shape, sources(i)%ud%umin, &
									sources(i)%ud%umax, &
									sources(i)%ejection_angle_distr, &
									sources(i)%sd
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
						
						call jet_direction(sources(i)%betaM, sources(i)%zeta, &
						     sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)
						     
						call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, &
						               sources(i)%sd, sources(i)%ud)
												
					enddo
				close(100)
			
			end subroutine read_sources_params
			
			
			
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
