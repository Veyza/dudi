! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! File: distributions_fun.f95
! Description: The functions describing the ejection process and auxilary
!              functions used by them

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module distributions_fun
	implicit none
		contains
		
		
			! The axisymmetric distribution of ejection direction
			! distribution_shape is the parameter used to select
			! the expression for the PDF
			! wpsi is the polar angle in the coordinate system where
			! the distribution is axisymmetrical
			! psi is the polar angle in the horizontal coordinate system
			! lambdaM is the azimuth in the horizontal CS
			! zeta and eta are respectively zenith angle and azimuth
			! of the distribution symmetry axis in the horizontal CS
			function ejection_direction_distribution(distribution_shape, &
			                    wpsi, psi, lambdaM, zeta, eta) result(fpsi)
				use const
				implicit none
				integer N, i
				real(8), parameter :: normconst1 = 7.5960829056967811d-3
				real(8), parameter :: normconst3 = 8.5760756217641998d-1
				real(8), parameter :: psimax0 = 0d0
				real(8), parameter :: psimax45 = 0.7853982d0
				real(8), parameter :: omega3 = 0.05235988d0
				real(8), parameter :: omega5 = 0.08726646d0
				real(8), parameter :: omega10 = 0.1745329d0
				real(8), parameter :: omega45 = 0.7853982d0
				real(8) fpsi, Jpsi
				integer, intent(in) :: distribution_shape
				real(8), intent(in) :: psi, wpsi, lambdaM, zeta, eta
				
				select case(distribution_shape)
					case(1)
					! pseudo Gaussian distribution of polar angle, uniform distribution of azimuth
						if(psi < halfpi * 0.99 .and. wpsi < halfpi * 0.99) then
							fpsi = Exp(-(wpsi-psimax0)**2 / 2d0 / omega5 / omega5)
							! this factor is normalization due to the fact that fpsi
							! domain is from 0 to pi/2 and not from -infinity to +infinity
							fpsi = fpsi / normconst1
							fpsi = fpsi / twopi
						else
							fpsi = 0d0
						endif
					case(2)
					! Uniform distribution of polar angle inside a cone,
					! uniform distribution of azimuth
						if(wpsi <= omega10) then
							fpsi = 1d0 / (1d0 - cos(omega10)) / twopi
						else
							fpsi = 0d0
						endif
					case(3)
					! pseudo Gaussian distribution of polar angle,
					! uniform distribution of azimuth
						if(wpsi < halfpi * 0.99) then
							fpsi = Exp(-(wpsi-psimax45)**2 / 2d0 / omega45 / omega45)
							! this factor is normalization due to the fact that fpsi
							! domain is from 0 to pi/2 and not from -infinity to +infinity
							fpsi = fpsi / normconst3
							fpsi = fpsi / twopi
						else
							fpsi = 0d0
						endif
					case(4)
					! HERE IS THE PLACE FOR WRITING YOUR OWN PDF
					fpsi = 0d0
				endselect
				if(zeta /= 0d0) then
					Jpsi = Jacobian_tilt(psi, lambdaM, zeta, eta)
					fpsi = fpsi * Jpsi
				endif
				
				fpsi = fpsi * sin(wpsi)
					
			end function ejection_direction_distribution

			
						
			! Jacobian of coordinate transformation from vertical CS
			! to the CS with z-axis coinciding with the jet axis of symmetry
			! zeta and A are respectively zenith angle ang azimuth
			! of the distribution symmetry axis in the horizontal CS
			! psi and lambdaM are respectively polar angle
			! and azimuth of ejection in the horizontal CS
			function Jacobian_tilt(psi, lambdaM, zeta, A) result(J)
				use const
				implicit none
				real(8), intent(in) :: psi, lambdaM, zeta, A
				real(8) sinpsi, cospsi, sinzeta, coszeta
				real(8) J
				
				if(psi < 0.12d0 .and. zeta < 0.12d0) then
					J = psi / sqrt(psi*psi + zeta*zeta &
					               - 2d0 * psi * zeta * cos(lambdaM - A))
				else
					sinpsi = sin(psi) ; cospsi = cos(psi)
					sinzeta = sin(zeta) ; coszeta = cos(zeta)

					J = 4d0 * sinpsi / sqrt(10d0 - 2d0 * (cospsi - sinpsi) * (cospsi + sinpsi) &
						- 3d0 * cos(2d0 * (psi - zeta)) - 2d0 * (coszeta - sinzeta) * (coszeta + sinzeta) &
						- 3d0 * cos(2d0 * (psi + zeta)) &
						- 8d0 * cos(2d0 * (lambdaM - A)) * sinpsi * sinpsi * sinzeta * sinzeta &
						- 32d0 * cos(lambdaM - A) * sinpsi * cospsi * sinzeta * coszeta)
				endif
				
			end function Jacobian_tilt
			
			
			! This function represents the ejection speed distribution
			! (possibly time-dependent)
			! ud is  the parameter used to select the expression for the distribution
			! u is the ejection speed
			! R is the particle size
			function ejection_speed_distribution(ud, u, R) result(fu)
				use const
				use define_types
				implicit none
				real(8), parameter :: Rc = 0.5
				real(8) urel, Rrel, fu
				type(ejection_speed_properties) ud
				real(8) R, u
								
				select case(ud%ud_shape)
					case(1)
						Rrel = R / Rc
						urel = u / ud%umax
						fu = Rrel * (1d0 + Rrel) * (1d0 - urel)**(Rrel - 1.0) * urel / ud%umax
					case(2)
						fu = 0d0
						if(u < ud%umax .and. u > ud%umin) fu = 1d0 / (ud%umax - ud%umin)
					case(3)
					! HERE IS THE PLACE FOR WRITING YOUR OWN PDF
						fu = 0d0
				endselect
			
			end function ejection_speed_distribution
			
			
			
			
			! This function represents the size distribution of the dust particles.
			! It can be used also to obtain the mean radius, cross section
			! or volume of the dust particles.
			! R is a particle radius
			! sd is a parameter used to select the expression for the distribution
			! mom defines the obtained quantity: 0 -- number density,
			! 1 -- mean radius, 2 -- cross section, 3 -- volume
			function size_distribution(R, sd, mom) result(fR)
				use const
				implicit none
				real(8), parameter :: mu = -1d0
				real(8), parameter :: sigma = 1.5d0
				real(8)  q
				real(8), parameter :: r1 = 0.2d0
				real(8), parameter :: r2 = 20d0
				real(8), intent(in) :: R
				integer, intent(in) :: sd, mom
				real(8) fR, C_size_distr
				
				select case(sd)
					case(1)
						fR = exp(-(log(R) - mu)**2 / 2d0 / sigma**2) / R
						C_size_distr = sigma * sqrtpi * sqrt2d0
					case(2)
						q = 3d0
						if(r1 <= R .and. R <= r2) then
							fR = R**(-3)
						else
							fR = 0d0
						endif
						C_size_distr = (r2**(1d0-q) - r1**(1d0-q)) / (1d0 - q)
					case(3)
						q = 5d0
						if(r1 <= R .and. R <= r2) then
							fR = R**(-3)
						else
							fR = 0d0
						endif
						C_size_distr = (r2**(1d0-q) - r1**(1d0-q)) / (1d0 - q)
					case(4)
						! HERE IS THE PLACE FOR WRITING YOUR OWN PDF
						fR = 0d0
				endselect
				fR =  R**mom * fR / C_size_distr
				
				
			end function size_distribution
			
			
			! This function represents the factor gamma(t) in Formula 43
			! t is the moment of ejection
			! gamma0 is a parameter which can be used in the definition of
			! the function. In the implementet examples gamma0 is the production
			! rate at maximum.
			! ratefun is the parameter used to choose the expression for the gammarate
			function production_rate(t, gamma0, ratefun) result(gammarate)
				implicit none
				real(8) tmax
				integer ratefun
				real(8) t, gamma0, gammarate
				
				select case(ratefun)
					case(1)
						gammarate = 0d0
						if(t < tmax .and. t >= 0d0) gammarate = gamma0
					case(2)
						gammarate = 0d0
						tmax = 5d2
						if(t > 0d0 .and. t < 2d0 * tmax) then
							gammarate = gamma0 * (-t**2 + 2d0 * t * tmax) / tmax**2
						endif
					case(3)
						! HERE IS THE PLACE TO WRITE YOUR OWN FUNCTION FOR THE PRODUCTION RATE
						gammarate = 0d0
				endselect
				
			end function production_rate
			

end module distributions_fun
