! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! File: twobody_fun.f95
! Description: The subroutines used to compute the integrand

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module TwoBody_fun
	implicit none
	contains

		! Compute various quantities needed for evaluation of the integrand
		subroutine Apu_u_angles_ddphidtheta(point, s, v, theta, e, dbeta, dphi, &
						u, psi, wpsi, lambdaM, lambda, ddphidtheta, &
						sindphi, dphi_is_large)
			use const
			use define_types
			use help
			implicit none
			type(position_in_space), intent(in) :: point
			type(source_properties), intent(in) :: s
			real(8), intent(in) :: v, u, theta, e
			real(8), intent(in) :: dbeta, dphi
			logical, intent(in) :: dphi_is_large
			real(8), intent(out) :: psi, wpsi, lambdaM, lambda, ddphidtheta
			real(8), intent(out) :: sindphi
			real(8) cosal, cosalM, sinal, sinalM
			real(8) sinlambdaM, coslambdaM, sinlambda, coslambda
			real(8) dee, pp, dpp, wrr, wvv, hh
			real(8) cosphi, cosphim, cosdphi, sintheta
			real(8) dphi1, dphi2, dphi3, dphi4, delta, numder
			integer i
			
			sinal = sin(point%alpha) ; cosal = cos(point%alpha)
			sinalM = sin(s%alphaM) ; cosalm = cos(s%alphaM)
			sindphi = sin(dphi) ; cosdphi = cos(dphi)
			sintheta = sin(theta)
							
			!	 angular momentum (eq 27)		
			hh = point%r * v * sintheta		
				!  psi
			psi = asin(hh / s%r / u)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(psi /= psi .or. abs(psi-halfpi) < 1d-8) then
				if(psi /= psi) then
					write(666,*) 'Apu << sin(psi) =', hh / s%r / u, 'corrections applied'
				endif
				if(abs(psi-halfpi) < 1d-8) then
					write(666,*) 'psi is close to pi/2, corrections applied'
				endif
				psi = halfpi - 1d-5
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			! lambda and lambdaM are to be found from a spherical triangle
			! the direction and the length of arc along which a particle
			! traveled from rM to r matters for exact geometry of
			! the spherical triangle namely, sings of lambda and lambdaM
			! depend on it
			if(s%betaM > point%beta) then
				if(s%betaM - point%beta < pi) then
					sinlambda = -sinalM * sin(dbeta) / sindphi
					coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
				else
					sinlambda = sinalM * sin(dbeta) / sindphi
					coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
				endif				
			else
				if(point%beta - s%betaM < pi) then
					sinlambda = sinalM * sin(dbeta) / sindphi
					coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
				else
					sinlambda = -sinalM * sin(dbeta) / sindphi
					coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
				endif
			endif
			sinlambdaM = sinal * sinlambda / sinalM
			coslambdaM = (cosal - cosalM * cosdphi) / sinalM / sindphi

			if(dphi_is_large) then
				sinlambda = - sinlambda
				coslambda = - coslambda
				sinlambdaM = - sinlambdaM
				coslambdaM = - coslambdaM
			endif
			lambda = myatan1(coslambda, sinlambda)
			
			lambdaM = myatan1(coslambdaM, sinlambdaM)
			
			wpsi = acos(cos(psi) * cos(s%zeta) &
					+ cos(lambdaM - s%eta) * sin(psi) * sin(s%zeta))

			
			wrr = point%r_scaled ; wvv = v / vesc

			pp = 2d0 * wrr * wrr * wvv * wvv * sintheta**2
			dpp = 2d0 * pp / tan(theta)
			dee = (wvv * wvv - 1d0 / wrr) * dpp / e
			
			cosphim = (pp - 1d0) / e
			cosphi = (pp / wrr - 1d0) / e

			ddphidtheta = (((1d0 + pp) * (wrr * wvv * wvv - 1d0) + wrr) * 2d0 * cos(theta)) &
				/ (wvv * sqrt(wrr * (-1d0 + wrr + wrr * wvv**2 - wrr**3 * wvv**2 * sintheta**2))) &
				- (2d0 * sintheta**2 * (-2d0 + 2d0 * wrr * wvv**2 + 1d0 / sintheta**2))

			ddphidtheta = ddphidtheta * wvv * wvv * wrr / e / e
			if(ddphidtheta /= ddphidtheta) then
				delta = 1d-3
				dphi1 = deltaphi(theta+2d0*delta, wrr, wvv)
				dphi2 = deltaphi(theta+delta, wrr, wvv)
				dphi3 = deltaphi(theta-delta, wrr, wvv)
				dphi4 = deltaphi(theta-2d0*delta, wrr, wvv)
				
				numder = (-dphi1 + 8d0 * dphi2 - 8d0 * dphi3 + dphi4) / 12d0 / delta 
				numder = (dphi2 - dphi3) / 2d0 / delta
				write(666,*) 'the derivative d\Delta\phi/d\theta was &
								obtained numerically because the analytical &
								expression contains numerically difficult parts'
				ddphidtheta = numder
			endif			

			
		end subroutine Apu_u_angles_ddphidtheta



		! compute \Delta\phi from \theta, velocity and spacecraft position
		function deltaphi(theta, wrr, wvv)
			use const
			implicit none
			real(8), intent(in) :: theta
			real(8), intent(in) :: wrr, wvv
			real(8) deltaphi, hh, e, pp, cosphi, cosphim
			
			pp = 2d0 * wrr * wrr * wvv * wvv * sin(theta)**2
			e = sqrt(1d0 + 2d0 * pp * (wvv * wvv - 1d0 / wrr))
			
			cosphim = (pp - 1d0) / e
			if(abs(cosphim) > 1d0) cosphim = sign(1d0, cosphim)
			cosphi = (pp / wrr - 1d0) / e
			if(abs(cosphi) > 1d0) cosphi = sign(1d0, cosphi)
			if(theta < halfpi) then
				deltaphi = acos(cosphi) - acos(cosphim)
			else
				deltaphi = 2d0 * pi - acos(cosphi) - acos(cosphim)
			endif
		
		end function deltaphi



		! Integrand_number_density performs integration over theta and lambda,
		! returns the expression standing under the integral over v
		subroutine Integrand_number_density(Integrand, velocity, amin, &
			                                point, dphi, dbeta, s, tnow)
			use const
			use define_types
			use help
			use distributions_fun
			implicit none
			integer i
			real(8), intent(in) :: velocity, dphi, dbeta, amin, tnow
			type(source_properties), intent(in) :: s
			type(position_in_space), intent(in) :: point
			real(8) lambdaM, sinlambdaM, coslambdaM, lambda
			real(8) uu, Ekep, ee(2), psi(2), wpsi(2), semi_major_axis
			real(8), intent(out) :: Integrand
			real(8) theta(2), ddphidtheta(2), tmp, deltat(2)
			real(8) fac1, fac2, Jpsi, tmpIntegrand, rate(2), cf, alt
			real(8) sindphi
			logical dphi_is_large(2)

			semi_major_axis = (2d0 / point%r - velocity**2 / gm)**(-1)
			theta = -999d0
			! find solutions for theta
			if(semi_major_axis < 0d0 .or. semi_major_axis >= amin) then
				if(semi_major_axis > 0d0) then
					call theta_geometry_ellipse(point%r, s%r, velocity, &
					                  dphi, semi_major_axis, ee, theta, &
					                  deltat, dphi_is_large, s%production_fun > 0)
				else
					call theta_geometry_hyperbola(point%r, s%r, velocity, &
					               dphi, abs(semi_major_axis), ee, theta, &
					            deltat, dphi_is_large, s%production_fun > 0)
!~ 			write(*,*) 'intergand << theta', theta
				endif
			endif
			uu = sqrt(vesc * vesc &
						+ 2.0 * (velocity * velocity / 2d0 - gm / point%r))
			if(uu > s%ud%umax .or. uu < s%ud%umin) then
				fac1 = 0d0
			else
		!	 interpolate Gu from a precalculated table
				if(uu / s%ud%umax > s%ui(GRN)) then
					fac1 = velocity * s%Gu_precalc(GRN) / uu / uu
				else
					fac1 = LiNTERPOL(GRN, s%Gu_precalc, s%ui, uu) 
					fac1 = velocity * fac1 / uu / uu
				endif
			endif
			
			Integrand = 0d0
			wpsi = halfpi	
			do i = 1, 2
				if(s%production_fun > 0) then
					rate(i) = production_rate(tnow - deltat(i), s%production_rate, &
				                                               s%production_fun)
				else
					rate(i) = s%production_rate
				endif
				if(theta(i) >= 0d0 .and. rate(i) > 0) then

					call Apu_u_angles_ddphidtheta(point, s, velocity, &
					                    theta(i), ee(i), dbeta, dphi, &
					                             uu, psi(i), wpsi(i), &
					                 lambdaM, lambda, ddphidtheta(i), &
					                       sindphi, dphi_is_large(i))
					
			! the distribution of ejection angle is defined in coordinates (wpsi, wlambdaM) 
			! where wpsi is an angle between the jet main axis of symmetry and the direction of ejection
			! wlambdaM is a longitude in the plane perpendicular to the jet's main axis
			! However, the factor 1/cos(psi) comes from the Jacobian of transformation 
			! (alphaM, betaM, u, psi, lambdaM) -> (alpha, beta, v, theta, lambda)
			! and here psi is an angle between the direction of ejection and the normal to surface
					fac2 = ejection_direction_distribution(s%ejection_angle_distr, wpsi(i), &
															psi(i), lambdaM, s%zeta, s%eta)
					fac2 = fac2 / cos(psi(i))
											
					tmpIntegrand = fac1 * fac2 / abs(ddphidtheta(i)) * rate(i)
					! to calculate the flux
					if(flux) then
						tmpIntegrand = tmpIntegrand * velocity * abs(cos(theta(i)))
					endif
					
					Integrand = Integrand + tmpIntegrand
					if(tmpIntegrand /= tmpIntegrand) then
						write(*,*) 'NaN is obtained for an integrand value'
						write(*,*) 'factor related to ejection speed distribution: fac1 =', fac1
						write(*,*) 'factor related to ejection direction distribution: fac2 =', fac2
						write(*,*) 'the partial derivative of delta phi by theta =', ddphidtheta(i)
						write(*,*) 'theta =', theta(i), 'psi =', psi(i), 'lambdaM =', lambdaM
						write(*,*) 'dust production rate =', rate(i)
						stop
					endif
				endif
				tmpIntegrand = 0d0
			enddo
			
		end subroutine Integrand_number_density	


		
		
		! subroutine theta_geometry_hyperbola recieves vectors' absolute
		! values and the angle between these vectors
		! it's assumed that the point (0, 0, 0) is a focus of an hyperbola
		! value of a semi major axis of the hyperbola, the points lay on is
		! also an input parameter of the subroutine
		! the subroutine returns theta the angle between radius-vector r
		! and a tangent to the hyperbola in the point r
		! the choice between two possible values of this angle is made
		! in the way that movement from point rm to the point r along
		! the hyperbola is possible
		! (means: r and rm lay on the same quadrant in the CS in which
		! the equation of the hyperbola is in canonical form)
		! In general case there are two possible hyperbolae
		! => two possible values of theta are to be investigated
		! Returns array of 2 vallues of theta.
		! If theta = large negative number,
		! it means "no physically plausible solution can be found"
		subroutine theta_geometry_hyperbola(r0, rm0, vv, phi, a0, &
										ee, theta, deltat, dphi_is_large, timeDependence)
			use help
			use const
			implicit none
			real(8), intent(in) :: r0, rm0, a0, phi,  vv
			real(8), intent(out) :: ee(2), theta(2), deltat(2)
			real(8) r, rmoon, a
			real(8) c, x(2), y(2), xm, ym, b
			real(8) shift(2), r2d(2), rm2d(2), angle, tangent(2), tmp(2)
			real(8) cosf1, eanm, ean, psi, f1, f2, discr
			real(8) one_plus_e, one_minus_e, one_minus_e2, aux, tmp1, tmp2
			integer i
			logical solved
			logical, intent(out) :: dphi_is_large(2)
			logical, intent(in) :: timeDependence
			
			solved = .FALSE.
			theta = -555d0
			dphi_is_large = .FALSE.
															
			r = r0 / rm0 ; rmoon = rm0 / rm0; a = a0 / rm0							
			! define x-axis in the same direction as r-vector
			r2d(1) = r ; r2d(2) = 0d0					
			
			! we don't have enough information to define the sign of r
			! vector in the right-handed coordinate system
			! but the value of theta that we are looking for are the same
			! in both cases
			
			rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
			
			! (x(1),y(1)) and (x(2),y(2)) are coordinates
			! of 2 possible position of the hyperbola's second focus
			call circle_intersection(rm2d(1), rm2d(2), 2d0 * a + rmoon, &
										r2d(1), 2d0 * a + r, x, y)		
				
			do i = 1, 2
				! shift is coordinates of the ellipse's center
				! in the CS centered at the focus
				shift(1) = x(i) / 2d0 ; shift(2) = y(i)	/ 2d0		
				
				! coords of vector r in CS centered at center of the ellipse
				r2d = r2d - shift	
				! coords of vector rm in CS centered at center of the ellipse								
				rm2d = rm2d - shift								
				
				! angle between major axis of the ellipse and the current x-axis
				angle = atan(y(i) / x(i))
				! vector r in the CS with its center at the center of the ellipse
				! and the x-axis along major axis of the ellipse									
				r2d = rot2d(r2d, -angle)	
				! vector rm in the CS with its center at the center of the ellipse
				! and the x-axis along major axis of the ellipse								
				rm2d = rot2d(rm2d, -angle)	
				! vector shift in the CS with its center at the center of the ellipse
				! and the x-axis along major axis of the ellipse							
				shift = rot2d(shift, -angle)
				! the hyperbola solves our problem only if r and rm lay
				! on the same branch and the trajectory doesn't intersect
				! the moon's surface (the particle doesn't pass the pericenter
				if(r2d(1) / rm2d(1) > 0d0 &								
				.and. r2d(2) / rm2d(2) > 0d0) then					
					c = 0.5d0 * sqrt(x(i)**2 + y(i)**2)						
					b = sqrt(c**2 - a**2)									
					ee(i) = c / a
					one_plus_e = 1d0 + ee(i)
					one_minus_e = 1d0 - ee(i)
					one_minus_e2 = one_plus_e * one_minus_e
					aux = -a * one_minus_e2
					cosf1 = (aux - 1d0) / ee(i)
					f1 = acos(cosf1)
					f2 = f1 + phi
					
					theta(i) = halfpi - atan((ee(i) * sin(f2)) / (1d0 + ee(i) * cos(f2)))
											
					call control(theta(i), ee(i), phi, vv, r0, rm0, dphi_is_large(i), solved, discr)
					if(timeDependence) then
						ean = 2d0 * atanh(tan(f2/2d0) * sqrt(-one_minus_e / one_plus_e))
						eanm = 2d0 * atanh(tan(f1/2d0) * sqrt(-one_minus_e / one_plus_e))
						tmp1 = a0 * sqrt(a0 / gm) * (ee(i) * sinh(eanm) - eanm)
						tmp2 = a0 * sqrt(a0 / gm) * (ee(i) * sinh(ean) - ean)
						deltat(i) = tmp2 - tmp1
					else
						deltat(i) = 0d0
					endif
					if(.not. solved .and. theta(i) > 0d0) then
						write(666,*) '   '
						write(666,*) 'Theta was found with an insufficient accuracy of', &
										discr, 'from the geometry of hyperbola'
						write(666,*) 'for r = ', r0
						write(666,*) 'dphi = ', phi
						write(666,*) 'the obtained value of eccentricity = ', ee(i)
						write(666,*) 'and the obtained value of theta = ', theta(i)
						if(dphi_is_large(i)) then
							write(666,*) 'the case of \Delta\phi > pi has been encoutered'
						endif
					endif
						
				else
					theta(i) = -444d0
					ee(i) = -444d0
					dphi_is_large(i) = .FALSE.
				endif
				! we have changed the vectors r and rmoon we started from
				! so we need to go back to the beginning
				! to find the second value of theta
				r = r0 / rm0 ; rmoon = rm0 / rm0; a = a0 / rm0
				rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
				r2d(1) = r ; r2d(2) = 0d0
				
			enddo
		
		end subroutine theta_geometry_hyperbola



		! subroutine theta_geometry_ellipse recieves absolute values
		! of two vectors and the angle between these vectors
		! it's assumed that the point (0, 0, 0) is a focus of an ellipse
		! value of a semi major axis of an ellipse, the points lay on
		! is also an input parameter of the subroutine
		! the subroutine returns theta that is the angle between
		! radius-vector r and a tangent to the ellipse in the point r
		! the choice between two possible values of this angle is made
		! in the way that movement happens from point rm to the point r 
		! In general case there are two possible ellipses
		! => two possible values of theta are to be found
		! Returns array of 2 vallues of theta.
		! If theta = large negative number
		! it means "no physically plausible solution can be found"
		subroutine theta_geometry_ellipse(r0, rm0, vv, phi, a0, &
								ee, theta, deltat, dphi_is_large, timeDependence)
			use help
			use const
			implicit none
			real(8), intent(in) :: r0, rm0, a0, phi,  vv
			real(8), intent(out) :: ee(2), theta(2), deltat(2)
			real(8) r, rmoon, a, discr
			real(8) cc(2), x(2), y(2)
			real(8) r2d(2), rm2d(2), tmp(2), tmp1, tmp2
			real(8) E1(1), f1, f2, sinf1, cosf1, ean, rtest
			real(8) one_plus_e, one_minus_e, one_minus_e2, eanm
			real(8) aux
			integer i
			logical solved
			logical, intent(out) :: dphi_is_large(2)
			logical, intent(in) :: timeDependence
			
			solved = .FALSE.
			dphi_is_large = .FALSE.
			theta = -888d0
			
			r = r0 / rm0  ; rmoon = rm0 / rm0; a = a0 / rm0
			! define x-axis in the same direction as rmoonvector
			r2d(1) = r; r2d(2) = 0d0					
			! we don't have enough information to define the sign of r
			! vector in the right-handed coordinate system
			! but the value of theta that we are looking for are the same
			! in both cases
			rm2d(1) = rmoon * cos(phi) ; rm2d(2) = rmoon * sin(phi)
			
			! (x(1),y(1)) and (x(2),y(2)) are coordinates of 2 possible
			! position of the ellsipse's second focus
			call circle_intersection(rm2d(1), rm2d(2), 2d0 * a - rmoon, &
									r2d(1), 2d0 * a - r, x, y) 
			! distance between the foci of the ellipse
			cc = sqrt(x**2 + y**2)
			
			do i = 1, 2
				if(cc(i) == cc(i)) then
					ee(i) = cc(i) / 2d0 / a

					one_plus_e = 1d0 + ee(i)
					one_minus_e = 1d0 - ee(i)
					one_minus_e2 = one_plus_e * one_minus_e
					cosf1 = (-x(i) * rm2d(1) - y(i) * rm2d(2)) / rmoon / cc(i)
					f1 = acos(cosf1)
					f2 = f1 + phi
					
					if(r < rmoon) then
						rtest = a * one_minus_e2 / (1d0 + ee(i) * cos(f2))
						! if r and rm are both located after apocenter
						if(abs(1d0 - rtest / r ) > 1d-6) then				
							f2 = (twopi - f1) + phi
							rtest = a * one_minus_e2 / (1d0 + ee(i) * cos(f2))
							! rm is before apocenter, the case of large dphi encountered	
							! phi is the angle between vectors r and rm.
							! it can be that between r and rm the particle
							! traveled the arc of 2pi - phi
							! if the direction of movement along the ellipse
							! is chosen incorrect rtest /= r
							if(abs(1d0 - rtest / r ) > 1d-6) then			
								f2 = f1 + (2d0 * pi - phi)
								dphi_is_large(i) = .TRUE.
								rtest = a * one_minus_e2 / (1d0 + ee(i) * cos(f2))
								! rm is after apocenter, the case of large dphi encountered
								if(abs(1d0 - rtest / r ) > 1d-6) then		
									f1 = twopi -f1
									f2 = f1 + (2d0 * pi - phi)
									dphi_is_large(i) = .TRUE.
								endif
							else
								f1 = twopi - f1
							endif
						endif
						if(f2 > twopi) f2 = f2 - twopi
					else						
						rtest = a * one_minus_e2 / (1d0 + ee(i) * cos(f2))
						if(abs(1d0 - rtest / r ) > 1d-6) then
							f2 = f1 + (2d0 * pi - phi)
							dphi_is_large(i) = .TRUE.
						endif
					endif
					! No ejection downward even if rM > r
					if(f1 < pi - 1d-4) then			
						theta(i) = halfpi - atan((ee(i) * sin(f2)) &
												/ (1d0 + ee(i) * cos(f2)))
						
						call control(theta(i), ee(i), phi, vv, r0, rm0, &
									dphi_is_large(i), solved, discr)
					
						if(timeDependence) then
							ean = 2d0 * atan(tan(f2/2d0) &
									* sqrt(one_minus_e / one_plus_e))
							if(ean < 0d0) ean = twopi + ean
							eanm = 2d0 * atan(tan(f1/2d0) &
									* sqrt(one_minus_e / one_plus_e))
							tmp1 = a0 * sqrt(a0 / gm) * (eanm - ee(i) * sin(eanm))
							tmp2 = a0 * sqrt(a0 / gm) * (ean - ee(i) * sin(ean))
							deltat(i) = tmp2 - tmp1
						else
							deltat(i) = 0d0
						endif
					else
						theta = -777d0
					endif
											
					if(.not. solved .and. theta(i) > 0d0) then
						write(666,*) '   '
						write(666,*) 'Theta was found with an insufficient accuracy of', &
										discr, 'from the geometry of ellipse'
						write(666,*) 'for r =  ', r0
						write(666,*) 'rt =     ', rtest * rm0
						write(666,*) 'for rm = ', rm0
						write(666,*) 'f1 = ', f1
						write(666,*) 'f2 = ', f2
						write(666,*) 'dphi = ', phi
						write(666,*) 'the obtained value of eccentricity = ', ee(i)
						write(666,*) 'and the obtained value of theta = ', theta(i)
						if(dphi_is_large(i)) then
							write(666,*) 'the case of \Delta\phi > pi has been encoutered'
						endif
					endif
				else
					theta(i) = -888d0
					ee(i) = -888d0
					dphi_is_large(i) = .FALSE.
				endif
									
			enddo
		
		
		end subroutine theta_geometry_ellipse


		

		
		! subroutine control tests if the obtained value of theta is correct
		! the criterion is: using  the obtained value of theta
		! one gets the same value of dphi which was used to calculate the theta
		! also the eccentricity values are compared
		! input parameters: r0 and rm0 - lengths of two vectors - positions on the orbit,
		! vv - speed at position r0, phi - angle between r0 and rm0 used
		! in subroutine theta_geometry_...to calculate theta
		! theta is the angle between r0 and velosity at the position r0
		! ee is eccentricity obtained in theta_geometry...
		! dphi_is_large tells if it is the case when the particle traveled
		! from rm to r over an arc of 2pi - dphi
		! the subroutine control uses formulae for energy
		! and angular momentum to obtain values of ee1 and dphi
		! they must differ in less then eps, in this case the function returns TRUE
		! otherwise it returns FALSE in the variable solved
		! it also returns the estimated accuracy stored in the variable discr
		subroutine control(theta, ee, phi, vv, r0, rm0, dphi_is_large, solved, discr)
			use const
			implicit none
			real(8), parameter :: eps = 1d-4
			real(8) hh, hh2, ee1, cosp, cospm, phi1, phi1m, dphi, Ekep, darccosdr
			real(8), intent(in) :: vv, phi, ee, r0, rm0, theta
			logical, intent(out) :: solved
			real(8), intent(out) :: discr
			logical, intent(in) :: dphi_is_large
			
			Ekep = vv * vv / 2d0 - gm / r0
			hh = r0 * vv * sin(theta)
			hh2 = hh * hh
		!  eccentricity (eq 31)
			ee1 = sqrt(1d0 + 2d0 * Ekep * (hh / gm) * (hh / gm))
			if(abs(ee1 - ee) > eps) then
				write(666,*) 'eccentricity is incorrect:', ee, 'instead of', ee1
			endif

		!   delta phi (equation 32)
			cosp = (hh2 / r0 / gm - 1d0) / ee1
			cospm = (hh2 / rm0 / gm - 1d0) / ee1
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(cosp > 1d0) then
				cosp = 1d0
				write(666,*) 'cos(phi) > 1 obtained, corrections applied'
			endif
			if(cosp < -1d0) then
				cosp = -1d0
				write(666,*) 'cos(phi) < -1 obtained, corrections applied'
			endif
			
			if(cospm > 1d0) then
				cospm = 1d0
				write(666,*) 'cos(phiM) > 1 obtained, corrections applied'
			endif
			if(cospm < -1d0) then
				cospm = -1d0
				write(666,*) 'cos(phiM) < -1 obtained, corrections applied'
			endif
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			phi1 = acos(cosp)
			phi1m = acos(cospm)
			if(theta < halfpi) then
				dphi = phi1 - phi1m
			else
				dphi = (2d0 * pi - phi1) - phi1m
			endif
			
			if(dphi_is_large) then
				discr = abs(2d0 * pi - dphi - phi) / abs(phi)
				solved = discr < eps
			else
				discr = abs(dphi - phi) / abs(phi)
				solved = discr < eps
			endif
			
		end subroutine control
			
			



end module twobody_fun
