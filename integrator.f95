! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License 
! (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! File: integrator.f95
! Description: The subroutines that manage the numerical integration

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module integrator
	implicit none
	! change of veps is not recommended
	! difference of the actual minimal velocity value
	! from the left boundary used in integration
	real(8), parameter :: veps = 1d-10
	real(8), parameter :: trpower = 4
		contains


		! The function rules the integration determining the parameters and the limits
		subroutine DUDI(density, point, source, tnow)
			use const
			use define_types
			use help
			use TwoBody_fun
			implicit none
			integer Nprestep
			! fraction which the interval of the pole integration constitutes
			! to the total interval of integration
			real(8), parameter :: prestep_relative_size = 1d-4
			integer i
			real, intent(out) :: density(2)
			real(8), intent(in) :: tnow
			real(8) vc, vmax, v_limits(2), amin, vmin, vmin2, amin2
			type(position_in_space), intent(in) :: point
			type(source_properties), intent(in) :: source
			real(8) ldif, lsum, term, vi, f1, f2, vinterval, viprev
			real(8) angle, dphi, dbeta
			real xel(order_v_el), wel(order_v_el)
			real xhy(order_v_hy), why(order_v_hy)
			logical pole
			
			call ApuTrajectory(point, dphi, dbeta, angle, source)
			! escape velosity at distance rr
			vc = sqrt(2d0 * gm / point%r)
			
			amin = (point%r + source%r) / 4d0 &
			    + 0.5d0 * sqrt((point%r**2 + source%r**2) &
			    / 4d0 - point%r * source%r * cos(dphi) / 2d0)

			amin2 = (2d0 / source%r - source%ud%umin**2 / gm)**(-1)
			if(amin2 /= amin2) amin2 = 0d0
			! minimal velocity possible at given position (defines minimal energy
			! or "size" of the orbit on which a particle can go from rm to rr
			vmin = sqrt(gm * (2d0 / point%r - 1d0 / amin))
			vmin2 = sqrt(gm * (2d0 / point%r - 1d0 / amin2))
			! v_max is then the maximal *possible* speed at radius r,
			! assuming that the ejection velocity is limited by gas velocity
			vmax = sqrt(source%ud%umax * source%ud%umax &
			        + 2d0 * gm * (1d0 / point%r - 1d0 / source%r))
			
			pole = amin > amin2
			if(.not. pole) vmin = vmin2
			
			density = 0.0
			! if the maximal possible velocity is enough to get from rm to rr
			if(vmax > vmin) then
				
				if(pole) then
					call estimate_N_steps_pole_integration(Nprestep, &
								source%zeta * rad2deg, point%r_scaled, &
								angle, source%is_jet)
					v_limits(1) = vmin + veps
					v_limits(2) = (vmax - vmin) * prestep_relative_size + vmin
					vinterval = (v_limits(2) - v_limits(1))
					
					call Integrand_number_density(f1, v_limits(1), amin, &
					                             point, dphi, dbeta, source, tnow)
					viprev = v_limits(1)
					do i = 2, Nprestep
						vi = (dble(i-1) / dble(Nprestep))**trpower * vinterval + v_limits(1)
						call Integrand_number_density(f2, vi, amin, &
						                           point, dphi, dbeta, source, tnow)
						density(1) = density(1) + (vi - viprev) * 0.5d0 * (f1 + f2)
						f1 = f2
						viprev = vi
					enddo
				else
					v_limits(2) = vmin
				endif
				
				if(vc > v_limits(2)) then			
					v_limits(1) = v_limits(2)
					v_limits(2) = min(vc, vmax)
					
				! the particles on the elliptic orbits 
					call GaussLegendreQuadra(xel, wel, order_v_el)
					ldif = v_limits(2) - v_limits(1)
					ldif = ldif * 0.5d0
					lsum = v_limits(2) + v_limits(1)
					lsum = lsum * 0.5d0
					do i = 1, order_v_el
						call Integrand_number_density(term, ldif * xel(i) + lsum, amin, &
									point, dphi, dbeta, source, tnow)
						density(1) = density(1) + ldif * wel(i) * term
					enddo
				endif
					
				if(vc < vmax) then
				! the particles on the escaping trajectories
					call GaussLegendreQuadra(xhy, why, order_v_hy)
					v_limits(1) = v_limits(2)
					v_limits(2) = vmax
					ldif = v_limits(2) - v_limits(1)
					ldif = ldif * 0.5d0
					lsum = v_limits(2) + v_limits(1)
					lsum = lsum * 0.5d0
					do i = 1, order_v_hy
						call Integrand_number_density(term, ldif * xhy(i) + lsum, amin, &
									point, dphi, dbeta, source, tnow)
						density(2) = density(2) + ldif * why(i) * term
					enddo
				endif
				! factor independent on velocity
				density = density / point%r / source%r / sin(dphi)
			endif
			
		end subroutine DUDI

				
		
		
	subroutine estimate_N_steps_pole_integration(Nprestep, z, r, xi, isjet)
			implicit none
			real(8), parameter :: ximin = 0.1745329d0		! 10 degree in radians
			real(8), parameter :: ximax = 0.7853982d0		! 45 degree in radians
			integer, intent(out) :: Nprestep
			real(8), intent(in) :: xi, r, z
			logical, intent(in) :: isjet
			
			if(isjet) then
				if(ximin < xi .and. xi < ximax) then
					if(r < 2d0) then
						Nprestep = 15 + 10 * int(z)
						return
					else
						Nprestep = 10 + 5 * int(z)
						return
					endif
				else
					if(r < 1.05) then
						Nprestep = 80
						return
					else
						Nprestep = 0
						return
					endif
				endif
			else
				Nprestep = 15
			endif
			
			return
		
		end subroutine estimate_N_steps_pole_integration
		

		
    
 end module integrator

