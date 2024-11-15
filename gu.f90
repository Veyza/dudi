! This file is a part of DUDI, the Fortran-90 implementation 
! of the two-body model for dust dynamics
! Version 1.1.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 

! File: gu.f90
! Description: The subroutines that performe integration over the observed
!              interval of particle sizes

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com
module Gu
	implicit none
		contains
		
		! If rlim1 + 1d-16  < rlim2, Gu_integral integrates the product
		! of size distribution and ejection speed distribution (possibly
		! size-dependent) over the interval (rlim1, rlim2). This integral is 
		! denoted G^p_u in the paper (Ershova & Schmidt, 2021).
		! The function G^p_u is evaluated for the array of u-values (stored in ui)
		! there are GRN values of ui uniformly distributed between umin and umax
		! If rlim1 + 1d-16 >= rlim2, Gu_integral computes an array with 
		! product of size distribution and ejection speed distribution
		! evaluated at rlim2 and an array of u-values.
			subroutine Gu_integral(ui, S, sd, ud, rlim1, rlim2)
				use const
				use define_types
				use help
				use distributions_fun
				real(8) x, cc, ldif, lsum
				real xi(order_R), wi(order_R)
				real(8), intent(in) :: rlim1, rlim2
				integer, intent(in) :: sd
				type(ejection_speed_properties), intent(in) :: ud
				real(8), intent(out) :: ui(GRN), S(GRN)
				real(8) f1, f2
				integer i, ii
				real(8) up, low, tmp, integral_limits(4), r_i, C_size_distr
				
				forall(i = 0:(GRN-1)) ui(i+1) = ud%umin &
				                     + dble(i) / dble(GRN-1) * (ud%umax - ud%umin)
				if(rlim2 - rlim1 < 0d0) then
					write(*,*) 'size limits are likely wrong:', rlim1, rlim2 
					stop
				endif
				! if a range of particle radii is considered
				if(rlim2 - rlim1 > 1e-16) then			
					up = rlim2
					low = rlim1
					S = 0d0
					
					call GaussLegendreQuadra(xi, wi, order_R)
					do i = 1, GRN
						ldif = up - low
						ldif = ldif * 0.5d0
						lsum = up + low
						lsum = lsum * 0.5d0					
						do ii = 1, order_R
							r_i = ldif * xi(ii) + lsum
							S(i) = S(i) + ldif * wi(ii) &
								* size_distribution(r_i, sd) * r_i**p &
								* ejection_speed_distribution(ud, ui(i), r_i)
						enddo
						if(S(i) < 0 .or. S(i) /= S(i)) then
							write(*,*) 'Gu has an incorrect value of ', S(i), sd, i
							stop
						endif
					enddo
				! if no integration over grain size is required
				else
					do i = 1, GRN
						S(i) = size_distribution(rlim1, sd) * rlim1**p &
								* ejection_speed_distribution(ud, ui(i), rlim1)
					enddo
				endif
				! S * ice density [kg/m^3] * 1d-18 m^3/micron^3 * 4/3 pi, this is to obtain mass in kilogramms
				if(p == 3) S = S * rho * 1d-18 * 4d0 / 3d0 * pi
				! S * pi * 1d-12 [m^2 / micron^2] to obtain area of crossection in m^2
				if(p == 2) S = S * pi * 1d-12
								
			end subroutine Gu_integral
			
			
			
			! computes the mass confined between r1 and r2
			! resuls is in kg
			subroutine mass_production(mass, sd, r1, r2)
				use const
				use distributions_fun
				use help
				integer, intent(in) :: sd
				real(8), intent(in) :: r1, r2
				integer i
				real(8) mass, ldif, lsum, r_i
				real xi(order_R), wi(order_R)
				
				mass = 0d0
				call GaussLegendreQuadra(xi, wi, order_R)
				ldif = r2 - r1
				ldif = ldif * 0.5d0
				lsum = r2 + r1
				lsum = lsum * 0.5d0	
				
				do i = 1, order_R
					r_i = ldif * xi(i) + lsum
					mass = mass + ldif * wi(i) * size_distribution(r_i, sd) &
					       * r_i**3
				enddo
				
				! upon integration we obtained the volume of particles as if they were cubes
				! we multiply it by 4pi/3 to obtain the volume of spherical particles
				! 1d-18 is the factor for the conversion from microns^3 to meters^3
				! then, by multiplying by density rho, we obtain the mass
				mass = mass * rho * 1d-18 * 4d0 / 3d0 * pi
			
			end subroutine mass_production
			
			

end module Gu
