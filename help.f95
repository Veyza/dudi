! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! File: help.f95
! Description: General auxiliary functions and subroutines

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module help
	implicit none
		contains

			! rotation of a vector in the plain
			pure function rot2d(r, phi) result(x)
				implicit none
				real(8), intent(in) :: r(2), phi
				real(8) x(2)
				
				x(1) = r(1) * cos(phi) - r(2) * sin(phi)
				x(2) = r(1) * sin(phi) + r(2) * cos(phi)
				
			end function rot2d
			
			
			
			
			! solves a system of 2 quadratic equations
			! (x - x0)**2 + (y - y0)**2 = R0**2, x0**2 + y0**2 = 1
			!(x - x1)**2 + y**2 = R1**2
			! returns 2-d vectors x and y where (x(1), y(1))
			! and (x(2), y(2)) are solutions of the SoE
			pure subroutine circle_intersection(x0, y0, R0, x1, R1, x, y)
				implicit none
				real(8), intent(in) :: x0, y0, R0, x1, R1
				real(8), intent(out) :: x(2), y(2)
				real(8) sR, dx, sx, dx2, dy2, R02, R12, x1x0
				real(8) sumdifs2, sqrtshort, ybracket
				real(8) tmp(2), eps, eps1, eps0
				
				sR = R0 + R1
				dx = x0 - x1 ; sx = x0 + x1
				x1x0 = x0 * x1
				dx2 = dx * dx ; dy2 = y0 * y0
				R02 = R0 * R0 ; R12 = R1 * R1
				
				sumdifs2 = dx2 + dy2
				
				sqrtshort = sqrt(-2d0 * (x1x0 - abs(x1)) * (sR**2 - sumdifs2) * dy2)
				
				ybracket = (R12 + sumdifs2) * dy2

				x(1) = ((R12 - R02) * dx + sx * sumdifs2 &
					- sqrtshort) / (2d0 * sumdifs2)


				y(1) = (-R02 * dy2 + (x0 - x1) * sqrtshort + ybracket) &
						/ (2d0 * sumdifs2 * y0)

				x(2) = ((R12 - R02) * dx + sx * sumdifs2 &
					+ sqrtshort) / (2d0 * sumdifs2)


				y(2) = (-R02 * dy2 + (x1 - x0) * sqrtshort + ybracket) &
						/ (2d0 * sumdifs2 * y0)
				
			end subroutine circle_intersection
			
			

		
			
			pure function vector_product(x,y) result(z)
				implicit none
				real(8), intent(in) :: x(3), y(3)
				real(8) z(3)
				
				z(1) = x(2) * y(3) - x(3) * y(2)
				z(2) = x(3) * y(1) - x(1) * y(3)
				z(3) = x(1) * y(2) - x(2) * y(1)
			
			end function vector_product
		
	
	
		
			! Linear interpolation. Returns a value of y at xout
			! x must be in ascending order
			! if xout < x(1) or xout > x(N) then 
			! yout = y(1) or yout = y(N) respectively
			pure function LiNTERPOL(N, y, x, xout) result(yout)
			  implicit none
			  real(8), intent(in) ::  xout, y(N), x(N)
			  real(8) yout, x1, x2
			  integer i, i1, i2
			  integer, intent(in) :: N
			  
			  if(xout < x(1)) then
					yout = y(1)
					return
			  endif
			  
			  if(xout > x(N)) then
					yout = y(N)
					return
			  endif
			  
			  i = 1
			  i1 = 0
			  do while(i1 == 0 .and. i .lt. N) 
			    if(xout .lt. x(i+1) .and. xout .ge. x(i)) then
			      i1 = i;       i2 = i+1;
			      x1 = x(i1);   x2 = x(i2);
			    endif
			    i = i + 1
			  enddo

			  yout = y(i1) + (y(i2) - y(i1)) / (x2 - x1) * (xout - x1)
			  return
			  
			end function LiNTERPOL	
	
	
	
			
			pure function norma3d(v)
				implicit none
				real(8), intent(in) :: v(3)
				real(8) norma3d
				
					norma3d = sqrt(dot_product(v,v))
					
			end function norma3d
			
			
			
			pure function norma2d(v)
				implicit none
				real(8), intent(in) :: v(2)
				real(8) norma2d
				
					norma2d = sqrt(dot_product(v,v))
					
			end function norma2d
				
	
			
			
		! modified inverse tanget function 
		! returns values in [0, 2pi]			
			pure function myatan(N, re0, im)
				use const
				implicit none
				integer i
				integer, intent(in):: N
				real(8), intent(in) :: re0(N), im(N)
				real(8) myatan(N), re(N)
				
				re = re0

			! avoid problems with zero re
				do i = 1, N
					if(re(i) < 1d-12 .and. re(i) > 0d0) re(i) = 1d-12 
					if(re(i) > 1d-12 .and. re(i) < 0d0) re(i) = -1d-12 
				enddo

				myatan = atan(im / re)

				do i = 1, N
					if(re(i) < 0d0) myatan(i) = myatan(i) + pi
					if(re(i) > 0d0 .and. im(i) < 0d0) myatan(i) = myatan(i) + twopi
				enddo

			end function myatan
			
			
			
			
			
		! modified inverse tanget function 
		! returns values in [0, 2pi]			
			pure function myatan1(re0, im)
				use const
				implicit none
				real(8), intent(in) :: re0, im
				real(8) re, myatan1

				re = re0

				! avoid problems with zero re
				if(re < 1d-12 .and. re > 0d0) re = 1d-12 
				if(re > 1d-12 .and. re < 0d0) re = -1d-12 

				myatan1 = atan(im / re)

				if(re < 0d0) myatan1 = myatan1 + pi
				if(re > 0d0 .and. im < 0d0) myatan1 = myatan1 + twopi

			end function myatan1
			
			
			
			! Euler's rotation
			pure subroutine eulrot(phiE, thetaE, psiE, xin, yin, zin, xout, yout, zout, inverse)
				implicit none
				integer i
				logical, intent(in) :: inverse
				real(8), intent(in) :: phiE, thetaE, psiE, xin, yin, zin
				real(8), intent(out) :: xout, yout, zout
				real(8) cp, sp, ct, st, cps, sps
				
				cp = cos(phiE)
				sp = sin(phiE)
				ct = cos(thetaE)
				st = sin(thetaE)
				cps = cos(psiE)
				sps = sin(psiE)

				if (.not. inverse) then			! if NOT INVERSE case
				
					xout= (cps*cp-ct*sp*sps)*xin + (cps*sp+ct*cp*sps)*yin + sps*st*zin
					yout=(-sps*cp-ct*sp*cps)*xin + (-sps*sp+ct*cp*cps)*yin + cps*st*zin
					zout=st*sp*xin - st*cp*yin + ct*zin	

				else 

					xout=(cps*cp-ct*sp*sps)*xin + (-sps*cp-ct*sp*cps)*yin + st*sp*zin
					yout=(cps*sp+ct*cp*sps)*xin + (-sps*sp+ct*cp*cps)*yin - st*cp*zin
					zout=st*sps*xin + st*cps*yin + ct*zin

				endif

			end subroutine eulrot




			! for given coordinates of the source and the spacecraft compute:
			! dphi (\Delta\phi) the angle between
			! the moon-centered radius-vectors of the source and the spacecraft
			! dbeta (\Delta\beta) the angle between the projections of these vectors
			! to the moon equatorial plane
			! xi the angle between the vector pointind from the source to the spacecraft
			! and the source ejection axis of symmetry
			pure subroutine ApuTrajectory(point, dphi, dbeta, xi, source)
				use const
				use define_types
				implicit none
				integer i
				type(position_in_space), intent(in) :: point
				type(source_properties), intent(in) :: source
				real(8), intent(out) :: dphi, dbeta, xi
				real(8) Rsource(3)
				
				Rsource = source%rrM / source%r
				
				! dphi is an angle between a vector pointing to the source
				! and a vector pointing to the spacecraft
				dphi = acos(dot_product(Rsource, point%rvector) / point%r)
				
				! dbeta is an angle between proections of the same vectors
				! in the longitudinal plane
				dbeta = acos((Rsource(1) * point%rvector(1) &
									+ Rsource(2) * point%rvector(2)) &
									/ (norma2d(Rsource(1:2)) * norma2d(point%rvector(1:2))))
				
				! xi is an angle between the direction of source symmetry axis
				! and the direction from the source position to the spacecraft
				xi = acos(dot_product(source%symmetry_axis, point%rvector-source%rrM) &
							/ norma3d(point%rvector-source%rrM))
				
			end subroutine ApuTrajectory	
			
			
			
			
			


			pure subroutine invert_matrix3(A, B)
				implicit none
				real(8), intent(in) :: A(3,3)
				real(8), intent(out) :: B(3,3)
				real(8) detinv

				detinv = (A(1,1) * A(2,2) * A(3,3) - A(1,1) * A(2,3)  * A(3,2) &
						- A(1,2) * A(2,1) * A(3,3) + A(1,2) * A(2,3) * A(3,1) &
						+ A(1,3) * A(2,1) * A(3,2) - A(1,3) * A(2,2) * A(3,1))
				detinv = 1d0 / detinv

				B(1,1) = detinv * (A(2,2) * A(3,3) - A(2,3) * A(3,2))
				B(2,1) = - detinv * (A(2,1) * A(3,3) - A(2,3) * A(3,1))
				B(3,1) = detinv * (A(2,1) * A(3,2) - A(2,2) * A(3,1))
				B(1,2) = - detinv * (A(1,2) * A(3,3) - A(1,3) * A(3,2))
				B(2,2) = detinv * (A(1,1) * A(3,3) - A(1,3) * A(3,1))
				B(3,2) = - detinv * (A(1,1) * A(3,2) - A(1,2) * A(3,1))
				B(1,3) = detinv * (A(1,2) * A(2,3) - A(1,3) * A(2,2))
				B(2,3) = - detinv * (A(1,1) * A(2,3) - A(1,3) * A(2,1))
				B(3,3) = detinv * (A(1,1) * A(2,2) - A(1,2) * A(2,1))
				
			end subroutine invert_matrix3




			! returns the shortest distance (d) between 2 lines in 3d
			! and the closest point to the 2nd line on the 1st line
			! and the closest point to the 1st line on the 2nd line
			! M1 and M2 are points through which 1st and 2nd lines
			! go and s1 and s2 are pointing unit-vectors of the 1st
			! and the 2nd lines respectively
			subroutine dist_between_2lines(d, K1, K2, M1, s1, M2, s2)
				implicit none
				real(8), intent(out) :: d, K1(3), K2(3)
				real(8), intent(in) :: M1(3), s1(3), M2(3), s2(3)
				real(8) r(3), s(3), t1(3), t2(3)
				real(8) A(3,3), B(3), Ainv(3,3), A0(3,3), tmp(3)
				integer i
				
				s = vector_product(s1, s2)		! s is normal to both s1 and s2
				if(norma3d(s) < 1d-3) then
					! Lines are parallel to each other
					K2 = M2
					K1 = M1 + s1 * norma3d(M2 - M1)
					d = norma3d(vector_product(s1, M2 - M1))
					return
				endif
				
				t1 = vector_product(s1, s)				! t1 is normal to s and s1
				t2 = vector_product(s2, s)				! t2 is normal to s and s2
				! t1(1)*x + t1(2)*y + t1(3)*z = t1*M1 - equation of a plain (1)
				! which is parallel to s and contains line 1
				A(1,:) = t1
				B(1) = dot_product(t1, M1)
				! s2(2)*x - s2(1)*y = s2(2)*M2(1) - s2(1)*M2(2)  (2)
				A(2,:) = (/s2(2), -s2(1), 0d0/)
				B(2) = -(s2(1) * M2(2) - s2(2) * M2(1))
				! s2(3)*x - s2(3)*y = s2(3)*M2(1) - s2(1)*M2(3)  (3)
				A(3,:) = (/s2(3), 0d0, -s2(1)/)
				B(3) = -(s2(1) * M2(3) - s2(3) * M2(1))
				! together (2) and (3)  define the line 2 (s2, M2)
				! K2 is the point where line 2 intersects the plain containing
				! line 1 and parallel to common perpendicular of lines 1 and 2 (vector s)
				call invert_matrix3(A, Ainv)
				K2(1) = sum(Ainv(1,:) * B)
				K2(2) = sum(Ainv(2,:) * B)
				K2(3) = sum(Ainv(3,:) * B)
				
				A(1,:) = t2
				B(1) = dot_product(t2, M2)
				
				A(2,:) = -(/s1(2), -s1(1), 0d0/)
				B(2) = s1(1) * M1(2) - s1(2) * M1(1)
				A(3,:) = -(/s1(3), 0d0, -s1(1)/)
				B(3) = s1(1) * M1(3) - s1(3) * M1(1)
				! K1 is the point where line 1 intersects the plain containing
				! line 1 and parallel to common perpendicular of lines 1 and 2 (vector s)
				call invert_matrix3(A, Ainv)
				K1(1) = sum(Ainv(1,:) * B)
				K1(2) = sum(Ainv(2,:) * B)
				K1(3) = sum(Ainv(3,:) * B)
				
				! two plains parallel to vector s but not parallel to each other
				! intersect defining the line parallel to s,
				! the common perpendicular for lines 1 and 2
				
				d = norma3d(K1 - K2)

							
			end subroutine dist_between_2lines



			
			
			! returns nodes and weights of the Gauss-Legendre quadrature formula
			subroutine GaussLegendreQuadra(xi, wi, order)
				use const
				implicit none
				integer, intent(in) :: order 
				real, intent(out) :: xi(order), wi(order)
				
				if(order /= 5 .and. order /= 10 &
							.and. order /= 20 .and. order /= 30 .and. order /= 40 &
							.and. order /= 50) then
					write(*,*) 'order of integration must be equal to 5, 10, 20, 30, 40 or 50'
					write(*,*) 'check variables order_v and order_R'
					stop
				endif
				
				select case(order)
					case(50)
						 xi(1) = -0.9988664044;    wi(1) = 0.0029086226
						 xi(2) = -0.9940319694;    wi(2) = 0.0067597992
						 xi(3) = -0.9853540841;    wi(3) = 0.0105905484
						 xi(4) = -0.9728643851;    wi(4) = 0.0143808228
						 xi(5) = -0.9566109552;    wi(5) = 0.0181155607
						 xi(6) = -0.9366566189;    wi(6) = 0.0217802432
						 xi(7) = -0.9130785567;    wi(7) = 0.0253606736
						 xi(8) = -0.8859679795;    wi(8) = 0.0288429936
						 xi(9) = -0.8554297694;    wi(9) = 0.032213728
						 xi(10) = -0.8215820709;   wi(10) = 0.035459836
						 xi(11) = -0.7845558329;   wi(11) = 0.0385687566
						 xi(12) = -0.7444943022;   wi(12) = 0.0415284631
						 xi(13) = -0.7015524687;   wi(13) = 0.0443275043
						 xi(14) = -0.6558964657;   wi(14) = 0.0469550513
						 xi(15) = -0.6077029272;   wi(15) = 0.0494009384
						 xi(16) = -0.5571583045;   wi(16) = 0.0516557031
						 xi(17) = -0.5044581449;   wi(17) = 0.0537106219
						 xi(18) = -0.449806335;    wi(18) = 0.0555577448
						 xi(19) = -0.393414312;    wi(19) = 0.05718992565
						 xi(20) = -0.335500245;    wi(20) = 0.05860084981
						 xi(21) = -0.2762881938;   wi(21) = 0.0597850587
						 xi(22) = -0.2160072369;   wi(22) = 0.0607379708
						 xi(23) = -0.15489059;     wi(23) = 0.0614558996
						 xi(24) = -0.0931747016;   wi(24) = 0.0619360674
						 xi(25) = -0.0310983383;   wi(25) = 0.0621766167
						 xi(26) = 0.0310983383;    wi(26) = 0.0621766167
						 xi(27) = 0.0931747016;    wi(27) = 0.0619360674
						 xi(28) = 0.15489059 ;     wi(28) = 0.0614558996
						 xi(29) = 0.2160072369;    wi(29) = 0.0607379708
						 xi(30) = 0.2762881938;    wi(30) = 0.0597850587
						 xi(31) = 0.335500245;     wi(31) = 0.0586008498
						 xi(32) = 0.3934143119;    wi(32) = 0.057189926
						 xi(33) = 0.449806335;     wi(33) = 0.0555577448
						 xi(34) = 0.504458145;     wi(34) = 0.0537106219
						 xi(35) = 0.5571583045;    wi(35) = 0.051655703
						 xi(36) = 0.6077029272;    wi(36) = 0.04940093845
						 xi(37) = 0.6558964657;    wi(37) = 0.0469550513
						 xi(38) = 0.7015524687;    wi(38) = 0.044327504
						 xi(39) = 0.7444943022;    wi(39) = 0.0415284631
						 xi(40) = 0.7845558329;    wi(40) = 0.03856875661
						 xi(41) = 0.8215820709;    wi(41) = 0.0354598356
						 xi(42) = 0.8554297694;    wi(42) = 0.0322137282
						 xi(43) = 0.8859679795;    wi(43) = 0.0288429936
						 xi(44) = 0.9130785567;    wi(44) = 0.025360674
						 xi(45) = 0.9366566189;    wi(45) = 0.021780243
						 xi(46) = 0.9566109552;    wi(46) = 0.018115561
						 xi(47) = 0.9728643851;    wi(47) = 0.0143808228
						 xi(48) = 0.9853540841;    wi(48) = 0.0105905484
						 xi(49) = 0.9940319694;    wi(49) = 0.006759799
						 xi(50) = 0.9988664044;    wi(50) = 0.00290862255
					case(40)
						 xi(1) = -0.998237710;     wi(1) = 0.0045212771
						 xi(2) = -0.990726239;     wi(2) = 0.0104982845
						 xi(3) = -0.977259950;     wi(3) = 0.0164210584
						 xi(4) = -0.957916819;     wi(4) = 0.0222458492
						 xi(5) = -0.932812808;     wi(5) = 0.0279370070
						 xi(6) = -0.902098807;     wi(6) = 0.0334601953
						 xi(7) = -0.865959503;     wi(7) = 0.0387821680
						 xi(8) = -0.824612231;     wi(8) = 0.0438709082
						 xi(9) = -0.778305651;     wi(9) = 0.0486958076
						 xi(10) = -0.727318255;     wi(10) = 0.053227847
						 xi(11) = -0.671956685;     wi(11) = 0.057439769
						 xi(12) = -0.612553890;     wi(12) = 0.061306242
						 xi(13) = -0.549467125;     wi(13) = 0.064804013
						 xi(14) = -0.483075802;     wi(14) = 0.067912046
						 xi(15) = -0.413779204;     wi(15) = 0.070611647
						 xi(16) = -0.341994091;     wi(16) = 0.072886582
						 xi(17) = -0.268152185;     wi(17) = 0.074723169
						 xi(18) = -0.192697581;     wi(18) = 0.076110362
						 xi(19) = -0.116084071;     wi(19) = 0.077039818
						 xi(20) = -0.038772418;     wi(20) = 0.077505948
						 xi(21) = 0.0387724175;     wi(21) = 0.077505948
						 xi(22) = 0.1160840707;     wi(22) = 0.077039818
						 xi(23) = 0.1926975807;     wi(23) = 0.076110362
						 xi(24) = 0.2681521850;     wi(24) = 0.074723169
						 xi(25) = 0.3419940908;     wi(25) = 0.072886582
						 xi(26) = 0.4137792044;     wi(26) = 0.070611647
						 xi(27) = 0.4830758017;     wi(27) = 0.067912046
						 xi(28) = 0.5494671251;     wi(28) = 0.064804013
						 xi(29) = 0.6125538897;     wi(29) = 0.061306242
						 xi(30) = 0.6719566846;     wi(30) = 0.057439769
						 xi(31) = 0.7273182552;     wi(31) = 0.053227847
						 xi(32) = 0.7783056514;     wi(32) = 0.048695808
						 xi(33) = 0.8246122308;     wi(33) = 0.043870908
						 xi(34) = 0.8659595032;     wi(34) = 0.038782168
						 xi(35) = 0.9020988070;     wi(35) = 0.033460195
						 xi(36) = 0.9328128083;     wi(36) = 0.027937007
						 xi(37) = 0.9579168192;     wi(37) = 0.022245849
						 xi(38) = 0.9772599500;     wi(38) = 0.016421058
						 xi(39) = 0.9907262387;     wi(39) = 0.010498285
						 xi(40) = 0.9982377097;     wi(40) = 0.004521277
					case(30)
						 xi(1) = -0.996893466;     wi(1) = 7.96819292E-03
						 xi(2) = -0.983668149;     wi(2) = 1.84664689E-02
						 xi(3) = -0.960021853;     wi(3) = 2.87847072E-02
						 xi(4) = -0.926200032;     wi(4) = 3.87991928E-02
						 xi(5) = -0.882560551;     wi(5) = 4.84026745E-02
						 xi(6) = -0.829565763;     wi(6) = 5.74931577E-02
						 xi(7) = -0.767777443;     wi(7) = 6.59742281E-02
						 xi(8) = -0.697850466;     wi(8) = 7.37559721E-02
						 xi(9) = -0.620526195;     wi(9) = 8.07558969E-02
						 xi(10) = -0.536624134;    wi(10) = 8.68997872E-02
						 xi(11) = -0.447033763;    wi(11) = 9.21225250E-02
						 xi(12) = -0.352704734;    wi(12) = 9.63687375E-02
						 xi(13) = -0.254636914;    wi(13) = 9.95934233E-02
						 xi(14) = -0.153869912;    wi(14) = 0.101762392    
						 xi(15) = -5.14718443E-02; wi(15) = 0.102852650    
						 xi(16) = 5.14718443E-02;  wi(16) = 0.102852650    
						 xi(17) = 0.153869912;     wi(17) = 0.101762392    
						 xi(18) = 0.254636914;     wi(18) = 9.95934233E-02
						 xi(19) = 0.352704734;     wi(19) = 9.63687375E-02
						 xi(20) = 0.447033763;     wi(20) = 9.21225250E-02
						 xi(21) = 0.536624134;     wi(21) = 8.68997872E-02
						 xi(22) = 0.620526195;     wi(22) = 8.07558969E-02
						 xi(23) = 0.697850466;     wi(23) = 7.37559721E-02
						 xi(24) = 0.767777443;     wi(24) = 6.59742281E-02
						 xi(25) = 0.829565763;     wi(25) = 5.74931577E-02
						 xi(26) = 0.882560551;     wi(26) = 4.84026745E-02
						 xi(27) = 0.926200032;     wi(27) = 3.87991928E-02
						 xi(28) = 0.960021853;     wi(28) = 2.87847072E-02
						 xi(29) = 0.983668149;     wi(29) = 1.84664689E-02
						 xi(30) = 0.996893466;     wi(30) = 7.96819292E-03
					case(20)
						xi(1) = -0.993128598;     wi(1) = 1.76140070E-02
						xi(2) = -0.963971913;     wi(2) = 4.06014286E-02
						xi(3) = -0.912234426;     wi(3) = 6.26720488E-02
						xi(4) = -0.839116991;     wi(4) = 8.32767412E-02
						xi(5) = -0.746331930;     wi(5) = 0.101930119    
						xi(6) = -0.636053681;     wi(6) = 0.118194535    
						xi(7) = -0.510867000;     wi(7) = 0.131688640    
						xi(8) = -0.373706102;     wi(8) = 0.142096102    
						xi(9) = -0.227785856;     wi(9) = 0.149172992    
						xi(10) = -7.65265226E-02; wi(10) = 0.152753383    
						xi(11) = 7.65265226E-02;  wi(11) = 0.152753383    
						xi(12) = 0.227785856;     wi(12) = 0.149172992    
						xi(13) = 0.373706102;     wi(13) = 0.142096102    
						xi(14) = 0.510867000;     wi(14) = 0.131688640    
						xi(15) = 0.636053681;     wi(15) = 0.118194535    
						xi(16) = 0.746331930;     wi(16) = 0.101930119    
						xi(17) = 0.839116991;     wi(17) = 8.32767412E-02
						xi(18) = 0.912234426;     wi(18) = 6.26720488E-02
						xi(19) = 0.963971913;     wi(19) = 4.06014286E-02
						xi(20) = 0.993128598;     wi(20) = 1.76140070E-02
					case(10)
						xi(1) = -0.973906517;     wi(1) = 6.66713417E-02
						xi(2) = -0.865063369;     wi(2) = 0.149451345    
						xi(3) = -0.679409564;     wi(3) = 0.219086364    
						xi(4) = -0.433395386;     wi(4) = 0.269266725    
						xi(5) = -0.148874342;     wi(5) = 0.295524240    
						xi(6) = 0.148874342;      wi(6) = 0.295524240    
						xi(7) = 0.433395386;      wi(7) = 0.269266725    
						xi(8) = 0.679409564;      wi(8) = 0.219086364    
 						xi(9) = 0.865063369;      wi(9) = 0.149451345    
						xi(10) = 0.973906517;     wi(10) = 6.66713417E-02
					case(5)
						xi(1) = -0.906179845;    wi(1) = 0.236926883    
						xi(2) = -0.538469315;    wi(2) = 0.478628665    
						xi(3) = 0.00000000;      wi(3) = 0.568888903    
						xi(4) = 0.538469315;     wi(4) = 0.478628665    
						xi(5) = 0.906179845;     wi(5) = 0.236926883 
				
				endselect
			
			
			end subroutine GaussLegendreQuadra
			
			

end module help
