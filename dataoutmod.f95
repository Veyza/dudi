! This file is a part of DUDI, the Fortran-95 implementation 
! of the two-body model for dust dynamics
! Version 1.0.0
! This is free software. You can use and redistribute it 
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt, 
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186 
! File: dataoutmod.f95
! Description: The subroutines used to write the result to a text-file

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module dataoutmod
	implicit none
		contains
			
			subroutine result_out(density, nt, point)
				use const
				use define_types
				implicit none
				integer, intent(in) :: nt
				type(position_in_space), intent(in) :: point(nt)
				real, intent(in) :: density(nt,2)
				integer i
				
				open(111, file = './results/twobody_model_result.dat', &
				          status = 'replace')
					do i = 1, nt
						write(111,*) density(i,:), density(i,:), &
							point(i)%r, (halfpi - point(i)%alpha) * rad2deg, &
							point(i)%beta * rad2deg
					enddo
				close(111)
			
			end subroutine result_out
			
			
			
			subroutine cassini_flyby_out(density, ttab, bg, nt)
				implicit none
				integer, intent(in) :: nt
				real, intent(in) :: density(nt,2), ttab(nt), bg
				integer i
				
				open(111, file = './results/E2_profile.dat', status = 'replace')
					do i = 1, nt
						write(111,'(f9.3, x, e11.4)') ttab(i), sum(density(i,:)) + bg
					enddo
				close(111)
			
			end subroutine cassini_flyby_out	
			
			
			
			
			subroutine surface_deposition_out(num, deposition, nt, dphi)
				use const
				implicit none
				integer, intent(in) :: num, nt
				real, intent(in) :: deposition(nt)
				real(8), intent(in) :: dphi(nt)
				integer i
				character(len = 39) fname
				
				select case(num)
					case(1)
						fname = './results/narrow_jet_shallow_sd.dat'
					case(2)
						fname = './results/diffuse_source_steep_sd.dat'
					case(3)
						fname = './results/diffuse_source_shallow_sd.dat'
					case(4)
						fname = './results/narrow_jet_steep_sd.dat'
				endselect
				
				open(111, file = fname, status = 'replace')
					do i = 1, nt
						write(111,*) dphi(i) * rm * 1d-3, deposition(i)
					enddo
				close(111)
			
			end subroutine surface_deposition_out		
			
			
		
			

end module dataoutmod
