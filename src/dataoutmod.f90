! This file is a part of DUDI, the Fortran-90 implementation
! of the two-body model for dust dynamics
! Version 1.2.3
! This is free software. You can use and redistribute it
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt,
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186

! File: dataoutmod.f90
! Description: The subroutines used to write the result to a text-file

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module dataoutmod
	implicit none
	contains

		  ! A helper subroutine to remove leading and trailing blanks from formatted_value
		  subroutine trim_and_format(value)
			character(len=*), intent(inout) :: value
			integer :: start_pos, end_pos
			character(len=len(value)) :: trimmed_value

			! Find the first non-blank character
			start_pos = 1
			do while (value(start_pos:start_pos) == ' ')
			   start_pos = start_pos + 1
			end do

			! Find the last non-blank character
			end_pos = len_trim(value)

			! Extract the trimmed portion of the string
			trimmed_value = value(start_pos:end_pos)

			! Reassign the trimmed value back to the input
			value = trimmed_value
		  end subroutine trim_and_format



			subroutine vertical_slicematrix_out(resmat, nt, fnum, cellsize)
				real, intent(in) :: resmat(nt, nt)
				integer, intent(in) :: nt
				real, intent(in) :: fnum
				real(8), intent(in) :: cellsize
				integer i
				character(len = 80) fname

				if(fnum == 0.1) then
					if(cellsize < 2.9d3) then
						fname = './results/dens_in_plane_dust_a-panel.dat'
					else
						fname = './results/dens_in_plane_dust_b-panel.dat'
					endif
					open(111, file = fname, status = 'replace')
				endif
				if(fnum == 0.2) then
					if(cellsize < 2.9d3) then
						fname = './results/dens_in_plane_gas_a-panel.dat'
					else
						fname = './results/dens_in_plane_gas_b-panel.dat'
					endif
					open(111, file = fname, status = 'replace')
				endif
					do i = 1, nt
						write(111,*) resmat(i,:)
					enddo
				close(111)


			end subroutine vertical_slicematrix_out



			subroutine composition_matrix_out(salt_poor, salt_rich, nt, fnum, cellsize)
				real, intent(in) :: salt_poor(nt, nt), salt_rich(nt, nt)
				integer, intent(in) :: nt
				real, intent(in) :: fnum
				real(8), intent(in) :: cellsize
				integer i
				character(len = 80) fname_poor, fname_rich

				if(fnum == 0.4) then
					if(cellsize < 2.9d3) then
						fname_poor = './results/salt_poor_plane_a-panel.dat'
						fname_rich = './results/salt_rich_plane_a-panel.dat'
					else
						fname_poor = './results/salt_poor_plane_b-panel.dat'
						fname_rich = './results/salt_rich_plane_b-panel.dat'
					endif
					open(111, file = fname_poor, status = 'replace')
						do i = 1, nt
							write(111,*) salt_poor(i,:)
						enddo
					close(111)
					open(131, file = fname_rich, status = 'replace')
						do i = 1, nt
							write(131,*) salt_rich(i,:)
						enddo
					close(131)
				endif
			end subroutine composition_matrix_out


			subroutine cassini_flyby_out(type1, type2, type3, ttab, nt, fname)
				implicit none
				integer, intent(in) :: nt
				real, intent(in) :: type1(nt), type2(nt), type3(nt), ttab(nt)
				character(*), intent(in) :: fname
				real ndens(nt)
				integer i

				ndens = type1 + type2 + type3
				open(111, file = fname, status = 'replace')
					do i = 1, nt
						write(111,*) ttab(i), ndens(i), type1(i) / ndens(i), type2(i) / ndens(i), type3(i) / ndens(i)
					enddo
				close(111)
				write(*,*) 'result is in the file', fname

			end subroutine cassini_flyby_out






			subroutine matrix_out(fname, image, nt1, nt2)
				use const
				implicit none
				integer, intent(in) :: nt1, nt2
				real, intent(in) :: image(nt1,nt2)
				integer, parameter :: outchannel = 115
				integer i
				character(*), intent(in) :: fname
				character(len = 20) outformat
				character(len = 2) tmp1, tmp2
				character(len = 12) tmp0
				character(len = 18) tmpformat


				tmp1 = 'I4'; tmp2 = 'I4'
				if(nt1*2 < 1000) tmp1 = 'I3'
				if(nt1*2 < 100) tmp1 = 'I2'
				if(nt2*2 < 1003) tmp2 = 'I3'
				if(nt2*2 < 100) tmp2 = 'I2'
				tmpformat = '(' // tmp1 // ', A1,' // tmp2 // ')'
				write(outformat,'(A1,I3,A1,I3,A11)') '(', nt1, 'x',nt2, '(ES12.4E2))'


				open(outchannel, file = fname, status = 'replace')

					do i = 1, nt1
						write(outchannel,outformat) image(i,:)
					enddo

				close(outchannel)
				write(*,*) 'result is in the file ', fname

             end subroutine matrix_out


end module dataoutmod
