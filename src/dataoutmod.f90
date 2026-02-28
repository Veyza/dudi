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

		  subroutine alt_decrease_out(fname, density, rs, alphas, Nalts, Nlats, Nlons)
			use const
			implicit none
			character(*), intent(in) :: fname
			integer, intent(in) :: Nalts, Nlats, Nlons
			real, intent(in) :: density(Nalts, Nlats, Nlons)
			real(8), intent(in) :: rs(Nalts), alphas(Nlats)
			real aver_at_circle(Nalts, Nlats)
			integer i, ii

			forall(i = 1:Nalts)
			forall(ii = 1:Nlats)
				aver_at_circle(i,ii) = sum(density(i,ii,:)) / Nlons
			endforall
			endforall

			open(111, file = fname, status = 'replace')
				write(111,*) 999, (halfpi - alphas) * rad2deg
				do i = 1, Nalts
					write(111,*) rs(i), aver_at_circle(i,:)
				enddo
			close(111)


		  end subroutine alt_decrease_out


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
