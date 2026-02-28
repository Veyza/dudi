module BGmod
	implicit none
		contains



		! function integrateBG calculates number density of particles with sizes between r1 and r2 in the E ring
		! assuming that we know nnorm which is number density of particles with sizes between rnorm1 and rnorm2
		! and the size-distribution of all the particles in the E ring which is given as a table in the file 'background_size_distribution.dat'
		! r1, r1, rnorm1, rnorm2, r and rstep given in microns
			function integrateBG(r1, r2, rnorm1, rnorm2, nnorm) result(ndout)
				use help
				use const
				implicit none
				real(8) rstep
				real(8) r1, r2, rnorm1, rnorm2, nnorm
				character(len = 55), parameter :: fname = './input_data_files/background_size_distribution.dat'
				integer, parameter :: N = 36
				real(8) r(N), distr(N), ndout, ndtmp
				integer i, imin, imax
				character(len = 65) str

				rstep = (r2 - r1) / GRN
				open(700, file = fname, status = 'old')
					do i = 1, 3
						read(700,*) str
						write(*,*) str
					enddo
					do i = 1, N
						read(700,*) r(i), distr(i)				! reading values of a tabulated function
					enddo

				close(700)

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				r = r * 2d0
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				call findSubsetInd(r, N, rnorm1, rnorm2, imin, imax)

				ndtmp = 0d0
				do i = imin, imax
					ndtmp = ndtmp + distr(i) * rstep			! integration of the tabulated function by rectangles method
				enddo

				call findSubsetInd(r, N, r1, r2, imin, imax)

				ndout = 0d0
				do i = imin, imax
					ndout = ndout + distr(i) * rstep
				enddo

				ndout = ndout / ndtmp * nnorm
				write(*,*) 'BG = ', ndout

			end function integrateBG


			! finds the index of the first element in the array arr which is greater than bottom
			! and the index of the last element of the array arr which is smaller than top
			! arr is assumed sorted from smallest to largest number
			subroutine findSubsetInd(arr, N, bottom, top, imin, imax)
				implicit none
				integer, intent(in) :: N
				real(8), intent(in) :: arr(N), bottom, top
				integer, intent(out) :: imin, imax

				imin = 1 ; imax = N
				do while(arr(imin) < bottom)
					imin = imin + 1
				enddo
				if(imin > 1) imin = imin - 1
				do while(arr(imax) > top)
					imax = imax - 1
				enddo
				if(imax < N) imax = imax + 1

			end subroutine findSubsetInd

end module BGmod
