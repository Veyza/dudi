! This file is a part of DUDI, the Fortran-90 implementation
! of the two-body model for dust dynamics
! Version 1.2.3
! This is free software. You can use and redistribute it
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt,
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186

! File: comparison_utils.f90
! Description: Contains functions performing comparisons of real(8) variables.
! This is FORTRAN-95 way to do such operations, newer FORTRAN versions have
! inbuilt functions for them.

! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com

module comparison_utils
  implicit none
contains

  pure logical function is_nan_r8(x)
    implicit none
    real(8), intent(in) :: x
    ! Fortran 95 portable NaN check: NaN is not equal to itself
    is_nan_r8 = (x /= x)
  end function is_nan_r8

  pure logical function is_finite_r8(x)
    implicit none
    real(8), intent(in) :: x
    ! finite <=> not NaN and not infinite
    is_finite_r8 = .not.(x /= x) .and. (abs(x) < huge(0.0d0))
  end function is_finite_r8

  elemental pure logical function is_zero_r8(x, atol)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in), optional :: atol
    real(8) :: a0
    if (present(atol)) then
      a0 = atol
    else
      a0 = 1d-14
    end if
    is_zero_r8 = abs(x) <= a0
  end function is_zero_r8


  ! helper that works well across scales
  pure logical function nearly_equal_r8(a,b,rtol,atol) result(eq)
    implicit none
    real(8), intent(in) :: a, b
    real(8), intent(in), optional :: rtol, atol
    real(8) :: r, a0, m
    if (present(rtol)) then
      r = rtol
    else
      r = 1d-12
    end if
    if (present(atol)) then
      a0 = atol
    else
      a0 = 1d-14
    end if
    m  = max(1d0, abs(a), abs(b))
    eq = abs(a-b) <= max(a0, r*m)
  end function nearly_equal_r8


end module comparison_utils
