! This file is a part of DUDI, the Fortran-95 implementation
! of the two-body model for the dynamics of dust ejected from an
! atmosphereless body.
! Version 1.2.x
! This is free software. You can use and redistribute it
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt,
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186
!
! Author: Anastasiia Ershova
! E-mail: vveyzaa@gmail.com
!
! python_interface/fortran_bridge/py_bridge.f90
!
! C-interoperable wrappers around the DUDI integrator.
! - Only interoperable scalars/1D arrays cross the boundary.
! - We reconstruct derived types locally and call the real routines.
! - Fortran LOGICAL is replaced by INTEGER(C_INT) on the C boundary.
! - DUDI returns default REAL; we convert to REAL(8) on output.

module py_dudi_bridge
  use iso_c_binding,   only: c_int, c_double
  use iso_fortran_env, only: real64
  use const
  use define_types, only: &
       position_in_space, &
       source_properties, &
       ejection_speed_properties
  use integrator,  only: DUDI, DUDI_mean_velocity
  use inputdata,   only: jet_direction
  use Gu,          only: Gu_integral
  use batching,    only: DUDI_batch_points, DUDI_batch_sources, DUDI_batch_sources_points, &
                        DUDI_mean_velocity_batch_points, DUDI_mean_velocity_batch_sources, &
                        DUDI_mean_velocity_batch_sources_points
  implicit none

contains

  ! Helper: environment-controlled diagnostics (kept for parity with DUDI-hc)
  logical function env_enabled(name)
    implicit none
    character(len=*), intent(in) :: name
    character(len=8) :: val
    integer :: stat, lenval

    env_enabled = .false.
    call get_environment_variable(name, value=val, length=lenval, status=stat, trim_name=.true.)
    if (stat == 0) then
       if (lenval > 0) then
          ! treat any non-'0' first char as enabled
          env_enabled = (val(1:1) /= '0')
       end if
    end if
  end function env_enabled


  !--------------------------------------------------------------------
  !  Single-point density wrapper around DUDI
  !
  !  C interface (Python ctypes prototype):
  !
  !  void py_dudi_density(
  !      double point_r,
  !      double point_alpha,
  !      double point_beta,
  !      double src_r,
  !      double src_alphaM,
  !      double src_betaM,
  !      double src_zeta,
  !      double src_eta,
  !      int    src_eject_distr,
  !      int    src_ud_shape,
  !      double src_umin,
  !      double src_umax,
  !      int    src_sd,
  !      int    src_production_fun,
  !      double src_production_rate,
  !      int    src_is_jet,          ! 0/1
  !      double tnow,
  !      double density_out[2]       ! bound, unbound
  !  );
  !--------------------------------------------------------------------
  subroutine py_dudi_density( &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, &
       tnow, density_out) &
       bind(C, name="py_dudi_density")

    use iso_fortran_env, only: real64
    implicit none

    ! C boundary types
    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM
    real(c_double), value, intent(in) :: src_zeta, src_eta
    integer(c_int),  value, intent(in) :: src_eject_distr
    integer(c_int),  value, intent(in) :: src_ud_shape
    real(c_double),  value, intent(in) :: src_umin, src_umax
    integer(c_int),  value, intent(in) :: src_sd
    integer(c_int),  value, intent(in) :: src_production_fun
    real(c_double),  value, intent(in) :: src_production_rate
    integer(c_int),  value, intent(in) :: src_is_jet
    real(c_double),  value, intent(in) :: tnow

    real(c_double), intent(out) :: density_out(2)

    ! Local Fortran-side types
    type(position_in_space)        :: point
    type(source_properties)        :: source
    type(ejection_speed_properties):: ud
    real                           :: density_local(2)

    ! Build POINT (spherical coordinates in radians, radius in meters)
    point%r     = real(point_r,   kind=real64)
    point%alpha = real(point_alpha, kind=real64)
    point%beta  = real(point_beta,  kind=real64)
    point%rvector(1) = point%r * sin(point%alpha) * cos(point%beta)
    point%rvector(2) = point%r * sin(point%alpha) * sin(point%beta)
    point%rvector(3) = point%r * cos(point%alpha)
    point%r_scaled   = point%r / rm
    point%compute    = .true.

    ! Build SOURCE
    source%r      = real(src_r,      kind=real64)
    source%alphaM = real(src_alphaM, kind=real64)
    source%betaM  = real(src_betaM,  kind=real64)
    source%zeta   = real(src_zeta,   kind=real64)
    source%eta    = real(src_eta,    kind=real64)

    source%rrM(1) = source%r * sin(source%alphaM) * cos(source%betaM)
    source%rrM(2) = source%r * sin(source%alphaM) * sin(source%betaM)
    source%rrM(3) = source%r * cos(source%alphaM)

    source%ejection_angle_distr = int(src_eject_distr, kind=kind(source%ejection_angle_distr))
    ud%ud_shape   = int(src_ud_shape, kind=kind(ud%ud_shape))
    ud%umin       = real(src_umin, kind=real64)
    ud%umax       = real(src_umax, kind=real64)
    source%ud     = ud

    source%sd               = int(src_sd, kind=kind(source%sd))
    source%production_fun   = int(src_production_fun, kind=kind(source%production_fun))
    source%production_rate  = real(src_production_rate, kind=real64)
    source%is_jet           = (src_is_jet /= 0_c_int)

    ! Symmetry axis in moon-centred Cartesian coordinates
    call jet_direction(source%betaM, source%zeta, source%eta, &
                       source%rrM, source%symmetry_axis)

    ! Precompute Gu(ui) on [rmin, rmax] defined in const.f90
    call Gu_integral(source%ui, source%Gu_precalc, source%sd, &
                     source%ud, rmin, rmax)

    ! Call core integrator
    call DUDI(density_local, point, source, real(tnow, kind=real64))

    density_out(1) = real(density_local(1), kind=real64)
    density_out(2) = real(density_local(2), kind=real64)

  end subroutine py_dudi_density


  !--------------------------------------------------------------------
  !  Single-point mean-velocity wrapper around DUDI_mean_velocity
  !
  !  C interface:
  !
  !  void py_dudi_mean_velocity(
  !      ... same inputs as py_dudi_density ...,
  !      double tnow,
  !      double integral_out[8]
  !  );
  !--------------------------------------------------------------------
  subroutine py_dudi_mean_velocity( &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, &
       tnow, integral_out) &
       bind(C, name="py_dudi_mean_velocity")

    use iso_fortran_env, only: real64
    implicit none

    ! C boundary types
    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM
    real(c_double), value, intent(in) :: src_zeta, src_eta
    integer(c_int),  value, intent(in) :: src_eject_distr
    integer(c_int),  value, intent(in) :: src_ud_shape
    real(c_double),  value, intent(in) :: src_umin, src_umax
    integer(c_int),  value, intent(in) :: src_sd
    integer(c_int),  value, intent(in) :: src_production_fun
    real(c_double),  value, intent(in) :: src_production_rate
    integer(c_int),  value, intent(in) :: src_is_jet
    real(c_double),  value, intent(in) :: tnow

    real(c_double), intent(out) :: integral_out(8)

    ! Local Fortran-side types
    type(position_in_space)        :: point
    type(source_properties)        :: source
    type(ejection_speed_properties):: ud
    real                           :: integral_local(8)

    ! Build POINT
    point%r     = real(point_r,   kind=real64)
    point%alpha = real(point_alpha, kind=real64)
    point%beta  = real(point_beta,  kind=real64)
    point%rvector(1) = point%r * sin(point%alpha) * cos(point%beta)
    point%rvector(2) = point%r * sin(point%alpha) * sin(point%beta)
    point%rvector(3) = point%r * cos(point%alpha)
    point%r_scaled   = point%r / rm
    point%compute    = .true.

    ! Build SOURCE
    source%r      = real(src_r,      kind=real64)
    source%alphaM = real(src_alphaM, kind=real64)
    source%betaM  = real(src_betaM,  kind=real64)
    source%zeta   = real(src_zeta,   kind=real64)
    source%eta    = real(src_eta,    kind=real64)

    source%rrM(1) = source%r * sin(source%alphaM) * cos(source%betaM)
    source%rrM(2) = source%r * sin(source%alphaM) * sin(source%betaM)
    source%rrM(3) = source%r * cos(source%alphaM)

    source%ejection_angle_distr = int(src_eject_distr, kind=kind(source%ejection_angle_distr))
    ud%ud_shape   = int(src_ud_shape, kind=kind(ud%ud_shape))
    ud%umin       = real(src_umin, kind=real64)
    ud%umax       = real(src_umax, kind=real64)
    source%ud     = ud

    source%sd               = int(src_sd, kind=kind(source%sd))
    source%production_fun   = int(src_production_fun, kind=kind(source%production_fun))
    source%production_rate  = real(src_production_rate, kind=real64)
    source%is_jet           = (src_is_jet /= 0_c_int)

    call jet_direction(source%betaM, source%zeta, source%eta, &
                       source%rrM, source%symmetry_axis)

    call Gu_integral(source%ui, source%Gu_precalc, source%sd, &
                     source%ud, rmin, rmax)

    call DUDI_mean_velocity(integral_local, point, source, real(tnow, kind=real64))

    integral_out = real(integral_local, kind=real64)

  end subroutine py_dudi_mean_velocity


  !--------------------------------------------------------------------
  !  Batched wrapper over points for a single source.
  !
  !  C interface:
  !
  !  void py_dudi_batch_points(
  !      int    n_points,
  !      double density_out[n_points * 2],   // interleaved [b0,u0,b1,u1,...]
  !      const double point_r[n_points],
  !      const double point_alpha[n_points],
  !      const double point_beta[n_points],
  !      double src_r, src_alphaM, src_betaM,
  !      double src_zeta, src_eta,
  !      int    src_eject_distr,
  !      int    src_ud_shape,
  !      double src_umin, src_umax,
  !      int    src_sd,
  !      int    src_production_fun,
  !      double src_production_rate,
  !      int    src_is_jet,
  !      double tnow
  !  );
  !--------------------------------------------------------------------
  subroutine py_dudi_batch_points( &
       n_points, density_out, &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, &
       tnow) &
       bind(C, name="py_dudi_batch_points")

    use iso_fortran_env, only: real64
    implicit none

    integer(c_int),  value, intent(in) :: n_points
    real(c_double), intent(out)        :: density_out(2*n_points)

    real(c_double), intent(in) :: point_r(n_points)
    real(c_double), intent(in) :: point_alpha(n_points)
    real(c_double), intent(in) :: point_beta(n_points)

    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM
    real(c_double), value, intent(in) :: src_zeta, src_eta
    integer(c_int),  value, intent(in) :: src_eject_distr
    integer(c_int),  value, intent(in) :: src_ud_shape
    real(c_double),  value, intent(in) :: src_umin, src_umax
    integer(c_int),  value, intent(in) :: src_sd
    integer(c_int),  value, intent(in) :: src_production_fun
    real(c_double),  value, intent(in) :: src_production_rate
    integer(c_int),  value, intent(in) :: src_is_jet
    real(c_double),  value, intent(in) :: tnow

    type(position_in_space)         :: points(n_points)
    type(source_properties)         :: source
    type(ejection_speed_properties) :: ud
    real                            :: density_2d(2, n_points)
    integer                         :: i

    ! Build POINT array
    do i = 1, n_points
       points(i)%r     = real(point_r(i),     kind=real64)
       points(i)%alpha = real(point_alpha(i), kind=real64)
       points(i)%beta  = real(point_beta(i),  kind=real64)
       points(i)%rvector(1) = points(i)%r * sin(points(i)%alpha) * cos(points(i)%beta)
       points(i)%rvector(2) = points(i)%r * sin(points(i)%alpha) * sin(points(i)%beta)
       points(i)%rvector(3) = points(i)%r * cos(points(i)%alpha)
       points(i)%r_scaled   = points(i)%r / rm
       points(i)%compute    = .true.
    end do

    ! Build SOURCE once
    source%r      = real(src_r,      kind=real64)
    source%alphaM = real(src_alphaM, kind=real64)
    source%betaM  = real(src_betaM,  kind=real64)
    source%zeta   = real(src_zeta,   kind=real64)
    source%eta    = real(src_eta,    kind=real64)

    source%rrM(1) = source%r * sin(source%alphaM) * cos(source%betaM)
    source%rrM(2) = source%r * sin(source%alphaM) * sin(source%betaM)
    source%rrM(3) = source%r * cos(source%alphaM)

    source%ejection_angle_distr = int(src_eject_distr, kind=kind(source%ejection_angle_distr))
    ud%ud_shape   = int(src_ud_shape, kind=kind(ud%ud_shape))
    ud%umin       = real(src_umin, kind=real64)
    ud%umax       = real(src_umax, kind=real64)
    source%ud     = ud

    source%sd               = int(src_sd, kind=kind(source%sd))
    source%production_fun   = int(src_production_fun, kind=kind(source%production_fun))
    source%production_rate  = real(src_production_rate, kind=real64)
    source%is_jet           = (src_is_jet /= 0_c_int)

    call jet_direction(source%betaM, source%zeta, source%eta, &
                       source%rrM, source%symmetry_axis)

    call Gu_integral(source%ui, source%Gu_precalc, source%sd, &
                     source%ud, rmin, rmax)

    call DUDI_batch_points(int(n_points, kind=kind(1)), density_2d, points, source, real(tnow, kind=real64))
    do i = 1, n_points
       density_out(2*(i-1)+1) = real(density_2d(1, i), kind=real64)
       density_out(2*(i-1)+2) = real(density_2d(2, i), kind=real64)
    end do

  end subroutine py_dudi_batch_points


  !--------------------------------------------------------------------
  !  Batched over sources (single point). density_out(2*n_sources).
  !  One [bound, unbound] pair per source.
  !--------------------------------------------------------------------
  subroutine py_dudi_batch_sources( &
       n_sources, density_out, &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, tnow) &
       bind(C, name="py_dudi_batch_sources")
    use iso_fortran_env, only: real64
    implicit none
    integer(c_int),  value, intent(in) :: n_sources
    real(c_double), intent(out)        :: density_out(2*n_sources)
    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: src_r(n_sources), src_alphaM(n_sources), src_betaM(n_sources)
    real(c_double), intent(in) :: src_zeta(n_sources), src_eta(n_sources)
    integer(c_int), intent(in) :: src_eject_distr(n_sources), src_ud_shape(n_sources)
    real(c_double), intent(in) :: src_umin(n_sources), src_umax(n_sources)
    integer(c_int), intent(in) :: src_sd(n_sources), src_production_fun(n_sources)
    real(c_double), intent(in) :: src_production_rate(n_sources)
    integer(c_int), intent(in) :: src_is_jet(n_sources)
    real(c_double), value, intent(in) :: tnow
    type(position_in_space)         :: point
    type(source_properties)         :: sources(n_sources)
    type(ejection_speed_properties) :: ud
    real                            :: density_2d(2, n_sources)
    integer                         :: i
    point%r     = real(point_r,   kind=real64)
    point%alpha = real(point_alpha, kind=real64)
    point%beta  = real(point_beta,  kind=real64)
    point%rvector(1) = point%r * sin(point%alpha) * cos(point%beta)
    point%rvector(2) = point%r * sin(point%alpha) * sin(point%beta)
    point%rvector(3) = point%r * cos(point%alpha)
    point%r_scaled   = point%r / rm
    point%compute    = .true.
    do i = 1, n_sources
       sources(i)%r      = real(src_r(i),      kind=real64)
       sources(i)%alphaM = real(src_alphaM(i), kind=real64)
       sources(i)%betaM  = real(src_betaM(i),  kind=real64)
       sources(i)%zeta   = real(src_zeta(i),   kind=real64)
       sources(i)%eta    = real(src_eta(i),    kind=real64)
       sources(i)%rrM(1) = sources(i)%r * sin(sources(i)%alphaM) * cos(sources(i)%betaM)
       sources(i)%rrM(2) = sources(i)%r * sin(sources(i)%alphaM) * sin(sources(i)%betaM)
       sources(i)%rrM(3) = sources(i)%r * cos(sources(i)%alphaM)
       sources(i)%ejection_angle_distr = int(src_eject_distr(i), kind=kind(sources(i)%ejection_angle_distr))
       ud%ud_shape   = int(src_ud_shape(i), kind=kind(ud%ud_shape))
       ud%umin       = real(src_umin(i), kind=real64)
       ud%umax       = real(src_umax(i), kind=real64)
       sources(i)%ud = ud
       sources(i)%sd              = int(src_sd(i), kind=kind(sources(i)%sd))
       sources(i)%production_fun  = int(src_production_fun(i), kind=kind(sources(i)%production_fun))
       sources(i)%production_rate = real(src_production_rate(i), kind=real64)
       sources(i)%is_jet          = (src_is_jet(i) /= 0_c_int)
       call jet_direction(sources(i)%betaM, sources(i)%zeta, sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)
       call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, sources(i)%sd, sources(i)%ud, rmin, rmax)
    end do

    call DUDI_batch_sources(int(n_sources, kind=kind(1)), density_2d, point, sources, &
         real(tnow, kind=real64))

    do i = 1, n_sources
       density_out(2*(i-1)+1) = real(density_2d(1, i), kind=real64)
       density_out(2*(i-1)+2) = real(density_2d(2, i), kind=real64)
    end do
  end subroutine py_dudi_batch_sources

  !--------------------------------------------------------------------
  !  Batched over points and sources; density summed. density_out(2*n_points).
  !--------------------------------------------------------------------
  subroutine py_dudi_batch_sources_points( &
       n_points, n_sources, density_out, &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, tnow) &
       bind(C, name="py_dudi_batch_sources_points")
    use iso_fortran_env, only: real64
    implicit none
    integer(c_int),  value, intent(in) :: n_points, n_sources
    real(c_double), intent(out)       :: density_out(2*n_points)
    real(c_double), intent(in) :: point_r(n_points), point_alpha(n_points), point_beta(n_points)
    real(c_double), intent(in) :: src_r(n_sources), src_alphaM(n_sources), src_betaM(n_sources)
    real(c_double), intent(in) :: src_zeta(n_sources), src_eta(n_sources)
    integer(c_int), intent(in) :: src_eject_distr(n_sources), src_ud_shape(n_sources)
    real(c_double), intent(in) :: src_umin(n_sources), src_umax(n_sources)
    integer(c_int), intent(in) :: src_sd(n_sources), src_production_fun(n_sources)
    real(c_double), intent(in) :: src_production_rate(n_sources)
    integer(c_int), intent(in) :: src_is_jet(n_sources)
    real(c_double), value, intent(in) :: tnow
    type(position_in_space)         :: points(n_points)
    type(source_properties)         :: sources(n_sources)
    type(ejection_speed_properties) :: ud
    real                            :: density_2d(2, n_points)
    integer                         :: i
    do i = 1, n_points
       points(i)%r     = real(point_r(i),     kind=real64)
       points(i)%alpha = real(point_alpha(i), kind=real64)
       points(i)%beta  = real(point_beta(i),  kind=real64)
       points(i)%rvector(1) = points(i)%r * sin(points(i)%alpha) * cos(points(i)%beta)
       points(i)%rvector(2) = points(i)%r * sin(points(i)%alpha) * sin(points(i)%beta)
       points(i)%rvector(3) = points(i)%r * cos(points(i)%alpha)
       points(i)%r_scaled   = points(i)%r / rm
       points(i)%compute    = .true.
    end do
    do i = 1, n_sources
       sources(i)%r      = real(src_r(i),      kind=real64)
       sources(i)%alphaM = real(src_alphaM(i), kind=real64)
       sources(i)%betaM  = real(src_betaM(i),  kind=real64)
       sources(i)%zeta   = real(src_zeta(i),   kind=real64)
       sources(i)%eta    = real(src_eta(i),    kind=real64)
       sources(i)%rrM(1) = sources(i)%r * sin(sources(i)%alphaM) * cos(sources(i)%betaM)
       sources(i)%rrM(2) = sources(i)%r * sin(sources(i)%alphaM) * sin(sources(i)%betaM)
       sources(i)%rrM(3) = sources(i)%r * cos(sources(i)%alphaM)
       sources(i)%ejection_angle_distr = int(src_eject_distr(i), kind=kind(sources(i)%ejection_angle_distr))
       ud%ud_shape   = int(src_ud_shape(i), kind=kind(ud%ud_shape))
       ud%umin       = real(src_umin(i), kind=real64)
       ud%umax       = real(src_umax(i), kind=real64)
       sources(i)%ud = ud
       sources(i)%sd              = int(src_sd(i), kind=kind(sources(i)%sd))
       sources(i)%production_fun  = int(src_production_fun(i), kind=kind(sources(i)%production_fun))
       sources(i)%production_rate = real(src_production_rate(i), kind=real64)
       sources(i)%is_jet          = (src_is_jet(i) /= 0_c_int)
       call jet_direction(sources(i)%betaM, sources(i)%zeta, sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)
       call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, sources(i)%sd, sources(i)%ud, rmin, rmax)
    end do
    call DUDI_batch_sources_points(int(n_points, kind=kind(1)), int(n_sources, kind=kind(1)), &
                                   density_2d, points, sources, real(tnow, kind=real64))
    do i = 1, n_points
       density_out(2*(i-1)+1) = real(density_2d(1, i), kind=real64)
       density_out(2*(i-1)+2) = real(density_2d(2, i), kind=real64)
    end do
  end subroutine py_dudi_batch_sources_points

  !--------------------------------------------------------------------
  !  Mean-velocity batch over points (single source). integral_out(8*n_points).
  !--------------------------------------------------------------------
  subroutine py_dudi_mean_velocity_batch_points( &
       n_points, integral_out, &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, tnow) &
       bind(C, name="py_dudi_mean_velocity_batch_points")
    use iso_fortran_env, only: real64
    implicit none
    integer(c_int),  value, intent(in) :: n_points
    real(c_double), intent(out)       :: integral_out(8*n_points)
    real(c_double), intent(in) :: point_r(n_points), point_alpha(n_points), point_beta(n_points)
    real(c_double), value, intent(in) :: src_r, src_alphaM, src_betaM, src_zeta, src_eta
    integer(c_int),  value, intent(in) :: src_eject_distr, src_ud_shape
    real(c_double), value, intent(in) :: src_umin, src_umax
    integer(c_int),  value, intent(in) :: src_sd, src_production_fun
    real(c_double), value, intent(in) :: src_production_rate
    integer(c_int),  value, intent(in) :: src_is_jet
    real(c_double), value, intent(in) :: tnow
    type(position_in_space)         :: points(n_points)
    type(source_properties)         :: source
    type(ejection_speed_properties) :: ud
    real                            :: integral_2d(8, n_points)
    integer                         :: i, k
    do i = 1, n_points
       points(i)%r     = real(point_r(i),     kind=real64)
       points(i)%alpha = real(point_alpha(i), kind=real64)
       points(i)%beta  = real(point_beta(i),  kind=real64)
       points(i)%rvector(1) = points(i)%r * sin(points(i)%alpha) * cos(points(i)%beta)
       points(i)%rvector(2) = points(i)%r * sin(points(i)%alpha) * sin(points(i)%beta)
       points(i)%rvector(3) = points(i)%r * cos(points(i)%alpha)
       points(i)%r_scaled   = points(i)%r / rm
       points(i)%compute    = .true.
    end do
    source%r      = real(src_r,      kind=real64)
    source%alphaM = real(src_alphaM, kind=real64)
    source%betaM  = real(src_betaM,  kind=real64)
    source%zeta   = real(src_zeta,   kind=real64)
    source%eta    = real(src_eta,    kind=real64)
    source%rrM(1) = source%r * sin(source%alphaM) * cos(source%betaM)
    source%rrM(2) = source%r * sin(source%alphaM) * sin(source%betaM)
    source%rrM(3) = source%r * cos(source%alphaM)
    source%ejection_angle_distr = int(src_eject_distr, kind=kind(source%ejection_angle_distr))
    ud%ud_shape   = int(src_ud_shape, kind=kind(ud%ud_shape))
    ud%umin       = real(src_umin, kind=real64)
    ud%umax       = real(src_umax, kind=real64)
    source%ud     = ud
    source%sd              = int(src_sd, kind=kind(source%sd))
    source%production_fun  = int(src_production_fun, kind=kind(source%production_fun))
    source%production_rate = real(src_production_rate, kind=real64)
    source%is_jet          = (src_is_jet /= 0_c_int)
    call jet_direction(source%betaM, source%zeta, source%eta, source%rrM, source%symmetry_axis)
    call Gu_integral(source%ui, source%Gu_precalc, source%sd, source%ud, rmin, rmax)
    call DUDI_mean_velocity_batch_points(int(n_points, kind=kind(1)), integral_2d, points, source, real(tnow, kind=real64))
    do i = 1, n_points
       do k = 1, 8
          integral_out(8*(i-1)+k) = real(integral_2d(k, i), kind=real64)
       end do
    end do
  end subroutine py_dudi_mean_velocity_batch_points

  !--------------------------------------------------------------------
  !  Mean-velocity batch over sources (single point). integral_out(8*n_sources).
  !--------------------------------------------------------------------
  subroutine py_dudi_mean_velocity_batch_sources( &
       n_sources, integral_out, &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, tnow) &
       bind(C, name="py_dudi_mean_velocity_batch_sources")
    use iso_fortran_env, only: real64
    implicit none
    integer(c_int),  value, intent(in) :: n_sources
    real(c_double), intent(out)       :: integral_out(8*n_sources)
    real(c_double), value, intent(in) :: point_r, point_alpha, point_beta
    real(c_double), intent(in) :: src_r(n_sources), src_alphaM(n_sources), src_betaM(n_sources)
    real(c_double), intent(in) :: src_zeta(n_sources), src_eta(n_sources)
    integer(c_int), intent(in) :: src_eject_distr(n_sources), src_ud_shape(n_sources)
    real(c_double), intent(in) :: src_umin(n_sources), src_umax(n_sources)
    integer(c_int), intent(in) :: src_sd(n_sources), src_production_fun(n_sources)
    real(c_double), intent(in) :: src_production_rate(n_sources)
    integer(c_int), intent(in) :: src_is_jet(n_sources)
    real(c_double), value, intent(in) :: tnow
    type(position_in_space)         :: point
    type(source_properties)         :: sources(n_sources)
    type(ejection_speed_properties) :: ud
    real                            :: integral_2d(8, n_sources)
    integer                         :: i, k
    point%r     = real(point_r,   kind=real64)
    point%alpha = real(point_alpha, kind=real64)
    point%beta  = real(point_beta,  kind=real64)
    point%rvector(1) = point%r * sin(point%alpha) * cos(point%beta)
    point%rvector(2) = point%r * sin(point%alpha) * sin(point%beta)
    point%rvector(3) = point%r * cos(point%alpha)
    point%r_scaled   = point%r / rm
    point%compute    = .true.
    do i = 1, n_sources
       sources(i)%r      = real(src_r(i),      kind=real64)
       sources(i)%alphaM = real(src_alphaM(i), kind=real64)
       sources(i)%betaM  = real(src_betaM(i),  kind=real64)
       sources(i)%zeta   = real(src_zeta(i),   kind=real64)
       sources(i)%eta    = real(src_eta(i),    kind=real64)
       sources(i)%rrM(1) = sources(i)%r * sin(sources(i)%alphaM) * cos(sources(i)%betaM)
       sources(i)%rrM(2) = sources(i)%r * sin(sources(i)%alphaM) * sin(sources(i)%betaM)
       sources(i)%rrM(3) = sources(i)%r * cos(sources(i)%alphaM)
       sources(i)%ejection_angle_distr = int(src_eject_distr(i), kind=kind(sources(i)%ejection_angle_distr))
       ud%ud_shape   = int(src_ud_shape(i), kind=kind(ud%ud_shape))
       ud%umin       = real(src_umin(i), kind=real64)
       ud%umax       = real(src_umax(i), kind=real64)
       sources(i)%ud = ud
       sources(i)%sd              = int(src_sd(i), kind=kind(sources(i)%sd))
       sources(i)%production_fun  = int(src_production_fun(i), kind=kind(sources(i)%production_fun))
       sources(i)%production_rate = real(src_production_rate(i), kind=real64)
       sources(i)%is_jet          = (src_is_jet(i) /= 0_c_int)
       call jet_direction(sources(i)%betaM, sources(i)%zeta, sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)
       call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, sources(i)%sd, sources(i)%ud, rmin, rmax)
    end do
    call DUDI_mean_velocity_batch_sources(int(n_sources, kind=kind(1)), integral_2d, point, sources, real(tnow, kind=real64))
    do i = 1, n_sources
       do k = 1, 8
          integral_out(8*(i-1)+k) = real(integral_2d(k, i), kind=real64)
       end do
    end do
  end subroutine py_dudi_mean_velocity_batch_sources

  !--------------------------------------------------------------------
  !  Mean-velocity batch over points and sources (summed). integral_out(8*n_points).
  !--------------------------------------------------------------------
  subroutine py_dudi_mean_velocity_batch_sources_points( &
       n_points, n_sources, integral_out, &
       point_r, point_alpha, point_beta, &
       src_r, src_alphaM, src_betaM, src_zeta, src_eta, &
       src_eject_distr, src_ud_shape, src_umin, src_umax, &
       src_sd, src_production_fun, src_production_rate, &
       src_is_jet, tnow) &
       bind(C, name="py_dudi_mean_velocity_batch_sources_points")
    use iso_fortran_env, only: real64
    implicit none
    integer(c_int),  value, intent(in) :: n_points, n_sources
    real(c_double), intent(out)       :: integral_out(8*n_points)
    real(c_double), intent(in) :: point_r(n_points), point_alpha(n_points), point_beta(n_points)
    real(c_double), intent(in) :: src_r(n_sources), src_alphaM(n_sources), src_betaM(n_sources)
    real(c_double), intent(in) :: src_zeta(n_sources), src_eta(n_sources)
    integer(c_int), intent(in) :: src_eject_distr(n_sources), src_ud_shape(n_sources)
    real(c_double), intent(in) :: src_umin(n_sources), src_umax(n_sources)
    integer(c_int), intent(in) :: src_sd(n_sources), src_production_fun(n_sources)
    real(c_double), intent(in) :: src_production_rate(n_sources)
    integer(c_int), intent(in) :: src_is_jet(n_sources)
    real(c_double), value, intent(in) :: tnow
    type(position_in_space)         :: points(n_points)
    type(source_properties)         :: sources(n_sources)
    type(ejection_speed_properties) :: ud
    real                            :: integral_2d(8, n_points)
    integer                         :: i, k
    do i = 1, n_points
       points(i)%r     = real(point_r(i),     kind=real64)
       points(i)%alpha = real(point_alpha(i), kind=real64)
       points(i)%beta  = real(point_beta(i),  kind=real64)
       points(i)%rvector(1) = points(i)%r * sin(points(i)%alpha) * cos(points(i)%beta)
       points(i)%rvector(2) = points(i)%r * sin(points(i)%alpha) * sin(points(i)%beta)
       points(i)%rvector(3) = points(i)%r * cos(points(i)%alpha)
       points(i)%r_scaled   = points(i)%r / rm
       points(i)%compute    = .true.
    end do
    do i = 1, n_sources
       sources(i)%r      = real(src_r(i),      kind=real64)
       sources(i)%alphaM = real(src_alphaM(i), kind=real64)
       sources(i)%betaM  = real(src_betaM(i),  kind=real64)
       sources(i)%zeta   = real(src_zeta(i),   kind=real64)
       sources(i)%eta    = real(src_eta(i),    kind=real64)
       sources(i)%rrM(1) = sources(i)%r * sin(sources(i)%alphaM) * cos(sources(i)%betaM)
       sources(i)%rrM(2) = sources(i)%r * sin(sources(i)%alphaM) * sin(sources(i)%betaM)
       sources(i)%rrM(3) = sources(i)%r * cos(sources(i)%alphaM)
       sources(i)%ejection_angle_distr = int(src_eject_distr(i), kind=kind(sources(i)%ejection_angle_distr))
       ud%ud_shape   = int(src_ud_shape(i), kind=kind(ud%ud_shape))
       ud%umin       = real(src_umin(i), kind=real64)
       ud%umax       = real(src_umax(i), kind=real64)
       sources(i)%ud = ud
       sources(i)%sd              = int(src_sd(i), kind=kind(sources(i)%sd))
       sources(i)%production_fun  = int(src_production_fun(i), kind=kind(sources(i)%production_fun))
       sources(i)%production_rate = real(src_production_rate(i), kind=real64)
       sources(i)%is_jet          = (src_is_jet(i) /= 0_c_int)
       call jet_direction(sources(i)%betaM, sources(i)%zeta, sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)
       call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, sources(i)%sd, sources(i)%ud, rmin, rmax)
    end do
    call DUDI_mean_velocity_batch_sources_points(int(n_points, kind=kind(1)), &
         int(n_sources, kind=kind(1)), integral_2d, points, sources, real(tnow, kind=real64))
    do i = 1, n_points
       do k = 1, 8
          integral_out(8*(i-1)+k) = real(integral_2d(k, i), kind=real64)
       end do
    end do
  end subroutine py_dudi_mean_velocity_batch_sources_points


end module py_dudi_bridge
