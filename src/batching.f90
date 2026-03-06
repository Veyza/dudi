! This file is a part of DUDI, the Fortran-95 implementation
! of the two-body model for the dynamics of dust ejected from an
! atmosphereless body.
! This is free software. You can use and redistribute it
! under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
! If you do, please cite the following paper
! Anastasiia Ershova and Jürgen Schmidt,
! Two-body model for the spatial distribution of dust ejected from
! an atmosphereless body, 2021, A&A, 650, A186
!
! Module: batching.f90
! Description:
!   Batched wrappers over DUDI's main subroutines DUDI and DUDI_mean_velocity:
!     * batch over points (single source),
!     * batch over sources (single point),
!     * batch over points and sources (sum contributions over sources).
!   OpenMP is used to parallelize over the loop index.

module batching
    use define_types
    use integrator, only: DUDI, DUDI_mean_velocity
    implicit none

contains

    !==================================================================
    !> Batched computation over an array of points for a single source.
    !!  For each i = 1..n_points, computes density(1:2, i) at points(i)
    !!  due to the same source via DUDI.
    !!  density(1,:) = bound, density(2,:) = unbound.
    !==================================================================
    subroutine DUDI_batch_points(n_points, density, points, source, tnow)
        implicit none
        integer, intent(in) :: n_points
        real, intent(out) :: density(2, n_points)
        type(position_in_space), intent(in) :: points(n_points)
        type(source_properties), intent(in) :: source
        real(8), intent(in) :: tnow

        integer :: i
        real :: dens(2)

        density(:, :) = 0.0

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dens) SCHEDULE(static)
        do i = 1, n_points
            call DUDI(dens, points(i), source, tnow)
            density(1, i) = dens(1)
            density(2, i) = dens(2)
        end do
        !$OMP END PARALLEL DO

    end subroutine DUDI_batch_points


    !==================================================================
    !> Batched computation over an array of sources for a single point.
    !!  For each i = 1..n_sources, computes density(1:2, i) at point
    !!  due to sources(i) via DUDI.
    !!  density(1,:) = bound, density(2,:) = unbound.
    !==================================================================
    subroutine DUDI_batch_sources(n_sources, density, point, sources, tnow)
        implicit none
        integer, intent(in) :: n_sources
        real, intent(out) :: density(2, n_sources)
        type(position_in_space), intent(in) :: point
        type(source_properties), intent(in) :: sources(n_sources)
        real(8), intent(in) :: tnow

        integer :: i
        real :: dens(2)

        density(:, :) = 0.0

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dens) SCHEDULE(static)
        do i = 1, n_sources
            call DUDI(dens, point, sources(i), tnow)
            density(1, i) = dens(1)
            density(2, i) = dens(2)
        end do
        !$OMP END PARALLEL DO

    end subroutine DUDI_batch_sources


    !==================================================================
    !> Batched computation over points and sources: for each point,
    !! sum density from all sources via DUDI.
    !!  density(1, i) = sum over sources of bound at point i,
    !!  density(2, i) = sum over sources of unbound at point i.
    !==================================================================
    subroutine DUDI_batch_sources_points(n_points, n_sources, density, points, sources, tnow)
        implicit none
        integer, intent(in) :: n_points, n_sources
        real, intent(out) :: density(2, n_points)
        type(position_in_space), intent(in) :: points(n_points)
        type(source_properties), intent(in) :: sources(n_sources)
        real(8), intent(in) :: tnow

    integer :: i, i_s
    real :: dens(2)
    real :: tmp(2, n_points)

    density(:, :) = 0.0

    do i_s = 1, n_sources
        tmp(:, :) = 0.0
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dens) SCHEDULE(static)
        do i = 1, n_points
            call DUDI(dens, points(i), sources(i_s), tnow)
            tmp(1, i) = dens(1)
            tmp(2, i) = dens(2)
        end do
        !$OMP END PARALLEL DO
        do i = 1, n_points
            density(1, i) = density(1, i) + tmp(1, i)
            density(2, i) = density(2, i) + tmp(2, i)
        end do
    end do

    end subroutine DUDI_batch_sources_points


    !==================================================================
    !> Batched mean-velocity over points for a single source.
    !!  integral(1:8, i) = DUDI_mean_velocity at points(i), source, tnow.
    !==================================================================
    subroutine DUDI_mean_velocity_batch_points(n_points, integral, points, source, tnow)
        implicit none
        integer, intent(in) :: n_points
        real, intent(out) :: integral(8, n_points)
        type(position_in_space), intent(in) :: points(n_points)
        type(source_properties), intent(in) :: source
        real(8), intent(in) :: tnow

        integer :: i
        real :: int8(8)

        integral(:, :) = 0.0

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, int8) SCHEDULE(static)
        do i = 1, n_points
            call DUDI_mean_velocity(int8, points(i), source, tnow)
            integral(:, i) = int8(:)
        end do
        !$OMP END PARALLEL DO

    end subroutine DUDI_mean_velocity_batch_points


    !==================================================================
    !> Batched mean-velocity over sources for a single point.
    !!  integral(1:8, i) = DUDI_mean_velocity at point, sources(i), tnow.
    !==================================================================
    subroutine DUDI_mean_velocity_batch_sources(n_sources, integral, point, sources, tnow)
        implicit none
        integer, intent(in) :: n_sources
        real, intent(out) :: integral(8, n_sources)
        type(position_in_space), intent(in) :: point
        type(source_properties), intent(in) :: sources(n_sources)
        real(8), intent(in) :: tnow

        integer :: i
        real :: int8(8)

        integral(:, :) = 0.0

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, int8) SCHEDULE(static)
        do i = 1, n_sources
            call DUDI_mean_velocity(int8, point, sources(i), tnow)
            integral(:, i) = int8(:)
        end do
        !$OMP END PARALLEL DO

    end subroutine DUDI_mean_velocity_batch_sources


    !==================================================================
    !> Batched mean-velocity over points and sources: for each point,
    !! sum integral from all sources via DUDI_mean_velocity.
    !!  integral(1:8, i) = sum over sources at point i.
    !==================================================================
    subroutine DUDI_mean_velocity_batch_sources_points(n_points, n_sources, integral, points, sources, tnow)
        implicit none
        integer, intent(in) :: n_points, n_sources
        real, intent(out) :: integral(8, n_points)
        type(position_in_space), intent(in) :: points(n_points)
        type(source_properties), intent(in) :: sources(n_sources)
        real(8), intent(in) :: tnow

        integer :: i, i_s
        real :: int8(8)
        real :: tmp8(8, n_points)

        integral(:, :) = 0.0

    do i_s = 1, n_sources
        tmp8(:, :) = 0.0
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, int8) SCHEDULE(static)
        do i = 1, n_points
            call DUDI_mean_velocity(int8, points(i), sources(i_s), tnow)
            tmp8(:, i) = int8(:)
        end do
        !$OMP END PARALLEL DO
        do i = 1, n_points
            integral(:, i) = integral(:, i) + tmp8(:, i)
        end do
    end do

    end subroutine DUDI_mean_velocity_batch_sources_points

end module batching
