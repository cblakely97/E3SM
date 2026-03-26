! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_gofs_coupling
!
!> \brief MPAS-Ocean GOFS3.1 one-way coupling
!> \author Caleb Blakely
!> \date   March 2026
!> \details
!>  This module implements one-way coupling from GOFS3.1 (HYCOM reanalysis)
!>  to MPAS-Ocean. Two physical mechanisms are applied:
!>    1. External pressure gradient forcing: depth-integrated baroclinic PG
!>       from GOFS plus a surface-density-weighted SSH gradient, applied
!>       uniformly at all depth levels.
!>    2. Internal tide drag: GOFS-derived drag coefficient applied to tidal
!>       barotropic velocity isolated via a 25-hour sub-sampled boxcar filter.
!>
!>  See: https://doi.org/10.1029/2024MS004502/
!
!-----------------------------------------------------------------------

module ocn_gofs_coupling

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_forcing
   use mpas_stream_manager
   use mpas_timekeeping
   use ocn_framework_forcing
   use ocn_constants
   use ocn_mesh
   use ocn_diagnostics_variables
   use ocn_config
   use mpas_log

   implicit none
   private
   save

   public :: ocn_gofs_coupling_init, &
             ocn_gofs_coupling_forcing_init, &
             ocn_gofs_coupling_get_forcing, &
             ocn_gofs_coupling_tend_pgrad, &
             ocn_gofs_coupling_tend_drag, &
             ocn_gofs_coupling_compute_lowpass, &
             ocn_gofs_coupling_update_filter

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------

   logical :: gofsCouplingOn
   real (kind=RKIND) :: gofsDragCoeff
   real (kind=RKIND) :: gdensity0Inv_gofs
   integer :: gofsNFilterInterval, gofsNFilterSamples

   ! Pointers to time-varying input fields (from timeVaryingForcing pool)
   real (kind=RKIND), dimension(:), pointer :: gofsBaroclinicPGrad => null()
   real (kind=RKIND), dimension(:), pointer :: gofsSurfaceDensity => null()
   real (kind=RKIND), dimension(:), pointer :: gofsInternalTideDragCoeff => null()

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_gofs_coupling_init
!
!> \brief   Initialize GOFS coupling module flags and constants
!> \author  Caleb Blakely
!> \date    March 2026
!> \details
!>  Sets module-level flags and constants from namelist config.
!>  Called from ocn_vel_forcing_init (no domain access needed).
!
!-----------------------------------------------------------------------

   subroutine ocn_gofs_coupling_init(err)!{{{

      integer, intent(out) :: err

      err = 0
      gofsCouplingOn = .false.

      if (config_use_gofs_coupling) then
         gofsCouplingOn = .true.
         gofsDragCoeff = config_gofs_tidal_drag_coeff
         gofsNFilterInterval = config_gofs_filter_interval
         gofsNFilterSamples = config_gofs_filter_nsamples
         gdensity0Inv_gofs = gravity / rho_sw

         call mpas_log_write('')
         call mpas_log_write(' GOFS3.1 one-way coupling enabled')
         call mpas_log_write('   tidal drag coeff multiplier: $r', &
              realArgs=(/gofsDragCoeff/))
         call mpas_log_write('   filter sub-sample interval:  $i timesteps', &
              intArgs=(/gofsNFilterInterval/))
         call mpas_log_write('   filter buffer size:          $i samples', &
              intArgs=(/gofsNFilterSamples/))
         call mpas_log_write('')
      end if

   end subroutine ocn_gofs_coupling_init!}}}

!***********************************************************************
!
!  routine ocn_gofs_coupling_forcing_init
!
!> \brief   Initialize GOFS time-varying forcing via mpas_forcing framework
!> \author  Caleb Blakely
!> \date    March 2026
!> \details
!>  Sets up the mpas_forcing group and fields for reading time-varying
!>  GOFS input data with linear interpolation. Also retrieves pointers
!>  to the interpolated fields in the timeVaryingForcing pool.
!>  Called from mpas_ocn_forward_mode after ocn_time_varying_forcing_init.
!
!-----------------------------------------------------------------------

   subroutine ocn_gofs_coupling_forcing_init(domain)!{{{

      type (domain_type), intent(inout) :: domain

      type (mpas_pool_type), pointer :: tvfPool

      if (.not. gofsCouplingOn) return

      ! Create forcing group for GOFS coupling
      call mpas_forcing_init_group( &
           forcingGroupHead, &
           "ocn_gofs_coupling", &
           domain, &
           config_gofs_forcing_start_time, &
           config_gofs_forcing_cycle_start, &
           config_gofs_forcing_cycle_duration, &
           config_do_restart)

      ! Initialize each input field with linear interpolation
      call mpas_forcing_init_field( &
           domain % streamManager, &
           forcingGroupHead, &
           "ocn_gofs_coupling", &
           "gofsBaroclinicPGrad", &
           "gofs_coupling", &
           "timeVaryingForcing", &
           "gofsBaroclinicPGrad", &
           "linear", &
           config_gofs_forcing_reference_time, &
           config_gofs_forcing_interval)

      call mpas_forcing_init_field( &
           domain % streamManager, &
           forcingGroupHead, &
           "ocn_gofs_coupling", &
           "gofsSurfaceDensity", &
           "gofs_coupling", &
           "timeVaryingForcing", &
           "gofsSurfaceDensity", &
           "linear", &
           config_gofs_forcing_reference_time, &
           config_gofs_forcing_interval)

      call mpas_forcing_init_field( &
           domain % streamManager, &
           forcingGroupHead, &
           "ocn_gofs_coupling", &
           "gofsInternalTideDragCoeff", &
           "gofs_coupling", &
           "timeVaryingForcing", &
           "gofsInternalTideDragCoeff", &
           "linear", &
           config_gofs_forcing_reference_time, &
           config_gofs_forcing_interval)

      ! Initialize field data (reads initial time records)
      call mpas_forcing_init_field_data( &
           forcingGroupHead, &
           "ocn_gofs_coupling", &
           domain % streamManager, &
           config_do_restart, &
           .false.)

      ! Retrieve pointers to interpolated input fields
      call mpas_pool_get_subpool(domain % blocklist % structs, &
           'timeVaryingForcing', tvfPool)
      call mpas_pool_get_array(tvfPool, 'gofsBaroclinicPGrad', &
           gofsBaroclinicPGrad)
      call mpas_pool_get_array(tvfPool, 'gofsSurfaceDensity', &
           gofsSurfaceDensity)
      call mpas_pool_get_array(tvfPool, 'gofsInternalTideDragCoeff', &
           gofsInternalTideDragCoeff)

      ! Initialize filter state
      gofsFilterIndex = 1
      gofsFilterCount = 0
      gofsFilterStepCounter = 0

   end subroutine ocn_gofs_coupling_forcing_init!}}}

!***********************************************************************
!
!  routine ocn_gofs_coupling_get_forcing
!
!> \brief   Update time-varying GOFS forcing fields via interpolation
!> \author  Caleb Blakely
!> \date    March 2026
!> \details
!>  Calls the mpas_forcing framework to update the interpolated GOFS
!>  input fields for the current time. Called inside the RK4 loop
!>  alongside atmospheric forcing, following the same clock management
!>  pattern (advance by forcingTimeIncrement, then reverse).
!
!-----------------------------------------------------------------------

   subroutine ocn_gofs_coupling_get_forcing(streamManager, domain)!{{{

      type (MPAS_streamManager_type), intent(inout) :: streamManager
      type (domain_type), intent(inout) :: domain

      real (kind=RKIND), pointer :: forcingTimeIncrement
      type (mpas_pool_type), pointer :: forcingPool
      real (kind=RKIND) :: dtSim

      if (.not. gofsCouplingOn) return

      ! Get the forcing time increment for this RK substep
      call mpas_pool_get_subpool(domain % blocklist % structs, &
           'forcing', forcingPool)
      call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', &
           forcingTimeIncrement)
      dtSim = forcingTimeIncrement

      ! Get interpolated forcing at current time + dtSim
      call mpas_forcing_get_forcing( &
           forcingGroupHead, &
           "ocn_gofs_coupling", &
           streamManager, &
           dtSim)

      ! Reverse the clock advance (same pattern as atmospheric forcing)
      call mpas_advance_forcing_clock(forcingGroupHead, -dtSim)

   end subroutine ocn_gofs_coupling_get_forcing!}}}

!***********************************************************************
!
!  routine ocn_gofs_coupling_tend_pgrad
!
!> \brief   Compute external pressure gradient tendency from GOFS coupling
!> \author  Caleb Blakely
!> \date    March 2026
!> \details
!>  Computes the external pressure gradient at each edge:
!>    externalPG = gofsBaroclinicPGrad(iEdge)
!>               + (g/rho_0) * [rho_s(cell2)*ssh(cell2)
!>                            - rho_s(cell1)*ssh(cell1)] / dcEdge(iEdge)
!>  Applied uniformly at all depth levels.
!>  Sign follows existing PG convention (tend -= PG).
!
!-----------------------------------------------------------------------

   subroutine ocn_gofs_coupling_tend_pgrad(ssh, tend, err)!{{{

      real (kind=RKIND), dimension(:), intent(in) :: ssh
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend
      integer, intent(out) :: err

      integer :: iEdge, k, cell1, cell2, kMin, kMax
      real (kind=RKIND) :: externalPG, invdcEdge

      err = 0
      if (.not. gofsCouplingOn) return

      call mpas_timer_start('vel gofs pgrad')

      !$omp do schedule(runtime) &
      !$omp    private(cell1, cell2, invdcEdge, kMin, kMax, externalPG, k)
      do iEdge = 1, nEdgesOwned
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)
         kMin = minLevelEdgeBot(iEdge)
         kMax = maxLevelEdgeTop(iEdge)

         if (kMax >= kMin) then
            ! External PG: baroclinic PG + surface-density-weighted SSH gradient
            externalPG = gofsBaroclinicPGrad(iEdge) &
                       + gdensity0Inv_gofs * invdcEdge &
                       * (gofsSurfaceDensity(cell2) * ssh(cell2) &
                        - gofsSurfaceDensity(cell1) * ssh(cell1))

            ! Store diagnostic
            gofsExternalPGrad(iEdge) = externalPG

            ! Apply uniformly at all depth levels
            do k = kMin, kMax
               tend(k, iEdge) = tend(k, iEdge) &
                               - edgeMask(k, iEdge) * externalPG
            end do
         end if
      end do
      !$omp end do

      call mpas_timer_stop('vel gofs pgrad')

   end subroutine ocn_gofs_coupling_tend_pgrad!}}}

!***********************************************************************
!
!  routine ocn_gofs_coupling_tend_drag
!
!> \brief   Compute internal tide drag tendency from GOFS coupling
!> \author  Caleb Blakely
!> \date    March 2026
!> \details
!>  Applies GOFS-derived drag at the bottom level only:
!>    tend(kBot,iEdge) -= dragCoeff * gofsInternalTideDragCoeff(iEdge)
!>                      * gofsTidalVelocity(iEdge)
!>  gofsTidalVelocity is pre-computed before the RK loop and held
!>  constant across substeps.
!
!-----------------------------------------------------------------------

   subroutine ocn_gofs_coupling_tend_drag(normalVelocity, tend, err)!{{{

      real (kind=RKIND), dimension(:,:), intent(in) :: normalVelocity
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend
      integer, intent(out) :: err

      integer :: iEdge, k
      real (kind=RKIND) :: dragTend

      err = 0
      if (.not. gofsCouplingOn) return

      call mpas_timer_start('vel gofs tidal drag')

      !$omp do schedule(runtime) private(k, dragTend)
      do iEdge = 1, nEdgesOwned
         k = maxLevelEdgeTop(iEdge)
         if (k > 0) then
            dragTend = gofsDragCoeff &
                     * gofsInternalTideDragCoeff(iEdge) &
                     * gofsTidalVelocity(iEdge)
            tend(k, iEdge) = tend(k, iEdge) - dragTend
            gofsTidalDragTend(iEdge) = dragTend
         end if
      end do
      !$omp end do

      call mpas_timer_stop('vel gofs tidal drag')

   end subroutine ocn_gofs_coupling_tend_drag!}}}

!***********************************************************************
!
!  routine ocn_gofs_coupling_compute_lowpass
!
!> \brief   Compute barotropic, lowpass, and tidal velocities
!> \author  Caleb Blakely
!> \date    March 2026
!> \details
!>  Called once per timestep, before the RK loop.
!>  1. Computes depth-averaged (barotropic) velocity at each edge.
!>  2. Computes lowpass-filtered barotropic velocity from the circular
!>     buffer (mean of valid entries).
!>  3. Computes tidal velocity = barotropic - lowpass.
!>  During spinup (buffer not full), lowpass is the mean of available
!>  samples. If buffer is empty, tidal velocity is zero.
!
!-----------------------------------------------------------------------

   subroutine ocn_gofs_coupling_compute_lowpass(normalVelocity, &!{{{
                     layerThickEdgeMean, err)

      real (kind=RKIND), dimension(:,:), intent(in) :: normalVelocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickEdgeMean
      integer, intent(out) :: err

      integer :: iEdge, k, kMin, kMax, iSample
      real (kind=RKIND) :: sumVelH, sumH, lowpass

      err = 0
      if (.not. gofsCouplingOn) return

      call mpas_timer_start('gofs compute lowpass')

      !$omp do schedule(runtime) &
      !$omp    private(kMin, kMax, sumVelH, sumH, k, lowpass, iSample)
      do iEdge = 1, nEdgesOwned
         kMin = minLevelEdgeBot(iEdge)
         kMax = maxLevelEdgeTop(iEdge)

         ! Step 1: Compute barotropic velocity
         sumVelH = 0.0_RKIND
         sumH = 0.0_RKIND
         do k = kMin, kMax
            sumVelH = sumVelH + normalVelocity(k, iEdge) &
                              * layerThickEdgeMean(k, iEdge)
            sumH = sumH + layerThickEdgeMean(k, iEdge)
         end do

         if (sumH > 0.0_RKIND) then
            gofsBarotropicVelocity(iEdge) = sumVelH / sumH
         else
            gofsBarotropicVelocity(iEdge) = 0.0_RKIND
         end if

         ! Step 2: Compute lowpass from circular buffer
         if (gofsFilterCount > 0) then
            lowpass = 0.0_RKIND
            do iSample = 1, gofsFilterCount
               lowpass = lowpass + gofsFilterBuffer(iSample, iEdge)
            end do
            lowpass = lowpass / real(gofsFilterCount, RKIND)
            gofsLowpassVelocity(iEdge) = lowpass
         else
            gofsLowpassVelocity(iEdge) = gofsBarotropicVelocity(iEdge)
         end if

         ! Step 3: Tidal velocity = barotropic - lowpass
         gofsTidalVelocity(iEdge) = gofsBarotropicVelocity(iEdge) &
                                  - gofsLowpassVelocity(iEdge)
      end do
      !$omp end do

      call mpas_timer_stop('gofs compute lowpass')

   end subroutine ocn_gofs_coupling_compute_lowpass!}}}

!***********************************************************************
!
!  routine ocn_gofs_coupling_update_filter
!
!> \brief   Update the 25-hour filter circular buffer
!> \author  Caleb Blakely
!> \date    March 2026
!> \details
!>  Called once per timestep, after the RK loop.
!>  Increments the step counter. When it reaches gofsNFilterInterval:
!>    - Recomputes barotropic velocity from the new state (time level 2)
!>    - Stores it in the circular buffer
!>    - Advances the circular index
!>    - Increments the valid sample count (capped at buffer size)
!>    - Resets the step counter
!
!-----------------------------------------------------------------------

   subroutine ocn_gofs_coupling_update_filter(normalVelocity, &!{{{
                     layerThickEdgeMean, err)

      real (kind=RKIND), dimension(:,:), intent(in) :: normalVelocity
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickEdgeMean
      integer, intent(out) :: err

      integer :: iEdge, k, kMin, kMax
      real (kind=RKIND) :: sumVelH, sumH

      err = 0
      if (.not. gofsCouplingOn) return

      ! Increment step counter
      gofsFilterStepCounter = gofsFilterStepCounter + 1

      ! Check if it's time to store a sample
      if (gofsFilterStepCounter >= gofsNFilterInterval) then

         call mpas_timer_start('gofs update filter')

         ! Compute barotropic velocity from the new state and store in buffer
         !$omp do schedule(runtime) private(kMin, kMax, sumVelH, sumH, k)
         do iEdge = 1, nEdgesOwned
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)

            sumVelH = 0.0_RKIND
            sumH = 0.0_RKIND
            do k = kMin, kMax
               sumVelH = sumVelH + normalVelocity(k, iEdge) &
                                 * layerThickEdgeMean(k, iEdge)
               sumH = sumH + layerThickEdgeMean(k, iEdge)
            end do

            if (sumH > 0.0_RKIND) then
               gofsFilterBuffer(gofsFilterIndex, iEdge) = sumVelH / sumH
            else
               gofsFilterBuffer(gofsFilterIndex, iEdge) = 0.0_RKIND
            end if
         end do
         !$omp end do

         ! Advance circular index
         gofsFilterIndex = mod(gofsFilterIndex, gofsNFilterSamples) + 1

         ! Increment valid sample count (capped at buffer size)
         if (gofsFilterCount < gofsNFilterSamples) then
            gofsFilterCount = gofsFilterCount + 1
         end if

         ! Reset step counter
         gofsFilterStepCounter = 0

         call mpas_timer_stop('gofs update filter')

      end if

   end subroutine ocn_gofs_coupling_update_filter!}}}

!***********************************************************************

end module ocn_gofs_coupling

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
