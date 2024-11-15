! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_remap
!
!> \brief MPAS ocean initialize case -- remap
!> \author FILL_IN_AUTHOR
!> \date   FILL_IN_DATE
!> \details
!>
!
!-----------------------------------------------------------------------

module ocn_init_remap
   
   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_config
   use ocn_init_vertical_grids
   use ocn_init_cell_markers
   ! CPB
   use ocn_spatial_mode_splitting
   ! end CPB
   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_setup_remap, &
             ocn_init_validate_remap

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_remap
!
!> \brief   Setup for this initial condition
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details
!>  This routine sets up the initial conditions for this case.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_remap(domain, iErr)!{{{
  implicit none

    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    type (domain_type), pointer :: domain_btr
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    
    ! local variables
    integer :: iCell, iEdge, k, idx
    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity
    integer :: nEdges

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, dcEdge
    integer, dimension(:,:), pointer :: cellsOnEdge
    real (kind=RKIND), dimension(:), pointer :: angleEdge
    
    real (kind=RKIND), dimension(:), pointer :: barotropicForcingZonal_bc,&
                                                barotropicForcingMeridional_bc
    real (kind=RKIND), dimension(:), pointer :: barotropicThicknessFluxZonal_btr,&
                                                barotropicThicknessFluxMeridional_btr
    integer, pointer :: nCellsSolve_bc, nCellsAll_bc, nCellsSolve_btr, nCellsAll_btr
    real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge
    iErr = 0

    if(config_init_configuration .ne. trim('remap')) return

    ! cpb
    domain_btr => domain % next
    ! end cpb

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
    ! CPB question: Since in this case we are using a depth-averaged quantity we don't really care
    ! about the vertical grid unless I am crazy. I think I can just set nVertLevelsP1 = 1
    nVertLevelsP1=2
    ! end CPB

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

    !--------------------------------------------------------------------
    ! Use this section to make boundaries non-periodic
    !--------------------------------------------------------------------

    ! Initalize min/max values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    xMin = 1.0E10_RKIND
    xMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max values.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

       yMin = min( yMin, minval(yCell(1:nCellsSolve)))
       yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
       xMin = min( xMin, minval(xCell(1:nCellsSolve)))
       xMax = max( xMax, maxval(xCell(1:nCellsSolve)))
       dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

       block_ptr => block_ptr % next
    end do

    ! Determine global min and max values.
    call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

    ! zero out halos to test halo exchanges
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

       call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
       call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

       ! get input barotropicForcing on baroclinic mesh and zero out halos
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_bc)
       call mpas_pool_get_dimension(meshPool, 'nCellsOwned', nCellsAll_bc)
       call mpas_pool_get_array(diagnosticsPool,&
                                'barotropicForcingZonal',&
                                barotropicForcingZonal_bc)
       call mpas_pool_get_array(diagnosticsPool,&
                                'barotropicForcingMeridional',&
                                barotropicForcingMeridional_bc)

!       do iCell = nCellsSolve_bc,nCellsAll_bc
!          barotropicForcingZonal_bc(iCell) = 0d0
!          barotropicForcingMeridional_bc(iCell) = 0d0
!       enddo

       block_ptr => block_ptr % next
    end do

    ! zero out halos to test halo exchanges
    block_ptr => domain_btr % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

       call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
       call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

       ! get input barotropicForcing on baroclinic mesh and zero out halos
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_btr)
       call mpas_pool_get_dimension(meshPool, 'nCellsAll', nCellsAll_btr)
       call mpas_pool_get_array(diagnosticsPool,&
                                'barotropicThicknessFluxZonal',&
                                barotropicThicknessFluxZonal_btr)
       call mpas_pool_get_array(diagnosticsPool,&
                                'barotropicThicknessFluxMeridional',&
                                barotropicThicknessFluxMeridional_btr)

!       do iCell = nCellsSolve_btr,nCellsAll_btr
!          barotropicThicknessFluxZonal_btr(iCell) = 0d0
!          barotropicThicknessFluxMeridional_btr(iCell) = 0d0
!       enddo

       block_ptr => block_ptr % next
    end do


    !--------------------------------------------------------------------
    ! initialize moab for spatial mode splitting
    !--------------------------------------------------------------------

    call ocn_spatial_mode_splitting_init_moab_domain(domain)
    call ocn_spatial_mode_splitting_init_moab_domain(domain % next)
    ! remap bathymetry from fine to coarse mesh so that we have consistent bathy. We will just use
    ! the first half of bc2btrVar_btr to send to the first half of bc2btrVar rather than create new
    ! tags that will only be used once
    !call ocn_spatial_mode_splitting_bottomDepth_remapping(domain)

    !--------------------------------------------------------------------
    ! baroclinic => barotropic
    !--------------------------------------------------------------------

    call ocn_spatial_mode_splitting_bc2btr_remapping(domain)

    !--------------------------------------------------------------------
    ! barotropic => baroclinic
    !--------------------------------------------------------------------

    call ocn_spatial_mode_splitting_btr2bc_remapping(domain)
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_remap!}}}


!***********************************************************************
!
!  routine ocn_init_validate_remap
!
!> \brief   Validation for this initial condition
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_remap(configPool, packagePool, iocontext, iErr)!{{{
   implicit none

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_remap_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('remap')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_remap_vert_levels', config_remap_vert_levels)

      if(config_vert_levels <= 0 .and. config_remap_vert_levels > 0) then
         config_vert_levels = config_remap_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for remap. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_remap!}}}


!***********************************************************************
end module ocn_init_remap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
