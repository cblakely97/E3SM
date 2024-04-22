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
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    
    ! local variables
    integer :: iCell, iEdge, k, idx
    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations
    ! cpb
    ! barotropic forcing temp arrays
    ! bc domain
    real (kind=RKIND), dimension(:), allocatable :: btForcingMeridional, btForcingZonal
    real (kind=RKIND), dimension(:), allocatable :: btForcingX, btForcingY, btForcingZ
    real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors
    ! bt domain
    real (kind=RKIND), dimension(:), allocatable :: btForcingMeridional_btr, btForcingZonal_btr
    real (kind=RKIND), dimension(:), allocatable :: btForcingX_btr, btForcingY_btr, btForcingZ_btr
    real (kind=RKIND) :: zonal_bt_edge, merid_bt_edge ! for cell => edge
    integer :: cell1, cell2
    ! end cpb

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
    
    ! cpb
    real (kind=RKIND), dimension(:), pointer :: barotropicForcing
    real (kind=RKIND), dimension(:), pointer :: barotropicForcing_btr
    real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge
    ! iMOAB variables
    character*100 :: outfn, writeopts, appname
    character(:), allocatable :: weights_identifier1
    character(:), allocatable :: disc_methods1, disc_methods2, dof_tag_names1, dof_tag_names2
    integer :: disc_orders1, disc_orders2
    integer :: fNoBubble, fMonotoneTypeID, fVolumetric, fNoConserve, fValidate, fInverseDistanceMap
    integer, dimension(2) :: tagIndices, tagTypes
    character(:), allocatable :: bcFields, btrFields
    integer, dimension(3) :: nverts, nelem, nblocks, nsbc, ndbc
    real (kind=RKIND), dimension(:), allocatable :: varbc, varbtr ! for sending back and forth
    ! end cpb
    iErr = 0
    tagTypes(1) = 1 ! DENSE_DOUBLE
    tagTypes(2) = 1 ! DENSE_DOUBLE


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

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
       call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

       block_ptr => block_ptr % next
    end do

    call ocn_spatial_mode_splitting_init_moab_domain(domain)
    call ocn_spatial_mode_splitting_init_moab_domain(domain % next)

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
