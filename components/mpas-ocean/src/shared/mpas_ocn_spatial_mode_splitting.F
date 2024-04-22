
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_spatial_mode_splitting
!
!> \brief MPAS ocean spatial mode splitting 
!> \author Coleman Blakely
!> \date   April 2024
!> \details
!>  This module contains the subroutines used in the spatial mode 
!>  splitting wherein the barotropic and baroclinic modes are run on 
!>  separate meshes. 
!
!-----------------------------------------------------------------------

module ocn_spatial_mode_splitting

    use mpas_kind_types
    use mpas_io_units
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_constants

    use ocn_constants
    use ocn_config
    use ocn_init_vertical_grids
    use ocn_init_cell_markers
    use mpas_geometry_utils
#ifdef HAVE_MOAB
    use mpas_log     
    use mpas_field_routines
    use mpas_sort
    use mpas_stream_manager
    use mpas_moabmesh
    use mpas_dmpar
    use mpas_vector_reconstruction
    use iMOAB
   
    use ocn_diagnostics_variables
#include "moab/MOABConfig.h"

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

    public :: ocn_spatial_mode_splitting_init_moab_domain,&
              ocn_spatial_mode_splitting_bc2btr_remapping,&
              ocn_spatial_mode_splitting_btr2bc_remapping

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------
   
    ! identifiers for moab meshes
    ! id's on compute core 
    integer :: cmpbc    = 314
    integer :: cmpbtr   = 314
    integer :: bc2btrid = 316
    integer :: btr2bcid = 317
    ! id's on coupler *NOTE: we don't presently have separate coupler cores
    ! but iMOAB requires us to register a coupler anyways.
    integer :: cplbc, cplbtr
    ! PID's for MOAB
    integer,public :: cmpbcPID, cmpbtrPID, bc2btrPID, btr2bcPID
    integer,public :: cplbcPID, cplbtrPID
    ! iMOAB communication variables
    character*100 :: appname
    character(:), allocatable :: bc2btrFields, bc2btrFields_btr,&
                                 btr2bcFields, btr2bcFields_btr
    integer, dimension(2) :: bc2btrIndices, btr2bcIndices,&
                             bc2btrTypes, btr2bcTypes
    real (kind=RKIND), dimension(:), allocatable :: bc2btrVar, btr2bcVar, &
                                                    bc2btrVar_btr, btr2bcVar_btr
    ! iMOAB weights variables
    character(:), allocatable :: bc2btr_weights_identifier,&
                                 btr2bc_weights_identifier,&
                                 disc_method,&
                                 dof_tag_names,&
                                 fv_method
    integer :: disc_order, fNoBubble, fMonotoneTypeID, fVolumetric, fNoConserve,&
               fValidate, fInverseDistanceMap
    ! The following variables are the cartesian coordinate components of 
    ! the fields that are used for remapping. They are outputs of 
    ! mpas_reconstruct but are not needed so we allocate them as private
    ! module variables but they are never used.
    real (kind=RKIND), dimension(:), allocatable :: barotropicForcingX,&
                                                    barotropicForcingY,&
                                                    barotropicForcingZ,&
                                                    barotropicThicknessFluxX_btr,&
                                                    barotropicThicknessFluxY_btr,&
                                                    barotropicThicknessFluxZ_btr

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_init_moab_domain
!
!> \brief   Initializes the moab mesh for bc/btr domains
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!>  Based on the domainID, initializes the baroclinic/barotropic domains
!>  within moab for the spatial mode splitting. Only does so if we have
!>  configured the run to do so
!
!-----------------------------------------------------------------------
    subroutine ocn_spatial_mode_splitting_init_moab_domain(domain)!{{{
    type (domain_type), intent(inout) :: domain
    type (mpas_pool_type), pointer :: meshPool
    integer, pointer :: nCellsSolve
    logical, pointer :: config_spatial_mode_splitting
    integer :: ierr=0
    integer, dimension(3) :: nverts, nelem, nblocks, nsbc, ndbc
    ! number of vertices on meshes
    integer, pointer :: nCellsSolve, nCellsSolve_btr
    ! edge tangent vectors for calculating cell-centered quantities
    real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors
    real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors_btr

    ! get configuration to see if we are using the mode splitting, if not return
    call mpas_pool_get_config(domain % configs, 'config_spatial_mode_splitting',&
       config_spatial_mode_splitting)
    if (.not. config_spatial_mode_splitting) then
       return
    endif
    ! if this is baroclinic domain, call mpas_moabmesh and then return. If barotropic
    ! compute overlap/weights, etc.
    if ( domain % domainID .eq. 0 ) then
       ! initialize moab for baroclinic domain
       call init_moab_mpas(domain, cmpbc, cmpbcPID)
       ! allocate barotropicForcingX, barotropicForcingY, barotropicForcingZ for baroclinic
       ! domain. They are not used in calculations but are outputs from mpas_reconstruct
       ! so we need to have them hanging around.
       call mpas_pool_get_subpool(domain % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       allocate ( barotropicForcingX(nCellsSolve),&
                  barotropicForcingY(nCellsSolve),&
                  barotropicForcingZ(nCellsSolve) )
       ! initialize reconstruction routines to go from edge based quantities to cell 
       ! centered quantities
       call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)
       call mpas_rbf_interp_initialize(meshPool)
       call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors)
       call mpas_init_reconstruct(meshPool)
       return
    elseif ( domain % domainID .eq. 1 ) then
       call init_moab_mpas(domain, cmpbtr, cmpbtrPID)
       ! allocate barotropicThicknessFluxX, etc. for barotropic
       ! domain. They are not used in calculations but are outputs from mpas_reconstruct
       ! so we need to have them hanging around.
       call mpas_pool_get_subpool(domain % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       allocate ( barotropicThicknessFluxX_btr(nCellsSolve),&
                  barotropicThicknessFluxY_btr(nCellsSolve),&
                  barotropicThicknessFluxZ_btr(nCellsSolve) )
       ! initialize reconstruction routines to go from edge based quantities to cell 
       ! centered quantities
       call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors_btr)
       call mpas_rbf_interp_initialize(meshPool)
       call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors_btr)
       call mpas_init_reconstruct(meshPool)
    endif
    ! once both the baroclinic and barotropic meshes are registered in moab, we need to calculate
    ! the overlap meshes and compute the weights for interpolation as well as define our tags
    !
    ! barotropic2baroclinic
    !
    ! register baroclinic to barotropic application
    appname = 'bc2btr_cpl'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, &
               domain % dminfo % comm, &
               bc2btrid, &
               bc2btrPID)
    call errorout(ierr, "failed to register baroclinic to barotropic application")
    ! compute baroclinic to barotropic intersection
    ierr = iMOAB_ComputeMeshIntersectionOnSphere(cmpbcPID,&
                                          cmpbtrPID, &
                                          bc2btrPID)
    call errorout(ierr,"failed to compute intersection of baroclinic to barotropic mesh")
    ! get projection weights for baroclinic to barotropic. Note that both domains are finite
    ! volume discretizations with the same discretization order. Due to this we just use the
    ! same variables for both source and target mesh
    bc2btr_weights_identifier = 'bc2btr'//C_NULL_CHAR
    disc_method = 'fv'//C_NULL_CHAR
    fv_method = ''//C_NULL_CHAR ! not sure what this is for, figure it out later
    disc_order = 1
    dof_tag_names = 'GLOBAL_ID'//C_NULL_CHAR
    fNoBubble = 1
    fMonotoneTypeID = 0
    fVolumetric = 0
    fNoConserve = 0
    fValidate = 0
    fInverseDistanceMap = 0

    ierr = iMOAB_ComputeScalarProjectionWeights( &
                     bc2btrPID, &
                     bc2btr_weights_identifier, &
                     disc_method, &
                     disc_order, &
                     disc_method, &
                     disc_order, &
                     fv_method, &
                     fNoBubble, &
                     fMonotoneTypeID, &
                     fVolumetric, &
                     fInverseDistanceMap, &
                     fNoConserve, &
                     fValidate, &
                     dof_tag_names, &
                     dof_tag_names)
    call errorout(ierr, "falied to compute weights for baroclinic to barotropic remapping")
    !
    ! baroclinic2barotropic
    !
    ! register barotropic to baroclinic application
    appname = 'btr2bc_cpl'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, &
               domain % dminfo % comm, &
               btr2bcid, &
               btr2bcPID)
    call errorout(ierr, "failed to register barotropic to baroclinic application")
    ! compute baroclinic to barotropic intersection
    ierr = iMOAB_ComputeMeshIntersectionOnSphere(cmpbtrPID,&
                                          cmpbcPID, &
                                          btr2bcPID)
    call errorout(ierr,"failed to compute intersection of barotropic to baroclinic mesh")
    ! get projection weights for baroclinic to barotropic. Note that both domains are finite
    ! volume discretizations with the same discretization order. Due to this we just use the
    ! same variables for both source and target mesh
    btr2bc_weights_identifier = 'btr2bc'//C_NULL_CHAR

    ierr = iMOAB_ComputeScalarProjectionWeights( &
                     btr2bcPID, &
                     btr2bc_weights_identifier, &
                     disc_method, &
                     disc_order, &
                     disc_method, &
                     disc_order, &
                     fv_method, &
                     fNoBubble, &
                     fMonotoneTypeID, &
                     fVolumetric, &
                     fInverseDistanceMap, &
                     fNoConserve, &
                     fValidate, &
                     dof_tag_names, &
                     dof_tag_names)
    call errorout(ierr, "failed to compute weights for barotropic to baroclinic remapping")
    !
    ! define tags on meshes
    !
    ! define baroclinic to barotropic variables on baroclinic mesh
    bc2btrTypes(:) = 1 ! DENSE_DOUBLE
    bc2btrFields = 'barotropicForcingZonal:barotropicForcingMeridional'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(cmpbcPID, &
                     bc2btrFields, &
                     bc2btrTypes(1), &
                     1, &
                     bc2btrIndices(1))
    call errorout(ierr, "failed to define bc2btr tags on bc mesh")
    ! define barotropic to baroclinic variables on baroclinic mesh
    btr2bcTypes(:) = 1 ! DENSE_DOUBLE
    btr2bcFields = 'barotropicThicknessFluxZonal:barotropicThicknessFluxMeridional'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(cmpbcPID, &
                     btr2bcFields, &
                     btr2bcTypes(1), &
                     1, &
                     btr2bcIndices(1))
    call errorout(ierr, "failed to define btr2bc tags on bc mesh")
    ! initialize buffers for transmitting the data (both baroclinic to barotropic and
    ! baroclinic to barotropic while we have the mesh size)
    ierr = iMOAB_GetMeshInfo(cmpbcPID, nverts, nele, nblocks, nsbc, ndbc)
    call errorout(ierr, "failed to retrieve mesh info from baroclinic mesh")
    allocate ( bc2btrVar(2*nelem(3)), btr2bcVar(2*nelem(3)) )
    bc2btrVar(:) = 0d0
    ! define baroclinic to barotropic variables on barotropic mesh
    bc2btrFields_btr = 'barotropicForcingZonal_btr:barotropicForcingMeridional_btr'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(cmpbtrPID, &
                  bc2btrFields_btr, &
                  bc2btrTypes(2),&
                  1, &
                  bc2btrIndices(2))
    call errorout(ierr, "failed to define bc2btr tags on btr mesh")
    ! define barotropic to baroclinic variables on barotropic mesh
    btr2bcFields_btr = &
      'barotropicThicknessFluxZonal_btr:barotropicThicknessFluxMeridional_btr'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(cmpbtrPID, &
                  btr2bcFields_btr, &
                  btr2bcTypes(2), &
                  1, &
                  btr2bcIndices(2))
    call errorout(ierr, "failed to define btr2bc tags on btr mesh")
    ! initialize buffers for transmitting data from barotropic to baroclinic mesh
    ierr = iMOAB_GetMeshInfo(cmpbtrPID, nverts, nele, nblocks, nsbc, ndbc)
    call errorout(ierr, "failed to retrieve mesh info from barotropic mesh")
    allocate ( bc2btrVar_btr(2*nelem(3)), btr2bcVar_btr(2*nelem(3)) )

    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_init_moab_domain!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_bc2btr_remapping
!
!> \brief   Remaps barotropicForcing from baroclinic to barotropic domain
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!>  This routine remaps barotropicForcing from the baroclinic domain to
!>  the barotropic domain when using spatial mode splitting
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_bc2btr_remapping(domain)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type(domain_type), intent(inout) :: domain

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------


    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------
    type (domain_type), pointer :: domain_btr
    ! 
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    ! mesh information
    integer, pointer :: nCellsSolve, nCellsSolve_btr,&
                        nEdgesSolve, nEdgesSolve_btr
    ! arrays for storing data
    real (kind=RKIND), dimension(:), pointer :: barotropicForcing,&
                                                barotropicForcing_btr,&
                                                barotropicForcingZonal,&
                                                barotropicForcingZonal_btr,&
                                                barotropicForcingMeridional,&
                                                barotropicForcingMeridional_btr
    ! for going from cell centered to edge based quantities
    integer :: cell1, cell2
    real (kind=RKIND) :: zonal_edge, meridional_edge
    real (kind=RKIND), dimension(:), pointer :: angleEdge
    integer, dimension(:,:), pointer :: cellsOnEdge
    ! looping variables
    integer :: iCell, iEdge
    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    ! "next" domain is barotropic domain
    domain_btr => domain % next
    ! get info from baroclinic domain
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
        ! get barotropicForcing from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingZonal',&
                  barotropicForcingZonal)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingMeridional',&
                  barotropicForcingMeridional)
        ! reconstruct from edge based to cell-centered quantities
        call mpas_reconstruct(meshPool, barotropicForcing, &
                  barotropicForcingX, barotropicForcingY, barotropicForcingZ, &
                  barotropicForcingZonal, barotropicForcingMeridional)
        ! pack buffers for iMOAB
        do iCell = 1,nCellsSolve
            bc2btrVar(iCell) = barotropicForcingZonal(iCell)
            bc2btrVar(iCell+nCellsSolve) = barotropicForcingMeridional(iCell)
        enddo
        ! put buffers into moab structure
        ierr=iMOAB_SetDoubleTagStorage(cmpbcPID,&
                           bc2btrFields,&
                           2*nCellsSolve,&
                           1,&
                           bc2btrVar)
        call errorout(ierr,"failed to assign bc2btr data to baroclinic mesh")
        ! do remapping
        ierr = iMOAB_ApplyScalarProjectionWeights(bc2btrPID,&
                           bc2btr_weights_identifier,&
                           bc2btrFields,&
                           bc2btrFields_btr)
        call errorout(ierr,"failed to reproject bc2btr data to barotropic mesh")
        block_ptr => block_ptr % next
    enddo
    
    ! send info to barotropic comain
    block_ptr => domain_btr % blocklist
    do while(associated(block_ptr))
        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_btr)
        ! get barotropicForcing from pools
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingZonal',&
                  barotropicForcingZonal_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingMeridional',&
                  barotropicForcingMeridional_btr)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
        ! get remapped data from iMOAB
        ierr = iMOAB_GetDoubleTagStorage(cmpbtrPID,&
                  bc2btrFields_btr,&
                  2*nCellsSolve_btr,&
                  1,&
                  bc2btrVar_btr)
        call errorout(ierr,"failed to get bc2btrVar_btr")
        ! unpack buffers
        do iCell = 1,nCellsSolve_btr
            barotropicForcingZonal_btr(iCell) = bc2btrVar_btr(iCell)
            barotropicForcingMeridional_btr(iCell) = bc2btrVar_btr(iCell+nCellsSolve_btr)
        enddo
        ! put back into edge based arrays
        barotropicForcing_btr(:) = 0d0
        do iEdge = 1,nEdgesSolve_btr
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (cell1.gt.nCellsSolve_btr.OR.cell1.lt.1.OR. &
                cell2.gt.nCellsSolve_btr.OR.cell2.lt.1) then
                cycle
             endif
             zonal_edge = 0.5d0 * (barotropicForcingZonal_btr(cell1) + &
                                      barotropicForcingZonal_btr(cell2))
             meridional_edge = 0.5d0 * (barotropicForcingMeridional_btr(cell1) + &
                                      barotropicForcingMeridional_btr(cell2))
             barotropicForcing_btr(iedge) = barotropicForcing_btr(iEdge) + &
                              cos(angleEdge(iEdge))*zonal_bt_edge + &
                              sin(angleEdge(iEdge))*merid_bt_edge
        enddo
        block_ptr => block_ptr % next
    enddo
    end subroutine ocn_spatial_mode_splitting_bc2btr_remapping!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_btr2bc_remapping
!
!> \brief   Remaps barotropicThickness Flux from barotropic to baroclinic domain
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!>  This routine remaps barotropicThicknessFlux from the barotropic
!>  to the baroclinic domain when using spatial mode spligging
!
!-----------------------------------------------------------------------


    subroutine ocn_spatial_mode_splitting_btr2bc_remapping(domain)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type(domain_type), intent(inout) :: domain

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------


    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------
    type (domain_type), pointer :: domain_btr
    ! 
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    ! mesh information
    integer, pointer :: nCellsSolve, nCellsSolve_btr,&
                        nEdgesSolve, nEdgesSolve_btr
    ! arrays for storing data
    real (kind=RKIND), dimension(:), pointer :: barotropicThicknessFlux,&
                                                barotropicThicknessFlux_btr,&
                                                barotropicThicknessFluxZonal,&
                                                barotropicThicknessFluxZonal_btr,&
                                                barotropicThicknessFluxMeridional,&
                                                barotropicThicknessFluxMeridional_btr
    ! for going from cell centered to edge based quantities
    integer :: cell1, cell2
    real (kind=RKIND) :: zonal_edge, meridional_edge
    real (kind=RKIND), dimension(:), pointer :: angleEdge
    integer, dimension(:,:), pointer :: cellsOnEdge
    ! looping variables
    integer :: iCell, iEdge
    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    ! "next" domain is barotropic domain
    domain_btr => domain % next
    ! get info from barotropic domain
    block_ptr => domain_btr % blocklist
    do while(associated(block_ptr))
        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_btr)
        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                  barotropicThicknessFlux_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxZonal',&
                  barotropicThicknessFluxZonal_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxMeridional',&
                  barotropicThicknessFluxMeridional_btr)
        ! reconstruct from edge based to cell-centered quantities
        call mpas_reconstruct(meshPool, barotropicThicknessFlux_btr,&
                  barotropicThicknessFluxX_btr,&
                  barotropicThicknessFluxY_btr,&
                  barotropicThicknessFluxZ_btr,&
                  barotropicThicknessFluxZonal_btr,&
                  barotropicThicknessFluxMeridional_btr)
        ! pack buffers
        do iCell = 1,nCellsSolve_btr
            btr2bcVar_btr(iCell) = barotropicThicknessFluxZonal(iCell)
            btr2bcVar_btr(iCell+nCellsSolve_btr) = barotropicThicknessFluxMeridional(iCell)
        enddo
        ! put buffers into moab structure
        ierr = iMOAB_SetDoubleTagStorage(cmpbtrPID,&
                  btr2bcFields_btr,&
                  2*nCellsSolve_btr,&
                  1,&
                  btr2bcVar_btr)
        call errorout(ierr,"failed to assign btr2bc data to barotropic mesh")
        ! do remapping
        ierr = iMOAB_ApplyScalarProjectionWeights(btr2bcPID,&
                  btr2bc_weights_identifier,&
                  btr2bcFields_btr,&
                  btr2bcFields)
        block_ptr => block_ptr % next
    enddo
    ! put data on baroclinic domain
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                  barotropicThicknessFlux)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxZonal',&
                  barotropicThicknessFluxZonal)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxMeridional',&
                  barotropicThicknessFluxMeridional)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
        ! get remapped data from iMOAB
        ierr = iMOAB_GetDoubleTagStorage(cmpbcPID,&
                  btr2bcFields,&
                  2*nCellsSolve,&
                  1,&
                  btr2bcVar)
        call errorout(ierr,"failed to get btr2bcVar from iMOAB")
        ! unpack buffers
        do iCell = 1,nCellsSolve
            barotropicThicknessFluxZonal(iCell) = btr2bcVar(iCell)
            barotropicThicknessFluxMeridional(iCell) = btr2bcVar(iCell+nCellsSolve)
        enddo
        ! put into edge-based arrays
        barotropicThicknessFlux(:)=0d0
        do iEdge = 1,nEdgesSolve
           cell1 = cellsOnEdge(1, iEdge)
           cell2 = cellsOnEdge(2,iEdge)
            if (cell1.gt.nCellsSolve.OR.cell1.lt.1.OR. &
                cell2.gt.nCellsSolve.OR.cell2.lt.1) then
                cycle
             endif
             zonal_edge = 0.5d0 * (barotropicThicknessFluxZonal(cell1) + &
                                      barotropicThicknessFluxZonal(cell2))
             meridional_edge = 0.5d0 * (barotropicThicknessFluxMeridional(cell1) + &
                                      barotropicThicknessFluxMeridional(cell2))
             barotropicThicknessFlux(iEdge) = barotropicThicknessFlux(iEdge) + &
                              cos(angleEdge(iEdge))*zonal_bt_edge + &
                              sin(angleEdge(iEdge))*merid_bt_edge
        enddo
        block_ptr => block_ptr % next
    enddo
    
    end subroutine ocn_spatial_mode_splitting_btr2bc_remapping!}}}
#endif
end module ocn_spatial_mode_splitting
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
