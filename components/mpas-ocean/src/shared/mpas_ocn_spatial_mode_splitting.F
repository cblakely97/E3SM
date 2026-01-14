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

    use mpas_geometry_utils
    use mpas_log     
    use mpas_field_routines
    use mpas_sort
    use mpas_stream_manager
    use mpas_moabmesh
    use mpas_dmpar
    use mpas_vector_reconstruction
    use mpas_vector_operations

    use ocn_diagnostics_variables
    use ocn_diagnostics
    use ocn_mesh
    use ocn_constants
    use ocn_config
    use ocn_vel_self_attraction_loading, only : check, interpolate
    use ocn_vel_tidal_potential

    use netcdf

#ifdef HAVE_MOAB
    use iMOAB
#include "moab/MOABConfig.h"
#endif

    use mpi

    implicit none
    private
    save

    !--------------------------------------------------------------------
    !
    ! Public parameters
    !
    !--------------------------------------------------------------------

    ! Multiplicative masks and vectors for various conditions
!    real(kind=RKIND), public, target, allocatable :: edgeMask_spatial_split_btr(:,:),&
!                                                     edgeMask_spatial_split_bc(:,:)
    integer, public, target, allocatable :: maxLevelEdgeTop_spatial_split(:)
    real(kind=RKIND), public, pointer :: edgeMask_spatial_split_btr(:,:),&
                                                     edgeMask_spatial_split_bc(:,:)

    !--------------------------------------------------------------------
    !
    ! Public member functions
    !
    !--------------------------------------------------------------------

    public :: ocn_spatial_mode_splitting_init_moab_domain,&
              ocn_spatial_mode_splitting_setup_domain,&
              ocn_spatial_mode_splitting_bc2btr_remapping,&
              ocn_spatial_mode_splitting_btr2bc_remapping,&
              ocn_spatial_mode_splitting_time_int_split_setup_btr,&
              ocn_spatial_mode_splitting_dmpar_init_bc,&
              ocn_spatial_mode_splitting_dmpar_init_btr,&
              ocn_spatial_mode_splitting_set_edgeMask_btr

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------

    ! mpi comm pointer and group
    integer :: moabComm
    integer :: moabGroup
   
    ! identifiers for moab meshes
    ! barotropic mesh identifiers
    integer :: btrGID
    integer :: btrXbcGID

    ! baroclinic mesh identifiers
    integer :: bcGID
    integer :: bcXbtrGID

    ! PID's for MOAB
    integer,public :: btrPID, btrXbcPID
    integer,public :: bcPID,  bcXbtrPID

    ! iMOAB communication variables
    character*100 :: appname
    character(:), allocatable :: bc2btrFields, bc2btrFields_btr,&
                                 btr2bcFields, btr2bcFields_btr
    integer, dimension(2) :: bc2btrIndices, btr2bcIndices,&
                             bc2btrTypes, btr2bcTypes
    real (kind=RKIND), dimension(:), allocatable :: bc2btrVar, btr2bcVar, &
                                                    bc2btrVar_btr, btr2bcVar_btr
    integer :: partScheme = 0 ! trivial partitioner

    ! iMOAB weights variables
    character(:), allocatable :: bc2btr_weights_identifier,&
                                 btr2bc_weights_identifier,&
                                 disc_method,&
                                 dof_tag_names,&
                                 fv_method

    integer :: disc_order, fNoBubble, fMonotoneTypeID, fVolumetric, fNoConserve,&
               fValidate, fInverseDistanceMap, filter_type

    ! These edge-based quantities are used in remapping but are local variables only.
   real (kind=RKIND), dimension(:), allocatable :: barotropicForcingOnEdgeZonal_bc,&
                                                   barotropicForcingOnEdgeMeridional_bc,&
                                                   barotropicForcingOnEdgeZonal_btr,&
                                                   barotropicForcingonEdgeMeridional_btr,&
                                                   barotropicThicknessFluxOnEdgeZonal_bc,&
                                                   barotropicThicknessFluxOnEdgeMeridional_bc,&
                                                   barotropicThicknessFluxOnEdgeZonal_btr,&
                                                   barotropicThicknessFluxOnEdgeMeridional_btr


   ! The following variables are the cartesian coordinate components of 
   ! the fields that are used for remapping. They are outputs of 
   ! mpas_reconstruct but are not needed so we allocate them as private
   ! module variables but they are never used.
   real (kind=RKIND), dimension(:), allocatable :: barotropicForcingX_bc,&
                                                   barotropicForcingY_bc,&
                                                   barotropicForcingZ_bc,&
                                                   barotropicThicknessFluxX_btr,&
                                                   barotropicThicknessFluxY_btr,&
                                                   barotropicThicknessFluxZ_btr

   ! variables for setting edgeMask for the boundary between the baroclinic and barotropic domains
   integer, pointer :: nVertLevels_bc, nEdgesAll_bc

   ! The following variables are for performing the edge-based remapping with the mpas_dmpar
   ! routines. This is a temporary measure until we get moab fully functioning but we want a
   ! functioning baseline.

   ! weights stored as sparse array
   integer, dimension(:), allocatable :: rowValues_bc2btr, colValues_bc2btr,&
                                         rowValues_btr2bc, colValues_btr2bc
   real(kind=RKIND), dimension(:), allocatable :: sValues_bc2btr, sValues_btr2bc

   ! for MPI communication
   integer :: curProc, nProcs
   character(:), allocatable :: weights_filename_btr2bc
   integer :: nEdgesGlobal_bc, nEdgesGlobal_btr
   integer, dimension(:), allocatable :: nEdgesDisplacement_bc,&
                                         nEdgesDisplacement_btr,&
                                         indexToEdgeIDGathered_bc,&
                                         indexToEdgeIDGathered_btr,&
                                         nEdgesPerProc_bc,&
                                         nEdgesPerProc_btr

   ! for debugging (temporary)
   logical :: carryOverVals=.FALSE.

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

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (domain_type), intent(inout) :: domain
    type (mpas_pool_type), pointer :: meshPool
    logical, pointer :: config_spatial_mode_splitting
    integer :: ierr=0
    integer, dimension(3) :: nverts, nelem, nblocks, nsbc, ndbc

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------
    ! number of vertices on meshes
    integer, pointer :: nCellsSolve_bc, nCellsSolve_btr

    ! number of edges in mesh
    integer, pointer :: nEdgesSolve_bc, nEdgesSolve_btr

    ! edge tangent vectors for calculating cell-centered quantities
    real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors_bc
    real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors_btr

    ! iMOAB variables
    integer :: nghlay = 0
    integer :: context_id
    integer :: nTags

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code
    ! identifiers for moab meshes
    ! barotropic mesh identifiers
    btrGID    = 314
    btrXbcGID  = 316

    ! baroclinic mesh identifiers
    bcGID      = 414
    bcXbtrGID = 416

    ! number of tags we define in this subroutine
    nTags=2

    ! get configuration to see if we are using the mode splitting, if not return
    call mpas_pool_get_config(domain % configs, 'config_spatial_mode_splitting',&
                              config_spatial_mode_splitting)
    if (.not. config_spatial_mode_splitting) then
       return
    endif
#ifndef HAVE_MOAB
    call mpas_log_write('MPAS-Ocean must be compiled with HAVE_MOAB=true for ',//
                         'spatial mode splitting.', MPAS_LOG_CRIT)
#else

    ! if this is baroclinic domain, call mpas_moabmesh and then return. If barotropic
    ! compute overlap/weights, etc.
    if ( domain % domainID .eq. 0 ) then

       !--------------------------------------------------------------------
       ! baroclinic domain
       !--------------------------------------------------------------------

       ! get mpi communicator for baroclinic domain
       moabComm = MPI_COMM_NULL
       call mpi_comm_dup(domain % dminfo % comm, &
                         moabComm, ierr)
       call mpi_comm_group(moabComm, moabGroup, ierr)

       ! initialize moab domain
       call init_moab_mpas(domain, bcGID, bcPID, moabComm)

       ! get mesh size for allocating local variables
       call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_bc)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_bc)

       ! allocate edge-based meridional and zonal component arrays
       allocate ( barotropicForcingOnEdgeZonal_bc(nEdgesSolve_bc),&
                  barotropicForcingOnEdgeMeridional_bc(nEdgesSolve_bc),&
                  barotropicThicknessFluxOnEdgeZonal_bc(nEdgesSolve_bc),&
                  barotropicThicknessFluxOnEdgeMeridional_bc(nEdgesSolve_bc) )

       ! allocate barotropicForcingX_bc, barotropicForcingY_bc, barotropicForcingZ_bc for baroclinic
       ! domain. They are not used in calculations but are outputs from mpas_reconstruct
       ! so we need to have them hanging around.
       allocate ( barotropicForcingX_bc(nCellsSolve_bc),&
                  barotropicForcingY_bc(nCellsSolve_bc),&
                  barotropicForcingZ_bc(nCellsSolve_bc) )

       ! initialize reconstruction routines to go from edge based quantities to cell 
       ! centered quantities
       call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors_bc)
       call mpas_rbf_interp_initialize(meshPool)
       call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors_bc)
       call mpas_init_reconstruct(meshPool,.true.)

       return
    elseif ( domain % domainID .eq. 1 ) then

       !--------------------------------------------------------------------
       ! barotropic domain
       !--------------------------------------------------------------------

       ! initialize moab domain
       call init_moab_mpas(domain, btrGID, btrPID, moabComm)

       ! get mesh size for allocating local variables
       call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_btr)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_btr)

       ! allocate edge-based meridional and zonal component arrays
       allocate ( barotropicForcingOnEdgeZonal_btr(nEdgesSolve_btr),&
                  barotropicForcingOnEdgeMeridional_btr(nEdgesSolve_btr),&
                  barotropicThicknessFluxOnEdgeZonal_btr(nEdgesSolve_btr),&
                  barotropicThicknessFluxOnEdgeMeridional_btr(nEdgesSolve_btr) )

       ! allocate barotropicThicknessFluxX, etc. for barotropic
       ! domain. They are not used in calculations but are outputs from mpas_reconstruct
       ! so we need to have them hanging around.
       allocate ( barotropicThicknessFluxX_btr(nCellsSolve_btr),&
                  barotropicThicknessFluxY_btr(nCellsSolve_btr),&
                  barotropicThicknessFluxZ_btr(nCellsSolve_btr) )

       ! initialize reconstruction routines to go from edge based quantities to cell 
       ! centered quantities
       call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors_btr)
       call mpas_rbf_interp_initialize(meshPool)
       call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors_btr)
       call mpas_init_reconstruct(meshPool)

    endif

    !--------------------------------------------------------------------
    ! Finish MOAB initialization
    !--------------------------------------------------------------------

    ! register applications other than btrGID and bcGID

    ! register baroclinic to barotropic intersection
    appname = 'bcXbtr'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, &
                                     moabComm, &
                                     bcXbtrGID, &
                                     bcXbtrPID)
    call errorout(ierr, "failed to register bcXbtr application")

    ! barotropic to baroclinic intersection on baroclinic side
    appname = 'btrXbc'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname,&
                                     moabComm, &
                                     btrXbcGID, &
                                     btrXbcPID)
    call errorout(ierr, "failed to register btrXbc application")

    !--------------------------------------------------------------------
    ! calculate the mesh intersections
    !--------------------------------------------------------------------

    ! baroclinic => barotropic
    call ocn_spatial_mode_splitting_calc_mesh_intersection(bcPID,&
                                                           btrPID,&
                                                           bcXbtrPID,&
                                                           bcGID,&
                                                           btrGID,&
                                                           bcXbtrGID)

    ! barotropic => baroclinic
    call ocn_spatial_mode_splitting_calc_mesh_intersection(btrPID,&
                                                           bcPID,&
                                                           btrXbcPID,&
                                                           btrGID,&
                                                           bcGID,&
                                                           btrXbcGID)

    !--------------------------------------------------------------------
    ! calculate the projection weights
    !--------------------------------------------------------------------

    ! baroclinic => barotropic
    bc2btr_weights_identifier='bc2btr'//C_NULL_CHAR
    fv_method='bilin'//C_NULL_CHAR
    call ocn_spatial_mode_splitting_calc_weights(bcXbtrPID,&
                                                 bc2btr_weights_identifier,&
                                                 fv_method)

    ! barotropic => baroclinic
    btr2bc_weights_identifier='btr2bc'//C_NULL_CHAR
    fv_method=''//C_NULL_CHAR ! conservative
    call ocn_spatial_mode_splitting_calc_weights(btrXbcPID,&
                                                 btr2bc_weights_identifier,&
                                                 fv_method)

    !--------------------------------------------------------------------
    ! define tags
    !--------------------------------------------------------------------

    ! all field types are dense double
    bc2btrTypes(:)    = 1 ! DENSE_DOUBLE
    btr2bcTypes(:)    = 1 ! DENSE_DOUBLE

    ! baroclinic to barotropic fields
    bc2btrFields     = &
      'barotropicForcingZonal_bc:barotropicForcingMeridional_bc'//C_NULL_CHAR
    bc2btrFields_btr = &
      'barotropicForcingZonal_btr:barotropicForcingMeridional_btr'//C_NULL_CHAR

    ! barotropic to baroclinic fields
    btr2bcFields     = &
      'barotropicThicknessFluxZonal_bc:barotropicThicknessFluxMeridional_bc'//C_NULL_CHAR
    btr2bcFields_btr = &
      'barotropicThicknessFluxZonal_btr:barotropicThicknessFluxMeridional_btr'//C_NULL_CHAR

    ! baroclinic => barotropic
    call ocn_spatial_mode_splitting_define_tags(bcPID,&
                                                btrPID,&
                                                bc2btrFields,&
                                                bc2btrFields_btr,&
                                                bc2btrTypes,&
                                                bc2btrIndices)

    ! barotropic => baroclinic
    call ocn_spatial_mode_splitting_define_tags(btrPID,&
                                                bcPID,&
                                                btr2bcFields_btr,&
                                                btr2bcFields,&
                                                btr2bcTypes,&
                                                btr2bcIndices)


    !--------------------------------------------------------------------
    ! allocate buffers
    !--------------------------------------------------------------------

    ! start from here

    ! baroclinic to barotropic
    call ocn_spatial_mode_splitting_allocate_buffers(bcPID,&
                                                     bc2btrVar,&
                                                     btr2bcVar,&
                                                     nTags)

    ! barotropic to baroclinic
    call ocn_spatial_mode_splitting_allocate_buffers(btrPID,&
                                                     btr2bcVar_btr,&
                                                     bc2btrVar_btr,&
                                                     nTags)

    ! destroy mesh structure (NOTE: THIS IS DONE IN update_vertical_grid EVENTUALLY)
    call ocn_meshDestroy(iErr)

#endif
    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_init_moab_domain!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_calc_mesh_intersection
!
!> \brief   Calculates the weights for spatial mode splitting
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!>  Sends source mesh to target "side", calculates mesh intersection, 
!>  calculates coverage graph.
!
!-----------------------------------------------------------------------
    subroutine ocn_spatial_mode_splitting_calc_mesh_intersection(sourcePID,&
                                                       targetPID,&
                                                       overlapPID,&
                                                       sourceGID,&
                                                       targetGID,&
                                                       overlapGID)

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! moab PID's
    integer,intent(inout) :: sourcePID,&
                             targetPID,&
                             overlapPID

    ! global id's
    integer,intent(inout) :: sourceGID,&
                             targetGID,&
                             overlapGID

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! moab variables
    integer :: meshType=3 ! finite volume

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    ! calculate mesh intersection
    ierr=iMOAB_ComputeMeshIntersectionOnSphere(sourcePID,&
                                               targetPID,&
                                               overlapPID)
    call errorout(ierr,'failed to compute mesh intersection')

    ! calculate coverage graph
    ierr=iMOAB_ComputeCommGraph(sourcePID,&
                                overlapPID,&
                                moabComm,&
                                moabGroup,&
                                moabGroup,&
                                meshType,&
                                meshType,&
                                sourceGID,&
                                overlapGID)
    call errorout(ierr,'failed to compute coverage graph')

!    ierr=iMOAB_DumpCommGraph(sourcePID,&
!                             overlapGID,&
!                             1,&
!                             'commgraph'//C_NULL_CHAR)
!    call errorout(ierr,'failed to dump coverage graph')

#endif

    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_calc_mesh_intersection!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_calc_weights
!
!> \brief   Calculates the weights for spatial mode splitting
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!>  Sends source mesh to target "side", calculates mesh intersection, 
!>  calculates coverage graph, calculates remapping weights.
!
!-----------------------------------------------------------------------
    subroutine ocn_spatial_mode_splitting_calc_weights(overlapPID,&
                                                       weights_identifier,&
                                                       fv_method)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! moab PID's
    integer,intent(inout) :: overlapPID

    ! remapping tags and method
    character(:),allocatable,intent(inout) :: weights_identifier,fv_method

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! iMOAB variables
    character(:),allocatable :: discMethod
    character(:),allocatable :: dofTagNames
    integer :: discOrder,fNoBubble,fMonotoneTypeID,fVolumetric,fNoConserve,&
               fValidate,fInverseDistanceMap

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    ! set some of the iMOAB parameters
    discMethod='fv'//C_NULL_CHAR
    discOrder=1
    dofTagNames='GLOBAL_ID'//C_NULL_CHAR
    fNoBubble=1
    fMonotoneTypeID=0
    fVolumetric=0
    fNoConserve=0
    fValidate=0
    fInverseDistanceMap=0

    ! calculate projection weights
    ierr=iMOAB_ComputeScalarProjectionWeights(overlapPID,&
                                              weights_identifier,&
                                              discMethod,&
                                              discOrder,&
                                              discMethod,&
                                              discOrder,&
                                              fv_method,&
                                              fNoBubble,&
                                              fMonotoneTypeID,&
                                              fVolumetric,&
                                              fInverseDistanceMap,&
                                              fNoConserve,&
                                              fValidate,&
                                              dofTagNames,&
                                              dofTagNames)
    call errorout(ierr,'unable to compute weights')

#endif

    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_calc_weights!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_define_tags
!
!> \brief   Adds tags for moab remapping
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!
!-----------------------------------------------------------------------
    subroutine ocn_spatial_mode_splitting_define_tags(sourcePID,&
                                                      targetPID,&
                                                      fieldTags,&
                                                      fieldTagsProjected,&
                                                      tagTypes,&
                                                      tagIndex)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! mesh PIDs
    integer,intent(inout) :: sourcePID,&
                             targetPID


    ! tag types and indices
    integer,dimension(2),intent(inout) :: tagTypes
    integer,dimension(2),intent(inout) :: tagIndex

    ! names of fields to be projected
    character(:),allocatable,intent(inout) :: fieldTags,fieldTagsProjected

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! number of degrees of freedom (always 1 for us)
    integer :: nDOFs

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    nDOFS=1 ! 

    ! define tags on source
    ierr=iMOAB_DefineTagStorage(sourcePID,&
                                fieldTags,&
                                tagTypes(1),&
                                nDOFS,&
                                tagIndex(1))
    call errorout(ierr,'failed to define tags on source')

    ! target fields
    ierr=iMOAB_DefineTagStorage(targetPID,&
                                fieldTagsProjected,&
                                tagTypes(2),&
                                nDOFS,&
                                tagIndex(2))
    call errorout(ierr,'failed to define tags on target')

#endif

    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_define_tags!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_allocate_buffers
!
!> \brief   Allocates buffers for iMOAB
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!
!-----------------------------------------------------------------------
    subroutine ocn_spatial_mode_splitting_allocate_buffers(sourcePID,&
                                                           source2target,&
                                                           target2source,&
                                                           nTags)!{{{

    !-----------------------------------------------------------------
    ! in/out variables
    !-----------------------------------------------------------------

    ! buffers
    real(kind=RKIND),dimension(:),allocatable,intent(inout) :: source2target,&
                                                               target2source

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! mesh PIDs
    integer,intent(inout) :: sourcePID

    ! number of tags
    integer,intent(in) :: nTags

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! mesh dimensions
    integer, dimension(3) :: nverts, nelem, nblocks, nsbc, ndbc

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB
    ! get mesh dimensions from moab
    ierr=iMOAB_GetMeshInfo(sourcePID,nverts,nelem,nblocks,nsbc,ndbc)
    call errorout(ierr,'failed to retrieve mesh info (source on source)')

    ! initialize buffers
    if ( .NOT. allocated( source2target ) ) then
       allocate ( source2target(nTags*nelem(3)) )
    endif
    if ( .NOT. allocated( target2source ) ) then
       allocate ( target2source(nTags*nelem(3)) )
    endif

    source2target(:)=0.0_RKIND
    target2source(:)=0.0_RKIND

#endif

    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_allocate_buffers!}}}
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

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------
    type(domain_type), intent(inout) :: domain

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------
    ! domain pointer for barotropic domain
    type (domain_type), pointer :: domain_btr

    ! configuration info
    logical, pointer :: config_spatial_mode_splitting

    ! pools
    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: statePool

    ! mesh information
    integer, pointer :: nCellsSolve_bc, nCellsSolve_btr,&
                        nEdgesSolve_bc, nEdgesSolve_btr,&
                        vertexDegree_btr, nEdgesAll_btr,&
                        nCellsAll_btr

    ! size info for moab
    integer :: totalNumVals, entity_type

    ! arrays for storing data
    real (kind=RKIND), dimension(:), pointer :: barotropicForcing_bc, &
       barotropicForcing_btr

    real (kind=RKIND), dimension(:,:), pointer :: normalBaroclinicVelocity_bc, &
       normalBaroclinicVelocity_btr

    real (kind=RKIND), dimension(:), allocatable :: normalBaroclinicVelocityBottom_bc,&
                                               normalBaroclinicVelocityBottom_btr

    integer, dimension(:), pointer :: maxLevelEdgeTop_bc, maxLevelEdgeTop_btr

    ! for going from cell centered to edge based quantities
    integer :: cell1, cell2
    real (kind=RKIND) :: zonal_edge, meridional_edge
    real (kind=RKIND), dimension(:), pointer :: angleEdge_btr
    integer, dimension(:,:), pointer :: cellsOnEdge_btr,&
                                        verticesOnEdge_btr,&
                                        cellsOnVertex_btr
    real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex_btr

    ! looping variables
    integer :: iCell, iEdge, k
    ! error handling
    integer :: ierr


    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    ! get configuration to see if we are using the mode splitting, if not return
    call mpas_pool_get_config(domain % configs, 'config_spatial_mode_splitting',&
                              config_spatial_mode_splitting)
    if (.not. config_spatial_mode_splitting) then
       return
    endif

#ifdef HAVE_MOAB

    ! "next" domain is barotropic domain
    domain_btr => domain % next

    !--------------------------------------------------------------------
    ! get barotropicForcing from baroclinic domain
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))

        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_bc)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_bc)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesAll_bc)

        ! get maxLevelEdgeTop
        call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop_bc)

        ! get variables on baroclinic domain
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing_bc)

        call mpas_dmpar_field_halo_exch(domain_btr, 'barotropicForcing')

        if ( config_additional_barotropic_terms ) then
           call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
              normalBaroclinicVelocity_bc, 1)

           call mpas_dmpar_field_halo_exch(domain_btr,'normalBaroclinicVelocity', &
               timeLevel=1)

           ! get baroclinic velocities only on bottom active layer
           allocate( normalBaroclinicVelocityBottom_bc(nEdgesAll_bc+1) )
           normalBaroclinicVelocityBottom_bc(:) = 0.0_RKIND
           
           do iEdge = 1,nEdgesAll_bc
              if (maxLevelEdgeTop_bc(iEdge).EQ.0) then
                 cycle
              endif
              normalBaroclinicVelocityBottom_bc(iEdge) = &
                  normalBaroclinicVelocity_bc(maxLevelEdgeTop(iEdge),iEdge)
           end do
        endif

        block_ptr => block_ptr % next
    enddo

    !--------------------------------------------------------------------
    ! remap barotropicForcing to barotropic domain
    !--------------------------------------------------------------------

    block_ptr => domain_btr % blocklist
    do while(associated(block_ptr))

       ! setup pools
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

       ! get num cells/edges
       call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesSolve_btr)
       call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesAll_btr)

       ! get variables on barotropic domain
       call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing_btr)

       call ocn_spatial_mode_splitting_dmpar_interpolate(domain, nEdgesSolve_bc, &
          nEdgesGlobal_bc, nEdgesPerProc_bc, nEdgesDisplacement_bc, indexToEdgeIDGathered_bc, &
          barotropicForcing_bc, nEdgesSolve_btr, nEdgesGlobal_btr, nEdgesPerProc_btr, &
          nEdgesDisplacement_btr, indexToEdgeIDGathered_btr, barotropicForcing_btr, &
          colValues_bc2btr, rowValues_bc2btr, sValues_bc2btr, nEdgesAll_btr, nEdgesAll_bc)

       call mpas_dmpar_field_halo_exch(domain_btr, 'barotropicForcing')

       if ( config_additional_barotropic_terms ) then
          call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
             normalBaroclinicVelocity_btr, 1)

          call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop_btr)

          allocate (normalBaroclinicVelocityBottom_btr(nEdgesAll_btr+1) )
          normalBaroclinicVelocityBottom_btr(:) = 0.0_RKIND

          call ocn_spatial_mode_splitting_dmpar_interpolate(domain, nEdgesSolve_bc, &
             nEdgesGlobal_bc, nEdgesPerProc_bc, nEdgesDisplacement_bc, indexToEdgeIDGathered_bc, &
             normalBaroclinicVelocityBottom_bc, nEdgesSolve_btr, nEdgesGlobal_btr, nEdgesPerProc_btr, &
             nEdgesDisplacement_btr, indexToEdgeIDGathered_btr, normalBaroclinicVelocityBottom_btr, &
             colValues_bc2btr, rowValues_bc2btr, sValues_bc2btr, nEdgesAll_btr, nEdgesAll_bc)

          do iEdge = 1,nEdgesAll_btr
             if (maxLevelEdgeTop_btr(iEdge).EQ.0) then
                cycle
             endif
             normalBaroclinicVelocity_btr(maxLevelEdgeTop_btr(iEdge),iEdge) = &
                normalBaroclinicVelocityBottom_btr(iEdge)
          enddo

          call mpas_dmpar_field_halo_exch(domain_btr, 'normalBaroclinicVelocity', &
               timeLevel=1)

          deallocate( normalBaroclinicVelocityBottom_bc, normalBaroclinicVelocityBottom_btr )
       endif

       block_ptr => block_ptr % next
    enddo
#endif
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
    ! domain pointer for barotropic domain
    type (domain_type), pointer :: domain_btr
    ! configuration info
    logical, pointer :: config_spatial_mode_splitting
    ! 
    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: statePool

    ! mesh information
    integer, pointer :: nCellsSolve_bc, nCellsSolve_btr,&
                        nEdgesSolve_bc, nEdgesSolve_btr,&
                        nEdgesAll_bc, nEdgesAll_btr

    ! arrays for storing data
    real (kind=RKIND), dimension(:), pointer :: barotropicThicknessFlux_bc, &
       barotropicThicknessFlux_btr, accumulatedDragOverDepth_bc, &
       accumulatedDragOverDepth_btr, accumulatedDrag_bc, accumulatedDrag_btr

    ! for going from cell centered to edge based quantities
    integer :: cell1, cell2
    real (kind=RKIND) :: zonal_edge, meridional_edge
    real (kind=RKIND), dimension(:), pointer :: angleEdge_bc
    integer, dimension(:,:), pointer :: cellsOnEdge_bc

    ! looping variables
    integer :: iCell, iEdge

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    ! get configuration to see if we are using the mode splitting, if not return
    call mpas_pool_get_config(domain % configs, 'config_spatial_mode_splitting',&
                              config_spatial_mode_splitting)
    if (.not. config_spatial_mode_splitting) then
       return
    endif

#ifdef HAVE_MOAB

    ! "next" domain is barotropic domain
    domain_btr => domain % next

    !--------------------------------------------------------------------
    ! get barotropicThicknessFlux from barotropic domain
    !--------------------------------------------------------------------

    block_ptr => domain_btr % blocklist
    do while(associated(block_ptr))

        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesAll_btr)

        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                                 barotropicThicknessFlux_btr)
        call mpas_dmpar_field_halo_exch(domain_btr,'barotropicThicknessFlux')

        if ( config_additional_barotropic_terms ) then
           call mpas_pool_get_array(diagnosticsPool, 'accumulatedDragOverDepth',&
                                    accumulatedDragOverDepth_btr)
           call mpas_dmpar_field_halo_exch(domain_btr,'accumulatedDragOverDepth')

           call mpas_pool_get_array(diagnosticsPool, 'accumulatedDrag',&
                                    accumulatedDrag_btr)
           call mpas_dmpar_field_halo_exch(domain_btr,'accumulatedDrag')

        endif

        block_ptr => block_ptr % next
    enddo

    !--------------------------------------------------------------------
    ! remap barotropicThicknessFlux to baroclinic domain
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))

        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_bc)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_bc)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesAll_bc)

        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                                 barotropicThicknessFlux_bc)

        call ocn_spatial_mode_splitting_dmpar_interpolate(domain, nEdgesSolve_btr, &
           nEdgesGlobal_btr, nEdgesPerProc_btr, nEdgesDisplacement_btr, indexToEdgeIDGathered_btr, &
           barotropicThicknessFlux_btr, nEdgesSolve_bc, nEdgesGlobal_bc, nEdgesPerProc_bc, &
           nEdgesDisplacement_bc, indexToEdgeIDGathered_bc, barotropicThicknessFlux_bc, &
           colValues_btr2bc, rowValues_btr2bc, sValues_btr2bc, nEdgesAll_btr, nEdgesAll_bc)

        call mpas_dmpar_field_halo_exch(domain, 'barotropicThicknessFlux')

        if ( config_additional_barotropic_terms ) then
           call mpas_pool_get_array(diagnosticsPool, 'accumulatedDragOverDepth',&
                                    accumulatedDragOverDepth_bc)

           call ocn_spatial_mode_splitting_dmpar_interpolate(domain, nEdgesSolve_btr, &
              nEdgesGlobal_btr, nEdgesPerProc_btr, nEdgesDisplacement_btr, indexToEdgeIDGathered_btr, &
              accumulatedDragOverDepth_btr, nEdgesSolve_bc, nEdgesGlobal_bc, nEdgesPerProc_bc, &
              nEdgesDisplacement_bc, indexToEdgeIDGathered_bc, accumulatedDragOverDepth_bc, &
              colValues_btr2bc, rowValues_btr2bc, sValues_btr2bc,nEdgesAll_btr, nEdgesAll_bc)

           call mpas_dmpar_field_halo_exch(domain,'accumulatedDragOverDepth')

           call mpas_pool_get_array(diagnosticsPool, 'accumulatedDrag',&
                                    accumulatedDrag_bc)

           call ocn_spatial_mode_splitting_dmpar_interpolate(domain, nEdgesSolve_btr, &
              nEdgesGlobal_btr, nEdgesPerProc_btr, nEdgesDisplacement_btr, indexToEdgeIDGathered_btr, &
              accumulatedDrag_btr, nEdgesSolve_bc, nEdgesGlobal_bc, nEdgesPerProc_bc, &
              nEdgesDisplacement_bc, indexToEdgeIDGathered_bc, accumulatedDrag_bc, &
              colValues_btr2bc, rowValues_btr2bc, sValues_btr2bc, nEdgesAll_btr, nEdgesAll_bc)

           call mpas_dmpar_field_halo_exch(domain,'accumulatedDrag')
        endif

        block_ptr => block_ptr % next
    enddo
#endif
    end subroutine ocn_spatial_mode_splitting_btr2bc_remapping!}}}

!***********************************************************************
!
!  ocn_spatial_mode_splitting_mesh_redirect
!
!> \brief redirects pointers for barotropic subcycling
!> \author Coleman Blakely
!> \date   April, 2024
!> \details
!> This subroutine redirects the pointers to various arrays, values, etc.
!> between the barotropic and baroclinic domains for spatial mode 
!> splitting before and after the barotropic subcycling occurs.
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_setup_domain(domain,&
                                                       statePool,&
                                                       forcingPool,&
                                                       tendPool,&
                                                       sshCur,&
                                                       sshNew,&
                                                       normalBarotropicVelocityCur,&
                                                       normalBarotropicVelocityNew,&
                                                       sshTend,&
                                                       normalBaroclinicVelocityCur,&
                                                       accumulatedDragOverDepth,&
                                                       accumulatedDrag) !{{{

    !-----------------------------------------------------------------
    ! IMPORTANT NOTE:
    !-----------------------------------------------------------------
    ! In addition to redirecting the in/out variables defined below,
    ! this subroutine also redirects the following public variables
    ! for use in the barotropic subcycle when spatial mode splitting 
    ! is used:
    !
    !    diagnostics pool:
    !       barotropicForcing
    !       barotropicThicknessFlux
    !
    ! This subroutine is called immediately before and after the 
    ! barotropic subcycling and thus resets to the baroclinic
    ! domain after each barotropic subcycle is complete.
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------
    type(domain_type), intent(inout) :: domain
    type(mpas_pool_type), pointer, intent(inout) :: statePool,&
                                           forcingPool,&
                                           tendPool
    
    ! statePool
    real (kind=RKIND), pointer, intent(inout), dimension(:) :: sshCur,&
                                                      sshNew,&
                                                      normalBarotropicVelocityCur,&
                                                      normalBarotropicVelocityNew
    real (kind=RKIND), pointer, intent(inout), dimension(:,:) :: normalBaroclinicVelocityCur
                                                       
    ! tendPool
    real (kind=RKIND), pointer, intent(inout), dimension(:) :: sshTend

    ! diagnosticsPool
    real (kind=RKIND), pointer, intent(inout), dimension(:) :: accumulatedDragOverDepth,&
       accumulatedDrag

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    type(block_type), pointer :: block  ! variables in current subblock

    type(mpas_pool_type), pointer :: diagnosticsPool

    ! loop indices
    integer :: i, k, n, iEdge
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    ! perform remapping: **NOTE** if the argument to this subroutine has 
    ! domainID = 1, that means we are switching from the baroclinic to
    ! the barotropic domain and thus have to remap barotropicForcing. 
    ! This is opposite what one might think initially.
    if ( domain % domainID .eq. 0 ) then

       ! baroclinic domain was passed. end of barotropic subcycle. 
       ! send barotropicThicknessFlux from barotropic to baroclinic
       !call ocn_spatial_mode_splitting_btr2bc_remapping(domain)
    elseif ( domain % domainID .eq. 1 ) then

       ! barotropic domain was passed. beginning of barotropic subcycle.
       ! send barotropicForcing from baroclinic to barotropic.
       !call ocn_spatial_mode_splitting_bc2btr_remapping(domain)

    endif

    ! destroy current mesh data
    call ocn_meshDestroy(ierr)

    ! create new mesh
    call ocn_meshCreate(domain)

    ! redirect pools
    block => domain%blocklist

    call mpas_pool_get_subpool(block%structs, 'state', statePool)

    call mpas_pool_get_subpool(block%structs, 'forcing', forcingPool)

    call mpas_pool_get_subpool(block%structs, 'tend', tendPool)

    ! redirect arrays from statePool
    call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)

    call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

    call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                             normalBarotropicVelocityCur, 1)

    call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                             normalBarotropicVelocityNew, 2)

    ! redirect arrays from tendency pool
    call mpas_pool_get_array(tendPool, 'ssh', sshTend)

    ! redirect barotropicForcing and barotropicThicknessFlux
    call mpas_pool_get_subpool(block%structs, 'diagnostics', diagnosticsPool)

    call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

    call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

    ! if we are using additional barotropic terms also redirect those
    if ( config_additional_barotropic_terms ) then 
       call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
          normalBaroclinicVelocityCur, 1)
       call mpas_pool_get_array(diagnosticsPool, 'accumulatedDragOverDepth', accumulatedDragOverDepth)
       call mpas_pool_get_array(diagnosticsPool, 'accumulatedDrag', accumulatedDrag)
    endif

    ! reset edgeMask appropriately if on baroclinic domain
    if (domain % domainID .eq. 0) then
       ! populate edgeMask with edgeMask_spatial_split_bc
       do iEdge=1,nEdgesAll
          !maxLevelEdgeTop(iEdge)=maxLevelEdgeTop_spatial_split(iEdge)
          do k=1,nVertLevels
             edgeMask(k,iEdge)=edgeMask_spatial_split_bc(k,iEdge)
          enddo
       enddo
    endif

!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_setup_domain !}}}

!***********************************************************************
!
!  ocn_spatial_mode_splitting_time_int_spit_setup_btr
!
!> \brief initializes split-explicit variables on barotropic domain
!> \author Coleman Blakely
!> \date   April, 2024
!> \details
!> This subroutine initializes the variables needed in the barotropic
!> subcycling of the split-explicit time integration on the barotropic
!> domain when spatial mode splitting is active. At present, that
!> only involves setting sshNew_btr equal to sshCur_btr because 
!> normalBarotropicVelocity is handled within the split-explicit time
!> stepping routine
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_time_int_split_setup_btr(domain) !{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type(domain_type), intent(inout) :: domain

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    type(block_type), pointer :: block_ptr  ! variables in current subblock

    type(mpas_pool_type), pointer :: meshPool, statePool

    ! meshPool
    integer, pointer :: nCellsAll_btr
    
    ! statePool
    real (kind=RKIND), pointer, dimension(:) :: sshCur_btr,&
                                                sshNew_btr

    ! loop indices
    integer :: iCell

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    block_ptr => domain % blocklist
    do while(associated(block_ptr))

        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

        ! get nCells and sea surface heights
        call mpas_pool_get_dimension(meshPool, 'nCells', nCellsAll_btr)
        call mpas_pool_get_array(statePool, 'ssh', sshCur_btr, 1)
        call mpas_pool_get_array(statePool, 'ssh', sshNew_btr, 2)

        do iCell = 1, nCellsAll_btr
           sshNew_btr(iCell) = sshCur_btr(iCell)
        enddo

        block_ptr => block_ptr % next
    enddo
!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_time_int_split_setup_btr !}}}

!***********************************************************************
!
!  ocn_spatial_mode_splitting_dmpar_remapping_init
!
!> \brief initializes dmpar remapping
!> \author Coleman Blakely
!> \date   March 2025
!> \details
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_dmpar_remapping_init(domain,&
                                                               weights_filename,&
                                                               rowValues,&
                                                               colValues,&
                                                               sValues,&
                                                               nEdgesGlobal,&
                                                               nEdgesDisplacement,&
                                                               indexToEdgeIDGathered,&
                                                               nEdgesPerProc) !{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type(domain_type), intent(in) :: domain
    character(:), allocatable, intent(in) :: weights_filename

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------

    ! weights stored as sparse array
    integer, dimension(:), allocatable, intent(inout) :: rowValues, colValues
    real(kind=RKIND), dimension(:), allocatable, intent(inout) :: sValues

    ! for communication
    integer, intent(inout) :: nEdgesGlobal
    integer, dimension(:), allocatable, intent(inout) :: nEdgesDisplacement,&
                                                         indexToEdgeIDGathered,&
                                                         nEdgesPerProc

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    type(block_type), pointer :: block_ptr  ! variables in current subblock
    type(mpas_pool_type), pointer :: meshPool

    ! meshPool
    integer, dimension(:), pointer :: nEdgesArrayLocal
    integer :: nEdgesOwnedLocal
    integer, dimension(:), pointer :: indexToEdgeIDLocal

    ! netcdf variables
    integer :: ncid ! file id
    integer :: ns_dim_id, na_dim_id ! dimension id's
    integer :: row_id, col_id, s_id ! variable id's
    character(len=NF90_MAX_NAME) :: dummy_char
    integer :: n_s, n_a

    ! loop indices
    integer :: iCell, iProc

    ! for error handling
    integer :: err,err_tmp

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    ! Begin MPI portion
    call MPI_COMM_RANK( domain % dminfo % comm, curProc, err_tmp)
    err = ior(err, err_tmp)
    call MPI_COMM_SIZE( domain % dminfo % comm, nProcs, err_tmp)
    err = ior(err, err_tmp)
 
    if (curProc.eq.0) then
        allocate(nEdgesPerProc(nProcs))
        allocate(nEdgesDisplacement(nProcs))
    else
        allocate(nEdgesPerProc(1))
        allocate(nEdgesDisplacement(1))
    endif

    ! Gather nEdgesOwned
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArrayLocal)
    call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeIDLocal)
    nEdgesOwnedLocal=nEdgesArrayLocal(1)
    call MPI_GATHER( nEdgesOwnedLocal, 1, MPI_INTEGER, nEdgesPerProc, 1, MPI_INTEGER, &
                     0, domain % dminfo % comm, err_tmp)
    err = ior(err, err_tmp)
 
    ! Set Displacement variable for GATHERV command
    if (curProc.eq.0) then
        nEdgesGlobal = sum(nEdgesPerProc)
        allocate(indexToEdgeIDGathered(nEdgesGlobal))
        nEdgesDisplacement(1) = 0
        if (nProcs > 1) then
            do iProc=2,nProcs
                nEdgesDisplacement(iProc) = nEdgesDisplacement(iProc-1) + nEdgesPerProc(iProc-1)
            enddo
        endif
    else
        allocate(indexToEdgeIDGathered(1))
    endif
 
    ! Gather indexToEdgeID
    call MPI_GATHERV( indexToEdgeIDLocal, nEdgesOwnedLocal, MPI_INTEGER, indexToEdgeIDGathered, &
            nEdgesPerProc, nEdgesDisplacement, MPI_INTEGER, 0, domain % dminfo % comm, err_tmp)
    err = ior(err, err_tmp)
 
    ! initialize map only on proc 0

    ! don't do this on baroclinic domain since we don't have the weights yet
!    if (domain % domainID .eq. 0) then
!       return
!    endif
    if (curProc.eq.0) then
        ! Open netcdf weights files
        call check( nf90_open(path = trim(weights_filename), mode = nf90_nowrite, ncid = ncid) ,err)
 
        ! Get dimension ID
        call check( nf90_inq_dimid(ncid, "n_s", ns_dim_id) ,err)
        call check( nf90_inq_dimid(ncid, "n_a", na_dim_id) ,err)
 
        ! Get Variable IDs
        call check( nf90_inq_varid(ncid, "row", row_id) ,err)
        call check( nf90_inq_varid(ncid, "col", col_id) ,err)
        call check( nf90_inq_varid(ncid, "S", s_id) ,err)
 
        ! Get Dimension Length 
        call check( nf90_inquire_dimension(ncid, ns_dim_id, dummy_char, n_s) ,err)
        call check( nf90_inquire_dimension(ncid, na_dim_id, dummy_char, n_a) ,err)
 
        ! Allocate matrices to read data into
        allocate ( rowValues (n_s) )
        allocate ( colValues (n_s) )
        allocate ( sValues (n_s) )
 
        ! Retrieve data
        call check( nf90_get_var(ncid, col_id, colValues(:) ) ,err)
        call check( nf90_get_var(ncid, row_id, rowValues(:) ) ,err)
        call check( nf90_get_var(ncid, s_id, sValues(:) ) ,err)
    endif
!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_dmpar_remapping_init !}}}

!***********************************************************************
!
!  ocn_spatial_mode_splitting_dmpar_interpolate
!
!> \brief gathers and interpolates using dmpar routines
!> \author Coleman Blakely
!> \date   March, 2025
!> \details
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_dmpar_interpolate(domain,&
                                                            nEdgesLocal_source,&
                                                            nEdgesGlobal_source,&
                                                            nEdgesPerProc_source,&
                                                            nEdgesDisplacement_source,&
                                                            indexToEdgeIDGathered_source,&
                                                            localArray_source,&
                                                            nEdgesLocal_target,&
                                                            nEdgesGlobal_target,&
                                                            nEdgesPerProc_target,&
                                                            nEdgesDisplacement_target,&
                                                            indexToEdgeIDGathered_target,&
                                                            localArray_target,&
                                                            colValues,&
                                                            rowValues,&
                                                            sValues,&
                                                            nEdgesAll_source,&
                                                            nEdgesAll_target) !{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (domain_type), intent(in) :: domain

    ! source
    integer, intent(in) :: nEdgesLocal_source,nEdgesGlobal_source
    integer, dimension(:), intent(in) :: nEdgesPerProc_source,&
                                         nEdgesDisplacement_source,&
                                         indexToEdgeIDGathered_source

    ! target
    integer, intent(in) :: nEdgesLocal_target,nEdgesGlobal_target
    integer, dimension(:), intent(in) :: nEdgesPerProc_target,&
                                         nEdgesDisplacement_target,&
                                         indexToEdgeIDGathered_target

    ! weights
    integer, dimension(:), intent(in) :: colValues, rowValues
    real(kind=RKIND), dimension(:), intent(in) :: sValues

    ! TEMP: for carry over
    integer, intent(in) :: nEdgesAll_source, nEdgesAll_target

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------
    real(kind=RKIND), dimension(:), intent(inout) :: localArray_source,&
                                                     localArray_target

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    type(block_type), pointer :: block_ptr  ! variables in current subblock
    type(mpas_pool_type), pointer :: meshPool

    ! send buffers
    real(kind=RKIND), dimension(:), allocatable :: sendBuffer_source, sendBuffer_target

    ! meshPool
    integer :: nEdges

    ! loop indices
    integer :: iEdge,idx

    ! arrays for storing
    real(kind=RKIND), dimension(:), allocatable :: gatheredArray_source,&
                                                   globalArray_source,&
                                                   gatheredArray_target,&
                                                   globalArray_target

    ! for error handling
    integer :: err

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    if ( carryOverVals ) then
       if (nEdgesAll_source .NE. nEdgesAll_target) then
          call mpas_log_write( 'ERROR: You are trying to carry over values but your &
             &meshes are not the same size',MPAS_LOG_CRIT)
       endif
       do iEdge = 1,nEdgesAll_target
          localArray_target(iEdge) = localArray_source(iEdge)
       enddo
       return
    endif


    ! allocate global arrays
    if (curProc.eq.0) then
       allocate(globalArray_source(nEdgesGlobal_source),&
                gatheredArray_source(nEdgesGlobal_source),&
                globalArray_target(nEdgesGlobal_target),&
                gatheredArray_target(nEdgesGlobal_target))
       globalArray_source(:)=0.0_RKIND
       gatheredArray_source(:)=0.0_RKIND
       globalArray_target(:)=0.0_RKIND
       gatheredArray_target(:)=0.0_RKIND
             
    else
       allocate(globalArray_source(1),&
                gatheredArray_source(1),&
                globalArray_target(1),&
                gatheredArray_target(1))
    endif

    ! gather values on proc 0
    call MPI_GATHERV(localArray_source,&
                     nEdgesLocal_source,&
                     MPI_DOUBLE,&
                     gatheredArray_source,&
                     nEdgesPerProc_source,&
                     nEdgesDisplacement_source,&
                     MPI_DOUBLE,&
                     0,&
                     domain % dminfo % comm,&
                     err)

    ! interpolate
    if (curProc.eq.0) then

       ! put values in correct order
       do iEdge = 1,nEdgesGlobal_source
          idx=indexToEdgeIDGathered_source(iEdge)
          globalArray_source(idx)=gatheredArray_source(iEdge)
       enddo

       ! interpolate
       call interpolate(colValues, rowValues, sValues,&
         globalArray_source, globalArray_target)

       ! put back in order to scatter
       do iEdge = 1,nEdgesGlobal_target
          idx=indexToEdgeIDGathered_target(iEdge)
          gatheredArray_target(iEdge)=globalArray_target(idx)
       enddo
    endif

    ! scatter
    localArray_target(:)=0.0_RKIND
    call MPI_SCATTERV(gatheredArray_target,&
                      nEdgesPerProc_target,&
                      nEdgesDisplacement_target,&
                      MPI_DOUBLE,&
                      localArray_target,&
                      nEdgesLocal_target,&
                      MPI_DOUBLE,&
                      0,&
                      domain % dminfo % comm,&
                      err)

    ! deallocate arrays
    deallocate(globalArray_source, globalArray_target,&
               gatheredArray_source, gatheredArray_target)

!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_dmpar_interpolate !}}}

!***********************************************************************
!
!  ocn_spatial_mode_splitting_dmpar_init_bc
!
!> \brief gathers and interpolates using dmpar routines
!> \author Coleman Blakely
!> \date   March, 2025
!> \details
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_dmpar_init_bc(domain) !{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (domain_type), intent(in) :: domain

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    weights_filename_btr2bc='soma_weights_32km_to_16km.nc'
    call ocn_spatial_mode_splitting_dmpar_remapping_init(domain,&
                                                         weights_filename_btr2bc,&
                                                         rowValues_bc2btr,&
                                                         colValues_bc2btr,&
                                                         sValues_bc2btr,&
                                                         nEdgesGlobal_bc,&
                                                         nEdgesDisplacement_bc,&
                                                         indexToEdgeIDGathered_bc,&
                                                         nEdgesPerProc_bc)


!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_dmpar_init_bc !}}}

!***********************************************************************
!
!  ocn_spatial_mode_splitting_dmpar_init_btr
!
!> \brief gathers and interpolates using dmpar routines
!> \author Coleman Blakely
!> \date   March, 2025
!> \details
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_dmpar_init_btr(domain) !{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (domain_type), intent(in) :: domain

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    weights_filename_btr2bc='soma_weights_16km_to_32km.nc'
    call ocn_spatial_mode_splitting_dmpar_remapping_init(domain,&
                                                         weights_filename_btr2bc,&
                                                         rowValues_btr2bc,&
                                                         colValues_btr2bc,&
                                                         sValues_btr2bc,&
                                                         nEdgesGlobal_btr,&
                                                         nEdgesDisplacement_btr,&
                                                         indexToEdgeIDGathered_btr,&
                                                         nEdgesPerProc_btr)


    end subroutine ocn_spatial_mode_splitting_dmpar_init_btr !}}}

!***********************************************************************
!
!  ocn_spatial_mode_splitting_set_edgeMask_btr
!
!> \brief sets edgeMask_btr to allow for barotropic but not baroclinic 
!>        flow and fluxes to propagate through the edges of the
!>        baroclinic domain.
!> \author Coleman Blakely
!> \date   March, 2025
!> \details
!
!-----------------------------------------------------------------------

    subroutine ocn_spatial_mode_splitting_set_edgeMask_btr(domain) !{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type (domain_type), intent(in) :: domain

    !-----------------------------------------------------------------
    ! input/output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! pools
    type (mpas_pool_type), pointer :: meshPool

    ! mesh dimensions
    integer, pointer :: nEdgesAll_bc, nVertLevels_bc

    ! for finding the bottom of the cell
    real (kind=RKIND), dimension(:), pointer :: bottomDepth_bc
    integer, dimension(:,:), pointer :: cellsOnEdge_bc
    integer, dimension(:), pointer :: maxLevelEdgeTop_bc, maxLevelCell_bc

    ! for loops
    integer :: iCell, iVertex, iEdge, k, kmin, kmax

    !
    integer :: cell1,cell2,vertex1,vertex2

    !
    integer :: edge_index, vertex_index, cell_index

    logical :: is_boundary_edge, is_boundary_cell

    integer :: nVerticesOnEdge=2

    ! for identifying both cells which are on the boundary and cells which border those on the
    ! boundary
    integer, dimension(:), allocatable :: isBoundaryCell, isSecondRowCell

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

    ! preallocate our arrays to tag boundary cells and boundary edges
    allocate( isBoundaryCell(nCellsAll+1), isSecondRowCell(nCellsAll+1))
    isBoundaryCell(:)=0
    isSecondRowCell(:)=0

    ! loop over all cells and identify all boundary cells
    do iCell=1,nCellsAll
       do iEdge=1,nEdgesOnCell(iCell)
          edge_index=edgesOnCell(iEdge,iCell)
          ! only have to check if cell2=nCellsAll+1 since cell1 never will
          cell2=cellsOnEdge(2,edge_index)
          if (cell2.EQ.nCellsAll+1) then
             isBoundaryCell(iCell)=1
          endif
       enddo
    enddo

    ! loop over again and find all second row cells
    do iCell=1,nCellsAll
       do iEdge=1,nEdgesOnCell(iCell)
          edge_index=edgesOnCell(iEdge,iCell)
          cell1=cellsOnEdge(1,edge_index)
          cell2=cellsOnEdge(2,edge_index)
          if ((isBoundaryCell(cell1).EQ.1).OR.&
               (isBoundaryCell(cell2).EQ.1)) then
             isSecondRowCell(iCell)=1
          endif
       enddo
    enddo


    ! get bottom depth of baroclinic mesh
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_array(meshPool, 'edgeMask_spatial_split_bc',edgeMask_spatial_split_bc)
    call mpas_pool_get_array(meshPool, 'edgeMask_spatial_split_btr',edgeMask_spatial_split_btr)

    ! we want edgeMask_btr to essentially be the same as edge mask but with boundary cells having
    ! values of 1 down to maxLevelCell
    edgeMask_spatial_split_btr(:,:)=0.0_RKIND

    ! recreate edgeMask
    do iEdge=1,nEdgesAll
       do k=minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
          edgeMask_spatial_split_btr(k,iEdge)=1d0
       enddo
    enddo

    ! check boundary cells and fill them out
    do iCell=1,nCellsAll
       if (isBoundaryCell(iCell).EQ.1) then
          ! identify kmax 
          kmax=maxLevelCell(iCell)
          do k=1,maxEdges
             cell1=cellsOnCell(k,iCell)
             if (maxLevelCell(cell1).EQ.0) then
                cycle
             endif
             kmax=MIN(kmax,maxLevelCell(cell1))
          enddo
          do iEdge=1,nEdgesOnCell(iCell)
             do k=1,kmax
                edgeMask_spatial_split_btr(k,edgesOnCell(iEdge,iCell))=1d0
             enddo
          enddo
       endif
    enddo

    ! we want edgeMask_bc to have zeroes on all edges which are only part of boundary cells. This
    ! happens when at least one vertex of an edge has a value of nCellsAll+1 in cellsOnVertex
    edgeMask_spatial_split_bc(:,:)=0.0_RKIND

    ! recreate edgeMask
    do iEdge=1,nEdgesAll
       do k=minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
          edgeMask_spatial_split_bc(k,iEdge)=1d0
       enddo
    enddo

    do iEdge=1,nEdgesAll
       cell1=cellsOnEdge(1,iEdge)
       cell2=cellsOnEdge(2,iEdge)
       if ((isBoundaryCell(cell1).EQ.1).OR.(isBoundaryCell(cell2).EQ.1)) then
          edgeMask_spatial_split_bc(:,iEdge)=0.0_RKIND
       elseif ((isBoundaryCell(cell1).EQ.1).AND.(isSecondRowCell(cell2).EQ.1)) then
          edgeMask_spatial_split_bc(:,iEdge)=0.0_RKIND
       elseif ((isBoundaryCell(cell2).EQ.1).AND.(isSecondRowCell(cell1).EQ.1)) then
          edgeMask_spatial_split_bc(:,iEdge)=0.0_RKIND
       elseif ((isSecondRowCell(cell1).EQ.1).AND.(isSecondRowCell(cell2).EQ.1)) then
          edgeMask_spatial_split_bc(:,iEdge)=0.0_RKIND
       endif
    enddo

    ! populate edgeMask with edgeMask_spatial_split_bc
    do iEdge=1,nEdgesAll
       do k=1,nVertLevels
          edgeMask(k,iEdge)=edgeMask_spatial_split_bc(k,iEdge)
       enddo
    enddo

!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_set_edgeMask_btr !}}}

end module ocn_spatial_mode_splitting
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
