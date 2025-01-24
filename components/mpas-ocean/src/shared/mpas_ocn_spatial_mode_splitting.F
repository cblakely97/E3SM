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
              ocn_spatial_mode_splitting_bottomDepth_remapping

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------

    ! mpi groups/comms. All of these span the same processors
    integer :: bcComm, btrComm
    integer :: bcGroup, btrGroup
   
    ! identifiers for moab meshes
    ! barotropic mesh identifiers
    integer :: btrOnBtr
    integer :: btrOnBc
    integer :: btrXbcOnBc

    ! baroclinic mesh identifiers
    integer :: bcOnBc
    integer :: bcOnBtr
    integer :: bcXbtrOnBtr

    ! PID's for MOAB
    integer,public :: btrOnBtrPID, btrOnBcPID, btrXbcOnBcPID
    integer,public :: bcOnBcPID, bcOnBtrPID, bcXbtrOnBtrPID

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

   ! The following variables are for conservatively remapping the bathymetry from barotropic to
   ! baroclinic domain.
   character(:), allocatable :: bottomDepthTag,bottomDepthTag_btr
   integer, dimension(2) :: bottomDepthIndex, bottomDepthType
   character(:), allocatable :: bottomDepth_weights_identifier

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
    btrOnBtr    = 314
    btrOnBc     = 315
    btrXbcOnBc  = 316

    ! baroclinic mesh identifiers
    bcOnBc      = 414
    bcOnBtr     = 415
    bcXbtrOnBtr = 416

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
       bcComm = MPI_COMM_NULL
       call mpi_comm_dup(domain % dminfo % comm, &
                         bcComm, ierr)
       call mpi_comm_group(bcComm, bcGroup, ierr)

       ! initialize moab domain
       call init_moab_mpas(domain, bcOnBc, bcOnBcPID, bcComm)

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

       ! get mpi communicator for barotropic domain
       btrComm = MPI_COMM_NULL
       call mpi_comm_dup(domain % dminfo % comm, &
                         btrComm, ierr)
       call mpi_comm_group(btrComm, btrGroup, ierr)

       ! initialize moab domain
       call init_moab_mpas(domain, btrOnBtr, btrOnBtrPID, btrComm)

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

    ! register applications other than btrOnBtr and bcOnBc

    ! register baroclinic on barotropic side
    appname = 'bcOnBtr'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, &
                                     btrComm, &
                                     bcOnBtr, &
                                     bcOnBtrPID)
    call errorout(ierr, "failed to register bcOnBtr application")

    ! register baroclinic to barotropic intersection on baroclinic side
    appname = 'bcXbtrOnBtr'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, &
                                     btrComm, &
                                     bcXbtrOnBtr, &
                                     bcXbtrOnBtrPID)
    call errorout(ierr, "failed to register bcXbtrOnBtr application")

    ! barotropic on baroclinic side
    appname = 'btrOnBc'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname,&
                                     bcComm, &
                                     btrOnBc, &
                                     btrOnBcPID)
    call errorout(ierr, "failed to register btrOnBc application")

    ! barotropic to baroclinic intersection on baroclinic side
    appname = 'btrXbcOnBc'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname,&
                                     bcComm, &
                                     btrXbcOnBc, &
                                     btrXbcOnBcPID)
    call errorout(ierr, "failed to register btrXbcOnBc application")

    !--------------------------------------------------------------------
    ! calculate the mesh intersections
    !--------------------------------------------------------------------

    ! baroclinic => barotropic
    call ocn_spatial_mode_splitting_calc_mesh_intersection(bcOnBcPID,&
                                                           bcOnBtrPID,&
                                                           btrOnBtrPID,&
                                                           bcXbtrOnBtrPID,&
                                                           bcOnBc,&
                                                           bcOnBtr,&
                                                           btrOnBtr,&
                                                           bcXbtrOnBtr,&
                                                           bcComm,&
                                                           btrComm,&
                                                           bcGroup,&
                                                           btrGroup)

    ! barotropic => baroclinic
    call ocn_spatial_mode_splitting_calc_mesh_intersection(btrOnBtrPID,&
                                                           btrOnBcPID,&
                                                           bcOnBcPID,&
                                                           btrXbcOnBcPID,&
                                                           btrOnBtr,&
                                                           btrOnBc,&
                                                           bcOnBc,&
                                                           btrXbcOnBc,&
                                                           btrComm,&
                                                           bcComm,&
                                                           btrGroup,&
                                                           bcGroup)

    !--------------------------------------------------------------------
    ! calculate the projection weights
    !--------------------------------------------------------------------

    ! baroclinic => barotropic
    bc2btr_weights_identifier='bc2btr'//C_NULL_CHAR
    fv_method='bilin'//C_NULL_CHAR
    call ocn_spatial_mode_splitting_calc_weights(bcXbtrOnBtrPID,&
                                                 bcComm,&
                                                 bc2btr_weights_identifier,&
                                                 fv_method)

    ! barotropic => baroclinic
    btr2bc_weights_identifier='btr2bc'//C_NULL_CHAR
    fv_method=''//C_NULL_CHAR ! conservative
    call ocn_spatial_mode_splitting_calc_weights(btrXbcOnBcPID,&
                                                 bcComm,&
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
    call ocn_spatial_mode_splitting_define_tags(bcOnBcPID,&
                                                bcOnBtrPID,&
                                                btrOnBtrPID,&
                                                bcComm,&
                                                btrComm,&
                                                bc2btrFields,&
                                                bc2btrFields_btr,&
                                                bc2btrTypes,&
                                                bc2btrIndices)

    ! barotropic => baroclinic
    call ocn_spatial_mode_splitting_define_tags(btrOnBtrPID,&
                                                btrOnBcPID,&
                                                bcOnBcPID,&
                                                btrComm,&
                                                bcComm,&
                                                btr2bcFields_btr,&
                                                btr2bcFields,&
                                                btr2bcTypes,&
                                                btr2bcIndices)


    !--------------------------------------------------------------------
    ! allocate buffers
    !--------------------------------------------------------------------

    ! baroclinic to barotropic
    call ocn_spatial_mode_splitting_allocate_buffers(bcOnBcPID,&
                                                     bcOnBtrPID,&
                                                     bcComm,&
                                                     btrComm,&
                                                     bc2btrVar,&
                                                     btr2bcVar,&
                                                     nTags)

    ! barotropic to baroclinic
    call ocn_spatial_mode_splitting_allocate_buffers(btrOnBtrPID,&
                                                     btrOnBcPID,&
                                                     btrComm,&
                                                     bcComm,&
                                                     btr2bcVar_btr,&
                                                     bc2btrVar_btr,&
                                                     nTags)

    !--------------------------------------------------------------------
    ! remap bottomDepth from barotropic to baroclinic domain
    !--------------------------------------------------------------------

    call ocn_spatial_mode_splitting_bottomDepth_remapping(domain)

    !--------------------------------------------------------------------
    ! update the vertical grid on the baroclinic domain
    !--------------------------------------------------------------------

    call ocn_spatial_mode_splitting_update_vertical_grid(domain)

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
    subroutine ocn_spatial_mode_splitting_calc_mesh_intersection(sourceOnSourcePID,&
                                                       sourceOnTargetPID,&
                                                       targetOnTargetPID,&
                                                       overlapPID,&
                                                       sourceOnSourceGID,&
                                                       sourceOnTargetGID,&
                                                       targetOnTargetGID,&
                                                       overlapGID,&
                                                       sourceComm,&
                                                       targetComm,&
                                                       sourceGroup,&
                                                       targetGroup)

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! global id's
    integer,intent(inout) :: sourceOnSourceGID, sourceOnTargetGID,&
                             targetOnTargetGID, overlapGID

    ! moab PID's
    integer,intent(inout) :: sourceOnSourcePID, SourceOnTargetPID,&
                             targetOnTargetPID, overlapPID

    ! communicators and groups
    integer,intent(inout) :: sourceComm, targetComm, sourceGroup, targetGroup

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    ! send source mesh to target domain
    if ( sourceComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_SendMesh(sourceOnSourcePID,&
                           targetComm,&
                           targetGroup,&
                           sourceOnTargetGID,&
                           partScheme)
       call errorout(ierr,'failed to send mesh')
    endif

    ! receive mesh on target
    if ( targetComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_ReceiveMesh(sourceOnTargetPID,&
                              targetComm,&
                              sourceGroup,&
                              sourceOnSourceGID)
       call errorout(ierr,'failed to receive mesh')
    endif

    ! free sender buffers
    if ( sourceComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_FreeSenderBuffers(sourceOnSourcePID,&
                                    sourceOnTargetGID)
       call errorout(ierr,'failed to free sender buffers')
    endif

    ! calculate mesh intersection
    if ( targetComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_ComputeMeshIntersectionOnSphere(sourceOnTargetPID,&
                                                  targetOnTargetPID,&
                                                  overlapPID)
       call errorout(ierr,'failed to compute mesh intersection')
    endif

    ! calculate coverage graph
    if ( sourceComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_CoverageGraph(sourceComm,&
                                sourceOnSourcePID,&
                                sourceOnTargetPID,&
                                overlapPID,&
                                sourceOnSourceGID,&
                                sourceOnTargetGID,&
                                sourceOnTargetGID)
       call errorout(ierr,'failed to compute coverage graph')
    endif

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
                                                       targetComm,&
                                                       weights_identifier,&
                                                       fv_method)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! moab PID's
    integer,intent(inout) :: overlapPID

    ! target communicator
    integer,intent(inout) :: targetComm

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
    if ( targetComm .NE. MPI_COMM_NULL ) then
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
    endif

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
    subroutine ocn_spatial_mode_splitting_define_tags(sourceOnSourcePID,&
                                                      sourceOnTargetPID,&
                                                      targetOnTargetPID,&
                                                      sourceComm,&
                                                      targetComm,&
                                                      fieldTags,&
                                                      fieldTagsProjected,&
                                                      tagTypes,&
                                                      tagIndex)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! mesh PIDs
    integer,intent(inout) :: sourceOnSourcePID,sourceOnTargetPID,&
                             targetOnTargetPID


    ! tag types and indices
    integer,dimension(2),intent(inout) :: tagTypes
    integer,dimension(2),intent(inout) :: tagIndex

    ! mpi communicators
    integer,intent(inout) :: sourceComm,targetComm

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

    ! arrays to set all tags to 0
    real (kind=RKIND), dimension(:), allocatable :: vals

    ! mesh dimensions
    integer, dimension(3) :: nverts, nelem, nblocks, nsbc, ndbc

    ! number of values
    integer :: totalNumVals

    ! data type
    integer :: dtype

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    nDOFS=1 ! 
    dtype=1

    ! define source field tags on source
    if ( sourceComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_DefineTagStorage(sourceOnSourcePID,&
                                   fieldTags,&
                                   tagTypes(1),&
                                   nDOFS,&
                                   tagIndex(1))
       call errorout(ierr,'failed to define tags on sourceOnSource')

       ! set tags to 0
       ierr=iMOAB_GetMeshInfo(sourceOnSourcePID,nverts,nelem,nblocks,nsbc,ndbc)
       call errorout(ierr,'failed to retrieve mesh info (source on source)')
    endif

    ! define tags on target "side"
    if ( targetComm .NE. MPI_COMM_NULL ) then

      ! source fields on source
      ierr=iMOAB_DefineTagStorage(sourceOnTargetPID,&
                                  fieldTags,&
                                  tagTypes(1),&
                                  nDOFS,&
                                  tagIndex(1))
       call errorout(ierr,'failed to define tags on sourceOnTarget')

       ! set tags to 0
       ierr=iMOAB_GetMeshInfo(sourceOnTargetPID,nverts,nelem,nblocks,nsbc,ndbc)
       call errorout(ierr,'failed to retrieve mesh info (source on source)')

       ! target fields
       ierr=iMOAB_DefineTagStorage(targetOnTargetPID,&
                                   fieldTagsProjected,&
                                   tagTypes(2),&
                                   nDOFS,&
                                   tagIndex(2))
       call errorout(ierr,'failed to define tags on targetOnTarget')

       ! set tags to 0
       ierr=iMOAB_GetMeshInfo(targetOnTargetPID,nverts,nelem,nblocks,nsbc,ndbc)
       call errorout(ierr,'failed to retrieve mesh info (source on source)')

    endif

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
    subroutine ocn_spatial_mode_splitting_allocate_buffers(sourceOnSourcePID,&
                                                           sourceOnTargetPID,&
                                                           sourceComm,&
                                                           targetComm,&
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
    integer,intent(inout) :: sourceOnSourcePID,sourceOnTargetPID

    ! mpi communicators
    integer,intent(inout) :: sourceComm,targetComm

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

    ! on source
    if ( sourceComm .NE. MPI_COMM_NULL ) then

       ! get mesh dimensions from moab
       ierr=iMOAB_GetMeshInfo(sourceOnSourcePID,nverts,nelem,nblocks,nsbc,ndbc)
       call errorout(ierr,'failed to retrieve mesh info (source on source)')

       ! initialize buffers
       if ( .NOT. allocated( source2target ) ) then
          allocate ( source2target(nTags*nelem(3)) )
       endif
       if ( .NOT. allocated( target2source ) ) then
          allocate ( target2source(nTags*nelem(3)) )
       endif
    endif

    ! on target
    if ( targetComm .NE. MPI_COMM_NULL ) then

       ! get mesh dimensions from moab
       ierr=iMOAB_GetMeshInfo(sourceOntargetPID,nverts,nelem,nblocks,nsbc,ndbc)
       call errorout(ierr,'failed to retrieve mesh info (source on source)')

       ! initialize buffers
       if ( .NOT. allocated( source2target ) ) then
          allocate ( source2target(nTags*nelem(3)) )
       endif
       if ( .NOT. allocated( target2source ) ) then
          allocate ( target2source(nTags*nelem(3)) )
       endif
    endif

    source2target(:)=0d0
    target2source(:)=0d0

#endif

    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_allocate_buffers!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_bottomDepth_remap
!
!> \brief   Remaps bottomDepth from barotropic to baroclinic domain
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!>  This routine conservatively remaps bottomDepth from the barotropic 
!>  to baroclinic domain.
!
!-----------------------------------------------------------------------


    subroutine ocn_spatial_mode_splitting_bottomDepth_remapping(domain)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------
    type(domain_type), intent(inout) :: domain ! will alwas be barotropic

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------
    ! domain pointer
    type (domain_type), pointer :: domain_bc,domain_btr

    ! configuration info
    logical, pointer :: config_spatial_mode_splitting

    ! pool pointers
    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool

    ! number of cells on each mesh
    integer, pointer :: nCellsSolve_btr, nCellsSolve_bc

    ! arrays for storing data
    real (kind=RKIND), dimension(:), pointer :: bottomDepth_bc,&
                                                bottomDepth_btr

    ! other parameters moab needs
    integer :: totalNumVals, entity_type

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

    !--------------------------------------------------------------------
    ! add bottomDepth to MOAB
    !--------------------------------------------------------------------

    ! set filter type (no filter)
    filter_type=0

    ! create weights for remapping bottomDepth
    bottomDepth_weights_identifier='bottomDepthWeights'//C_NULL_CHAR
    fv_method=''//C_NULL_CHAR
    call ocn_spatial_mode_splitting_calc_weights(btrXbcOnBcPID,&
                                                 bcComm,&
                                                 bottomDepth_weights_identifier,&
                                                 fv_method)

    ! define tags for bottomDepth
    bottomDepthTag='bottomDepth_bc'//C_NULL_CHAR
    bottomDepthTag_btr='bottomDepth_btr'//C_NULL_CHAR
    bottomDepthType(:)=1
    call ocn_spatial_mode_splitting_define_tags(btrOnBtrPID,&
                                                btrOnBcPID,&
                                                bcOnBCPID,&
                                                btrComm,&
                                                bcComm,&
                                                bottomDepthTag_btr,&
                                                bottomDepthTag,&
                                                bottomDepthType,&
                                                bottomDepthIndex)


    !--------------------------------------------------------------------
    ! remap bottomDepth from barotropic to baroclinic domains
    !--------------------------------------------------------------------

    ! get bottomDepth array from baroclinic domain
    domain_bc => domain % core % domainlist
    block_ptr => domain_bc % blocklist
    do while(associated(block_ptr))

       ! setup pools
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       ! get nCells
       call mpas_pool_get_dimension(meshPool,'nCellsSolve',nCellsSolve_bc)

       ! get bottomDepth from baroclinic domain
       call mpas_pool_get_array(meshPool,'bottomDepth',bottomDepth_bc)

       block_ptr => block_ptr % next
    end do 

    ! get bottomDepth array on barotropic domain
    domain_btr => domain_bc % next
    block_ptr => domain_btr % blocklist
    do while(associated(block_ptr))

       ! setup pools
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       ! get nCells
       call mpas_pool_get_dimension(meshPool,'nCellsSolve',nCellsSolve_btr)

       ! get bottomDepth from baroclinic domain
       call mpas_pool_get_array(meshPool,'bottomDepth',bottomDepth_btr)

       block_ptr => block_ptr % next
    end do 

    ! remap bathymetry
    totalNumVals=nCellsSolve_btr
    entity_type=1
    call ocn_spatial_mode_splitting_remap_data(btrOnBtrPID,&
                                               btrOnBcPID,&
                                               bcOnBcPID,&
                                               btrXbcOnBcPID,&
                                               btrOnBtr,&
                                               btrOnBc,&
                                               bcOnBc,&
                                               btrXbcOnBc,&
                                               btrComm,&
                                               bcComm,&
                                               btrGroup,&
                                               bcGroup,&
                                               bottomDepthTag_btr,&
                                               bottomDepthTag,&
                                               totalNumVals,&
                                               entity_type,&
                                               bottomDepth_btr,&
                                               bottomDepth_weights_identifier)

    ! store remapped bathymetry in bottomDepth_bc
    bottomDepth_bc(:)=-9999d0
    totalNumVals=nCellsSolve_bc
    entity_type=1
    ierr = iMOAB_GetDoubleTagStorage(bcOnBcPID,&
                                     bottomDepthTag,&
                                     totalNumVals,&
                                     entity_type,&
                                     bottomDepth_bc)
    call errorout(ierr,'failed to get bottomDepth_bc')

    ! halo exchange new bathy
    call mpas_dmpar_field_halo_exch(domain_bc,'bottomDepth')

#endif
    end subroutine ocn_spatial_mode_splitting_bottomDepth_remapping!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_update_vertical_grid
!
!> \brief   Updates vertical grid following remapping of bathymetry
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!>  This routine updates the vertical grid following conservative remap-
!>  ping of bathymetry
!
!-----------------------------------------------------------------------


    subroutine ocn_spatial_mode_splitting_update_vertical_grid(domain)!{{{

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
    ! domain pointer
    type (domain_type), pointer :: domain_bc

    ! block pointer
    type (block_type), pointer :: block_ptr

    ! pools
    type(mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool

    ! dimensions
    integer, pointer :: nCellsSolve_bc, nEdgesSolve_bc

    ! variable pointers from verticalMesh pool
    real (kind=RKIND), dimension(:,:), pointer :: restingThickness

    ! variable pointers from state pool
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness_bc
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocity_bc
    real (kind=RKIND), dimension(:), pointer :: ssh_bc

    ! for resetting the 

    ! local work variables
    real (kind=RKIND) :: maxBottomDepth, minBottomDepth, wc_thickness, delta_z
    integer :: idz

    ! for loops
    integer :: iCell, iEdge, iVertex, iVertLevel, kk

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    ! ensure we are on the baroclinic domain
    domain_bc => domain % core % domainlist
    block_ptr => domain_bc % blocklist
    do while(associated(block_ptr))

      ! get pools
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_bc)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_bc)

      ! get verticalMesh pool variables
      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
      
      ! get state pool variables
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness_bc, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity_bc, 1)
      call mpas_pool_get_array(statePool, 'ssh', ssh_bc, 1)

      ! get new maximum bottomDepth and make sure it matches refBottomDepth. By this I mean that we
      ! need refBottomDepth(n-1) to be greater than or equal to the maximum bottom depth
      maxBottomDepth=MAXVAL(bottomDepth)
      if (maxBottomDepth .GT. refBottomDepth(nVertLevels)) then
         call mpas_log_write('new max(bottomDepth) > refBotomDepth(nVertLevels)',MPAS_LOG_CRIT)
      endif

      ! loop over each cell and update as needed
      do iCell = 1,nCellsSolve_bc

         ! TODO: figure out what to do when we aren't starting from a cold start (i.e.,
         ! sum(layerThickness) should equal bottomDepth+ssh

         ! find the sum of the restingThicknesses
         wc_thickness = restingThickness(minLevelCell(iCell),iCell)
         do iVertLevel = minLevelCell(iCell)+1,maxLevelCell(iCell)
            wc_thickness = wc_thickness + restingThickness(iVertLevel,iCell)
         end do

         ! check if the sum of the layer thicknesses is consistent with the bottom depth
         if ( wc_thickness .EQ. bottomDepth(iCell) ) then
            cycle
         end if

         ! if we are inconsistent, then we need to adjust. First, calculate the max and min bottom
         ! depths we allow for this layer (+/- 50% of refZMid)

         ! find which level our bottom depth resides in
         ! NOTE: because I do this, we will only ever need to remove a layer. Adding a layer will
         ! never be required
         do iVertLevel = 1,nVertLevels
            if ( refBottomDepth(iVertLevel) .GE. bottomDepth(iCell) ) then
               idz = iVertLevel
               exit
            end if
         end do

         ! find min depth (refBottomDepth - 50% of layer thickness)
         minBottomDepth = 0.5d0 * ( refBottomDepth(idz-1) + refBottomDepth(idz) )

         ! if bottomDepth <= minBottomDepth, subtract a layer
         if ( bottomDepth(iCell) .LT. minBottomDepth ) then
            idz = idz - 1
         end if

         maxLevelCell(iCell)=idz

         ! find the sum of the restingThicknesses
         wc_thickness = restingThickness(minLevelCell(iCell),iCell)
         do iVertLevel = minLevelCell(iCell)+1,maxLevelCell(iCell)
            wc_thickness = wc_thickness + restingThickness(iVertLevel,iCell)
         end do

         delta_z = wc_thickness - bottomDepth(iCell)
         restingThickness(idz,iCell) = restingThickness(idz,iCell) - delta_z

         ! TODO: figure out what to do when we aren't starting from a cold start (i.e.,
         ! sum(layerThickness) should equal bottomDepth+ssh
         do iVertLevel = 1,nVertLevels
            layerThickness_bc(iVertLevel,iCell) = restingThickness(iVertLevel,iCell)
            normalVelocity_bc(iVertLevel,iCell) = 0d0
         end do

      end do

      ! perform halo exchange for maxLevelCell
      call mpas_dmpar_field_halo_exch(domain_bc,'maxLevelCell')
      call mpas_dmpar_field_halo_exch(domain_bc,'restingThickness')
      call mpas_dmpar_field_halo_exch(domain_bc,'layerThickness',timeLevel=1)
      call mpas_dmpar_field_halo_exch(domain_bc,'normalVelocity',timeLevel=1)

      ! call meshDestroy (will be recreated in the init routines)
      call ocn_meshDestroy(ierr)

      block_ptr => block_ptr % next
    end do

#endif
    end subroutine ocn_spatial_mode_splitting_update_vertical_grid!}}}


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
    real (kind=RKIND), dimension(:), pointer :: barotropicForcing_bc,&
                                                barotropicForcing_btr,&
                                                barotropicForcingZonal_bc,&
                                                barotropicForcingZonal_btr,&
                                                barotropicForcingMeridional_bc,&
                                                barotropicForcingMeridional_btr

    ! for going from cell centered to edge based quantities
    integer :: cell1, cell2
    real (kind=RKIND) :: zonal_edge, meridional_edge
    real (kind=RKIND), dimension(:), pointer :: angleEdge_btr
    integer, dimension(:,:), pointer :: cellsOnEdge_btr,&
                                        verticesOnEdge_btr,&
                                        cellsOnVertex_btr
    real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex_btr

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
    ! get barotropicForcing from baroclinic domain
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))

        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_bc)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_bc)

        ! get barotropicForcing from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingZonal',&
                                 barotropicForcingZonal_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingMeridional',&
                                 barotropicForcingMeridional_bc)

        ! exchange halos for barotropicForcing
        ! TODO: test if this is actually necessary
        call mpas_dmpar_field_halo_exch(domain,'barotropicForcing')

        ! reconstruct from edge based to cell-centered quantities
        call mpas_reconstruct(meshPool, barotropicForcing_bc, &
                              barotropicForcingX_bc, &
                              barotropicForcingY_bc, &
                              barotropicForcingZ_bc, &
                              barotropicForcingZonal_bc, &
                              barotropicForcingMeridional_bc)

        ! exchange halos for barotropicForcing
        ! TODO: test if this is actually necessary
         call mpas_dmpar_field_halo_exch(domain,'barotropicForcingMeridional')
         call mpas_dmpar_field_halo_exch(domain,'barotropicForcingZonal')

        ! pack buffers for iMOAB
        do iCell = 1,nCellsSolve_bc
            bc2btrVar(iCell) = barotropicForcingZonal_bc(iCell)
            bc2btrVar(iCell+nCellsSolve_bc) = barotropicForcingMeridional_bc(iCell)
        end do

        ! perform remapping
        totalNumVals=2*nCellsSolve_bc
        entity_type=1
        call ocn_spatial_mode_splitting_remap_data(bcOnBcPID,&
                                                   bcOnBtrPID,&
                                                   btrOnBtrPID,&
                                                   bcXbtrOnBtrPID,&
                                                   bcOnBc,&
                                                   bcOnBtr,&
                                                   btrOnBtr,&
                                                   bcXbtrOnBtr,&
                                                   bcComm,&
                                                   btrComm,&
                                                   bcGroup,&
                                                   btrGroup,&
                                                   bc2btrFields,&
                                                   bc2btrFields_btr,&
                                                   totalNumVals,&
                                                   entity_type,&
                                                   bc2btrVar,&
                                                   bc2btr_weights_identifier)

        block_ptr => block_ptr % next
    end do

    !--------------------------------------------------------------------
    ! remap barotropicForcing to barotropic domain
    !--------------------------------------------------------------------

    block_ptr => domain_btr % blocklist
    do while(associated(block_ptr))

        ! setup pools
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree_btr)

        ! get connectivity/etc.
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge_btr)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge_btr)
        call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge_btr)
        call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex_btr)
        call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex_btr)

        ! get barotropicForcing from pools
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingZonal',&
                                 barotropicForcingZonal_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingMeridional',&
                                 barotropicForcingMeridional_btr)

        ! get remapped data from iMOAB
        totalNumVals=2*nCellsSolve_btr
        entity_type=1
        ierr = iMOAB_GetDoubleTagStorage(btrOnBtrPID,&
                                         bc2btrFields_btr,&
                                         totalNumVals,&
                                         entity_type,&
                                         bc2btrVar_btr)
        call errorout(ierr,"failed to get bc2btrVar_btr")

        ! unpack buffers
        do iCell = 1,nCellsSolve_btr
           barotropicForcingZonal_btr(iCell) = bc2btrVar_btr(iCell)
           barotropicForcingMeridional_btr(iCell) = bc2btrVar_btr(iCell+nCellsSolve_btr)
        end do

        ! perform halo exchange
        call mpas_dmpar_field_halo_exch(domain_btr, 'barotropicForcingZonal')
        call mpas_dmpar_field_halo_exch(domain_btr, 'barotropicForcingMeridional')

        ! put back into edge based arrays
        barotropicForcing_btr(:) = 0d0
        do iEdge = 1,nEdgesSolve_btr
            cell1 = cellsOnEdge_btr(1, iEdge)
            cell2 = cellsOnEdge_btr(2, iEdge)
            zonal_edge = 0.5d0 * (barotropicForcingZonal_btr(cell1) + &
                                  barotropicForcingZonal_btr(cell2))
            meridional_edge = 0.5d0 * (barotropicForcingMeridional_btr(cell1) + &
                                       barotropicForcingMeridional_btr(cell2))
            barotropicForcing_btr(iedge) = barotropicForcing_btr(iEdge) + &
                                           cos(angleEdge_btr(iEdge))*zonal_edge + &
                                           sin(angleEdge_btr(iEdge))*meridional_edge
        end do

        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesAll_btr)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCellsAll_btr)
        do iEdge = 1,nEdgesAll_btr
           cell1=cellsOnEdge_btr(1,iEdge)
           cell2=cellsOnEdge_btr(2,iEdge)
           if (cell1.lt.1.OR.cell1.gt.nCellsAll_btr) then
              barotropicForcing_btr(iEdge)=0d0
           elseif (cell2.lt.1.OR.cell2.gt.nCellsAll_btr) then
              barotropicForcing_btr(iEdge)=0d0
           endif
        end do
        
        call mpas_dmpar_field_halo_exch(domain_btr, 'barotropicForcing')

        block_ptr => block_ptr % next
    end do
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

    ! mesh information
    integer, pointer :: nCellsSolve_bc, nCellsSolve_btr,&
                        nEdgesSolve_bc, nEdgesSolve_btr,&
                        vertexDegree_bc, nEdgesAll_bc,&
                        nCellsAll_bc

    ! size info for moab
    integer :: totalNumVals, entity_type

    ! arrays for storing data
    real (kind=RKIND), dimension(:), pointer :: barotropicThicknessFlux_bc,&
                                                barotropicThicknessFlux_btr,&
                                                barotropicThicknessFluxZonal_bc,&
                                                barotropicThicknessFluxZonal_btr,&
                                                barotropicThicknessFluxMeridional_bc,&
                                                barotropicThicknessFluxMeridional_btr

    ! for going from cell centered to edge based quantities
    integer :: cell1, cell2
    real (kind=RKIND) :: zonal_edge, meridional_edge
    real (kind=RKIND), dimension(:), pointer :: angleEdge_bc
    integer, dimension(:,:), pointer :: cellsOnEdge_bc,&
                                        verticesOnEdge_bc,&
                                        cellsOnVertex_bc
    real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex_bc

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

        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                                 barotropicThicknessFlux_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxZonal',&
                                 barotropicThicknessFluxZonal_btr)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxMeridional',&
                                 barotropicThicknessFluxMeridional_btr)

        ! exchange halos for barotropicThicknessFlux
        ! TODO: test if this is actually necessary
         call mpas_dmpar_field_halo_exch(domain_btr,'barotropicThicknessFlux')

        ! reconstruct from edge based to cell-centered quantities
        call mpas_reconstruct(meshPool, barotropicThicknessFlux_btr,&
                              barotropicThicknessFluxX_btr,&
                              barotropicThicknessFluxY_btr,&
                              barotropicThicknessFluxZ_btr,&
                              barotropicThicknessFluxZonal_btr,&
                              barotropicThicknessFluxMeridional_btr)
         call mpas_dmpar_field_halo_exch(domain_btr,'barotropicThicknessFluxZonal')
         call mpas_dmpar_field_halo_exch(domain_btr,'barotropicThicknessFluxMeridional')

        ! pack buffers
        do iCell = 1,nCellsSolve_btr
            btr2bcVar_btr(iCell) = barotropicThicknessFluxZonal_btr(iCell)
            btr2bcVar_btr(iCell+nCellsSolve_btr) = barotropicThicknessFluxMeridional_btr(iCell)
        end do

        ! perform remapping
        totalNumVals=2*nCellsSolve_btr
        entity_type=1
        call ocn_spatial_mode_splitting_remap_data(btrOnBtrPID,&
                                                   btrOnBcPID,&
                                                   bcOnBcPID,&
                                                   btrXbcOnBcPID,&
                                                   btrOnBtr,&
                                                   btrOnBc,&
                                                   bcOnBc,&
                                                   btrXbcOnBc,&
                                                   btrComm,&
                                                   bcComm,&
                                                   btrGroup,&
                                                   bcGroup,&
                                                   btr2bcFields_btr,&
                                                   btr2bcFields,&
                                                   totalNumVals,&
                                                   entity_type,&
                                                   btr2bcVar_btr,&
                                                   btr2bc_weights_identifier)

        block_ptr => block_ptr % next
    end do

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
        call mpas_pool_get_dimension(meshpool, 'vertexDegree', vertexDegree_bc)

        ! get connectivity/etc.
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge_bc)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge_bc)
        call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge_bc)
        call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex_bc)
        call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex_bc)

        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                                 barotropicThicknessFlux_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxZonal',&
                                 barotropicThicknessFluxZonal_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxMeridional',&
                                 barotropicThicknessFluxMeridional_bc)

        ! get remapped data from iMOAB
        totalNumVals=2*nCellsSolve_bc
        entity_type=1
        ierr = iMOAB_GetDoubleTagStorage(bcOnBcPID,&
                                         btr2bcFields,&
                                         totalNumVals,&
                                         entity_type,&
                                         btr2bcVar)
        call errorout(ierr,"failed to get btr2bcVar from iMOAB")

        ! unpack buffers
        do iCell = 1,nCellsSolve_bc
            barotropicThicknessFluxZonal_bc(iCell) = btr2bcVar(iCell)
            barotropicThicknessFluxMeridional_bc(iCell) = btr2bcVar(iCell+nCellsSolve_bc)
        end do

        ! exchange halos
        call mpas_dmpar_field_halo_exch(domain, 'barotropicThicknessFluxZonal')
        call mpas_dmpar_field_halo_exch(domain, 'barotropicThicknessFluxMeridional')

        ! put into edge-based arrays
        barotropicThicknessFlux_bc(:)=0d0
        do iEdge = 1,nEdgesSolve_bc
           cell1 = cellsOnEdge_bc(1, iEdge)
           cell2 = cellsOnEdge_bc(2,iEdge)
           zonal_edge = 0.5d0 * (barotropicThicknessFluxZonal_bc(cell1) + &
                                 barotropicThicknessFluxZonal_bc(cell2))
           meridional_edge = 0.5d0 * (barotropicThicknessFluxMeridional_bc(cell1) + &
                                      barotropicThicknessFluxMeridional_bc(cell2))
           barotropicThicknessFlux_bc(iEdge) = barotropicThicknessFlux_bc(iEdge) + &
                                            cos(angleEdge_bc(iEdge))*zonal_edge + &
                                            sin(angleEdge_bc(iEdge))*meridional_edge
        end do

        call mpas_dmpar_field_halo_exch(domain, 'barotropicThicknessFlux')

        block_ptr => block_ptr % next
    end do
#endif
    end subroutine ocn_spatial_mode_splitting_btr2bc_remapping!}}}

!***********************************************************************
!
!  routine ocn_spatial_mode_splitting_remap_data
!
!> \brief   
!> \author  Coleman Blakely
!> \date    April 2024
!> \details
!
!-----------------------------------------------------------------------
    subroutine ocn_spatial_mode_splitting_remap_data(sourceOnSourcePID,&
                                                     sourceOnTargetPID,&
                                                     targetOnTargetPID,&
                                                     overlapPID,&
                                                     sourceOnSourceGID,&
                                                     sourceOnTargetGID,&
                                                     targetOnTargetGID,&
                                                     overlapGID,&
                                                     sourceComm,&
                                                     targetComm,&
                                                     sourceGroup,&
                                                     targetGroup,&
                                                     fieldTags,&
                                                     fieldTagsProjected,&
                                                     totalNumVals,&
                                                     entity_type,&
                                                     buffer,&
                                                     weights_identifier)!{{{

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    ! moab pids
    integer,intent(inout) :: sourceOnSourcePID,sourceOnTargetPID,&
                          targetOnTargetPID,&
                          overlapPID

    ! moab global id's
    integer,intent(inout) :: sourceOnSourceGID,sourceOnTargetGID,&
                          targetOnTargetGID,&
                          overlapGID

    ! mpi comms/groups
    integer,intent(inout) :: sourceComm,targetComm
    integer,intent(inout) :: sourceGroup,targetGroup

    ! names of fields
    character(:),allocatable,intent(inout) :: fieldTags,fieldTagsProjected

    ! data for moab
    integer,intent(inout) :: totalNumVals ! number of values being sent, e.g., 
                                       ! 2*nCells for barotropicForcing
                                       ! (zonal and meridional components)
    integer,intent(inout) :: entity_type  ! 0 for vertex, 1 for element, (2 for edge?)

    ! data being projected
    real (kind=RKIND),dimension(:),intent(inout) :: buffer

    ! what weights we are using
    character(:),allocatable,intent(inout) :: weights_identifier

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    ! imoab filter type (we use no filter)
    integer :: filter_type

    ! error handling
    integer :: ierr

    ! filename
    character(:),allocatable :: filename

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    ! set filter type
    filter_type=0

    ! store values in moab on source side and send to target side
    if ( sourceComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_SetDoubleTagStorage(sourceOnSourcePID,&
                                      fieldTags,&
                                      totalNumVals,&
                                      entity_type,&
                                      buffer)
       call errorout(ierr,'unable to set storage on source')

       ! send to target "side"
       ierr=iMOAB_SendElementTag(sourceOnSourcePID,&
                                 fieldTags,&
                                 targetComm,&
                                 sourceOnTargetGID)
       call errorout(ierr,'unable to send data from source "side" to target "side"')
    endif

    ! receive values on target side
    if ( targetComm .NE. MPI_COMM_NULL ) then

       ierr=iMOAB_ReceiveElementTag(sourceOnTargetPID,&
                                    fieldTags,&
                                    targetComm,&
                                    sourceOnTargetGID)
       call errorout(ierr,'unable to receive data on target "side"')
    endif

    ! free send buffers
    if ( sourceComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_FreeSenderBuffers(sourceOnSourcePID,sourceOnTargetGID)
       call errorout(ierr,"unable to free sender buffers")
    endif

    ! perform remapping
    if ( targetComm .NE. MPI_COMM_NULL ) then
       ierr=iMOAB_ApplyScalarProjectionWeights(overlapPID,&
                                               filter_type,&
                                               weights_identifier,&
                                               fieldTags,&
                                               fieldTagsProjected)
       call errorout(ierr,'unable to reproject :(')
    endif

#endif

    !--------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_remap_data!}}}

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
                                                       sshTend) !{{{

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
                                                       
    ! tendPool
    real (kind=RKIND), pointer, intent(inout), dimension(:) :: sshTend

    !-----------------------------------------------------------------
    ! output variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------

    type(block_type), pointer :: block  ! variables in current subblock

    type(mpas_pool_type), pointer :: diagnosticsPool

    ! loop indices
    integer :: i, k, n
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

    call mpas_pool_get_array(tendPool, 'ssh', sshTend)

    ! redirect barotropicForcing and barotropicThicknessFlux
    call mpas_pool_get_subpool(block%structs, 'diagnostics', diagnosticsPool)

    call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

    call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)


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
        end do

        block_ptr => block_ptr % next
    end do
!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_time_int_split_setup_btr !}}}

end module ocn_spatial_mode_splitting
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
