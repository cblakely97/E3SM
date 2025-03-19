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
              ocn_spatial_mode_splitting_bottomDepth_remapping,&
              ocn_spatial_mode_splitting_dmpar_init_bc,&
              ocn_spatial_mode_splitting_dmpar_init_btr

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

   ! The following variables are for conservatively remapping the bathymetry from barotropic to
   ! baroclinic domain.
   character(:), allocatable :: bottomDepthTag,bottomDepthTag_btr
   integer, dimension(2) :: bottomDepthIndex, bottomDepthType
   character(:), allocatable :: bottomDepth_weights_identifier

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

       ! as a stopgap measure we are going to use scatter/gather
       ! initialize
!       write(6,*) "about to init dmpar for baroclinic"
!       weights_filename_btr2bc='16km_to_32km_weights.nc'
!       call ocn_spatial_mode_splitting_dmpar_remapping_init(domain,&
!                                                            weights_filename_btr2bc,&
!                                                            rowValues_bc2btr,&
!                                                            colValues_bc2btr,&
!                                                            sValues_bc2btr,&
!                                                            nEdgesGlobal_bc,&
!                                                            nEdgesDisplacement_bc,&
!                                                            indexToEdgeIDGathered_bc,&
!                                                            nEdgesPerProc_bc)
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

       ! as a stopgap measure we are going to use scatter/gather
       ! initialize
!       weights_filename_btr2bc='32km_to_!6km_weights.nc'
!       call ocn_spatial_mode_splitting_dmpar_remapping_init(domain,&
!                                                            weights_filename_btr2bc,&
!                                                            rowValues_btr2bc,&
!                                                            colValues_btr2bc,&
!                                                            sValues_btr2bc,&
!                                                            nEdgesGlobal_btr,&
!                                                            nEdgesDisplacement_btr,&
!                                                            indexToEdgeIDGathered_btr,&
!                                                            nEdgesPerProc_btr)
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

    ! buffer for sending the data
    real (kind=RKIND), dimension(:), allocatable :: buffer_btr2bc, buffer_bc2btr
    real (kind=RKIND), dimension(:), allocatable :: buffer_btr2bc_btr, buffer_bc2btr_btr

    ! other parameters moab needs
    integer :: totalNumVals, entity_type

    ! for loops
    integer :: iCell

    ! error handling
    integer :: ierr

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

#ifdef HAVE_MOAB

    !--------------------------------------------------------------------
    ! add bottomDepth to MOAB
    !--------------------------------------------------------------------

    ! set filter type (no filter)
    filter_type=0

    ! create weights for remapping bottomDepth
    bottomDepth_weights_identifier='bottomDepthWeights'//C_NULL_CHAR
    fv_method=''//C_NULL_CHAR
    call ocn_spatial_mode_splitting_calc_weights(btrXbcPID,&
                                                 bottomDepth_weights_identifier,&
                                                 fv_method)

    ! define tags for bottomDepth
    bottomDepthTag='bottomDepth_bc'//C_NULL_CHAR
    bottomDepthTag_btr='bottomDepth_btr'//C_NULL_CHAR
    bottomDepthType(:)=1
    call ocn_spatial_mode_splitting_define_tags(btrPID,&
                                                bcPID,&
                                                bottomDepthTag_btr,&
                                                bottomDepthTag,&
                                                bottomDepthType,&
                                                bottomDepthIndex)

    ! allocate our buffers
    call ocn_spatial_mode_splitting_allocate_buffers(btrPID,&
                                                     buffer_btr2bc_btr,&
                                                     buffer_bc2btr_btr,&
                                                     1)
    call ocn_spatial_mode_splitting_allocate_buffers(bcPID,&
                                                     buffer_btr2bc,&
                                                     buffer_bc2btr,&
                                                     1)

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
    buffer_btr2bc_btr(:)=0d0
    do iCell = 1,nCellsSolve_btr
      buffer_btr2bc_btr(iCell)=bottomDepth_btr(iCell)
    enddo

    totalNumVals=nCellsSolve_btr
    entity_type=1
    call ocn_spatial_mode_splitting_remap_data(btrPID,&
                                               bcPID,&
                                               btrXbcPID,&
                                               btrGID,&
                                               bcGID,&
                                               btrXbcGID,&
                                               bottomDepthTag_btr,&
                                               bottomDepthTag,&
                                               totalNumVals,&
                                               entity_type,&
                                               buffer_btr2bc_btr,&
                                               bottomDepth_weights_identifier)

    ! store remapped bathymetry in bottomDepth_bc
    bottomDepth_bc(:)=-9999d0
    totalNumVals=nCellsSolve_bc
    entity_type=1
    ierr = iMOAB_GetDoubleTagStorage(bcPID,&
                                     bottomDepthTag,&
                                     totalNumVals,&
                                     entity_type,&
                                     buffer_btr2bc)
    call errorout(ierr,'failed to get bottomDepth_bc')

    do iCell = 1,nCellsSolve_bc
       bottomDepth_bc(iCell)=buffer_btr2bc(iCell)
    enddo

    ! halo exchange new bathy
    call mpas_dmpar_field_halo_exch(domain_bc,'bottomDepth')

    deallocate(buffer_btr2bc_btr,buffer_bc2btr_btr,buffer_btr2bc,buffer_bc2btr)

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
    integer, pointer :: nVertLevels_bc

    ! variable pointers from mesh pool
    integer, dimension(:), pointer :: minLevelCell_bc, maxLevelCell_bc
    real (kind=RKIND), dimension(:), pointer :: bottomDepth_bc
    real (kind=RKIND), dimension(:), pointer :: refBottomDepth_bc

    ! variable pointers from verticalMesh pool
    real (kind=RKIND), dimension(:,:), pointer :: restingThickness_bc

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
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels_bc)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_bc)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve_bc)

      ! get mesh pool variables
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell_bc)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell_bc)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth_bc)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth_bc)

      ! get verticalMesh pool variables
      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness_bc)
      
      ! get state pool variables
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness_bc, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity_bc, 1)
      call mpas_pool_get_array(statePool, 'ssh', ssh_bc, 1)

      ! get new maximum bottomDepth and make sure it matches refBottomDepth. By this I mean that we
      ! need refBottomDepth_bc(n-1) to be greater than or equal to the maximum bottom depth
      maxBottomDepth=MAXVAL(bottomDepth_bc)
      if (maxBottomDepth .GT. refBottomDepth_bc(nVertLevels_bc)) then
         call mpas_log_write('new max(bottomDepth) > refBotomDepth(nVertLevels)',MPAS_LOG_CRIT)
      endif

      ! loop over each cell and update as needed
      do iCell = 1,nCellsSolve_bc

         ! TODO: figure out what to do when we aren't starting from a cold start (i.e.,
         ! sum(layerThickness) should equal bottomDepth+ssh

         ! find the sum of the restingThickness_bces
         wc_thickness = restingThickness_bc(minLevelCell_bc(iCell),iCell)
         do iVertLevel = minLevelCell_bc(iCell)+1,maxLevelCell_bc(iCell)
            wc_thickness = wc_thickness + restingThickness_bc(iVertLevel,iCell)
         end do

         ! check if the sum of the layer thicknesses is consistent with the bottom depth
         if ( wc_thickness .eq. bottomDepth_bc(iCell) ) then
            cycle
         end if

         ! if we are inconsistent, then we need to adjust. First, calculate the max and min bottom
         ! depths we allow for this layer (+/- 50% of refZMid)

         ! find which level our bottom depth resides in
         ! NOTE: because I do this, we will only ever need to remove a layer. Adding a layer will
         ! never be required
         do iVertLevel = 1,nVertLevels_bc
            if ( refBottomDepth_bc(iVertLevel) .GE. bottomDepth_bc(iCell) ) then
               idz = iVertLevel
               exit
            end if
         end do

         ! find min depth (refBottomDepth_bc - 50% of layer thickness)
         if (idz.eq.1) then
            call mpas_log_write('we should fail here',MPAS_LOG_CRIT)
         endif
         minBottomDepth = 0.5d0 * ( refBottomDepth_bc(idz-1) + refBottomDepth_bc(idz) )

         ! if bottomDepth_bc <= minBottomDepth, subtract a layer
         if ( bottomDepth_bc(iCell) .LT. minBottomDepth ) then
            idz = idz - 1
         end if

         maxLevelCell_bc(iCell)=idz

         ! find the sum of the restingThickness_bces
         wc_thickness = restingThickness_bc(minLevelCell_bc(iCell),iCell)
         do iVertLevel = minLevelCell_bc(iCell)+1,maxLevelCell_bc(iCell)
            wc_thickness = wc_thickness + restingThickness_bc(iVertLevel,iCell)
         end do

         delta_z = wc_thickness - bottomDepth_bc(iCell)
         restingThickness_bc(idz,iCell) = restingThickness_bc(idz,iCell) - delta_z

         ! TODO: figure out what to do when we aren't starting from a cold start (i.e.,
         ! sum(layerThickness) should equal bottomDepth_bc+ssh
         do iVertLevel = 1,nVertLevels_bc
            layerThickness_bc(iVertLevel,iCell) = restingThickness_bc(iVertLevel,iCell)
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
        call ocn_spatial_mode_splitting_remap_data(bcPID,&
                                                   btrPID,&
                                                   bcXbtrPID,&
                                                   bcGID,&
                                                   btrGID,&
                                                   bcXbtrGID,&
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
        ierr = iMOAB_GetDoubleTagStorage(btrPID,&
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
                        vertexDegree_bc
    integer, pointer :: nEdgesAll_bc, nEdgesAll_btr

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
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesAll_btr)

        ! exchange halos for barotropicThicknessFlux

        ! TODO: test if this is actually necessary
        call mpas_dmpar_field_halo_exch(domain_btr,'barotropicThicknessFlux')

        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                                 barotropicThicknessFlux_btr)

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

        ! get barotropicThicknessFlux from 
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesAll_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                                 barotropicThicknessFlux_bc)

        ! interpolate using dmpar routine for now
        call ocn_spatial_mode_splitting_dmpar_interpolate(domain,&
                                                          nEdgesAll_btr,&
                                                          nEdgesGlobal_btr,&
                                                          nEdgesPerProc_btr,&
                                                          nEdgesDisplacement_btr,&
                                                          indexToEdgeIDGathered_btr,&
                                                          barotropicThicknessFlux_btr,&
                                                          nEdgesAll_btr,&
                                                          nEdgesGlobal_bc,&
                                                          nEdgesPerProc_bc,&
                                                          nEdgesDisplacement_bc,&
                                                          indexToEdgeIDGathered_bc,&
                                                          barotropicThicknessFlux_bc,&
                                                          colValues_btr2bc,&
                                                          rowValues_btr2bc,&
                                                          sValues_btr2bc)

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
    subroutine ocn_spatial_mode_splitting_remap_data(sourcePID,&
                                                     targetPID,&
                                                     overlapPID,&
                                                     sourceGID,&
                                                     targetGID,&
                                                     overlapGID,&
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
    integer,intent(inout) :: sourcePID,&
                             targetPID,&
                             overlapPID

    ! moab global id's
    integer,intent(inout) :: sourceGID,&
                             targetGID,&
                             overlapGID

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

    ierr=iMOAB_SetDoubleTagStorage(sourcePID,&
                                   fieldTags,&
                                   totalNumVals,&
                                   entity_type,&
                                   buffer)
    call errorout(ierr,'unable to set storage on source')

    ! send to target "side"
    ierr=iMOAB_SendElementTag(sourcePID,&
                              fieldTags,&
                              moabComm,&
                              overlapGID)
    call errorout(ierr,'unable to send data from source "side" to target "side"')

    ierr=iMOAB_ReceiveElementTag(overlapPID,&
                                 fieldTags,&
                                 moabComm,&
                                 sourceGID)
    call errorout(ierr,'unable to receive data on target "side"')

    ierr=iMOAB_FreeSenderBuffers(sourcePID,overlapGID)
    call errorout(ierr,"unable to free sender buffers")

    ierr=iMOAB_ApplyScalarProjectionWeights(overlapPID,&
                                            filter_type,&
                                            weights_identifier,&
                                            fieldTags,&
                                            fieldTagsProjected)
    call errorout(ierr,'unable to reproject :(')

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
    if (domain % domainID .eq. 0) then
       return
    endif
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
                                                            sValues) !{{{

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

    ! allocate global arrays
    if (curProc.eq.0) then
       allocate(globalArray_source(nEdgesGlobal_source),&
                gatheredArray_source(nEdgesGlobal_source),&
                globalArray_target(nEdgesGlobal_target),&
                gatheredArray_target(nEdgesGlobal_target))
             
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
    localArray_target(:)=0d0
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

    weights_filename_btr2bc='16km_to_32km_weights.nc'
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

!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_dmpar_init_btr !}}}

!***********************************************************************

end module ocn_spatial_mode_splitting
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
