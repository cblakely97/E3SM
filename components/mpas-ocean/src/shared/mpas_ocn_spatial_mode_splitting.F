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
              ocn_spatial_mode_splitting_time_int_split_setup_btr

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------

    ! mpi groups/comms. All of these span the same processors
    integer :: bcComm, btrComm
    integer :: cplComm, bcCouComm, btrCouComm   ! coupler comms
    integer :: bcGroup, btrGroup, cplGroup      ! groups
   
    ! identifiers for moab meshes
    ! id's on compute core 
    integer :: cmpbc    = 314
    integer :: cmpbtr   = 315
    integer :: bc2btrid = 416
    integer :: btr2bcid = 417

    ! id's on coupler *NOTE: we don't presently have separate coupler cores
    ! but iMOAB requires us to register a coupler anyways.
    integer :: cplbc    = 414
    integer :: cplbtr   = 415

    ! PID's for MOAB
    integer,public :: cmpbcPID, cmpbtrPID, bc2btrPID, btr2bcPID
    integer,public :: cplbcPID, cplbtrPID

    ! iMOAB communication variables
    character*100 :: appname
    character(:), allocatable :: bc2btrFields, bc2btrFields_btr,&
                                 btr2bcFields, btr2bcFields_btr
    integer, dimension(2) :: bc2btrIndices, btr2bcIndices,&
                             bc2btrTypes, btr2bcTypes, &
                             bc2btrIndicesCou, btr2bcIndicesCou, &
                             bc2btrTypesCou, btr2bcTypesCou
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
    integer :: partScheme = 0 ! trivial partitioner
    integer :: context_id

    ! End preamble
    !-----------------------------------------------------------------
    ! Begin code

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
       call init_moab_mpas(domain, cmpbc, cmpbcPID, bcComm)

       ! get mesh size for allocating local variables
       call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCells', nCellsSolve_bc)
       call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesSolve_bc)

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
       call init_moab_mpas(domain, cmpbtr, cmpbtrPID, btrComm)

       ! get mesh size for allocating local variables
       call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCells', nCellsSolve_btr)
       call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesSolve_btr)

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
    ! define communicator groups
    !--------------------------------------------------------------------

    ! setup mpi communicators for coupler domains
    cplComm    = MPI_COMM_NULL
    bcCouComm  = MPI_COMM_NULL
    btrCouComm = MPI_COMM_NULL
    call mpi_comm_dup(domain % dminfo % comm, &
                      cplComm, ierr)
    call mpi_comm_dup(domain % dminfo % comm, &
                      bcCouComm, ierr)
    call mpi_comm_dup(domain % dminfo % comm, &
                      btrCouComm, ierr)

    ! coupler group
    call mpi_comm_group(cplComm, cplGroup, ierr)

    ! register coupler applications

    ! baroclinic coupler
    appname = 'bcCouple'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname,&
                                     cplComm, &
                                     cplbc, &
                                     cplbcPID)
    call errorout(ierr, "failed to register baroclinic coupler application")

    ! barotropic coupler
    appname = 'btCouple'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname,&
                                     cplComm, &
                                     cplbtr, &
                                     cplbtrPID)
    call errorout(ierr, "failed to register barotropic coupler application")

    ! register baroclinic to barotropic application on coupler
    appname = 'bc2btr_cpl'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, &
                                     cplComm, &
                                     bc2btrid, &
                                     bc2btrPID)
    call errorout(ierr, "failed to register baroclinic to barotropic application")

    ! register barotropic to baroclinic application on coupler
    appname = 'btr2bc_cpl'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, &
                                     cplComm, &
                                     btr2bcid, &
                                     btr2bcPID)
    call errorout(ierr, "failed to register barotropic to baroclinic application")

    !--------------------------------------------------------------------
    ! send baroclinic and barotropic meshes to the coupler
    !--------------------------------------------------------------------

    ! baroclinic
    if (bcComm .ne. MPI_COMM_NULL) then
       ierr = iMOAB_SendMesh(cmpbcPID, bcCouComm, cplGroup, cplbc, partScheme)
       call errorout(ierr, "failed to send baroclinic to coupler")
    endif
    if (cplComm .ne. MPI_COMM_NULL) then
       ierr = iMOAB_ReceiveMesh(cplbcPID, bcCouComm, bcGroup, cmpbc)
       call errorout(ierr,"failed to receive baroclinic on coupler")
    endif

    ! free send buffers
    if (bcComm .ne. MPI_COMM_NULL) then
       context_id = cplbc
       ierr = iMOAB_FreeSenderBuffers(cmpbcPID, context_id)
    endif

    ! barotropic
    if (btrComm .ne. MPI_COMM_NULL) then
       ierr = iMOAB_SendMesh(cmpbtrPID, btrCouComm, cplGroup, cplbtr, partScheme)
       call errorout(ierr, "failed to send barotropic to coupler")
    endif
    if (cplComm .ne. MPI_COMM_NULL) then
       ierr = iMOAB_ReceiveMesh(cplbtrPID, btrCouComm, btrGroup, cmpbtr)
       call errorout(ierr,"failed to receive baroclinic on coupler")
    endif

    ! free send buffers
    if (btrComm .ne. MPI_COMM_NULL) then
       context_id = cplbtr
       ierr = iMOAB_FreeSenderBuffers(cmpbtrPID, context_id)
    endif

    ! once both the baroclinic and barotropic meshes are registered in moab, we need to calculate
    ! the overlap meshes and compute the weights for interpolation as well as define our tags

    !--------------------------------------------------------------------
    ! compute mesh intersections on coupler
    !--------------------------------------------------------------------

    if (cplComm .ne. MPI_COMM_NULL) then
       ! baroclinic to barotropic
       ierr = iMOAB_ComputeMeshIntersectionOnSphere(cplbcPID,&
                                                    cplbtrPID, &
                                                    bc2btrPID)
       call errorout(ierr,"failed to compute intersection of baroclinic to barotropic mesh")

       ! baroclinic to barotropic
       ierr = iMOAB_ComputeMeshIntersectionOnSphere(cplbtrPID,&
                                                    cplbcPID, &
                                                    btr2bcPID)
       call errorout(ierr,"failed to compute intersection of barotropic to baroclinic mesh")
    endif

    !--------------------------------------------------------------------
    ! calculate coverage graphs on bcComm and btrComm
    !--------------------------------------------------------------------

    ! baroclinic coverage graph
    if (bcComm .ne. MPI_COMM_NULL) then
       ierr = iMOAB_CoverageGraph(bcCouComm, cmpbcPID, cplbcPID, bc2btrPID, &
                                  cmpbc, cplbc, cplbtr)
       call errorout(ierr, "failed to compute coverage graph for baroclinic")
    endif

    ! barotropic coverage graph
    if (btrComm .ne. MPI_COMM_NULL) then
       ierr = iMOAB_CoverageGraph(btrCouComm, cmpbtrPID, cplbtrPID, btr2bcPID, &
                                  cmpbtr, cplbtr, cplbc)
       call errorout(ierr, "failed to compute coverage graph for baroclinic")
    endif

    !--------------------------------------------------------------------
    ! compute scalar projection weights
    !--------------------------------------------------------------------

    ! get projection weights for baroclinic to barotropic. Note that both domains are finite
    ! volume discretizations with the same discretization order. Due to this we just use the
    ! same variables for both source and target mesh

    ! iMOAB variables for calculating the weights
    bc2btr_weights_identifier = 'bc2btr'//C_NULL_CHAR
    btr2bc_weights_identifier = 'btr2bc'//C_NULL_CHAR
    disc_method = 'fv'//C_NULL_CHAR
    disc_order = 1
    dof_tag_names = 'GLOBAL_ID'//C_NULL_CHAR
    fNoBubble = 1
    fMonotoneTypeID = 0
    fVolumetric = 0
    fNoConserve = 0
    fValidate = 0
    fInverseDistanceMap = 0

    if (cplComm .ne. MPI_COMM_NULL) then

       ! for coarse to fine (baroclinic to barotropic) use bilinear
       fv_method = 'bilin'//C_NULL_CHAR ! not sure what this is for, figure it out later

       ! baroclinic => barotropic
       ierr = iMOAB_ComputeScalarProjectionWeights(bc2btrPID, &
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

       ! for fine to coarse (barotropic to baroclinic) use conservative
       fv_method = ''//C_NULL_CHAR ! not sure what this is for, figure it out later

       ! barotropic => baroclinic
       ierr = iMOAB_ComputeScalarProjectionWeights(btr2bcPID, &
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
    endif

    !--------------------------------------------------------------------
    ! define tags on meshes
    !--------------------------------------------------------------------

    ! all field types are dense double
    bc2btrTypes(:)    = 1 ! DENSE_DOUBLE
    btr2bcTypes(:)    = 1 ! DENSE_DOUBLE
    bc2btrTypesCou(:) = 1
    btr2bcTypesCou(:) = 1

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

    ! define tags on bcComm
    if (bcComm .ne. MPI_COMM_NULL) then

       ! baroclinic to barotropic ON baroclinic
       ierr = iMOAB_DefineTagStorage(cmpbcPID, &
                                     bc2btrFields, &
                                     bc2btrTypes(1), &
                                     1, &
                                     bc2btrIndices(1))
       call errorout(ierr, "failed to define bc2btr tags on bc mesh")

       ! barotropic to baroclinic ON baroclinic
       ierr = iMOAB_DefineTagStorage(cmpbcPID, &
                                     btr2bcFields, &
                                     btr2bcTypes(1), &
                                     1, &
                                     btr2bcIndices(1))
       call errorout(ierr, "failed to define btr2bc tags on bc mesh")
    endif

    ! define tags on btrComm
    if (btrComm .ne. MPI_COMM_NULL) then

       ! baroclinic to barotropic ON barotropic
       ierr = iMOAB_DefineTagStorage(cmpbtrPID, &
                                     bc2btrFields_btr, &
                                     bc2btrTypes(2), &
                                     1, &
                                     bc2btrIndices(2))
       call errorout(ierr, "failed to define bc2btr tags on btr mesh")

       ! barotropic to baroclinic ON barotropic
       ierr = iMOAB_DefineTagStorage(cmpbtrPID, &
                                     btr2bcFields_btr, &
                                     btr2bcTypes(2), &
                                     1, &
                                     btr2bcIndices(2))
       call errorout(ierr, "failed to define btr2bc tags on btr mesh")
    endif

    ! define tags on cplComm
    if (cplComm .ne. MPI_COMM_NULL) then

       ! baroclinic to barotropic ON baroclinic
       ierr = iMOAB_DefineTagStorage(cplbcPID, &
                                     bc2btrFields, &
                                     bc2btrTypesCou(1), &
                                     1, &
                                     bc2btrIndicesCou(1))
       call errorout(ierr, "failed to define bc2btr tags on bc mesh (coupler)")

       ! barotropic to baroclinic ON baroclinic
       ierr = iMOAB_DefineTagStorage(cplbcPID, &
                                     btr2bcFields, &
                                     btr2bcTypesCou(1), &
                                     1, &
                                     btr2bcIndicesCou(1))
       call errorout(ierr, "failed to define btr2bc tags on bc mesh (coupler)")

       ! baroclinic to barotropic ON barotropic
       ierr = iMOAB_DefineTagStorage(cplbtrPID, &
                                     bc2btrFields_btr, &
                                     bc2btrTypesCou(2), &
                                     1, &
                                     bc2btrIndicesCou(2))
       call errorout(ierr, "failed to define bc2btr tags on btr mesh (coupler)")

       ! barotropic to baroclinic ON barotropic
       ierr = iMOAB_DefineTagStorage(cplbtrPID, &
                                     btr2bcFields_btr, &
                                     btr2bcTypesCou(2), &
                                     1, &
                                     btr2bcIndicesCou(2))
       call errorout(ierr, "failed to define btr2bc tags on btr mesh")
    endif

    ! initialize buffers for transmitting the data (both baroclinic to barotropic and
    ! baroclinic to barotropic while we have the mesh size)
    ierr = iMOAB_GetMeshInfo(cmpbcPID, nverts, nelem, nblocks, nsbc, ndbc)
    call errorout(ierr, "failed to retrieve mesh info from baroclinic mesh")
    allocate ( bc2btrVar(2*nelem(3)), btr2bcVar(2*nelem(3)) )
    bc2btrVar(:) = 0d0
    btr2bcVar(:) = 0d0

    ! initialize buffers for transmitting data from barotropic to baroclinic mesh
    ierr = iMOAB_GetMeshInfo(cmpbtrPID, nverts, nelem, nblocks, nsbc, ndbc)
    call errorout(ierr, "failed to retrieve mesh info from barotropic mesh")
    allocate ( bc2btrVar_btr(2*nelem(3)), btr2bcVar_btr(2*nelem(3)) )
    bc2btrVar_btr(:) = 0d0
    btr2bcVar_btr(:) = 0d0
#endif
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
    ! 
    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: statePool
    ! mesh information
    integer, pointer :: nCellsSolve_bc, nCellsSolve_btr,&
                        nEdgesSolve_bc, nEdgesSolve_btr,&
                        vertexDegree_btr
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
        call mpas_pool_get_dimension(meshPool, 'nCells', nCellsSolve_bc)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesSolve_bc)

        ! get barotropicForcing from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingZonal',&
                                 barotropicForcingZonal_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicForcingMeridional',&
                                 barotropicForcingMeridional_bc)

        ! reconstruct from edge based to cell-centered quantities
        call mpas_reconstruct(meshPool, barotropicForcing_bc, &
                              barotropicForcingX_bc, &
                              barotropicForcingY_bc, &
                              barotropicForcingZ_bc, &
                              barotropicForcingZonal_bc, &
                              barotropicForcingMeridional_bc)

        ! pack buffers for iMOAB
        !do iCell = 1,nCellsSolve_bc
        !    bc2btrVar(iCell) = barotropicForcingZonal_bc(iCell)
        !    bc2btrVar(iCell+nCellsSolve_bc) = barotropicForcingMeridional_bc(iCell)
        !enddo

        ! assign data to baroclinic mesh on compute processor 
        if (bcComm .ne. MPI_COMM_NULL) then

           ! assign values to moab
           ierr=iMOAB_SetDoubleTagStorage(cmpbcPID,&
                                          bc2btrFields,&
                                          2*nCellsSolve_bc,&
                                          1,&
                                          bc2btrVar)
           call errorout(ierr,"failed to assign bc2btr data to baroclinic mesh")

           ! send element tag
           ierr = iMOAB_SendElementTag(cmpbcPID,&
                                       bc2btrFields,&
                                       bcCouComm,&
                                       cplbtr)
           call errorout(ierr,"failed to send bc2btr data tag")
        endif

        ! receive data on coupler
        if (cplComm .ne. MPI_COMM_NULL) then
           
           ! receive data
           ierr = iMOAB_ReceiveElementTag(cplbcPID,&
                                          bc2btrFields,&
                                          bcCouComm,&
                                          cplbtr)
           call errorout(ierr,"failed to receive bc2btr data tag")
        endif

        ! free send buffers
        if (bcComm .ne. MPI_COMM_NULL) then
           ierr = iMOAB_FreeSenderBuffers(cmpbcPID, cplbtr)
           call errorout(ierr,"failed to free send buffers 698")
        endif

        ! do remapping on coupler
        if (cplComm .ne. MPI_COMM_NULL) then
           ierr = iMOAB_ApplyScalarProjectionWeights(bc2btrPID,&
                                                     0,&
                                                     bc2btr_weights_identifier,&
                                                     bc2btrFields,&
                                                     bc2btrFields_btr)
           call errorout(ierr,"failed to reproject bc2btr data to barotropic mesh")

           ! send data to cmpbtr
           ierr = iMOAB_SendElementTag(cplbtrPID, &
                                       bc2btrFields_btr, &
                                       btrCouComm, cmpbtr)
           call errorout(ierr,"failed to send bc2btr_btr tags")
        endif

        ! receive data from coupler on barotropic
        if (btrComm .ne. MPI_COMM_NULL) then
           ierr = iMOAB_ReceiveElementTag(cmpbtrPID, &
                                          bc2btrFields_btr, &
                                          btrCouComm, &
                                          cplbtr)
           call errorout(ierr,"failed to receive bc2btr_btr tags")
        endif

        ! free send buffers
        if (cplComm .ne. MPI_COMM_NULL) then
           ierr = iMOAB_FreeSenderBuffers(cplbtrPID, cmpbtr)
           call errorout(ierr,"failed to free send buffers 728")
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

        ! get num cells/edges
        call mpas_pool_get_dimension(meshPool, 'nCells', nCellsSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesSolve_btr)
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
        ierr = iMOAB_GetDoubleTagStorage(cmpbtrPID,&
                                         bc2btrFields_btr,&
                                         2*nCellsSolve_btr,&
                                         1,&
                                         bc2btrVar_btr)
        call errorout(ierr,"failed to get bc2btrVar_btr")

        ! unpack buffers
        do iCell = 1,nCellsSolve_btr
            !barotropicForcingZonal_btr(iCell) = bc2btrVar_btr(iCell)
            !barotropicForcingMeridional_btr(iCell) = bc2btrVar_btr(iCell+nCellsSolve_btr)
            !
!            barotropicForcingZonal_btr(iCell)=barotropicForcingZonal_bc(iCell)
!            barotropicForcingMeridional_btr(iCell)=barotropicForcingMeridional_bc(iCell)
        enddo

        ! perform halo exchange
        call mpas_dmpar_field_halo_exch(domain_btr, 'barotropicForcingZonal')
        call mpas_dmpar_field_halo_exch(domain_btr, 'barotropicForcingMeridional')

        ! reconstruct cell based zonal and meridional components to edge-based zonal and meridional
        ! components
        call mpas_vector_cell_to_edge_isotropic(barotropicForcingZonal_btr,&
                                                barotropicForcingMeridional_btr,&
                                                nEdgesSolve_btr,&
                                                vertexDegree_btr,&
                                                verticesOnEdge_btr,&
                                                cellsOnVertex_btr,&
                                                kiteAreasOnVertex_btr,&
                                                barotropicForcingOnEdgeZonal_btr,&
                                                barotropicForcingOnEdgeMeridional_btr)

        ! put back into edge based arrays
         barotropicForcing_btr(:) = 0_RKIND
         do iEdge = 1,nEdgesSolve_btr
            barotropicForcing_btr(iEdge) = barotropicForcing_btr(iEdge) + &
                          cos(angleEdge_btr(iEdge))*barotropicForcingOnEdgeZonal_btr(iEdge) + &
                          sin(angleEdge_btr(iEdge))*barotropicForcingOnEdgeMeridional_btr(iEdge)
           ! carry over edge-based values
           barotropicForcing_btr(iEdge)=barotropicForcing_bc(iEdge)
         enddo

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

    ! mesh information
    integer, pointer :: nCellsSolve_bc, nCellsSolve_btr,&
                        nEdgesSolve_bc, nEdgesSolve_btr

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
        call mpas_pool_get_dimension(meshPool, 'nCells', nCellsSolve_btr)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesSolve_btr)

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
        !do iCell = 1,nCellsSolve_btr
        !    btr2bcVar_btr(iCell) = barotropicThicknessFluxZonal_btr(iCell)
        !    btr2bcVar_btr(iCell+nCellsSolve_btr) = barotropicThicknessFluxMeridional_btr(iCell)
        !enddo

        ! assign data to barotropic mesh on compute processor
        if (btrComm .ne. MPI_COMM_NULL) then

           ! assign values to moab
           ierr = iMOAB_SetDoubleTagStorage(cmpbtrPID,&
                                            btr2bcFields_btr,&
                                            2*nCellsSolve_btr,&
                                            1,&
                                            btr2bcVar_btr)
           call errorout(ierr,"failed to assign btr2bcFields_btr to barotropic mesh")

           ! send element tag
           ierr = iMOAB_SendElementTag(cmpbtrPID,&
                                       btr2bcFields_btr,&
                                       btrCouComm,&
                                       cplbc)
           call errorout(ierr,"failed to send btr2bcFields_btr to coupler")
        endif

        ! receive data on coupler
        if (cplComm .ne. MPI_COMM_NULL) then
           
           ! receive data
           ierr = iMOAB_ReceiveElementTag(cplbtrPID,&
                                          btr2bcFields_btr,&
                                          btrCouComm,&
                                          cplbc)
           call errorout(ierr,"failed to receive btr2bcFields_btr on coupler")
        endif

        ! free send buffers
        if (btrComm .ne. MPI_COMM_NULL) then
           ierr = iMOAB_FreeSenderBuffers(cmpbtrPID, cplbc)
           call errorout(ierr,"failed to free send buffers 936")
        endif

        ! do remapping
        if (cplComm .ne. MPI_COMM_NULL) then
           ierr = iMOAB_ApplyScalarProjectionWeights(btr2bcPID,&
                                                     0,&
                                                     btr2bc_weights_identifier,&
                                                     btr2bcFields_btr,&
                                                     btr2bcFields)
           call errorout(ierr,"failed to reproject btr2bc data to baroclinic mesh")

           ! send data to cmpbc
           ierr = iMOAB_SendElementTag(cplbcPID,&
                                       btr2bcFields, &
                                       bcCouComm,cmpbc)
           call errorout(ierr,"failed to send btr2bcFields")
        endif

        ! receive data from coupler on baroclinic
        if (bcComm .ne. MPI_COMM_NULL) then
           ierr = iMOAB_ReceiveElementTag(cmpbcPID,&
                                          btr2bcFields,&
                                          bcCouComm,&
                                          cplbc)
           call errorout(ierr,"failed to receive btr2bcFields")
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
        call mpas_pool_get_dimension(meshPool, 'nCells', nCellsSolve_bc)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesSolve_bc)

        ! get barotropicThicknessFlux from 
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux',&
                                 barotropicThicknessFlux_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxZonal',&
                                 barotropicThicknessFluxZonal_bc)
        call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFluxMeridional',&
                                 barotropicThicknessFluxMeridional_bc)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge_bc)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge_bc)

        ! get remapped data from iMOAB
        ierr = iMOAB_GetDoubleTagStorage(cmpbcPID,&
                                         btr2bcFields,&
                                         2*nCellsSolve_bc,&
                                         1,&
                                         btr2bcVar)
        call errorout(ierr,"failed to get btr2bcVar from iMOAB")

        ! unpack buffers
        do iCell = 1,nCellsSolve_bc
         !   barotropicThicknessFluxZonal_bc(iCell) = btr2bcVar(iCell)
         !   barotropicThicknessFluxMeridional_bc(iCell) = btr2bcVar(iCell+nCellsSolve_bc)
            !
!            barotropicThicknessFluxZonal_bc(iCell)=barotropicThicknessFluxZonal_btr(iCell)
!            barotropicThicknessFluxMeridional_bc(iCell)=barotropicThicknessFluxMeridional_btr(iCell)
        enddo

        ! exchange halos
        call mpas_dmpar_field_halo_exch(domain, 'barotropicThicknessFluxZonal')
        call mpas_dmpar_field_halo_exch(domain, 'barotropicThicknessFluxMeridional')

        ! put into edge-based arrays
        barotropicThicknessFlux_bc(:)=0_RKIND
        do iEdge = 1,nEdgesSolve_bc
           cell1 = cellsOnEdge_bc(1, iEdge)
           cell2 = cellsOnEdge_bc(2,iEdge)
           zonal_edge = 0.5_RKIND * (barotropicThicknessFluxZonal_bc(cell1) + &
                                 barotropicThicknessFluxZonal_bc(cell2))
           meridional_edge = 0.5_RKIND * (barotropicThicknessFluxMeridional_bc(cell1) + &
                                      barotropicThicknessFluxMeridional_bc(cell2))
           barotropicThicknessFlux_bc(iEdge) = barotropicThicknessFlux_bc(iEdge) + &
                                            cos(angleEdge_bc(iEdge))*zonal_edge + &
                                            sin(angleEdge_bc(iEdge))*meridional_edge
          ! direct carry over
          barotropicThicknessFlux_bc(iEdge)=barotropicThicknessFlux_btr(iEdge)
        enddo
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
        enddo

        block_ptr => block_ptr % next
    enddo
!-------------------------------------------------------------------------------

    end subroutine ocn_spatial_mode_splitting_time_int_split_setup_btr !}}}

end module ocn_spatial_mode_splitting
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
