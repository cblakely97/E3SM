
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

   public :: ocn_spatial_mode_splitting_init_moab_domain

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
   ! edge tangent vectors for calculating cell-centered quantities
   real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors
   real (kind=RKIND), dimension(:,:), pointer :: edgeTangentVectors_btr
   real (kind=RKIND), dimension(:), pointer :: barotropicForcing, barotropicThicknessFlux, &
                                               barotropicForcing_btr, barotropicThicknessFlux_btr 
   ! number of vertices on meshes (for loops to avoid calling iMOAB_GetMeshInfo repeatedly)
   integer, pointer :: nCellsSolve, nCellsSolve_btr
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
    subroutine ocn_spatial_mode_splitting_init_moab_domain(domain)
    type (domain_type), intent(inout) :: domain
    type (mpas_pool_type), pointer :: meshPool
    integer, pointer :: nCellsSolve
    logical, pointer :: config_spatial_mode_splitting
    integer :: ierr=0

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
       allocate ( barotropicThicknessFluxX(nCellsSolve),&
                  barotropicThicknessFluxY(nCellsSolve),&
                  barotropicThicknessFluxZ(nCellsSolve) )
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







    end subroutine
#endif
end module ocn_spatial_mode_splitting
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
