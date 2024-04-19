
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
   ! id's on coupler
   integer :: cplbc, cplbtr
   ! PID's for MOAB
   integer,public :: cmpbcPID, cmpbtrPID, bc2btrPID, btr2bcPID
   integer,public :: cplbcPID, cplbtrPID
   ! iMOAB communication tags

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
    endif
    end subroutine
#endif
end module ocn_spatial_mode_splitting
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
