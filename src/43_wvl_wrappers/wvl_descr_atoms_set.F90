!!****f* defs_wvltypes/wvl_descr_atoms_set
!!
!! NAME
!! wvl_descr_atoms_set
!!
!! FUNCTION
!! Defines wvl%atoms% data structure
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=unit cell length scales (bohr)
!! dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>= wavelet type
!!                 | nat      =  number of atoms
!!                 | ntypes   =  number of species
!!                 | alat1    =  acell(1)
!!                 | alat2    =  acell(2)
!!                 | alat3    =  acell(3)
!!                 | iatype   =  types for atoms
!!                 | lfrztyp  =  flag for the movement of atoms.
!!                 | natpol   =  integer related to polarisation at the first step
!!
!! PARENTS
!!      gstate,wvl_memory
!!
!! CHILDREN
!!      allocate_atoms_nat,allocate_atoms_ntypes,astruct_set_n_atoms
!!      astruct_set_n_types,f_release_routine,f_routine
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_descr_atoms_set(acell, icoulomb, natom, ntypat, typat, wvl)

 use m_profiling_abi

  use defs_basis
  use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use  module_types, only: atoms_null
  use dynamic_memory, only: f_routine,f_release_routine
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_atoms_set'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in)                    :: icoulomb, natom, ntypat
  type(wvl_internal_type), intent(inout) :: wvl
  !arrays
  integer, intent(in)                    :: typat(natom)
  real(dp), intent(in)                   :: acell(3)

!Local variables-------------------------------
!scalars
#if defined HAVE_DFT_BIGDFT
  integer :: itype
#endif

! *********************************************************************

#if defined HAVE_DFT_BIGDFT
 call f_routine(ABI_FUNC)

 wvl%atoms=atoms_null()

!We create the atoms_data structure from this dataset
!to be used later in BigDFT routines.
 if (icoulomb == 0) then
   wvl%atoms%astruct%geocode = 'P'
 else if (icoulomb == 1) then
   wvl%atoms%astruct%geocode = 'F'
 else if (icoulomb == 2) then
   wvl%atoms%astruct%geocode = 'S'
 end if
 write(wvl%atoms%astruct%units, "(A)") "Bohr"

 call astruct_set_n_atoms(wvl%atoms, natom, ABI_FUNC)
 call astruct_set_n_types(wvl%atoms, ntypat, ABI_FUNC)

 do itype = 1, ntypat, 1
   write(wvl%atoms%astruct%atomnames(itype), "(A,I2)") "At. type", itype
 end do
 wvl%atoms%astruct%cell_dim(1)   =  acell(1)
 wvl%atoms%astruct%cell_dim(2)   =  acell(2)
 wvl%atoms%astruct%cell_dim(3)   =  acell(3)
 wvl%atoms%astruct%iatype   = typat

 wvl%atoms%astruct%sym%symObj = 0

 call allocate_atoms_nat(wvl%atoms, ABI_FUNC)
 call allocate_atoms_ntypes(wvl%atoms, ABI_FUNC)

 call f_release_routine()


#endif  
end subroutine wvl_descr_atoms_set
!!***
