!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_energies
!! NAME
!!  m_energies
!!
!! FUNCTION
!!  This module provides the definition of the energies used
!!  to store energies from GS calculations.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MT, DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


MODULE m_energies

 use defs_basis
 use m_profiling_abi
 use m_errors

 use defs_abitypes, only : dataset_type

 implicit none

 private

!public parameter
 integer, public, parameter :: n_energies=22
  ! Number of energies stored in energies datastructure

!!***

!!****t* m_energies/energies_type
!! NAME
!! energies_type
!!
!! FUNCTION
!! This structured datatype contains all parts of total energy. Not all
!! attributes may have a value, depending on the scheme used to
!! compute the total energy and several options read from dtset.
!!
!! SOURCE

 type, public :: energies_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  real(dp) :: e_localpsp
   ! Local psp energy (hartree)

  real(dp) :: e_eigenvalues
   ! Sum of the eigenvalues - Band energy (Hartree)
   ! (valid for double-counting scheme dtset%optene == 1)

  real(dp) :: e_ewald
   ! Ewald energy (hartree), store also the ion/ion energy for free boundary conditions.

  real(dp) :: e_hartree
   ! Hartree part of total energy (hartree units)

  real(dp) :: e_corepsp
   ! psp core-core energy

  real(dp) :: e_corepspdc
   ! psp core-core energy double-counting

  real(dp) :: e_kinetic
   ! Kinetic energy part of total energy.
   ! (valid for direct scheme, dtset%optene == 0)

  real(dp) :: e_nonlocalpsp
   ! Nonlocal pseudopotential part of total energy.

  real(dp) :: e_entropy
   ! Entropy energy due to the occupation number smearing (if metal)
   ! Value is multiplied by dtset%tsmear, see %entropy for the entropy alone.
   ! (valid for metals, dtset%occopt>=3 .and. dtset%occopt<=8)

  real(dp) :: entropy

  real(dp) :: e_xc
   ! Exchange-correlation energy (hartree)

  real(dp) :: e_vdw_dftd2
   ! Dispersion energy from DFT-D2 Van der Waals correction (hartree)

  real(dp) :: e_xcdc
   ! enxcdc=exchange-correlation double-counting energy (hartree)

  real(dp) :: e_paw
   ! PAW spherical part energy

  real(dp) :: e_pawdc
   ! PAW spherical part double-counting energy

  real(dp) :: e_elecfield
   ! Electric enthalpy, by adding both ionic and electronic contributions

  real(dp) :: e_magfield
   ! Orbital magnetic enthalpy, by adding orbital contribution

  real(dp) :: e_fermie
   ! Fermie energy

  real(dp) :: e_sicdc
   ! Self-interaction energy double-counting

  real(dp) :: e_exactX
   ! Fock exact-exchange energy (hartree)

  real(dp) :: h0
   ! h0=e_kinetic+e_localpsp+e_nonlocalpsp

  real(dp) :: e_electronpositron
   ! Electron-positron: electron-positron interaction energy

  real(dp) :: edc_electronpositron
   ! Electron-positron: double-counting electron-positron interaction energy

  real(dp) :: e0_electronpositron
   !  Electron-positron: energy only due to unchanged particles
   !                     (if calctype=1, energy due to electrons only)
   !                     (if calctype=2, energy due to positron only)

  real(dp) :: e_monopole
   ! Monopole correction to the total energy for charged supercells

  real(dp) :: e_xc_vdw
   ! vdW-DF correction to the XC energy

 end type energies_type

!public procedures.
 public :: energies_init
 public :: energies_copy
 public :: energies_to_array
 public :: energies_eval_eint
!!***


CONTAINS !===========================================================
!!***

!!****f* m_energies/energies_init
!!
!! NAME
!! energies_init
!!
!! FUNCTION
!! Set zero in all values of a type(energies_type) object
!!
!! INPUTS
!!
!! OUTPUT
!!   energies <type(energies_type)>=values to initialise
!!
!! PARENTS
!!      bethe_salpeter,gstate,m_electronpositron,m_results_gs,scfcv,screening
!!      setup_positron,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine energies_init(energies)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'energies_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(out) :: energies

! *************************************************************************

!@energies_type

 energies%entropy       = zero
 energies%e_entropy     = zero
 energies%e_fermie      = zero
 energies%e_paw         = zero
 energies%e_pawdc       = zero
 energies%e_kinetic     = zero
 energies%e_localpsp    = zero
 energies%e_nonlocalpsp = zero
 energies%e_eigenvalues = zero
 energies%e_hartree     = zero
 energies%e_ewald       = zero
 energies%e_xc          = zero
 energies%e_vdw_dftd2   = zero
 energies%e_xcdc        = zero
 energies%e_xc_vdw      = zero
 energies%e_elecfield   = zero
 energies%e_corepsp     = zero
 energies%e_corepspdc   = zero
 energies%e_exactX      = zero
 energies%e_electronpositron   = zero
 energies%edc_electronpositron = zero
 energies%e0_electronpositron  = zero
 energies%e_monopole    = zero

end subroutine energies_init
!!***

!----------------------------------------------------------------------

!!****f* m_energies/energies_copy
!!
!! NAME
!! energies_copy
!!
!! FUNCTION
!! Copy a type(energies_type) object into another
!!
!! INPUTS
!!   energies_in <type(energies_type)>=input values (to copy)
!!
!! OUTPUT
!!   energies_out <type(energies_type)>=output values
!!
!! PARENTS
!!      afterscfloop,m_electronpositron,m_results_gs,setup_positron
!!
!! CHILDREN
!!
!! SOURCE

 subroutine energies_copy(energies_in,energies_out)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'energies_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(in)  :: energies_in
 type(energies_type),intent(out) :: energies_out

!*************************************************************************

!@energies_type

 energies_out%entropy              = energies_in%entropy
 energies_out%e_entropy            = energies_in%e_entropy
 energies_out%e_fermie             = energies_in%e_fermie
 energies_out%e_paw                = energies_in%e_paw
 energies_out%e_pawdc              = energies_in%e_pawdc
 energies_out%e_kinetic            = energies_in%e_kinetic
 energies_out%e_localpsp           = energies_in%e_localpsp
 energies_out%e_nonlocalpsp        = energies_in%e_nonlocalpsp
 energies_out%e_eigenvalues        = energies_in%e_eigenvalues
 energies_out%e_hartree            = energies_in%e_hartree
 energies_out%e_ewald              = energies_in%e_ewald
 energies_out%e_xc                 = energies_in%e_xc
 energies_out%e_vdw_dftd2          = energies_in%e_vdw_dftd2
 energies_out%e_xcdc               = energies_in%e_xcdc
 energies_out%e_xc_vdw             = energies_in%e_xc_vdw
 energies_out%e_elecfield          = energies_in%e_elecfield
 energies_out%e_corepsp            = energies_in%e_corepsp
 energies_out%e_corepspdc          = energies_in%e_corepspdc
 energies_out%e_exactX             = energies_in%e_exactX
 energies_out%e_electronpositron   = energies_in%e_electronpositron
 energies_out%edc_electronpositron = energies_in%edc_electronpositron
 energies_out%e0_electronpositron  = energies_in%e0_electronpositron
 energies_out%e_monopole           = energies_in%e_monopole

end subroutine energies_copy
!!***

!----------------------------------------------------------------------

!!****f* m_energies/energies_to_array
!!
!! NAME
!! energies_to_array
!!
!! FUNCTION
!! Transfer a energies datastructure into a single array or
!! transfer an array into a energies datastructure
!!
!! INPUTS
!!   option= 1: copy energies datastructure into an array
!!   option=-1: copy an array into a energies datastructure
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   energies <type(energies_type)>=energies stored in a datastructure
!!   energies_array=energies stored in a single array
!!
!! PARENTS
!!      m_results_img
!!
!! CHILDREN
!!
!! SOURCE

 subroutine energies_to_array(energies,energies_array,option)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'energies_to_array'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
!arrays
 real(dp),intent(inout) :: energies_array(n_energies)
 type(energies_type),intent(inout)  :: energies

!*************************************************************************

!@energies_type

 ABI_CHECK(n_energies==22,"Bug in energies_to_array !")

 if (option==1) then
   energies_array(1) =energies%entropy
   energies_array(2) =energies%e_entropy
   energies_array(3) =energies%e_fermie
   energies_array(4) =energies%e_paw
   energies_array(5) =energies%e_pawdc
   energies_array(6) =energies%e_kinetic
   energies_array(7) =energies%e_localpsp
   energies_array(8) =energies%e_nonlocalpsp
   energies_array(9) =energies%e_eigenvalues
   energies_array(10)=energies%e_hartree
   energies_array(11)=energies%e_ewald
   energies_array(12)=energies%e_xc
   energies_array(13)=energies%e_vdw_dftd2
   energies_array(14)=energies%e_xcdc
   energies_array(15)=energies%e_elecfield
   energies_array(16)=energies%e_corepsp
   energies_array(17)=energies%e_electronpositron
   energies_array(18)=energies%edc_electronpositron
   energies_array(19)=energies%e0_electronpositron
   energies_array(20)=energies%e_monopole
   energies_array(21)=energies%e_corepspdc
   energies_array(22)=energies%e_exactX
 end if

 if (option==-1) then
   energies%entropy              = energies_array(1)
   energies%e_entropy            = energies_array(2)
   energies%e_fermie             = energies_array(3)
   energies%e_paw                = energies_array(4)
   energies%e_pawdc              = energies_array(5)
   energies%e_kinetic            = energies_array(6)
   energies%e_localpsp           = energies_array(7)
   energies%e_nonlocalpsp        = energies_array(8)
   energies%e_eigenvalues        = energies_array(9)
   energies%e_hartree            = energies_array(10)
   energies%e_ewald              = energies_array(11)
   energies%e_xc                 = energies_array(12)
   energies%e_vdw_dftd2          = energies_array(13)
   energies%e_xcdc               = energies_array(14)
   energies%e_elecfield          = energies_array(15)
   energies%e_corepsp            = energies_array(16)
   energies%e_electronpositron   = energies_array(17)
   energies%edc_electronpositron = energies_array(18)
   energies%e0_electronpositron  = energies_array(19)
   energies%e_monopole           = energies_array(20)
   energies%e_corepspdc          = energies_array(21)
   energies%e_exactX             = energies_array(22)
 end if

end subroutine energies_to_array
!!***

!----------------------------------------------------------------------

!!****f* m_energies/energies_eval_eint
!!
!! NAME
!! energies_eval_eint
!!
!! FUNCTION
!! Compute the internal energy (Direct and DC as it was in prtene)
!!
!! INPUTS
!!  energies <type(energies_type)>=values of parts of total energy
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryphase
!!   | kptopt
!!   | occopt
!!   | positron=option for electron-positron calculation
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!!
!! OUTPUT
!!  optdc=option for double counting scheme
!!  eint=internal energy with direct scheme
!!  eintdc=internal energy with double counting scheme
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      prtene
!!
!! CHILDREN
!!
!! SOURCE

 subroutine energies_eval_eint(energies,dtset,usepaw,optdc,eint,eintdc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'energies_eval_eint'
!End of the abilint section

  implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'energies_eval_eint'
!End of the abilint section

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(in) :: energies
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: usepaw
 integer , intent(out) :: optdc
 real(dp), intent(out) :: eint
 real(dp), intent(out) :: eintdc

!Local variables-------------------------------
! Do not modify the length of this string
!scalars
 integer :: ipositron
 logical :: wvlbigdft=.false.

! *************************************************************************

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 if(dtset%usewvl==1 .and. dtset%wvl_bigdft_comp==1) wvlbigdft=.true.

 optdc=-1;ipositron=0
 if (dtset%positron==0) then
   if (abs(energies%e_xcdc)<1.e-15_dp) optdc=0
   if (abs(energies%e_localpsp)<1.e-15_dp.and.abs(energies%e_xcdc)>1.e-15_dp) optdc=1
   if (abs(energies%e_localpsp)>1.e-15_dp.and.abs(energies%e_xcdc)>1.e-15_dp) optdc=2
   if (wvlbigdft .and. dtset%iscf > 0) optdc=1
 else
   ipositron=2
   if (abs(energies%e_ewald)<1.e-15_dp.and.abs(energies%e_hartree)<1.e-15_dp) ipositron=1
   if (abs(energies%edc_electronpositron)<1.e-15_dp) optdc=0
   if (abs(energies%e_electronpositron)<1.e-15_dp.and.abs(energies%edc_electronpositron)>1.e-15_dp) optdc=1
   if (abs(energies%e_electronpositron)>1.e-15_dp.and.abs(energies%edc_electronpositron)>1.e-15_dp) optdc=2
 end if

 eint  = zero
 eintdc = zero
!============= Evaluate some parts of the energy ===========

 if (optdc==0.or.optdc==2) then
   eint = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&   energies%e_localpsp + energies%e_corepsp 
   if (usepaw==0) eint = eint + energies%e_nonlocalpsp
   if (usepaw==1) eint = eint + energies%e_paw
   if (dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or.  &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17) eint=eint+energies%e_elecfield    !!HONG
   eint = eint + energies%e_ewald + energies%e_vdw_dftd2
   if (ipositron/=0) eint=eint+energies%e0_electronpositron+energies%e_electronpositron
 end if
 if (optdc>=1) then
   eintdc = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&   energies%e_xcdc + energies%e_corepsp - energies%e_corepspdc
   if (usepaw==1) eintdc = eintdc + energies%e_pawdc
   if (dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or.  &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17) eintdc = eintdc + energies%e_elecfield
   eintdc = eintdc + energies%e_ewald + energies%e_vdw_dftd2
   if (ipositron/=0) eintdc=eintdc-energies%edc_electronpositron &
&   +energies%e0_electronpositron+energies%e_electronpositron
 end if

end subroutine energies_eval_eint
!!***

END MODULE m_energies
!!***

