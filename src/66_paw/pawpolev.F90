!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpolev
!! NAME
!! pawpolev
!!
!! FUNCTION
!! Compute the PAW term for polarization, named expected value term
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (FJ, PH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  ntypat = number of atom types
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  pelev(3)= electronic polarisation. expectation value term (PAW only)
!!
!! PARENTS
!!      berryphase_new
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawpolev(my_natom,natom,ntypat,pawrhoij,pawtab,pelev,&
&                   mpi_comm_atom) ! optional argument (parallelism)

 use m_profiling_abi

 use defs_basis
 use m_pawtab, only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_errors
 use m_xmpi, only : xmpi_sum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawpolev'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat
 integer,optional,intent(in) :: mpi_comm_atom
!arrays
 real(dp),intent(out) :: pelev(3)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)


!Local variables ---------------------------------------
!scalars
 integer :: iatom,idir,ierr,irhoij,ispden,itypat,jrhoij,klmn
 logical :: paral_atom
 real(dp) :: c1,ro_dlt
!arrays
 integer,dimension(3) :: idirindx = (/4,2,3/)
 real(dp) :: tsec(2)


! *************************************************************************

 DBG_ENTER("COLL")

 call timab(560,1,tsec)

!Check for parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(my_natom/=natom))

!note that when vector r is expanded in real spherical harmonics, the factor
!sqrt(four_pi/three) appears, as in the following
!x = sqrt(four_pi/three)*r*S_{1,1}
!y = sqrt(four_pi/three)*r*S_{1,-1}
!z = sqrt(four_pi/three)*r*S_{1,0}
!
!the moments pawtab()%qijl do not include such normalization factors
!see pawinit.F90 for their definition and computation

 c1=sqrt(four_pi/three)

 pelev=zero
 do idir=1,3
   do iatom=1,my_natom
     itypat=pawrhoij(iatom)%itypat
     do ispden=1,pawrhoij(iatom)%nspden
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         ro_dlt=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         pelev(idir)=pelev(idir)+ro_dlt*c1*pawtab(itypat)%qijl(idirindx(idir),klmn)
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do
     end do
   end do
 end do

 if (paral_atom) then
   call xmpi_sum(pelev,mpi_comm_atom,ierr)
 end if

 call timab(560,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawpolev
!!***
