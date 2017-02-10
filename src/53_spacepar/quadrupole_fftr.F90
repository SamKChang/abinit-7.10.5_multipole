!{\src2tex{textfont=tt}}
!!****f* ABINIT/multipoles_fftr
!! NAME
!! multipoles_fftr
!!
!! FUNCTION
!!  Compute spatial multipole moments of input array on FFT grid
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!  call mean_fftr to deal with the averaging over several MPI OMP processors
!!
!! COPYRIGHT
!! Copyright (C) 2003-2014 ABINIT group (MJV, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
!!  nfft=number of FFT points stored by one proc
!!  nfftot=total number of FFT points
!!  ngfft =number of subdivisions along each lattice vector
!!  nspden=number of spin-density components
!!  rprimd = dimensionful lattice vectors
!!
!! OUTPUT
!!  dipole(nspden)=mean value of the dipole of input array, for each nspden component
!!
!! PARENTS
!!      multipoles_fftr
!!
!! CHILDREN
!!      multipoles_fftr,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!KYSC 20170209 quadrupole tensor
subroutine quadrupole_fftr(arraysp,dipole,mpi_enreg,nfft,ngfft,nspden,rprimd,neworigin)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_cgtools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'quadrupole_fftr'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ngfft(3),nspden
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: neworigin(3)
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: dipole(3,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft1, ifft2, ifft3, ifft, ispden, nfftot
 character(len=500) :: message
 real(dp) :: invn1, invn2, invn3
!arrays
 real(dp) :: meansp(nspden)
 real(dp),allocatable :: wrapfft(:)
 real(dp), allocatable :: tmpsp(:,:)

 print *, "yo working function of quadrupole_fftr from new file"

end subroutine quadrupole_fftr

subroutine quadrupole_tensor_out(arraysp,mpi_enreg,natom,nfft,ngfft,nspden,&
&  ntypat,rprimd,typat,ucvol,xred,ziontypat)

 use m_profiling_abi
 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'quadrupole_tensor_out'
 use interfaces_14_hidewrite
 use interfaces_53_spacepar, except_this_one => quadrupole_tensor_out
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom, ntypat
 integer,intent(in) :: nfft,ngfft(3),nspden
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp), intent(in) :: ucvol
!arrays
 integer, intent(in) :: typat(natom)
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(in) :: ziontypat(ntypat)

! local vars
 integer :: iatom, ispden
 character(len=500) :: message

 real(dp) :: center_of_charge(3)
 real(dp), allocatable :: dipole_el(:,:)
 real(dp) :: dipole_ions(3), ziontotal, dipole_tot(3)

! *************************************************************************

!nuclear part of dipole
 dipole_ions = zero
 ziontotal = zero
 do iatom = 1, natom
   dipole_ions = dipole_ions + xred(:,iatom)*ziontypat(typat(iatom))
   ziontotal = ziontotal + ziontypat(typat(iatom))
 end do

!find coordinates of center of charge on fft grid
!NOTE: wrt center of charge, dipole_ions is 0
 center_of_charge(:) = dipole_ions(:)/ziontotal * ngfft(:)
 write (message, '(a,3E20.10)') ' Center of charge for ionic distribution (red coordinates): ', dipole_ions(:)/ziontotal
 call wrtout(std_out, message, 'COLL')

!get electronic part of dipole with respect to center of charge of ions
 ABI_ALLOCATE(dipole_el,(3,nspden))

 !call multipoles_fftr(arraysp,dipole_el,mpi_enreg,nfft,ngfft,nspden,rprimd,center_of_charge)
 call quadrupole_fftr(arraysp,dipole_el,mpi_enreg,nfft,ngfft,nspden,rprimd,center_of_charge)
 dipole_el = dipole_el * ucvol

 dipole_tot(1) = -sum(dipole_el(1,:))
 dipole_tot(2) = -sum(dipole_el(2,:))
 dipole_tot(3) = -sum(dipole_el(3,:))

!output
 write (message, '(4a)') ch10,' Dipole in cartesian coord (atomic units charge*distance)',&
& ch10, ' with respect to the center of ionic charge'
 call wrtout(std_out, message, 'COLL')
 write (message, '(a)') ' Electronic part (absolute value): '
 call wrtout(std_out, message, 'COLL')
 do ispden = 1, nspden
   write (message, '(a,I5,a,3E18.6)') '   density component ', ispden, ' dipole = ',  dipole_el(:,ispden)
   call wrtout(std_out, message, 'COLL')
 end do

 write (message, '(a,3E18.6,a)') ' Total dipole =        ',  dipole_tot, ch10
 call wrtout(std_out, message, 'COLL')

 ABI_DEALLOCATE(dipole_el)
 
end subroutine quadrupole_tensor_out
!KYSC END
