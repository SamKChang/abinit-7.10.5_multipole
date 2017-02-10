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
subroutine quadrupole_fftr(arraysp,quadrupole,mpi_enreg,nfft,ngfft,nspden,rprimd,neworigin,ucvol)

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
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: neworigin(3)
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: quadrupole(6,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft1, ifft2, ifft3, ifft, ispden, nfftot
 real(dp) :: r_size, r_x, r_y, r_z
 character(len=500) :: message
 real(dp) :: invn1, invn2, invn3
!arrays
 real(dp),allocatable :: wrapfft1(:)
 real(dp),allocatable :: wrapfft2(:)
 real(dp),allocatable :: wrapfft3(:)
 real(dp) :: meansp_xx(nspden)
 real(dp) :: meansp_yy(nspden)
 real(dp) :: meansp_zz(nspden)
 real(dp) :: meansp_xy(nspden)
 real(dp) :: meansp_xz(nspden)
 real(dp) :: meansp_yz(nspden)
 real(dp), allocatable :: tmpsp_xx(:,:)
 real(dp), allocatable :: tmpsp_yy(:,:)
 real(dp), allocatable :: tmpsp_zz(:,:)
 real(dp), allocatable :: tmpsp_xy(:,:)
 real(dp), allocatable :: tmpsp_xz(:,:)
 real(dp), allocatable :: tmpsp_yz(:,:)

 print *, "Quadrupole tensor:"

 ABI_ALLOCATE(tmpsp_xx,(nfft,nspden))
 ABI_ALLOCATE(tmpsp_yy,(nfft,nspden))
 ABI_ALLOCATE(tmpsp_zz,(nfft,nspden))
 ABI_ALLOCATE(tmpsp_xy,(nfft,nspden))
 ABI_ALLOCATE(tmpsp_xz,(nfft,nspden))
 ABI_ALLOCATE(tmpsp_yz,(nfft,nspden))

 nfftot = ngfft(1)*ngfft(2)*ngfft(3)
 invn1 = one / ngfft(1) ! to get r_x
 invn2 = one / ngfft(2) ! to get r_y
 invn3 = one / ngfft(3) ! to get r_z

!for the moment impose no fft parallelization
!FIXME: needs to know somehow which fft points (i1i2i3 triplets) are on this processor...
 if (nfft /= nfftot) then
   write (message,'(3a)') 'Error: fft parallelization of multipoles_fftr is not coded yet.',ch10,&
&   ' return from routine and continue as if nothing happened'
   call wrtout(std_out, message, 'COLL')
   return
 end if

 ABI_ALLOCATE(wrapfft1,(ngfft(1)))
 ABI_ALLOCATE(wrapfft2,(ngfft(2)))
 ABI_ALLOCATE(wrapfft3,(ngfft(3)))
 do ifft1 = 1, ngfft(1)
   wrapfft1(ifft1) = mod((ifft1-1-neworigin(1)+half*ngfft(1)), dble(ngfft(1))) - half*ngfft(1)
 end do
 do ifft2 = 1, ngfft(2)
   wrapfft2(ifft2) = mod((ifft2-1-neworigin(2)+half*ngfft(2)), dble(ngfft(2))) - half*ngfft(2)
 end do
 do ifft3 = 1, ngfft(3)
   wrapfft3(ifft3) = mod((ifft3-1-neworigin(3)+half*ngfft(3)), dble(ngfft(3))) - half*ngfft(3)
 end do

!!!!! Qxx !!!!!
 ifft = 1
 do ifft3 = 1, ngfft(3)
   do ifft2 = 1, ngfft(2)
     do ifft1 = 1, ngfft(1)

       r_x = (wrapfft1(ifft1) * invn1)
       r_y = (wrapfft2(ifft2) * invn2)
       r_z = (wrapfft3(ifft3) * invn3)
       r_size = (r_x**2 + r_y**2 + r_z**2)

       tmpsp_xx(ifft,:) = arraysp(ifft,:) * (3 * r_x**2 - r_size)
       tmpsp_yy(ifft,:) = arraysp(ifft,:) * (3 * r_y**2 - r_size)
       tmpsp_zz(ifft,:) = arraysp(ifft,:) * (3 * r_z**2 - r_size)

       tmpsp_xy(ifft,:) = arraysp(ifft,:) * (3 * r_x * r_y)
       tmpsp_xz(ifft,:) = arraysp(ifft,:) * (3 * r_x * r_z)
       tmpsp_yz(ifft,:) = arraysp(ifft,:) * (3 * r_y * r_z)

       ifft = ifft + 1
     end do
   end do
 end do

 call mean_fftr(tmpsp_xx,meansp_xx,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 call mean_fftr(tmpsp_yy,meansp_yy,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 call mean_fftr(tmpsp_zz,meansp_zz,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 call mean_fftr(tmpsp_xy,meansp_xy,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 call mean_fftr(tmpsp_xz,meansp_xz,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 call mean_fftr(tmpsp_yz,meansp_yz,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)

 meansp_xx = meansp_xx * ucvol 
 meansp_yy = meansp_yy * ucvol
 meansp_zz = meansp_zz * ucvol
 meansp_xy = meansp_xy * ucvol
 meansp_xz = meansp_xz * ucvol
 meansp_yz = meansp_yz * ucvol

 print *, "quadrupole tensor component Qxx: ", meansp_xx
 print *, "quadrupole tensor component Qyy: ", meansp_yy
 print *, "quadrupole tensor component Qzz: ", meansp_zz
 print *, "quadrupole tensor component Qxy: ", meansp_xy
 print *, "quadrupole tensor component Qxz: ", meansp_xz
 print *, "quadrupole tensor component Qyz: ", meansp_yz
 quadrupole(1,:) = meansp_xx
 quadrupole(2,:) = meansp_yy
 quadrupole(3,:) = meansp_zz
 quadrupole(4,:) = meansp_xy
 quadrupole(5,:) = meansp_xz
 quadrupole(6,:) = meansp_yz

 ABI_DEALLOCATE(wrapfft1)
 ABI_DEALLOCATE(wrapfft2)
 ABI_DEALLOCATE(wrapfft3)

 ABI_DEALLOCATE(tmpsp_xx)
 ABI_DEALLOCATE(tmpsp_yy)
 ABI_DEALLOCATE(tmpsp_zz)
 ABI_DEALLOCATE(tmpsp_xy)
 ABI_DEALLOCATE(tmpsp_xz)
 ABI_DEALLOCATE(tmpsp_yz)
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
 real(dp), allocatable :: quadrupole(:,:)
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
 ABI_ALLOCATE(quadrupole,(6,nspden))

 call quadrupole_fftr(arraysp,quadrupole,mpi_enreg,nfft,ngfft,nspden,rprimd,center_of_charge,ucvol)
 !dipole_el = dipole_el * ucvol
 quadrupole = quadrupole * ucvol
 print *, "multiplicative factor", ucvol

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
