!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkorthsy
!! NAME chkorthsy
!! chkorthsy
!!
!! FUNCTION
!! Check the orthogonality of the symmetry operations
!! (lengths and absolute values of scalar products should be preserved)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gprimd(3,3)=dimensional primitive transl. for reciprocal space (bohr**-1)
!! rmet=Real space metric.
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symrel(3,3,1:nsym)=symmetry operations in real space in terms of primitive translations
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      chkinp,ingeo
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkorthsy(gprimd,nsym,rmet,rprimd,symrel)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkorthsy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rmet(3,3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,jj
 real(dp),parameter :: tol=1.0d-11
 real(dp) :: residual
 character(len=500) :: message
!arrays
 real(dp) :: prods(3,3),rmet_sym(3,3),rprimd_sym(3,3)

! *************************************************************************

!Loop over all symmetry operations
 do isym=1,nsym

!  Compute symmetric of primitive vectors under point symmetry operations
!  The symrel array might have to be transposed ???
   do ii=1,3
     rprimd_sym(:,ii)=symrel(1,ii,isym)*rprimd(:,1)+&
&     symrel(2,ii,isym)*rprimd(:,2)+&
&     symrel(3,ii,isym)*rprimd(:,3)
   end do

!  DEBUG
!  write(std_out,*)' chkorthsy : isym=',isym
!  write(std_out,*)rprimd(1:3,1)
!  write(std_out,*)rprimd(1:3,2)
!  write(std_out,*)rprimd(1:3,3)
!  write(std_out,*)rprimd_sym(1:3,1)
!  write(std_out,*)rprimd_sym(1:3,2)
!  write(std_out,*)rprimd_sym(1:3,3)
!  ENDDEBUG

!  If the new lattice is the same as the original one,
!  the lengths and angles are preserved
   do ii=1,3
     rmet_sym(ii,:)=rprimd_sym(1,ii)*rprimd_sym(1,:)+&
&     rprimd_sym(2,ii)*rprimd_sym(2,:)+&
&     rprimd_sym(3,ii)*rprimd_sym(3,:)
   end do

!  DEBUG
!  write(std_out,*)' rmet :'
!  write(std_out,*)rmet(1:3,1)
!  write(std_out,*)rmet(1:3,2)
!  write(std_out,*)rmet(1:3,3)
!  write(std_out,*)rmet_sym(1:3,1)
!  write(std_out,*)rmet_sym(1:3,2)
!  write(std_out,*)rmet_sym(1:3,3)
!  ENDDEBUG

   residual=0.0d0
   do ii=1,3
     do jj=1,3
       residual=residual+(abs(rmet_sym(ii,jj))-abs(rmet(ii,jj)))**2
     end do
   end do
   if(residual>tol)then
     write(message, '(a,i5,a,a,a,a,a,es12.4,a,a,a,a,a,a,a)' )&
&     'The symmetry operation number',isym,' does not preserve',ch10,&
&     'vector lengths and angles.',ch10,&
&     'The value of the residual is',residual,'.',ch10,&
&     'Action : modify rprim, acell and/or symrel so that',ch10,&
&     'vector lengths and angles are preserved.',ch10,&
&     'Beware, the tolerance on symmetry operations is very small.'
     MSG_ERROR(message)
!    Should use ierr
   end if

!  Also, the scalar product of rprimd_sym and gprimd must give integer numbers
   do ii=1,3
     prods(ii,:)=rprimd_sym(1,ii)*gprimd(1,:)+ &
&     rprimd_sym(2,ii)*gprimd(2,:)+ &
&     rprimd_sym(3,ii)*gprimd(3,:)
   end do

   do ii=1,3
     do jj=1,3
       residual=prods(ii,jj)-anint(prods(ii,jj))
       if(abs(residual)>tol)then
         write(message, '(a,i0,a,a,a,a,a,a,a)' )&
&         'The symmetry operation number',isym,' generates',ch10,&
&         'a different lattice.',ch10,&
&         'Action : modify rprim, acell and/or symrel so that',ch10,&
&         'the lattice is preserved.'
         MSG_ERROR(message)
!        Should use ierr
       end if
     end do
   end do
 end do

end subroutine chkorthsy
!!***
