!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_clib
!! NAME
!! m_clib
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2014 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_clib

#ifdef HAVE_FC_ISO_C_BINDING
#define USE_MODULE 
#else
#define USE_MODULE use m_iso_c_binding
#endif

 use iso_c_binding

 implicit none

 private 

 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0_dp,1.0_dp))  ! Complex should not be used presently

 type,public :: Mallinfo_t
   integer(C_LONG) :: arena 
   integer(C_LONG) :: hblkhd 
   integer(C_LONG) :: usmblks 
   integer(C_LONG) :: fsmblks 
   integer(C_LONG) :: uordblks 
   integer(C_LONG) :: fordblks
 end type Mallinfo_t

!FIXME the interfaces below have been commented out since abilint
! crashes during the analysis of the file (maybe due to the macro USE_MODULE!)      

! ===================================================
! ==== Fortran-bindings declared in intrinsics.c ====
! ===================================================
! interface 
!   subroutine clib_fflush()
!     import 
!     implicit none
!   end subroutine clib_fflush
! end interface
!
! interface 
!   subroutine clib_getenv(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_getenv
! end interface
!
!! ===================================================
!! ==== Fortran-bindings declared in fsi_posix.c ====
!! ===================================================
! interface 
!   subroutine clib_mkdir(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_mkdir
! end interface
!
! interface 
!   subroutine clib_chdir(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_chdir
! end interface
!
! interface 
!   subroutine clib_rename(ierr,from_fname,to_fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: from_fname,to_fname
!   end subroutine clib_rename
! end interface
!
! interface 
!   subroutine clib_remove(ierr,fname)
!     import 
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_remove
! end interface
!
! interface 
!   subroutine clib_getcwd(ierr,fname)
!     import 
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_getcwd
! end interface
!
! interface 
!   subroutine clib_gethname(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_gethname
! end interface
!
!! =====================================================
!! ==== Fortran-bindings declared in progress_bar.c ====
!! =====================================================
! interface
!   subroutine clib_progress_bar(actual, max)
!     import
!     implicit none
!     integer(C_INT),intent(in) :: actual
!     integer(C_INT),intent(in) :: max
!   end subroutine clib_progress_bar
! end interface
!
!! =================================================
!! ==== Fortran-bindings declared in mallinfo.c ====
!! =================================================
! interface
!   subroutine clib_mallinfo(arena, hblkhd, usmblks, fsmblks, uordblks, fordblks)
!     import
!     implicit none
!     integer(C_LONG),intent(out) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
!   end subroutine clib_mallinfo
! end interface
!
!
!! ==================================================
!! ==== Fortran-bindings declared in gnu_tools.c ====
!! ==================================================
!
! interface
!   subroutine clib_mtrace(ierr)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_mtrace
! end interface
!
! interface
!   subroutine clib_muntrace(ierr)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_muntrace
! end interface
!
! interface
!   subroutine clib_mcheck(ierr)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_mcheck
! end interface

CONTAINS  !===========================================================
!!***

!!****f* m_clib/fmallinfo
!! NAME
!!   fmallinfo
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      clib_alignment_of
!!
!! SOURCE

subroutine fmallinfo(Minfo)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fmallinfo'
!End of the abilint section

 type(Mallinfo_t),intent(out) :: Minfo

!Local variables-------------------------------
 integer(C_LONG) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
! *********************************************************************

  call clib_mallinfo(arena,hblkhd,usmblks,fsmblks,uordblks,fordblks) 

  Minfo%arena    = arena
  Minfo%hblkhd   = hblkhd
  Minfo%usmblks  = usmblks
  Minfo%fsmblks  = fsmblks
  Minfo%uordblks = uordblks
  Minfo%fordblks = fordblks

end subroutine fmallinfo 
!!***

!!****f* m_clib/clib_print_mallinfo
!! NAME
!!   fmallinfo
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine clib_print_mallinfo(Minfo,unt)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clib_print_mallinfo'
!End of the abilint section

 integer,intent(in) :: unt
 type(Mallinfo_t),intent(in) :: Minfo
! *********************************************************************

 write(unt,*)' Total space in arena            : ',Minfo%arena
 write(unt,*)' Space in holding block headers  : ',Minfo%hblkhd
 write(unt,*)' Space in small blocks in use    : ',Minfo%usmblks
 write(unt,*)' Space in free small blocks      : ',Minfo%fsmblks
 write(unt,*)' Space in ordinary blocks in use : ',Minfo%uordblks
 write(unt,*)' Space in free ordinary blocks   : ',Minfo%fordblks
 write(unt,*)' End memory statistics '

end subroutine clib_print_mallinfo
!!***

!----------------------------------------------------------------------


!!****f* m_clib/clib_show_fc_alignment
!! NAME
!!   fmallinfo
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine clib_show_fc_alignment(unt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clib_show_fc_alignment'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unt

!Local variables-------------------------------
 integer,parameter :: sz=1024
 integer :: algn
!arrays
 real(dp),allocatable :: arr(:)

! *************************************************************************

 !allocate(arr(sz))
 call clib_alignment_of(arr, 16, algn)
 write(unt,*)" dp_arr: p % 16 = ",algn
 call clib_alignment_of(arr, 64, algn)
 write(unt,*)" dp_arr: p % 64 = ",algn
 !deallocate(arr)

end subroutine clib_show_fc_alignment
!!***

END MODULE m_clib
!!***

!!  !!****f* m_clib/rename
!!  !! NAME
!!  !!  rename
!!  !!
!!  !! FUNCTION
!!  !!
!!  !! INPUTS 
!!  !!
!!  !! OUTPUT
!!  !!
!!  !! PARENTS
!!  !!
!!  !! CHILDREN
!!  !!
!!  !! SOURCE
!!  
!!  !subroutine rename(old_fname, new_fname, ierr)
!!  !
!!  !!Arguments ------------------------------------
!!  ! character(len=*),intent(in) :: old_fname, new_fname
!!  !
!!  !!Local variables-------------------------------
!!  !
!!  !! *********************************************************************
!!  !
!!  !#if def HAVE_CLIB
!!  ! call clib_rename(old, new)
!!  !
!!  !#elif defined HAVE_FC_COMMAND_LINE
!!  ! call command_line("mv old new")
!!  !
!!  !#elif defined HAVE_FC_SYSTEM
!!  ! call system("mv old new")
!!  !
!!  !#else
!!  !#error "get a decent compiler (e.g. gfortran) or a decent computer!"
!!  !#endif
!!  !
!!  !end subroutine rename
!!  !!***
