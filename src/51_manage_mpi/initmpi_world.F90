!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_world
!! NAME
!!  initmpi_world
!!
!! FUNCTION
!!  Initializes the mpi information for world.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2014 ABINIT group (FJ, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!
!!  mpi_enreg=informations about MPI parallelization
!!
!!
!! SIDE EFFECTS
!!  xmpi_world is redifined for the number of processors on which ABINIT is launched
!!
!! PARENTS
!!      finddistrproc
!!
!! CHILDREN
!!      abi_io_redirect
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_world(mpi_enreg,nproc)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_world'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in)::nproc
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays
 integer,allocatable :: ranks(:)

! ***********************************************************************

 DBG_ENTER("COLL")

 if(nproc==mpi_enreg%nproc) return

 ABI_ALLOCATE(ranks,(0:nproc-1))
 ranks(0:nproc-1)=(/((ii),ii=0,nproc-1)/)
 mpi_enreg%comm_world=xmpi_subcomm(xmpi_world,nproc,ranks)
 ABI_DEALLOCATE(ranks)

 if(mpi_enreg%me<nproc)  then
   mpi_enreg%me=xcomm_rank(mpi_enreg%comm_world)
   mpi_enreg%nproc=xcomm_size(mpi_enreg%comm_world)
   call abi_io_redirect(new_io_comm=mpi_enreg%comm_world,new_leave_comm=mpi_enreg%comm_world)
 else
   mpi_enreg%me=-1
 end if

 DBG_EXIT("COLL")

end subroutine initmpi_world
!!***
