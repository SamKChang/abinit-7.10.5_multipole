!{\src2tex{textfont=tt}}
!!****f* ABINIT/create_nc_file
!! NAME
!! create_nc_file
!!
!! FUNCTION
!! Create an NetCDF file including a dimension one definition
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!      etsf_io_low_open_create,etsf_io_low_write_dim
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine create_nc_file (filename,ncid)

 use defs_basis
 use m_profiling_abi
 use m_errors

#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'create_nc_file'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ncid
!arrays
character(len=*),intent(in) :: filename

!Local variables-------------------------------
#if defined HAVE_TRIO_NETCDF
integer :: one_id
logical :: lstat
#if defined HAVE_TRIO_ETSF_IO
type(etsf_io_low_error) :: Error
#else
integer :: ncerr
#endif
#endif

! *************************************************************************

!Create the NetCDF file

#if defined HAVE_TRIO_ETSF_IO
 call etsf_io_low_open_create(ncid, filename, etsf_file_format_version, &
& lstat, error_data= Error, with_etsf_header=.FALSE., overwrite=.TRUE.)
 call etsf_io_low_write_dim(ncid,'one',1,lstat,one_id,Error)

#elif defined HAVE_TRIO_NETCDF
 ncerr=nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=ncid)
 NCF_CHECK(ncerr,'nf90_create')
 ncerr=nf90_def_dim(ncid,'one',1,one_id)
 NCF_CHECK(ncerr,'nf90_def_dim')
 lstat=.true.

#endif

 end subroutine create_nc_file
!!***
