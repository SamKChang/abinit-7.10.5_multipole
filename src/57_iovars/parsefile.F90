!!****f* ABINIT/parsefile
!! NAME
!! parsefile
!!
!! FUNCTION
!!  Glue function, to read the given file, put it into a string,
!!  change everything to uppercase, remove carriage returns and
!!  non significant blank characters. May also read a XML input
!!  CML file if specified. Finaly read ndtset input variable.
!!
!! INPUTS
!!  filnamin= the file to read
!!  comm=MPI communicator
!!
!! OUTPUT
!!  lenstr= the length of the resulting string.
!!  ndtset= the number of declared datasets.
!!  string= contains on output the content of the file, ready for parsing.
!!
!! PARENTS
!!      abinit,m_ab7_invars_f90,ujdet
!!
!! CHILDREN
!!      importcml,importxyz,instrng,intagm,inupper,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine parsefile(filnamin,lenstr,ndtset,string,comm)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'parsefile'
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_47_xml
 use interfaces_57_iovars, except_this_one => parsefile
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: filnamin
 integer,intent(in) :: comm
 integer,intent(out) :: ndtset,lenstr
 character(len=strlen),intent(out) :: string

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: option,marr,tread,lenstr_nocml, lenstr_noxyz,ierr
 character(len=strlen) :: string_raw
 character(len=500) :: message
!arrays
 integer :: intarr(1)
 real(dp) :: dprarr(1)

! *************************************************************************

 ! Read the input file, and store the information in a long string of characters
 ! Note: this is done only by me=0, and then string and other output vars are BCASTED

 if (xcomm_rank(comm) == master) then
   !strlen from defs_basis module
   option=1
   call instrng (filnamin,lenstr,option,strlen,string)

   ! Copy original file, without change of case
   string_raw=string

   ! To make case-insensitive, map characters of string to upper case:
   call inupper(string(1:lenstr))

   lenstr_nocml = lenstr
   ! Might import data from CML file(s) into string
   ! Need string_raw to deal properly with CML filenames
   call importcml(lenstr,string_raw,string,strlen)

   lenstr_noxyz = lenstr
   call importxyz(lenstr,string_raw,string,strlen)
   if (abs(lenstr_noxyz-lenstr_nocml) > 0 .and. abs(lenstr-lenstr_noxyz) > 0) then
     write(message, '(5a)' )&
&     'Input was taken from CML and xyz files in different datasets ',ch10,&
&     'This is not allowed.  ',ch10,&
&     'Action : only use one of cmlfile _or_ xyzfile in the input file.'
     MSG_ERROR(message)
   end if
   
   ! Ensure we do not mix cml and xyz, which conflict for determination of znucl and ntypat
   
   !6) Take ndtset from the input string
   ndtset=0; marr=1
   call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),"ndtset",tread,'INT')
   if (tread==1) ndtset=intarr(1)
   ! Check that ndtset is not negative
   if (ndtset<0 .or. ndtset>9999) then
     write(message, '(a,i0,a,a,a,a)' )&
&     'Input ndtset must be non-negative and < 10000, but was ',ndtset,ch10,&
&     'This is not allowed.  ',ch10,&
&     'Action : modify ndtset in the input file.'
     MSG_ERROR(message)
   end if
 end if ! master 

 if (xcomm_size(comm) > 1) then
   ! Broadcast data.
   call xmpi_bcast(lenstr,master,comm,ierr)
   call xmpi_bcast(ndtset,master,comm,ierr)
   call xmpi_bcast(string,master,comm,ierr)
 end if 

end subroutine parsefile
!!***
