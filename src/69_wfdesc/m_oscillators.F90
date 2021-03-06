!{\src2tex{textfont=tt}}
!****m* ABINIT/m_oscillators
!! NAME
!!  m_oscillators
!!
!! FUNCTION
!!  This module contains procedures to calculate the oscillator matrix elements used in the GW code.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2014 ABINIT group (MG)
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

MODULE m_oscillators

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_fft     

 use m_gwdefs,    only : czero_gw
 use m_fstrings,  only : toupper
 use m_blas,      only : xcopy
 use m_gsphere,   only : gsphere_t

 implicit none

 private 
!!***

!----------------------------------------------------------------------

 public :: rho_tw_g 
 public :: get_uug
 public :: calc_wfwfg         ! Calculate the Fourier transform of the product u_{bk}^*(r).u_{b"k}(r) at an arbitrary k in the BZ.
 public :: sym_rhotwgq0       ! Symmetrize the oscillator matrix elements in the BZ in the special case of q=0.
!!***

!----------------------------------------------------------------------

CONTAINS  !=========================================================================================================================

!!****f* m_oscillators/rho_tw_g
!! NAME
!! rho_tw_g
!!
!! FUNCTION
!! Calculate rhotwg(G)=<wfn1|exp(-i(q+G).r)|wfn2>
!!
!! INPUTS
!! dim_rtwg=Define the size of the output array rhotwg
!!   === for nspinor==1 ===
!!    dim_rtwg=1
!!   === for nspinor==2 ===
!!    dim_rtwg=2 if only <up|up>, <dwn|dwn> matrix elements are required
!!    dim_rtwg=4 for <up|up>, <dwn|dwn>, <up|dwn> and <dwn|up>.
!! map2sphere= 1 to retrieve Fourier components indexed according to igfftg0.
!!             0 to retrieve Fourier components indexed according to the FFT box.
!!               NOTE: If map2sphere==0 npwvec must be equal to nr
!! use_padfft= Only compatible with map2sphere 1. 
!!             1 if matrix elements are calculated via zero-padded FFT.
!!             0 R-->G Transform in done on the full FFT box.
!! igfftg0(npwvec*map2sphere)=index of G-G_o in the FFT array for each G in the sphere.
!! i1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! i2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau} 
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npwvec=number of plane waves (in the sphere if map2sphere==1, in the FFT box if map2sphere==1)
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! nspinor=number of spinorial components.
!! spinrot1(4),spinrot2(4)=components of the spinor rotation matrix :
!!  spinrot(1)=$\cos\phi/2$
!!  spinrot(2)=$\sin\phi/2 \times u_x$
!!  spinrot(3)=$\sin\phi/2 \times u_y$
!!  spinrot(4)=$\sin\phi/2 \times u_z$
!!   where $\phi$ is the angle of rotation, and
!!   $(u_x,u_y,u_z)$ is the normalized direction of the rotation axis
!! wfn1(nr*nspinor*ndat),wfn2(nr*nspinor*ndat)=the two wavefunctions (periodic part)
!! [nhat12(2,nr,nspinor**2*ndat)]=Compensation charge in real space to be added to \Psi_1^*\Psi_2 -- Only for PAW.
!!
!! OUTPUT
!! rhotwg(npwvec)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,cchi0q0_intraband
!!      check_completeness,cohsex_me,exc_build_block,exc_build_ham,m_eet
!!      m_fft_prof,prep_calc_ucrpa
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many
!!
!! SOURCE

subroutine rho_tw_g(nspinor,npwvec,nr,ndat,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
& wfn1,i1,ktabr1,ktabp1,spinrot1,&
& wfn2,i2,ktabr2,ktabp2,spinrot2,&
& dim_rtwg,rhotwg) !& nhat12)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rho_tw_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: i1,i2,npwvec,nr,nspinor,dim_rtwg,map2sphere,use_padfft,ndat
 complex(dpc),intent(in) :: ktabp1,ktabp2
!arrays
 integer,intent(in) :: gbound(:,:) !gbound(2*mgfft+8,2)
 integer,intent(in) :: igfftg0(npwvec*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 real(dp),intent(in) :: spinrot1(4),spinrot2(4)
 complex(gwpc),intent(in) :: wfn1(nr*nspinor*ndat),wfn2(nr*nspinor*ndat)
 complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg*ndat)
! real(dp),optional,intent(in) :: nhat12(2,nr,nspinor**2*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: nd1=1
 integer :: ig,ir,ir1,ir2,igfft,iab,spad1,spad2,spad0,ispinor,dat
 integer :: upt,rhopt
 integer :: nx,ny,nz,ldx,ldy,ldz,mgfft
 complex(gwpc) :: u1a,u1b,u2a,u2b
 type(fftbox_plan3_t) :: plan
!arrays
 integer :: spinor_pad(2,4) 
 complex(dpc) :: spinrot_mat1(2,2),spinrot_mat2(2,2)
 complex(gwpc),allocatable :: u12prod(:),cwavef1(:),cwavef2(:)

! *************************************************************************

 SELECT CASE (nspinor)

 CASE (1) ! Collinear case.

#if 1
   call ts_usug_kkp_bz(npwvec,nr,ndat,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
&    wfn1,i1,ktabr1,ktabp1,&
&    wfn2,i2,ktabr2,ktabp2,rhotwg)
#else
!$OMP PARALLEL DO PRIVATE(upt,rhopt) SCHEDULE(DYNAMIC) 
   do dat=1,ndat
     upt   = 1 + (dat-1)*nr
     rhopt = 1 + (dat-1)*npwvec
     call ts_usug_kkp_bz(npwvec,nr,nd1,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
&      wfn1(upt),i1,ktabr1,ktabp1,&
&      wfn2(upt),i2,ktabr2,ktabp2,rhotwg(rhopt))
   end do
#endif

 CASE (2) ! Spinorial case.

   MSG_ERROR("Add zero-padded FFT") !TODO
   ABI_CHECK(ndat==1,"ndat != 1 not coded")
   dat = 1
   upt = 1
   rhopt = 1

   ABI_MALLOC(cwavef1,(nr*nspinor))
   ABI_MALLOC(cwavef2,(nr*nspinor))

   ! === Apply Time-reversal if required ===
   ! \psi_{-k}^1 =  (\psi_k^2)^*
   ! \psi_{-k}^2 = -(\psi_k^1)^*
   if (i1==1) then 
     cwavef1(:)=wfn1(:) 
   else if (i1==2) then
     cwavef1(1:nr)     = GWPC_CONJG(wfn1(nr+1:2*nr))
     cwavef1(nr+1:2*nr)=-GWPC_CONJG(wfn1(1:nr))
   else 
     MSG_ERROR('Wrong i1 in spinor')
   end if

   if (i2==1) then 
     cwavef2(:)=wfn2(:) 
   else if (i2==2) then
     cwavef2(1:nr)     = GWPC_CONJG(wfn2(nr+1:2*nr))
     cwavef2(nr+1:2*nr)=-GWPC_CONJG(wfn2(1:nr))
   else 
     MSG_ERROR('Wrong i2 in spinor')
   end if

   ! === Rotate wavefunctions in r-space ===
   do ispinor=1,nspinor
     spad0=(ispinor-1)*nr
     do ir=1,nr
       ir1=ktabr1(ir); ir2=ktabr2(ir)
       cwavef1(ir+spad0) = cwavef1(ir1+spad0)*ktabp1
       cwavef2(ir+spad0) = cwavef2(ir2+spad0)*ktabp2
     end do 
   end do !ispinor

   ! === Rotation in spinor space ===
   !spinrots1=spinrot1(1) ; spinrots2=spinrot2(1)
   !spinrotx1=spinrot1(2) ; spinrotx2=spinrot2(2)
   !spinroty1=spinrot1(3) ; spinroty2=spinrot2(3)
   !spinrotz1=spinrot1(4) ; spinrotz2=spinrot2(4)
   spinrot_mat1(1,1)= spinrot1(1) + j_dpc*spinrot1(4)
   spinrot_mat1(1,2)= spinrot1(3) + j_dpc*spinrot1(2)
   spinrot_mat1(2,1)=-spinrot1(3) + j_dpc*spinrot1(2)
   spinrot_mat1(2,2)= spinrot1(1) - j_dpc*spinrot1(4)

   spinrot_mat2(1,1)= spinrot2(1) + j_dpc*spinrot2(4)
   spinrot_mat2(1,2)= spinrot2(3) + j_dpc*spinrot2(2)
   spinrot_mat2(2,1)=-spinrot2(3) + j_dpc*spinrot2(2)
   spinrot_mat2(2,2)= spinrot2(1) - j_dpc*spinrot2(4)

   do ir=1,nr
     !ar=wavefspinor(1,ir)
     !ai=wavefspinor(2,ir)
     !br=wavefspinor(1,npw2+ir)
     !bi=wavefspinor(2,npw2+ir)
     u1a=cwavef1(ir) 
     u1b=cwavef1(ir+nr) 
     cwavef1(ir)   =spinrot_mat1(1,1)*u1a+spinrot_mat1(1,2)*u1b
     cwavef1(ir+nr)=spinrot_mat1(2,1)*u1a+spinrot_mat1(2,2)*u1b
     u2a=cwavef2(ir) 
     u2b=cwavef2(ir+nr) 
     cwavef2(ir)   =spinrot_mat2(1,1)*u2a+spinrot_mat2(1,2)*u2b
     cwavef2(ir+nr)=spinrot_mat2(2,1)*u2a+spinrot_mat2(2,2)*u2b
     !wavefspinor(1,ir)     = spinrots*ar-spinrotz*ai +spinroty*br-spinrotx*bi
     !wavefspinor(2,ir)     = spinrots*ai+spinrotz*ar +spinroty*bi+spinrotx*br
     !wavefspinor(1,npw2+ir)=-spinroty*ar-spinrotx*ai +spinrots*br+spinrotz*bi
     !wavefspinor(2,npw2+ir)=-spinroty*ai+spinrotx*ar +spinrots*bi-spinrotz*br
   end do

   spinor_pad(:,:)=RESHAPE((/0,0,nr,nr,0,nr,nr,0/),(/2,4/))

   do iab=1,dim_rtwg
     spad1=spinor_pad(1,iab)
     spad2=spinor_pad(2,iab)

     u12prod = GWPC_CONJG(cwavef1(spad1+1:spad1+nr)) * cwavef2(spad2+1:spad2+nr)

     ! Add compensation charge.
     !if (PRESENT(nhat12)) then 
     !  u12prod = u12prod + CMPLX(nhat12(1,:,iab),nhat12(2,:,iab))
     !end if
     spad0=(iab-1)*npwvec

     SELECT CASE (map2sphere)
     
     CASE (0) ! Need results on the full FFT box thus cannot use zero-padded FFT.

       call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
       call fftbox_execute(plan,u12prod)
       rhotwg(spad0+1:spad0+npwvec)=u12prod 

     CASE (1) ! Need results on the G-sphere. Call zero-padded FFT routines if required.

       if (use_padfft==1) then
         nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
         ldx=nx      ; ldy=ny      ; ldz=nz
         call fftpad(u12prod,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gbound)
       else 
         call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
         call fftbox_execute(plan,u12prod)
       end if

       do ig=1,npwvec      ! Have to map FFT to G-sphere.
         igfft=igfftg0(ig)
         if (igfft/=0) then ! G-G0 belong to the FFT mesh.
           rhotwg(ig+spad0)=u12prod(igfft) 
         else               ! Set this component to zero
           rhotwg(ig+spad0)=czero_gw
         end if
       end do

     CASE DEFAULT
       MSG_BUG("Wrong map2sphere")
     END SELECT
   end do !iab

   ABI_FREE(cwavef1)
   ABI_FREE(cwavef2)

 CASE DEFAULT 
   MSG_BUG('Wrong nspinor')
 END SELECT

end subroutine rho_tw_g
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/get_uug
!! NAME
!! get_uug
!!
!! FUNCTION
!! Calculate usug(G)= \intl \dr u1(r) exp(-i(q+G).r) u2(r) for ndat pair of wavefunctions
!!
!! INPUTS
!! npw=number of plane waves in the sphere.
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! use_padfft= 
!!             1 if matrix elements are calculated via zero-padded FFT.
!!             0 R-->G Transform in done on the full FFT box.
!! igfftg0(npw)=index of G-G_o in the FFT array for each G in the sphere.
!! u1(nr*ndat),u2(nr*ndat)=the two wavefunctions (periodic part)
!! [trans1] = "C" if the complex conjugate of u1 should be used. Default is "N"
!! [trans2] = "C" if the complex conjugate of u2 should be used. Default is "N"
!!
!! OUTPUT
!! usug(npw*ndat)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,cohsex_me,exc_build_block
!!      exc_build_ham
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many
!!
!! SOURCE

subroutine get_uug(npw,nr,ndat,ngfft,use_padfft,igfftg0,gbound,u1,u2,usug,trans1,trans2) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_uug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nr,use_padfft,ndat
 character(len=*),optional,intent(in) :: trans1,trans2
!arrays
 integer,intent(in) :: gbound(:,:) !gbound(2*mgfft+8,2)
 integer,intent(in) :: igfftg0(npw),ngfft(18)
 complex(gwpc),intent(in) :: u1(nr*ndat),u2(nr*ndat)
 complex(gwpc),intent(out) :: usug(npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: nx,ny,nz,ldx,ldy,ldz,mgfft
 character(len=1) :: my_trans1,my_trans2
 type(fftbox_plan3_t) :: plan
!arrays
 complex(gwpc),allocatable :: u12prod(:) 

! *************************************************************************

 ABI_MALLOC(u12prod,(nr*ndat))

 my_trans1 = "N"; if (PRESENT(trans1)) my_trans1 = toupper(trans1(1:1))
 my_trans2 = "N"; if (PRESENT(trans2)) my_trans2 = toupper(trans2(1:1))

 if (my_trans1=="N" .and. my_trans2=="N") then
   u12prod = u1 * u2
 else if (my_trans1=="H" .and. my_trans2=="N") then
   u12prod = GWPC_CONJG(u1) * u2
 else if (my_trans1=="N" .and. my_trans2=="H") then
   u12prod = u1 * GWPC_CONJG(u2)
 else if (my_trans1=="H" .and. my_trans2=="H") then
   u12prod = GWPC_CONJG(u1) * GWPC_CONJG(u2)
 else 
   MSG_ERROR("Wrong combination of trans1, trans2")
 end if

 ! Call zero-padded FFT routines if required.
 if (use_padfft==1) then
   nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
   ldx=nx      ; ldy=ny      ; ldz=nz
   call fftpad(u12prod,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gbound)
 else
   call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
   call fftbox_execute(plan,u12prod)
 end if
 !
 ! From the FFT to the G-sphere.
 call gw_box2gsph(nr,ndat,npw,igfftg0,u12prod,usug)

 ABI_FREE(u12prod)

end subroutine get_uug
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/ts_usug_kkp_bz
!! NAME
!! ts_usug_kkp_bz
!!
!! FUNCTION
!! Calculate usug(G)=<u1|exp(-i(q+G).r)|u2> for ndat pair of wavefunctions
!! TODO: The routine is thread-safe hence it can be called within an OMP parallel region.
!!
!! INPUTS
!! map2sphere= 1 to retrieve Fourier components indexed according to igfftg0.
!!             0 to retrieve Fourier components indexed according to the FFT box.
!!               NOTE: If map2sphere==0 npw must be equal to nr
!! use_padfft= Only compatible with map2sphere 1. 
!!             1 if matrix elements are calculated via zero-padded FFT.
!!             0 R-->G Transform in done on the full FFT box.
!! igfftg0(npw*map2sphere)=index of G-G_o in the FFT array for each G in the sphere.
!! time1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! time2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau} 
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npw=number of plane waves (in the sphere if map2sphere==1, in the FFT box if map2sphere==1)
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! u1(nr*ndat),u2(nr*ndat)=the two wavefunctions (periodic part)
!!
!! OUTPUT
!! usug(npw*ndat)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many
!!
!! SOURCE

subroutine ts_usug_kkp_bz(npw,nr,ndat,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
& u1,time1,ktabr1,ktabp1,&
& u2,time2,ktabr2,ktabp2,usug) !& nhat12)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ts_usug_kkp_bz'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: time1,time2,npw,nr,map2sphere,use_padfft,ndat
 complex(dpc),intent(in) :: ktabp1,ktabp2
!arrays
 integer,intent(in) :: gbound(:,:) !gbound(2*mgfft+8,2)
 integer,intent(in) :: igfftg0(npw*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 complex(gwpc),intent(in) :: u1(nr*ndat),u2(nr*ndat)
 complex(gwpc),intent(out) :: usug(npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: nx,ny,nz,ldx,ldy,ldz,mgfft
 type(fftbox_plan3_t) :: plan
!arrays
 complex(gwpc),allocatable :: u12prod(:) 

! *************************************************************************

 ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2(r,b2,kbz2), to account for symmetries:
 ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz) 
 !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 !
 ABI_MALLOC(u12prod,(nr*ndat))

 call usur_kkp_bz(nr,ndat,time1,ktabr1,ktabp1,u1,time2,ktabr2,ktabp2,u2,u12prod)
 !
 ! Add compensation charge.
 !if (PRESENT(nhat12)) then 
 !  u12prod = u1prod + CMPLX(nhat12(1,:,1),nhat12(2,:,1))
 !end if

 SELECT CASE (map2sphere)
 CASE (0) 
   ! Need results on the full FFT box thus cannot use zero-padded FFT.
   call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
   call fftbox_execute(plan,u12prod)
   call xcopy(nr*ndat,u12prod,1,usug,1)

 CASE (1) 
   ! Need results on the G-sphere. Call zero-padded FFT routines if required.
   if (use_padfft==1) then
     nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
     ldx=nx      ; ldy=ny      ; ldz=nz
     call fftpad(u12prod,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gbound)
   else
     call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
     call fftbox_execute(plan,u12prod)
   end if
   !
   ! From the FFT to the G-sphere.
   call gw_box2gsph(nr,ndat,npw,igfftg0,u12prod,usug)
   !
 CASE DEFAULT
   MSG_BUG("Wrong map2sphere")
 END SELECT

 ABI_FREE(u12prod)

end subroutine ts_usug_kkp_bz
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/usur_kkp_bz
!! NAME
!! usur_kkp_bz
!!
!! FUNCTION
!! Calculate u1_kbz^*(r) u2_kbz(r) in real space from the symmetric images in the IBZ.
!! Does not support spinor wavefunctions.
!!
!! INPUTS
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! u1(nr*ndat),u2(nr*ndat)=the two wavefunctions iin the IBZ (periodic part)
!! time1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! time2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau} 
!!
!! OUTPUT
!!  u12prod(nr*dat) = u1_kbz^*(r) u2_kbz(r) for the ndat pairs.
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many
!!
!! SOURCE

subroutine usur_kkp_bz(nr,ndat,time1,ktabr1,ktabp1,u1,time2,ktabr2,ktabp2,u2,u12prod)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'usur_kkp_bz'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nr,ndat,time1,time2
 complex(dpc),intent(in) :: ktabp1,ktabp2
!arrays
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 complex(gwpc),intent(in) :: u1(nr*ndat),u2(nr*ndat)
 complex(gwpc),intent(out) :: u12prod(nr*ndat)

!Local variables-------------------------------
!scalars
 integer :: ir,dat,padat
 complex(gwpc) :: my_ktabp1,my_ktabp2
!arrays
 complex(gwpc),allocatable :: u1_bz(:),u2_bz(:) 

! *************************************************************************

 !
 ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2(r,b2,kbz2), to account for symmetries:
 ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz) 
 !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 !
 ABI_MALLOC(u1_bz,(nr*ndat))
 ABI_MALLOC(u2_bz,(nr*ndat))

 my_ktabp1 = ktabp1
 my_ktabp2 = ktabp2

 if (ndat==1) then
   do ir=1,nr
     u1_bz(ir) = u1(ktabr1(ir))*my_ktabp1
   end do
   do ir=1,nr
     u2_bz(ir) = u2(ktabr2(ir))*my_ktabp2
   end do
 else
!$OMP PARALLEL PRIVATE(padat)
!$OMP DO
   do dat=1,ndat
     padat = (dat-1)*nr
     do ir=1,nr
       u1_bz(ir+padat) = u1(ktabr1(ir)+padat)*my_ktabp1
     end do
   end do
!$OMP END DO NOWAIT
!$OMP DO
   do dat=1,ndat
     padat = (dat-1)*nr
     do ir=1,nr
       u2_bz(ir+padat) = u2(ktabr2(ir)+padat)*my_ktabp2
     end do
   end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
 end if
 !
 ! Treat time-reversal.
 SELECT CASE (time1)
 CASE (1)
   !
   if (ndat==1) then
     if (time2==1) then
       do ir=1,nr
         u12prod(ir) = GWPC_CONJG(u1_bz(ir)) * u2_bz(ir)
       end do
     else if (time2==2) then
       do ir=1,nr
         u12prod(ir) = GWPC_CONJG(u1_bz(ir)) * GWPC_CONJG(u2_bz(ir))
       end do
     else
       MSG_ERROR("Wrong time2")
     end if
   else 
     if (time2==1) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = GWPC_CONJG(u1_bz(ir+padat)) * u2_bz(ir+padat)
         end do
       end do
     else if (time2==2) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = GWPC_CONJG(u1_bz(ir+padat)) * GWPC_CONJG(u2_bz(ir+padat))
         end do
       end do
     else
       MSG_ERROR("Wrong time2")
     end if
   end if
   !
 CASE (2)
   !
   if (ndat==1) then
     if (time2==1) then
       do ir=1,nr
         u12prod(ir) = u1_bz(ir) * u2_bz(ir)
       end do
     else if (time2==2) then
       do ir=1,nr
         u12prod(ir) = u1_bz(ir) * GWPC_CONJG(u2_bz(ir))
       end do
     else
       MSG_ERROR("Wrong time2")
     end if
   else
     if (time2==1) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = u1_bz(ir+padat) * u2_bz(ir+padat)
         end do
       end do
     else if (time2==2) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = u1_bz(ir+padat) * GWPC_CONJG(u2_bz(ir+padat))
         end do
       end do
     else
       MSG_ERROR("Wrong time2")
     end if
   end if
     !
 CASE DEFAULT
   MSG_ERROR("Wrong time1")
 END SELECT

 ABI_FREE(u1_bz)
 ABI_FREE(u2_bz)

end subroutine usur_kkp_bz
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/gw_box2gsph
!! NAME
!! gw_box2gsph
!!
!! FUNCTION
!! Trasnfer data from the FFT box to the G-sphere.
!!
!! INPUTS
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! npw=number of plane waves in the sphere 
!! igfftg0(npw)=index of G-G_o in the FFT array for each G in the sphere.
!! iarrbox(nr*ndat)=Input array on the FFT mesh
!!
!! OUTPUT
!! oarrsph(npw*ndat)=output array on the sphere.
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many
!!
!! SOURCE

subroutine gw_box2gsph(nr,ndat,npw,igfftg0,iarrbox,oarrsph)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_box2gsph'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nr,ndat,npw
!arrays
 integer,intent(in) :: igfftg0(npw)
 complex(gwpc),intent(in) :: iarrbox(nr*ndat)
 complex(gwpc),intent(out) :: oarrsph(npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: ig,igfft,dat,pgsp,pfft

! *************************************************************************

 if (ndat==1) then
   do ig=1,npw
     igfft=igfftg0(ig)
     if (igfft/=0) then ! G-G0 belong to the FFT mesh.
       oarrsph(ig) = iarrbox(igfft) 
     else               ! Set this component to zero.
       oarrsph(ig) = czero_gw 
     end if
   end do
 else
!$OMP PARALLEL DO PRIVATE(pgsp,pfft,igfft)
   do dat=1,ndat
     pgsp = (dat-1)*npw
     pfft = (dat-1)*nr
     do ig=1,npw
       igfft=igfftg0(ig)
       if (igfft/=0) then ! G-G0 belong to the FFT mesh.
         oarrsph(ig+pgsp) = iarrbox(igfft+pfft) 
       else               ! Set this component to zero.
         oarrsph(ig+pgsp) = czero_gw 
       end if
     end do
   end do
 end if

end subroutine gw_box2gsph
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/calc_wfwfg
!! NAME
!! calc_wfwfg
!!
!! FUNCTION
!!  Calculate the Fourier transform of the product u_{bk}^*(r).u_{b"k}(r) 
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sigc_me,cchi0,cchi0q0,check_completeness,cohsex_me,m_eet
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3_many
!!
!! SOURCE

subroutine calc_wfwfg(ktabr_k,ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_wfwfg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nfftot
!arrays
 integer,intent(in) :: ktabr_k(nfftot),ngfft_gw(18)
 complex(gwpc),intent(in) :: wfr_jb(nfftot),wfr_kb(nfftot)
 complex(gwpc),intent(out) :: wfg2_jk(nfftot)

!Local variables-------------------------------
 integer,parameter :: ndat1=1
 type(fftbox_plan3_t) :: plan
!arrays
 complex(gwpc),allocatable :: wfr2_dpcplx(:)

! *************************************************************************

 ! There is no need to take into account phases arising from non-symmorphic
 ! operations since the wavefunctions are evaluated at the same k-point.
 ABI_MALLOC(wfr2_dpcplx,(nfftot))

 SELECT CASE (ktabi_k)

 CASE (1)
   wfr2_dpcplx = GWPC_CONJG(wfr_jb(ktabr_k)) * wfr_kb(ktabr_k)

 CASE (2) ! Conjugate the product if time-reversal is used to reconstruct this k-point
   wfr2_dpcplx = wfr_jb(ktabr_k) * GWPC_CONJG(wfr_kb(ktabr_k))

 CASE DEFAULT
   MSG_ERROR("Wrong ktabi_k")
 END SELECT

 ! Transform to Fourier space (result in wfg2_jk)
 call fftbox_plan3_many(plan,ndat1,ngfft_gw(1:3),ngfft_gw(1:3),ngfft_gw(7),-1)

 call fftbox_execute(plan,wfr2_dpcplx,wfg2_jk)
 ABI_FREE(wfr2_dpcplx)

end subroutine calc_wfwfg
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/sym_rhotwgq0   
!! NAME
!!  sym_rhotwgq0
!!
!! FUNCTION
!!  Symmetrization of the oscillator matrix elements <k-q,b1|exp(-i(q+G).r)|k,b2> in the special case of q=0.
!!  The matrix elements in the full BZ is obtained from the matrix elements in the IBZ by
!!  rotating the wavefunctions and taking into account time reversal symmetry.
!!  strictly speaking the symmetrization can be performed only for non-degenerate states.
!!
!! INPUTS
!!  Gsph<gsphere_t>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).  
!!  npw=Number of G-vectors
!!  dim_rtwg=Number of spin-spin combinations, 1 for collinear spin, 4 is nspinor==2 (TODO NOT CODED)
!!  itim_k=2 if time reversal is used to reconstruct the k in the BZ, 1 otherwise.
!!  isym_k=The index of the symmetry symrec rotains k_IBZ onto k_BZ.
!!  rhxtwg_in(dim_rtwg*npw)=The input matrix elements in the IBZ.
!!
!! OUTPUT
!!  rhxtwg_sym(dim_rtwg*npw)=The symmetrized matrix elements in the BZ.
!!
!! NOTES
!! Let M_{G}(k,q) =<k-q,b1|exp(-i(q+G).r)|k,b2>
!!  At q ==0, supposing non-degenerate bands, one obtains:
!!
!!  1) M_{ SG}( Sk) = e^{-iSG.t} M_{G}   (k)
!!  2) M_{-SG}(-Sk) = e^{+iSG.t} M_{G}^* (k)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npw,rhxtwg_in,Gsph) result(rhxtwg_sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sym_rhotwgq0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,dim_rtwg,itim_k,isym_k
 type(gsphere_t),intent(in) :: Gsph
!arrays
 complex(gwpc),intent(in) :: rhxtwg_in(dim_rtwg*npw) 
 complex(gwpc) :: rhxtwg_sym(dim_rtwg*npw) 

!Local variables ------------------------------
!scalars
 integer :: ig
 character(len=500) :: msg

!************************************************************************

 if (dim_rtwg/=1) then
  MSG_ERROR("dim_rtwg/=1 not coded")
 end if

 SELECT CASE (isym_k)

 CASE (1) ! Fractional translation associated to E is assumed to be (zero,zero,zero).

   SELECT CASE (itim_k) 
   CASE (1) ! Identity, no time-reversal. No symmetrization is needed.
     rhxtwg_sym(:) = rhxtwg_in(:)

   CASE (2) ! Identity + Time-reversal.
     do ig=1,npw
       rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = GWPC_CONJG(rhxtwg_in(ig))
     end do

   CASE DEFAULT
     write(msg,'(a,i0)')"Wrong value of itim_k: ",itim_k
     MSG_ERROR(msg)
   END SELECT

 CASE DEFAULT ! Rotate wavefunctions.

   SELECT CASE (itim_k) 

   CASE (1) ! no time-reversal, only rotation.
    do ig=1,npw
      rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = rhxtwg_in(ig) * Gsph%phmSGt(ig,isym_k) 
    end do

   CASE (2) ! time-reversal + spatial rotation.
    do ig=1,npw
      rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = GWPC_CONJG( rhxtwg_in(ig) * Gsph%phmSGt(ig,isym_k) )
    end do

   CASE DEFAULT
     write(msg,'(a,i0)')"Wrong value of itim_k: ",itim_k
     MSG_ERROR(msg)
   END SELECT

 END SELECT 

end function sym_rhotwgq0
!!***

END MODULE m_oscillators
!!***
