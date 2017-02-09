!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_getghc
!! NAME
!! prep_getghc
!!
!! FUNCTION
!! this routine prepares the data to the call of getghc.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2015 ABINIT group (FBottin,MT,GZ,MD,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blocksize= size of block for FFT
!!  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
!!  cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!         Typically lambda is the eigenvalue (or its guess)
!!  mpi_enreg=informations about mpi parallelization
!!  natom=number of atoms in cell.
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                       if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!  vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!   kinetic energy density, in real space, on the augmented fft grid. (optional argument)
!!   This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!
!! OUTPUT
!!  gwavef=(2,npw*nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                  or <G|H-lambda.S|C> (if sij_opt=-1).
!!  swavef=(2,npw*nspinor*ndat)=matrix elements <G|S|C>.
!!
!! SIDE EFFECTS
!!  ====== if gs_ham%usepaw==1
!!  cwaveprj(natom,nspinor*(1+cpopt)*ndat)= wave functions at k projected with nl projectors
!!
!! PARENTS
!!      chebfi,lobpcgwf
!!
!! CHILDREN
!!      dcopy,getghc,prep_index_wavef_bandpp,prep_sort_wavef_spin
!!      prep_wavef_sym_do,prep_wavef_sym_undo,timab,xmpi_alltoallv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_getghc(cwavef,dimffnl,dtfil,gs_hamk,gvnlc,gwavef,swavef,ikpt,&
&          lambda,blocksize,mpi_enreg,natom,npw_k,nspinor,paral_kgb,&
&          prtvol,sij_opt,vlocal,fock,cpopt,cwaveprj,&
&          already_transposed,vxctaulocal) ! optional argument

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_bandfft_kpt

 use m_pawcprj,     only : pawcprj_type
 use m_fock,        only : fock_type
 use m_hamiltonian, only : gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_getghc'
 use interfaces_18_timing
 use interfaces_66_wfs, except_this_one => prep_getghc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,dimffnl,ikpt,natom,npw_k,nspinor
 integer,intent(in) :: paral_kgb,prtvol,sij_opt
 real(dp),intent(in) :: lambda
 type(datafiles_type),intent(in) :: dtfil
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(mpi_type),intent(inout) :: mpi_enreg
 type(fock_type),pointer,intent(inout):: fock
!arrays
 real(dp),intent(in) :: cwavef(2,npw_k*nspinor*blocksize)
 real(dp),intent(inout) :: gvnlc (2,npw_k*nspinor*blocksize)
 real(dp),intent(inout) :: gwavef(2,npw_k*nspinor*blocksize)
 real(dp),intent(inout) :: swavef(2,npw_k*nspinor*blocksize)
 real(dp),intent(inout) :: vlocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc)
 integer, intent(in) :: cpopt
 type(pawcprj_type), intent(inout) :: cwaveprj(gs_hamk%natom,nspinor*blocksize*(1+cpopt))
 logical, intent(in),optional :: already_transposed
 real(dp), intent(inout), optional :: vxctaulocal(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc,4)

!Local variables-------------------------------
!local variables for mpialltoallv
!local variable for bandpp and inversion by symetry of time
!scalars
 integer,parameter :: tim_getghc=6
 integer :: bandpp,bandpp_sym,ier,ikpt_this_proc
 integer :: iscalc,nbval,nproc_band,nproc_fft
 integer :: old_me_g0,spaceComm=0
 logical :: flag_inv_sym, do_transpose
!arrays
 integer,pointer :: idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,ABI_CONTIGUOUS pointer :: kg_k_gather(:,:),kg_k_gather_sym(:,:)
 integer,allocatable :: index_wavef_band(:),index_wavef_send(:),index_wavef_spband(:)
 integer,ABI_CONTIGUOUS pointer :: rdispls(:),rdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: recvcounts(:),recvcounts_sym(:),recvcounts_sym_tot(:)
 integer,ABI_CONTIGUOUS pointer :: sdispls(:),sdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: sendcounts(:),sendcounts_sym(:),sendcounts_sym_all(:)
 integer,ABI_CONTIGUOUS pointer :: tab_proc(:)
 integer,allocatable :: rdisplsloc(:)
 integer,allocatable :: recvcountsloc(:),sdisplsloc(:)
 integer,allocatable :: sendcountsloc(:)
 real(dp) :: tsec(2)
 real(dp),ABI_CONTIGUOUS pointer :: ffnl_gather(:,:,:,:),kinpw_gather(:),ph3d_gather(:,:,:)
 real(dp),allocatable,target :: cwavef_alltoall1(:,:),gvnlc_alltoall1(:,:)
 real(dp),allocatable,target :: gwavef_alltoall1(:,:),swavef_alltoall1(:,:)
 real(dp),allocatable,target :: cwavef_alltoall2(:,:),gvnlc_alltoall2(:,:)
 real(dp),allocatable,target :: gwavef_alltoall2(:,:),swavef_alltoall2(:,:)
 real(dp),allocatable :: tmp_ffnl_gather(:,:,:,:),tmp_kinpw_gather(:),tmp_ph3d_gather(:,:,:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:),gvnlc_alltoall_sym(:,:)
 real(dp),pointer :: gwavef_alltoall_sym(:,:)
 real(dp),pointer :: swavef_alltoall_sym(:,:)

! *************************************************************************
 

 call timab(630,1,tsec)
 call timab(631,3,tsec)
 
 do_transpose = .true.
 if(present(already_transposed)) then
   if(already_transposed) do_transpose = .false.
 end if

 nproc_band = mpi_enreg%nproc_band
 nproc_fft  = mpi_enreg%nproc_fft
 bandpp     = mpi_enreg%bandpp

 flag_inv_sym = (gs_hamk%istwf_k==2 .and. any(gs_hamk%ngfft(7) == [401,402,312]))

 if (flag_inv_sym) then
   gs_hamk%istwf_k = 1
   if (modulo(bandpp,2)==0) then
     bandpp_sym   = bandpp/2
   else
     bandpp_sym   = bandpp
   end if
 end if
!====================================================================================

 spaceComm=mpi_enreg%comm_fft
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_band
 ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)

 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 recvcounts   =>bandfft_kpt(ikpt_this_proc)%recvcounts(:)
 sendcounts   =>bandfft_kpt(ikpt_this_proc)%sendcounts(:)
 rdispls      =>bandfft_kpt(ikpt_this_proc)%rdispls   (:)
 sdispls      =>bandfft_kpt(ikpt_this_proc)%sdispls   (:)
 ndatarecv    =>bandfft_kpt(ikpt_this_proc)%ndatarecv

 kg_k_gather           =>bandfft_kpt(ikpt_this_proc)%kg_k_gather(:,:)
 kinpw_gather          =>bandfft_kpt(ikpt_this_proc)%kinpw_gather(:)
 ffnl_gather           =>bandfft_kpt(ikpt_this_proc)%ffnl_gather(:,:,:,:)
 ph3d_gather           =>bandfft_kpt(ikpt_this_proc)%ph3d_gather(:,:,:)
 gs_hamk%gbound(:,:)   = bandfft_kpt(ikpt_this_proc)%gbound(:,:)

 if (flag_inv_sym ) then
   idatarecv0           =>bandfft_kpt(ikpt_this_proc)%idatarecv0
   ndatarecv_tot        =>bandfft_kpt(ikpt_this_proc)%ndatarecv_tot
   ndatasend_sym        =>bandfft_kpt(ikpt_this_proc)%ndatasend_sym
   kg_k_gather_sym      =>bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym(:,:)
   rdispls_sym          =>bandfft_kpt(ikpt_this_proc)%rdispls_sym(:)
   recvcounts_sym       =>bandfft_kpt(ikpt_this_proc)%recvcounts_sym(:)
   recvcounts_sym_tot   =>bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot(:)
   sdispls_sym          =>bandfft_kpt(ikpt_this_proc)%sdispls_sym(:)
   sendcounts_sym       =>bandfft_kpt(ikpt_this_proc)%sendcounts_sym(:)
   sendcounts_sym_all   =>bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all(:)
   tab_proc             =>bandfft_kpt(ikpt_this_proc)%tab_proc(:)
 end if
 iscalc=(sij_opt+1)/2  ! 0 if S not calculated, 1 otherwise
 nbval=(ndatarecv*nspinor*bandpp)*iscalc

 if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. nspinor==2 ).or. &
&      ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
   ABI_ALLOCATE(cwavef_alltoall1,(2,ndatarecv*nspinor*bandpp))
   ABI_ALLOCATE(gwavef_alltoall1,(2,ndatarecv*nspinor*bandpp))
   ABI_ALLOCATE(swavef_alltoall1,(2,ndatarecv*nspinor*bandpp))
   ABI_ALLOCATE(gvnlc_alltoall1,(2,ndatarecv*nspinor*bandpp))
   swavef_alltoall1(:,:)=zero
   gvnlc_alltoall1(:,:)=zero
   cwavef_alltoall1(:,:)=zero
   gwavef_alltoall1(:,:)=zero
 endif
 ABI_ALLOCATE(cwavef_alltoall2,(2,ndatarecv*nspinor*bandpp))
 ABI_ALLOCATE(gwavef_alltoall2,(2,ndatarecv*nspinor*bandpp))
 ABI_ALLOCATE(swavef_alltoall2,(2,ndatarecv*nspinor*bandpp))
 ABI_ALLOCATE(gvnlc_alltoall2,(2,ndatarecv*nspinor*bandpp))
 swavef_alltoall2(:,:)=zero
 gvnlc_alltoall2(:,:)=zero
 cwavef_alltoall2(:,:)=zero
 gwavef_alltoall2(:,:)=zero

 recvcountsloc(:)=recvcounts(:)*2*nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*nspinor
 sdisplsloc(:)=sdispls(:)*2*nspinor
 call timab(631,2,tsec)

 if(do_transpose) then
   call timab(545,3,tsec)
   if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. nspinor==2 ).or. &
&       ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall1,&
&      recvcountsloc,rdisplsloc,spaceComm,ier)
   else
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall2,&
&      recvcountsloc,rdisplsloc,spaceComm,ier)
   endif
   call timab(545,2,tsec)
 else
   ! Here, we cheat, and use DCOPY to bypass some compiler's overzealous bound-checking
   ! (ndatarecv*nspinor*bandpp might be greater than the declared size of cwavef)
   call DCOPY(2*ndatarecv*nspinor*bandpp, cwavef, 1, cwavef_alltoall2, 1)
 end if

 if(gs_hamk%istwf_k==2) then
   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
 end if

!====================================================================
 if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then
   if (do_transpose .and. mpi_enreg%paral_spinor==0.and.nspinor==2)then
     call timab(632,3,tsec)
!    Sort to have all nspinor=1 first, then all nspinor=2
     call prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef_spband)
     cwavef_alltoall2(:,:)=cwavef_alltoall1(:,index_wavef_spband)
     call timab(632,2,tsec)
   end if

   call timab(635,3,tsec)
   if(present(vxctaulocal))then
     call getghc(cpopt,cwavef_alltoall2,cwaveprj,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall2,swavef_alltoall2(:,1:nbval),gs_hamk,gvnlc_alltoall2,kg_k_gather,&
&     kinpw_gather,lambda,mpi_enreg,natom,1,ndatarecv,nspinor,&
&     paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal,fock,ikpt_this_proc,vxctaulocal=vxctaulocal)
   else
     call getghc(cpopt,cwavef_alltoall2,cwaveprj,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall2,swavef_alltoall2(:,1:nbval),gs_hamk,gvnlc_alltoall2,kg_k_gather,&
&     kinpw_gather,lambda,mpi_enreg,natom,1,ndatarecv,nspinor,&
&     paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal,fock,ikpt_this_proc)
   end if
   call timab(635,2,tsec)

   if (do_transpose .and. mpi_enreg%paral_spinor==0.and.nspinor==2)then
     call timab(634,3,tsec)
     gwavef_alltoall1(:,index_wavef_spband)=gwavef_alltoall2(:,:)
     if (sij_opt==1) swavef_alltoall1(:,index_wavef_spband)=swavef_alltoall2(:,:)
     gvnlc_alltoall1(:,index_wavef_spband)=gvnlc_alltoall2(:,:)
     ABI_DEALLOCATE(index_wavef_spband)
     call timab(634,2,tsec)
   end if

 else if ((.not.(flag_inv_sym)) .and. (bandpp>1)) then
!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------

   if(do_transpose) then
     call timab(632,3,tsec)
     call prep_index_wavef_bandpp(nproc_band,bandpp,&
&     nspinor,ndatarecv, recvcounts,rdispls, index_wavef_band)
!  -------------------------------------------------------
!  Sorting of the waves functions below bandpp
!  -------------------------------------------------------
     cwavef_alltoall2(:,:) = cwavef_alltoall1(:,index_wavef_band)
     call timab(632,2,tsec)
   end if

!  ----------------------
!  Fourier transformation
!  ----------------------
   call timab(636,3,tsec)
   if(present(vxctaulocal))then
     call getghc(cpopt,cwavef_alltoall2,cwaveprj,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall2,swavef_alltoall2,gs_hamk,gvnlc_alltoall2,kg_k_gather,&
&     kinpw_gather,lambda,mpi_enreg,natom,bandpp,ndatarecv,nspinor,&
&     paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal,fock,ikpt_this_proc,vxctaulocal=vxctaulocal)

   else
     call getghc(cpopt,cwavef_alltoall2,cwaveprj,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall2,swavef_alltoall2,gs_hamk,gvnlc_alltoall2,kg_k_gather,&
&     kinpw_gather,lambda,mpi_enreg,natom,bandpp,ndatarecv,nspinor,&
&     paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal,fock,ikpt_this_proc)
   end if
   call timab(636,2,tsec)

!  -----------------------------------------------------
!  Sorting of waves functions below the processors
!  -----------------------------------------------------
   if(do_transpose) then
     call timab(634,3,tsec)
!    cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)   ! NOT NEEDED
     gwavef_alltoall1(:,index_wavef_band) = gwavef_alltoall2(:,:)
     if (sij_opt==1) swavef_alltoall1(:,index_wavef_band) = swavef_alltoall2(:,:)
     gvnlc_alltoall1(:,index_wavef_band)  = gvnlc_alltoall2(:,:)

     ABI_DEALLOCATE(index_wavef_band)
     call timab(634,2,tsec)
   end if


 else if (flag_inv_sym) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   if(do_transpose) then
     call timab(632,3,tsec)
     call prep_index_wavef_bandpp(nproc_band,bandpp,&
&     nspinor,ndatarecv,&
&     recvcounts,rdispls,&
&     index_wavef_band)

!  -------------------------------------------------------
!  Sorting the wave functions below bandpp
!  -------------------------------------------------------
     cwavef_alltoall2(:,:) = cwavef_alltoall1(:,index_wavef_band)
   end if

!  ------------------------------------------------------------
!  We associate the waves functions by two
!  ------------------------------------------------------------
   call prep_wavef_sym_do(mpi_enreg,bandpp,nspinor,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,tab_proc,&
&   cwavef_alltoall2,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   ewavef_alltoall_sym,&
&   index_wavef_send)

!  ------------------------------------------------------------
!  Allocation
!  ------------------------------------------------------------
   ABI_ALLOCATE(gwavef_alltoall_sym,(2,ndatarecv_tot*bandpp_sym))
   ABI_ALLOCATE(swavef_alltoall_sym,(2,(ndatarecv_tot*bandpp_sym)*iscalc))
   ABI_ALLOCATE(gvnlc_alltoall_sym ,(2,ndatarecv_tot*bandpp_sym))

   gwavef_alltoall_sym(:,:)=zero
   swavef_alltoall_sym(:,:)=zero
   gvnlc_alltoall_sym(:,:)=zero

   ABI_ALLOCATE(tmp_ffnl_gather ,(ndatarecv_tot,dimffnl,gs_hamk%lmnmax,gs_hamk%ntypat))
   ABI_ALLOCATE(tmp_kinpw_gather,(ndatarecv_tot))
   ABI_ALLOCATE(tmp_ph3d_gather ,(2,ndatarecv_tot,gs_hamk%matblk))
   call timab(632,2,tsec)

!  ------------------------------------------------------------
!  Fourier calculation
!  ------------------------------------------------------------
   call timab(637,3,tsec)
   if(present(vxctaulocal))then
     call getghc(cpopt,ewavef_alltoall_sym,&
&     cwaveprj,dimffnl,tmp_ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall_sym, swavef_alltoall_sym,gs_hamk,gvnlc_alltoall_sym,&
&     kg_k_gather_sym,tmp_kinpw_gather,lambda,mpi_enreg,natom,&
&     bandpp_sym,ndatarecv_tot,nspinor,paral_kgb,tmp_ph3d_gather,prtvol,&
&     sij_opt,tim_getghc,1,vlocal,fock,ikpt_this_proc,vxctaulocal=vxctaulocal)
   else
     call getghc(cpopt,ewavef_alltoall_sym,&
&     cwaveprj,dimffnl,tmp_ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall_sym, swavef_alltoall_sym,gs_hamk,gvnlc_alltoall_sym,&
&     kg_k_gather_sym,tmp_kinpw_gather,lambda,mpi_enreg,natom,&
&     bandpp_sym,ndatarecv_tot,nspinor,paral_kgb,tmp_ph3d_gather,prtvol,&
&     sij_opt,tim_getghc,1,vlocal,fock,ikpt_this_proc)
   end if
   call timab(637,2,tsec)

   call timab(633,3,tsec)
   ABI_DEALLOCATE(tmp_ffnl_gather)
   ABI_DEALLOCATE(tmp_kinpw_gather)
   ABI_DEALLOCATE(tmp_ph3d_gather)

!  ------------------------------------------------------------
!  We dissociate each wave function in two waves functions
!  gwavef is classed below of bandpp
!  ------------------------------------------------------------
   call prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,idatarecv0,&
&   gwavef_alltoall2,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   gwavef_alltoall_sym,&
&   index_wavef_send)
  if (sij_opt==1)then
    call prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,idatarecv0,&
&     swavef_alltoall2,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     swavef_alltoall_sym,&
&     index_wavef_send)
  endif
  call prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,idatarecv0,&
&   gvnlc_alltoall2,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   gvnlc_alltoall_sym,&
&   index_wavef_send)

   ABI_DEALLOCATE(ewavef_alltoall_sym)
   ABI_DEALLOCATE(index_wavef_send)
   ABI_DEALLOCATE(gwavef_alltoall_sym)
   ABI_DEALLOCATE(swavef_alltoall_sym)
   ABI_DEALLOCATE(gvnlc_alltoall_sym)

!  -------------------------------------------
!  We call getghc to calculate the nl matrix elements.
!  --------------------------------------------
   gs_hamk%istwf_k=2
   !!write(std_out,*)"Setting iswfk_k to 2"

   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
   call timab(633,2,tsec)

   call timab(638,3,tsec)
   if(present(vxctaulocal))then
     call getghc(cpopt,cwavef_alltoall2,cwaveprj,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall2,swavef_alltoall2,gs_hamk,gvnlc_alltoall2,kg_k_gather,&
&     kinpw_gather,lambda,mpi_enreg,natom,bandpp,ndatarecv,nspinor,&
&     paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,2,vlocal,fock,ikpt_this_proc,vxctaulocal=vxctaulocal)
   else
     call getghc(cpopt,cwavef_alltoall2,cwaveprj,dimffnl,ffnl_gather,dtfil%filstat,&
&     gwavef_alltoall2,swavef_alltoall2,gs_hamk,gvnlc_alltoall2,kg_k_gather,&
&     kinpw_gather,lambda,mpi_enreg,natom,bandpp,ndatarecv,nspinor,&
&     paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,2,vlocal,fock,ikpt_this_proc)
   end if

   call timab(638,2,tsec)
   call timab(634,3,tsec)
   mpi_enreg%me_g0=old_me_g0

   gs_hamk%istwf_k=1

!  -------------------------------------------------------
!  Sorting the wave functions below the processors
!  -------------------------------------------------------
   if(do_transpose) then
!    cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)   ! NOT NEEDED
     gwavef_alltoall1(:,index_wavef_band) = gwavef_alltoall2(:,:)
     if (sij_opt==1) swavef_alltoall1(:,index_wavef_band) = swavef_alltoall2(:,:)
     gvnlc_alltoall1(:,index_wavef_band)  = gvnlc_alltoall2(:,:)
     ABI_DEALLOCATE(index_wavef_band)
     call timab(634,2,tsec)
   end if

 end if
!====================================================================
 
 if (gs_hamk%istwf_k==2) mpi_enreg%me_g0=old_me_g0

 if(do_transpose) then
   
   call timab(545,3,tsec)
   if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. nspinor==2 ).or. &
&       ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
     if (sij_opt==1) then
       call xmpi_alltoallv(swavef_alltoall1,recvcountsloc,rdisplsloc,swavef,&
&       sendcountsloc,sdisplsloc,spaceComm,ier)
     end if
     call xmpi_alltoallv(gvnlc_alltoall1,recvcountsloc,rdisplsloc,gvnlc,&
&      sendcountsloc,sdisplsloc,spaceComm,ier)
     call xmpi_alltoallv(gwavef_alltoall1,recvcountsloc,rdisplsloc,gwavef,&
&      sendcountsloc,sdisplsloc,spaceComm,ier)
   else
     if (sij_opt==1) then
       call xmpi_alltoallv(swavef_alltoall2,recvcountsloc,rdisplsloc,swavef,&
&       sendcountsloc,sdisplsloc,spaceComm,ier)
     end if
     call xmpi_alltoallv(gvnlc_alltoall2,recvcountsloc,rdisplsloc,gvnlc,&
&      sendcountsloc,sdisplsloc,spaceComm,ier)
     call xmpi_alltoallv(gwavef_alltoall2,recvcountsloc,rdisplsloc,gwavef,&
&      sendcountsloc,sdisplsloc,spaceComm,ier)
   endif

   call timab(545,2,tsec)
 else
   if(sij_opt == 1) then
     call DCOPY(2*ndatarecv*nspinor*bandpp, swavef_alltoall2, 1, swavef, 1)
   end if
   call DCOPY(2*ndatarecv*nspinor*bandpp, gvnlc_alltoall2, 1, gvnlc, 1)
   call DCOPY(2*ndatarecv*nspinor*bandpp, gwavef_alltoall2, 1, gwavef, 1)
 end if
 
!====================================================================
 if (flag_inv_sym) then
   gs_hamk%istwf_k = 2
 end if
!====================================================================
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall2)
 ABI_DEALLOCATE(gwavef_alltoall2)
 ABI_DEALLOCATE(gvnlc_alltoall2)
 ABI_DEALLOCATE(swavef_alltoall2)

 if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. nspinor==2 ).or. &
&      ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
   ABI_DEALLOCATE(cwavef_alltoall1)
   ABI_DEALLOCATE(gwavef_alltoall1)
   ABI_DEALLOCATE(gvnlc_alltoall1)
   ABI_DEALLOCATE(swavef_alltoall1)
 endif

 call timab(630,2,tsec)

end subroutine prep_getghc
!!***
