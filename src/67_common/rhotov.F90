!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhotov
!! NAME
!! rhotov
!!
!! FUNCTION
!! This routine is called to compute, from a given total density
!! the trial (local) potential and the residual potential.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (XG, GMR, MT, EB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | spinmagntarget=input variable that governs fixed moment calculation
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!   | ntypat=number of types of atoms in unit cell.
!!   | occopt=option for occupancies
!!   | typat(natom)=type (integer) for each atom
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=option for the computation of additional energies
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  optxc=option to be used for the call to rhohxc
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  ucvol = unit cell volume (Bohr**3)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp(nfft)=array for holding local psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  ==== if optres==0
!!    vtrial(nfft,nspden)= old value of trial potential
!!
!! OUTPUT
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_hartree=Hartree part of total energy (hartree units)
!!   | e_xc=exchange-correlation energy (hartree)
!!  ==== if optene==0.or.2
!!   | e_localpsp=local psp energy (hartree)
!!  ==== if optene==1.or.2
!!   | e_xcdc=exchange-correlation double-counting energy (hartree)
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if optxc==2.
!!  strsxc(6)=xc contribution to stress tensor (hartree/bohr^3)
!!  vxc(nfft,nspden)=Vxc(r) (already computed above; gets recomputed below too)
!!  vxcavg=mean of the vxc potential
!!  ==== if optres==0
!!    vresidnew(nfft,nspden)=potential residual
!!    vnew_mean(nspden)=mean of the potential formed from vpsp, vhartr and vxc, might be spin-dependent
!!    vres_mean(nspden)=mean of the potential residual, might be spin-dependent
!!    vres2=square of the norm of the residual
!!  vxctau(nfftf,dtset%nspden*dtset%usekden,4)=derivative of XC energy density with respect to
!!    kinetic energy density (metaGGA cases) (optional output)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  vhartr(nfft)=array for holding Hartree potential
!!  ==== if optres==1
!!    vtrial(nfft,nspden)= new value of trial potential
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      prctfvw1,prctfvw2,prctfw3,scfcv
!!
!! CHILDREN
!!      dotprod_vn,mag_constr,mean_fftr,psolver_rhohxc,rhohxc,rhohxcpositron
!!      sqnorm_v,timab,wvl_psitohpsi,wvl_vtrial_abi2big,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rhotov(dtset,energies,gprimd,gsqcut,istep,kxc,mpi_enreg,nfft,ngfft,&
&  nhat,nhatgr,nhatgrdim,nkxc,vresidnew,n3xccc,optene,optres,optxc,&
&  rhog,rhor,rprimd,strsxc,ucvol,usepaw,usexcnhat,&
&  vhartr,vnew_mean,vpsp,vres_mean,vres2,vtrial,vxcavg,vxc,wvl,xccc3d,xred,&
&  electronpositron,taug,taur,vxctau) ! optional argument


 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_errors
 use m_profiling_abi
 use m_ab7_mixing
 use m_abi2big
 use m_xmpi
 use m_cgtools

 use m_energies,         only : energies_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhotov'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_62_poisson
 use interfaces_62_wvl_wfs
 use interfaces_67_common, except_this_one => rhotov
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n3xccc,nfft,nhatgrdim,nkxc,optene,optres,optxc,usepaw,istep
 integer,intent(in) :: usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: vres2,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(wvl_data), intent(inout) :: wvl
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft),vpsp(nfft)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vnew_mean(dtset%nspden)
 real(dp),intent(out) :: vres_mean(dtset%nspden),vresidnew(nfft,dtset%nspden)
 real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
 real(dp),intent(in),optional :: taur(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden*dtset%usekden,4)

!Local variables-------------------------------
!scalars
 integer :: nk3xc,ifft,ipositron,ispden,nfftot,offset
 integer :: mpi_comm_sphgrid
 real(dp) :: doti,e_xcdc_vxctau
 logical :: calc_xcdc,with_vxctau
 logical :: wvlbigdft=.false.
!arrays
 real(dp) :: evxc,tsec(2),vmean(dtset%nspden)
 real(dp) :: vzeeman(dtset%nspden)
 real(dp),allocatable :: Vmagconstr(:,:)
 real(dp),allocatable :: rhowk(:,:),vnew(:,:), xcart(:,:)
! *********************************************************************

 DBG_ENTER("COLL")

 call timab(940,1,tsec)

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = .false.
 if (present(vxctau) .and. present(taur) .and. dtset%usekden /= 0) with_vxctau = .true.

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 if(dtset%usewvl==1 .and. dtset%wvl_bigdft_comp==1) wvlbigdft=.true.

!mpi communicator for spherical grid
 mpi_comm_sphgrid=mpi_enreg%comm_fft
 if(dtset%usewvl==1) mpi_comm_sphgrid=mpi_enreg%comm_wvl

!Get size of FFT grid
 nfftot=PRODUCT(ngfft(1:3))

 ipositron=0
 if (present(electronpositron)) ipositron=electronpositron_calctype(electronpositron)

!------Compute Hartree and xc potentials----------------------------------

!allocate vnew here. 
!In wvl: vnew is used at call to wvl_psitohpsi
 if (optres==0) then
   ABI_ALLOCATE(vnew,(nfft,dtset%nspden))
   vmean(:)=zero ; vnew_mean(:)=zero
 end if

 if (ipositron/=1) then
!  Compute xc potential (separate up and down if spin-polarized)
   if (dtset%icoulomb == 0 .and. dtset%usewvl == 0) then
!    Use the periodic solver to compute Hxc.
     nk3xc=1
     call timab(941,1,tsec)
     if (ipositron==0) then
       if(with_vxctau)then
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur,vxctau=vxctau)
       else
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur)
       end if
     else
       if(with_vxctau)then
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,&
&         electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau)
       else
         call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&         rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,&
&         electronpositron=electronpositron,taug=taug,taur=taur)
       end if
     end if
     call timab(941,2,tsec)
   elseif (.not. wvlbigdft) then
!    Use the free boundary solver.
     call timab(943,1,tsec)
     call psolver_rhohxc(energies%e_hartree, energies%e_xc, evxc, &
&     dtset%icoulomb, dtset%ixc, &
&     mpi_enreg, nfft, &
&     ngfft, nhat,usepaw,&
&     dtset%nscforder, dtset%nspden, n3xccc, rhor,rprimd,&
&     usexcnhat,dtset%usepaw,dtset%usewvl,vhartr, vxc, vxcavg,&
&     wvl%descr,wvl%den,wvl%e,&
&     xccc3d,dtset%xclevel,dtset%xc_denpos)
!    
     call timab(943,2,tsec)
!    end if
   end if
!  For icoulomb==0 and usewvl Ehartree is calculated in psolver_rhohxc().
!  For PAW we recalculate this since nhat was not taken into account 
!  in psolver_rhohxc: E_H= int v_H (n+nhat) dr
   if(.not. wvlbigdft .and. (dtset%icoulomb==0 .or. dtset%usepaw==1 ) ) then
     call timab(942,1,tsec)
     call dotprod_vn(1,rhor,energies%e_hartree,doti,nfft,nfftot,1,1,vhartr,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
     energies%e_hartree=half*energies%e_hartree
     call timab(942,2,tsec)
   end if
 else
   call timab(944,1,tsec)
   energies%e_hartree=zero;energies%e_xc=zero
   call rhohxcpositron(electronpositron,gprimd,kxc,mpi_enreg,nfft,ngfft,nhat,nkxc,dtset%nspden,n3xccc,&
&   dtset%paral_kgb,rhor,strsxc,ucvol,usexcnhat,usepaw,vhartr,vxc,vxcavg,xccc3d,dtset%xc_denpos)
   call timab(944,2,tsec)
 end if

 call timab(945,1,tsec)
 if (ipositron/=0) then
   call dotprod_vn(1,rhor,electronpositron%e_hartree,doti,&
&   nfft,nfftot,1,1,electronpositron%vha_ep,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
   vhartr=vhartr+electronpositron%vha_ep
 end if


!------Compute parts of total energy depending on potentials--------
!usewvl: the next call is not needed, e_localpsp has already been calculated.
!write(*,*)'rhotov erase me: energies%e_localpsp',energies%e_localpsp
 if ( (optene==0.or.optene==2 ).and. .not. wvlbigdft) then
!  Compute local psp energy energies%e_localpsp
   call dotprod_vn(1,rhor,energies%e_localpsp,doti,nfft,nfftot,1,1,vpsp,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
 end if
!write(*,*)'energies%e_localpsp',energies%e_localpsp

 calc_xcdc=.false.
 if (optene==1.or.optene==2) calc_xcdc=.true.
 if (dtset%usewvl==1 .and. dtset%nnsclo>0) calc_xcdc=.true.
 if (wvlbigdft) calc_xcdc=.false.
 if (dtset%usefock==1) calc_xcdc=.true.

 if (calc_xcdc) then

!  Compute double-counting XC energy energies%e_xcdc
   if (ipositron/=1) then
     if (usepaw==0.or.usexcnhat/=0) then
       call dotprod_vn(1,rhor,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
       if(with_vxctau)then
         call dotprod_vn(1,taur,e_xcdc_vxctau,doti,nfft,nfftot,dtset%nspden,1,vxctau(:,:,1),ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
         energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
       end if
     else
       ABI_ALLOCATE(rhowk,(nfft,dtset%nspden))
       rhowk=rhor-nhat
       call dotprod_vn(1,rhowk,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
       ABI_DEALLOCATE(rhowk)
     end if
     if (ipositron==2) energies%e_xcdc=energies%e_xcdc-electronpositron%e_xcdc
   else
     energies%e_xcdc=zero
   end if

 end if


!------Produce residual vector and square norm of it-------------
!(only if requested ; if optres==0)

!Set up array for Zeeman field
!EB vzeeman(:) = factor*( Hz, Hx+iHy; Hx-iHy, -Hz)
!EB factor = -g/2 * mu_B * mu_0 = -1/2*B in a.u.
!EB <-- vzeeman might have to be allocated correctly --> to be checked
 vzeeman(:) = zero
 if(dtset%nspden==2)then
!  EB The collinear case has to be checked : 
!  EB Is it vzeeman(1) or (2) that has to be added here? to be checked in setvtr and energy as well
   vzeeman(1) =-half*dtset%zeemanfield(3)  ! For collinear ispden=2 is rho_up only
 end if
 if(dtset%nspden==4)then
!  EB  write (*,*) 'rhotov --> initialization of vzeeman '
!  EB  write (*,*) 'rhotov --> vzeeman before: ',vzeeman(:)
   vzeeman(1)=-half*dtset%zeemanfield(3)
   vzeeman(2)=half*dtset%zeemanfield(3)
   vzeeman(3)=-half*dtset%zeemanfield(1)
   vzeeman(4)=half*dtset%zeemanfield(2)
!  EB  write  (*,*) 'rhotov --> vzeeman after: ', vzeeman(:)
 end if

 ABI_ALLOCATE(Vmagconstr, (nfft,dtset%nspden))
 Vmagconstr = zero
 if (dtset%magconon==1.or.dtset%magconon==2) then
!  Compute the constrained potential for the magnetic moments
   call mag_constr(dtset%natom, dtset%spinat, dtset%nspden, dtset%magconon, dtset%magcon_lambda, rprimd, &
&   mpi_enreg, nfft, ngfft, dtset%ntypat, dtset%ratsph, rhor, &
&   dtset%typat, Vmagconstr, xred)
 else
!  NOTE: mjv 25 May 2013 added this for ibm6 - otherwise gives NaN in vnew after
!  the addition below, in case magconon==0 and for certain libxc
!  functionals!! May need to copy this to setvtr.F90 if the same effect appears.
!  should check libxc test with a proper memory checker (valgrind).
   do ispden=1,dtset%nspden
     do ifft=1,nfft
       Vmagconstr(ifft, ispden) = zero
     end do
   end do
 end if

 if (optres==0) then
!  Compute potential residual
   if (.not. wvlbigdft) then
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=1,min(dtset%nspden,2)
       do ifft=1,nfft
         vnew(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)+vzeeman(ispden)+Vmagconstr(ifft,ispden)
         vresidnew(ifft,ispden)=vnew(ifft,ispden)-vtrial(ifft,ispden)
       end do
     end do
     if(dtset%nspden==4)then
!$OMP PARALLEL DO COLLAPSE(2)
       do ispden=3,4
         do ifft=1,nfft
           vnew(ifft,ispden)=vxc(ifft,ispden)+vzeeman(ispden)+Vmagconstr(ifft,ispden)
!          EB ! might be a bug:
!          EB ! before the following line was: vresidnew(ifft,ispden)=vxc(ifft,ispden)-vtrial(ifft,ispden)
!          MT ! was not a bug as vtrial=vxc for ispden=4,3
           vresidnew(ifft,ispden)=vnew(ifft,ispden)-vtrial(ifft,ispden)
         end do
       end do
     end if
     offset   = 0

!    Pass vtrial to BigDFT object
     if(dtset%usewvl==1) then
       call wvl_vtrial_abi2big(nfft,dtset%nspden,1,vtrial,wvl%den)
     end if

   else
!    Compute with covering comms the different part of the potential.
!    only forwvlbigdft
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xred2xcart(dtset%natom, rprimd, xcart, xred)
     call wvl_psitohpsi(dtset%diemix,energies%e_exactX, energies%e_xc, energies%e_hartree, &
&     energies%e_kinetic, energies%e_localpsp, energies%e_nonlocalpsp, energies%e_sicdc, &
&     istep + 1, istep, dtset%iscf, dtset%ixc, mpi_enreg%me_wvl, dtset%natom, dtset%nfft,&
&     mpi_enreg%nproc_wvl, dtset%nspden, &
&     vres2, .true., energies%e_xcdc, wvl,&
&     wvlbigdft, xcart, strsxc,vnew, vxc)
     ABI_DEALLOCATE(xcart)

     vresidnew = vnew - vtrial
     vtrial = vnew

     call mean_fftr(vxc, vmean(1:1),  nfft, nfftot, dtset%nspden,&
&     mpi_comm_sphgrid=mpi_comm_sphgrid)
     vxcavg = vmean(1)
     offset = 0
   end if

!  Compute mean values of potential and residual
   call mean_fftr(vnew(1+offset, 1),vnew_mean,nfft,nfftot,dtset%nspden,&
&   mpi_comm_sphgrid=mpi_comm_sphgrid)
   call mean_fftr(vresidnew(1+offset, 1),vmean,nfft,nfftot,dtset%nspden,&
&   mpi_comm_sphgrid=mpi_comm_sphgrid)

   ABI_DEALLOCATE(vnew)

!  Subtract the mean of the residual
!  Must take into account fixed occupation number in case of spin-polarized
   do ispden=1,dtset%nspden
     if (dtset%nspden==2.and.dtset%occopt>=3.and. abs(dtset%spinmagntarget+99.99_dp)<1.0d-10)then
       vres_mean(ispden)=(vmean(1)+vmean(2))*half
     else
       vres_mean(ispden)=vmean(ispden)
     end if
   end do

!$OMP PARALLEL DO COLLAPSE(2)
   do ispden=1,dtset%nspden
     do ifft=1,nfft
       vresidnew(ifft,ispden)=vresidnew(ifft,ispden)-vres_mean(ispden)
     end do
   end do

!  Compute square norm vres2 of potential residual vresid
   call sqnorm_v(1,nfft,vres2,dtset%nspden,optres,vresidnew(1+offset, 1),mpi_comm_sphgrid=mpi_comm_sphgrid)

 else    ! optres==0
!  ------Produce new value of trial potential-------------
!  (only if requested ; if optres==1)

   if (.not. wvlbigdft) then
!    EB   write (*,*) 'rhotov  --> optres!=0 vzeeman is :',vzeeman(:)
!    EB   write (*,*) 'rhotov --> vtrial before zeeman:',vtrial(1,1),vtrial(1,2),vtrial(1,3),vtrial(1,4)
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=1,min(dtset%nspden,2)
       do ifft=1,nfft
         vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)+vzeeman(ispden)+Vmagconstr(ifft,ispden)
       end do
     end do
     if(dtset%nspden==4) then
!$OMP PARALLEL DO
       do ifft=1,nfft
         vtrial(ifft,3:4)=vxc(ifft,3:4)+vzeeman(3:4)+Vmagconstr(ifft,3:4)
       end do
     end if
!    EB   write (*,*) 'rhotov --> vtrial after zeeman:',vtrial(1,1),vtrial(1,2),vtrial(1,3),vtrial(1,4)
!    Pass vtrial to BigDFT object
     if(dtset%usewvl==1) then
       call wvl_vtrial_abi2big(nfft,dtset%nspden,1,vtrial,wvl%den)
     end if
   else
!    Compute with covering comms the different part of the potential.
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xred2xcart(dtset%natom, rprimd, xcart, xred)
     call wvl_psitohpsi(dtset%diemix,energies%e_exactX, energies%e_xc, energies%e_hartree, &
&     energies%e_kinetic, energies%e_localpsp, energies%e_nonlocalpsp, energies%e_sicdc, &
&     istep + 1, istep, dtset%iscf, dtset%ixc, mpi_enreg%me_wvl, &
&     dtset%natom, dtset%nfft, mpi_enreg%nproc_wvl,&
&     dtset%nspden,vres2, .true.,energies%e_xcdc,  wvl,&
&     wvlbigdft, xcart, strsxc, vtrial, vxc)
     ABI_DEALLOCATE(xcart)
!    Compute vxcavg
     call mean_fftr(vxc, vmean(1:1), nfft, nfftot, dtset%nspden,&
&     mpi_comm_sphgrid=mpi_comm_sphgrid)
     vxcavg = vmean(1)
   end if

 end if

 ABI_DEALLOCATE(Vmagconstr)

 call timab(945,2,tsec)
 call timab(940,2,tsec)

!debug:
!write(*,*)'rhotov: erase me l458'
!write(*,*)'energies%e_xcdc',energies%e_xcdc
!write(*,*)'energies%e_hartree',energies%e_hartree
!end debug

 DBG_EXIT("COLL")

end subroutine rhotov
!!***
