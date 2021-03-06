!!****m* ABINIT/interfaces_65_psp
!! NAME
!! interfaces_65_psp
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/65_psp
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_65_psp

 implicit none

interface
 subroutine atm2fft(atindx1,atmrho,atmvloc,dyfrn,dyfrv,eei,eltfrn,gauss,gmet,gprimd,&  
  &  grn,grv,gsqcut,mgfft,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&  
  &  optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,&  
  &  pawtab,ph1d,qgrid,qprtrb,rhog,strn,strv,ucvol,usepaw,vg,vg1,vg1_core,vprtrb,vspl,&  
  &  is2_in,mpi_comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments
  use defs_basis
  use m_distribfft
  use m_pawtab
  implicit none
  integer,optional,intent(in) :: is2_in
  integer,optional,intent(in) :: me_g0
  integer,intent(in) :: mgfft
  integer,optional,intent(in) :: mpi_comm_fft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,intent(in) :: optatm
  integer,intent(in) :: optdyfr
  integer,intent(in) :: opteltfr
  integer,intent(in) :: optgr
  integer,intent(in) :: optn
  integer,intent(in) :: optn2
  integer,intent(in) :: optstr
  integer,intent(in) :: optv
  integer,optional,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  type(distribfft_type),optional,intent(in),target :: distribfft
  real(dp),intent(in) :: eei
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: qprtrb(3)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(out) :: atmrho(nfft*optn)
  real(dp),intent(inout) :: atmvloc(nfft*optv)
  real(dp),intent(out) :: dyfrn(3,3,natom*optn*optdyfr)
  real(dp),intent(out) :: dyfrv(3,3,natom*optv*optdyfr)
  real(dp),intent(out) :: eltfrn(6+3*natom,6)
  real(dp),intent(in) :: gauss(2,ntypat*(optn2/3))
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(inout) :: grn(3,natom*optn*optgr)
  real(dp),intent(out) :: grv(3,natom*optv*optgr)
  integer,intent(in) :: nattyp(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rhog(2,nfft*optv*max(optgr,optstr,optdyfr,opteltfr))
  real(dp),intent(out) :: strn(6*optn*optstr)
  real(dp),intent(out) :: strv(6*optv*optstr)
  real(dp),intent(in) :: vg(2,nfft*optn*max(optgr,optstr,optdyfr,opteltfr))
  real(dp),intent(in) :: vg1(2,nfft*optn*opteltfr)
  real(dp),intent(in) :: vg1_core(2,nfft*optn*opteltfr)
  real(dp),intent(in) :: vprtrb(2)
  real(dp),intent(in) :: vspl(mqgrid,2,ntypat*optv)
 end subroutine atm2fft
end interface

interface
 subroutine atm2fft3(atindx,cplex,gmet,gprimd,gsqcut,idir,ipert,&  
  &  mgfft,mqgrid,natom,ndir,nfft,ngfft,ntypat,&  
  &  ph1d,qgrid,qphon,typat,ucvol,usepaw,xred,&  
  &  atmrho1,atmvloc1,distribfft,gauss,mpi_comm_fft,me_g0,optn_in,&  
  &  optn2_in,optv_in,pawtab,paral_kgb,vspl) ! optional arguments
  use defs_basis
  use m_distribfft
  use m_pawtab
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: idir
  integer,intent(in) :: ipert
  integer,optional,intent(in) :: me_g0
  integer,intent(in) :: mgfft
  integer,optional,intent(in) :: mpi_comm_fft
  integer,intent(in) :: mqgrid
  integer,intent(in) :: natom
  integer,intent(in) :: ndir
  integer,intent(in) :: nfft
  integer,intent(in) :: ntypat
  integer,optional,intent(in) :: optn2_in
  integer,optional,intent(in) :: optn_in
  integer,optional,intent(in) :: optv_in
  integer,optional,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  type(distribfft_type),optional,intent(in),target :: distribfft
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx(natom)
  real(dp),optional,intent(out) :: atmrho1(cplex*nfft,ndir)
  real(dp),optional,intent(out) :: atmvloc1(cplex*nfft,ndir)
  real(dp),optional,intent(in) :: gauss(2,ntypat)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  type(pawtab_type),optional,intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: qphon(3)
  integer,intent(in) :: typat(natom)
  real(dp),optional,intent(in) :: vspl(mqgrid,2,ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine atm2fft3
end interface

interface
 subroutine cc_derivatives(rad,ff,ff1,ff2,mmax,n1xccc,rchrg,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(in) :: ff(mmax)
  real(dp),intent(in) :: ff1(mmax)
  real(dp),intent(in) :: ff2(mmax)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine cc_derivatives
end interface

interface
 subroutine der_int(ff,df,rr,dr,nlast,smf)
  use defs_basis
  implicit none
  integer,parameter :: nmax=2000
  integer,intent(in) :: nlast
  real(dp),intent(out) :: smf
  real(dp), intent(out) :: df(0:nmax)
  real(dp), intent(in) :: dr(0:nmax)
  real(dp), intent(in) :: ff(0:nmax)
  real(dp), intent(in) :: rr(0:nmax)
 end subroutine der_int
end interface

interface
 subroutine gg1cc(gg1cc_xx,xx)
  use defs_basis
  implicit none
  real(dp),intent(out) :: gg1cc_xx
  real(dp),intent(in) :: xx
 end subroutine gg1cc
end interface

interface
 subroutine gp1cc(gp1cc_xx,xx)
  use defs_basis
  implicit none
  real(dp),intent(out) :: gp1cc_xx
  real(dp),intent(in) :: xx
 end subroutine gp1cc
end interface

interface
 subroutine gpp1cc(gpp1cc_xx,xx)
  use defs_basis
  implicit none
  real(dp),intent(out) :: gpp1cc_xx
  real(dp),intent(in) :: xx
 end subroutine gpp1cc
end interface

interface
 subroutine psden(ilog,ff,mesh,nc,rc,rad,ff1,ff2)
  use defs_basis
  implicit none
  integer,intent(in) :: ilog
  integer,intent(in) :: mesh
  real(dp),intent(in) :: rc
  real(dp),intent(out) :: ff(mesh)
  real(dp),intent(inout),optional :: ff1(mesh)
  real(dp),intent(inout),optional :: ff2(mesh)
  real(dp),intent(in) :: nc(mesh)
  real(dp),intent(in) :: rad(mesh)
 end subroutine psden
end interface

interface
 subroutine psp10in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax, nproj, psps, pspso, vlspl, zion)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipsp
  integer,intent(inout) :: lmax
  integer,intent(in) :: pspso
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: ekb(psps%lnmax)
  real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  integer,intent(out) :: nproj(psps%mpssoang)
  real(dp),intent(out) :: vlspl(psps%mqgrid_ff,2)
 end subroutine psp10in
end interface

interface
 subroutine psp10nl(ekb,ffspl,hij,lmax,mproj,mpsang,mqgrid,nproj,qgrid,rr)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: mproj
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(out) :: ekb(mpsang,mproj)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang,mproj)
  real(dp),intent(in) :: hij(0:lmax,3,3)
  integer,intent(in) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rr(0:lmax)
 end subroutine psp10nl
end interface

interface
 subroutine psp11lo(drdi,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&  
  &  vloc,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: drdi(mmax)
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
 end subroutine psp11lo
end interface

interface
 subroutine psp11nl(ffspl,indlmn,mmax,lnmax,lmnmax,mqgrid,n_proj,&  
  &  proj, proj_l, proj_np, qgrid, r, drdi, useylm)
  use defs_basis
  implicit none
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n_proj
  integer,intent(in) :: useylm
  real(dp),intent(in) :: drdi(mmax)
  real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)
  integer, intent(out) :: indlmn(6,lmnmax)
  real(dp),intent(in) :: proj(mmax,n_proj)
  integer, intent(in) :: proj_l(n_proj)
  integer, intent(in) :: proj_np(n_proj)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: r(mmax)
 end subroutine psp11nl
end interface

interface
 subroutine psp1cc(fchrg,n1xccc,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: fchrg
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp1cc
end interface

interface
 subroutine psp1in(dq,ekb,ekb1,ekb2,epsatm,epspsp,&  
  &  e990,e999,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mqgrid,nproj,n1xccc,pspcod,&  
  &  qchrg,qgrid,rcpsp,rms,useylm,vlspl,xcccrc,xccc1d,&  
  &  zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: pspcod
  integer,intent(in) :: useylm
  real(dp),intent(in) :: dq
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: e990(mpsang)
  real(dp),intent(out) :: e999(mpsang)
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ekb1(mpsang)
  real(dp),intent(out) :: ekb2(mpsang)
  real(dp),intent(out) :: epspsp(mpsang)
  real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: rcpsp(mpsang)
  real(dp),intent(out) :: rms(mpsang)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp1in
end interface

interface
 subroutine psp1lo(drad,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&  
  &  vloc,wksincos,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: drad(mmax)
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
  real(dp),intent(inout) :: wksincos(mmax,2,2)
 end subroutine psp1lo
end interface

interface
 subroutine psp1nl(dr,ekb,ffspl,lloc,lmax,mmax,mpsang,mqgrid,&  
  &  qgrid,rad,vloc,vpspll,wfll,wksincos)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: dr(mmax)
  real(dp),intent(out) :: ekb(mpsang)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
  real(dp),intent(in) :: vpspll(mmax,mpsang)
  real(dp),intent(in) :: wfll(mmax,mpsang)
  real(dp),intent(inout) :: wksincos(mmax,2,2)
 end subroutine psp1nl
end interface

interface
 subroutine psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipsp
  integer,intent(in) :: lmax
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
  real(dp),intent(out) :: ekb(psps%lnmax)
  real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  integer,intent(out) :: nproj(psps%mpsang)
  real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
 end subroutine psp2in
end interface

interface
 subroutine psp2lo(cc1,cc2,cc3,cc4,dvloc,epsatm,mqgrid,qgrid,q2vq,&  
  &  rloc,vlspl_recipSpace,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: cc1
  real(dp),intent(in) :: cc2
  real(dp),intent(in) :: cc3
  real(dp),intent(in) :: cc4
  real(dp),intent(out) :: epsatm
  real(dp),intent(in) :: rloc
  logical,intent(in) :: vlspl_recipSpace
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: dvloc(mqgrid)
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
 end subroutine psp2lo
end interface

interface
 subroutine psp2nl(ekb,ffspl,h1p,h1s,h2s,lnmax,mqgrid,qgrid,rrp,rrs)
  use defs_basis
  implicit none
  integer,intent(in) :: lnmax
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: h1p
  real(dp),intent(in) :: h1s
  real(dp),intent(in) :: h2s
  real(dp),intent(in) :: rrp
  real(dp),intent(in) :: rrs
  real(dp),intent(inout) :: ekb(lnmax)
  real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)
  real(dp),intent(in) :: qgrid(mqgrid)
 end subroutine psp2nl
end interface

interface
 subroutine psp3in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax, nproj, psps, pspso, vlspl, zion)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ipsp
  integer,intent(inout) :: lmax
  integer,intent(in) :: pspso
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(in) :: zion
  real(dp),intent(inout) :: ekb(psps%lnmax)
  real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  integer,intent(out) :: nproj(psps%mpssoang)
  real(dp),intent(out) :: vlspl(psps%mqgrid_ff,2)
 end subroutine psp3in
end interface

interface
 subroutine psp3nl(ekb,ffspl,h11s,h22s,h33s,h11p,h22p,h33p,h11d,h22d,&  
  &  h33d,h11f,mproj,mpsang,mqgrid,qgrid,rrd,rrf,rrp,rrs)
  use defs_basis
  implicit none
  integer,intent(in) :: mproj
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: h11d
  real(dp),intent(in) :: h11f
  real(dp),intent(in) :: h11p
  real(dp),intent(in) :: h11s
  real(dp),intent(in) :: h22d
  real(dp),intent(in) :: h22p
  real(dp),intent(in) :: h22s
  real(dp),intent(in) :: h33d
  real(dp),intent(in) :: h33p
  real(dp),intent(in) :: h33s
  real(dp),intent(in) :: rrd
  real(dp),intent(in) :: rrf
  real(dp),intent(in) :: rrp
  real(dp),intent(in) :: rrs
  real(dp),intent(out) :: ekb(mpsang,mproj)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang,mproj)
  real(dp),intent(in) :: qgrid(mqgrid)
 end subroutine psp3nl
end interface

interface
 subroutine psp4cc(fchrg,n1xccc,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: fchrg
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp4cc
end interface

interface
 subroutine psp5in(ekb,ekb1,ekb2,epsatm,epspsp,e990,e999,ffspl,indlmn,&  
  &  lloc,lmax,lmnmax,lnmax,mmax,mpsang,mpssoang,mqgrid,&  
  &  nproj,n1xccc,pspso,qchrg,qgrid,rcpsp,rms,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: pspso
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: e990(mpssoang)
  real(dp),intent(out) :: e999(mpssoang)
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ekb1(mpssoang)
  real(dp),intent(out) :: ekb2(mpssoang)
  real(dp),intent(out) :: epspsp(mpssoang)
  real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(inout) :: nproj(mpssoang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: rcpsp(mpssoang)
  real(dp),intent(out) :: rms(mpssoang)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp5in
end interface

interface
 subroutine psp5lo(al,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&  
  &  vloc,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: al
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
 end subroutine psp5lo
end interface

interface
 subroutine psp5nl(al,ekb,ffspl,lmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,vpspll,wfll)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: al
  real(dp),intent(out) :: ekb(mpsang)
  real(dp),intent(out) :: ffspl(mqgrid,2,mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
  real(dp),intent(in) :: vpspll(mmax,mpsang)
  real(dp),intent(in) :: wfll(mmax,mpsang)
 end subroutine psp5nl
end interface

interface
 subroutine psp6cc(mmax,n1xccc,rchrg,xccc1d,znucl,&  
  &  vh_tnzc) ! optional argument
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(in) :: znucl
  real(dp),intent(out),optional :: vh_tnzc(mmax)
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp6cc
end interface

interface
 subroutine psp6cc_drh(mmax,n1xccc,rchrg,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp6cc_drh
end interface

interface
 subroutine psp6in(ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mqgrid,nproj,n1xccc,optnlxccc,positron,qchrg,qgrid,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: optnlxccc
  integer,intent(in) :: positron
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpsang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp6in
end interface

interface
 subroutine psp8cc(mmax,n1xccc,rchrg,xccc1d)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: n1xccc
  real(dp),intent(in) :: rchrg
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp8cc
end interface

interface
 subroutine psp8in(ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mpssoang,mqgrid,nproj,n1xccc,pspso,qchrg,qgrid,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: pspso
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpssoang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp8in
end interface

interface
 subroutine psp8lo(amesh,epsatm,mmax,mqgrid,qgrid,q2vq,rad,vloc,yp1,ypn,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: amesh
  real(dp),intent(out) :: epsatm
  real(dp),intent(out) :: yp1
  real(dp),intent(out) :: ypn
  real(dp),intent(in) :: zion
  real(dp),intent(out) :: q2vq(mqgrid)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vloc(mmax)
 end subroutine psp8lo
end interface

interface
 subroutine psp8nl(amesh,ffspl,indlmn,lmax,lmnmax,lnmax,mmax,&  
  &  mqgrid,qgrid,rad,vpspll)
  use defs_basis
  implicit none
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mmax
  integer,intent(in) :: mqgrid
  real(dp),intent(in) :: amesh
  real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)
  integer,intent(in) :: indlmn(6,lmnmax)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(in) :: rad(mmax)
  real(dp),intent(in) :: vpspll(mmax,lnmax)
 end subroutine psp8nl
end interface

interface
 subroutine psp9in(filpsp,ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&  
  &  mmax,mpsang,mpssoang,mqgrid,nproj,n1xccc,pspso,qchrg,qgrid,&  
  &  useylm,vlspl,xcccrc,xccc1d,zion,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: lloc
  integer,intent(in) :: lmax
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(out) :: mmax
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: pspso
  integer,intent(in) :: useylm
  real(dp),intent(out) :: epsatm
  character(len=fnlen),intent(in) :: filpsp
  real(dp),intent(out) :: qchrg
  real(dp),intent(out) :: xcccrc
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
  real(dp),intent(out) :: ekb(lnmax)
  real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)
  integer,intent(out) :: indlmn(6,lmnmax)
  integer,intent(out) :: nproj(mpssoang)
  real(dp),intent(in) :: qgrid(mqgrid)
  real(dp),intent(out) :: vlspl(mqgrid,2)
  real(dp),intent(inout) :: xccc1d(n1xccc,6)
 end subroutine psp9in
end interface

interface
 subroutine pspatm(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&  
  &  psps,vlspl,dvlspl,xcccrc,xccc1d,comm_mpi)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawrad
  use m_pawtab
  implicit none
  integer, optional,intent(in) :: comm_mpi
  integer,intent(in) :: ipsp
  real(dp),intent(in) :: dq
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pawrad_type),intent(inout) :: pawrad
  type(pawtab_type),intent(inout) :: pawtab
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(out) :: xcccrc
  real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: ekb(psps%dimekb*(1-psps%usepaw))
  real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(inout) :: indlmn(6,psps%lmnmax)
  real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6)
 end subroutine pspatm
end interface

interface
 subroutine pspatm_abinit(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&  
  &  psps,vlspl,dvlspl,xcccrc,xccc1d,comm_mpi)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawrad
  use m_pawtab
  implicit none
  integer, optional,intent(in) :: comm_mpi
  integer,intent(in) :: ipsp
  real(dp),intent(in) :: dq
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pawrad_type),intent(inout) :: pawrad
  type(pawtab_type),intent(inout) :: pawtab
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(out) :: xcccrc
  real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: ekb(psps%dimekb*(1-psps%usepaw))
  real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6)
 end subroutine pspatm_abinit
end interface

interface
 subroutine pspatm_pspio(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&  
  &  psps,vlspl,dvlspl,xcccrc,xccc1d,comm_mpi)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawrad
  use m_pawtab
  implicit none
  integer, optional,intent(in) :: comm_mpi
  integer,intent(in) :: ipsp
  real(dp),intent(in) :: dq
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: epsatm
  type(pawrad_type),intent(out) :: pawrad
  type(pawtab_type),intent(out) :: pawtab
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(out) :: xcccrc
  real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: ekb(psps%dimekb*(1-psps%usepaw))
  real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer,intent(out) :: indlmn(6,psps%lmnmax)
  real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp),intent(inout) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6)
 end subroutine pspatm_pspio
end interface

interface
 subroutine pspcor(ecore,epsatm,natom,ntypat,typat,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(out) :: ecore
  real(dp),intent(in) :: epsatm(ntypat)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine pspcor
end interface

interface
 subroutine pspini(dtset,dtfil,ecore,gencond,gsqcut,gsqcutdg,level,pawrad,pawtab,psps,rprimd,comm_mpi)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_pawrad
  use m_pawtab
  implicit none
  integer, optional,intent(in) :: comm_mpi
  integer,intent(out) :: gencond
  integer,intent(in) :: level
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: ecore
  real(dp),intent(in) :: gsqcut
  real(dp),intent(in) :: gsqcutdg
  type(pseudopotential_type), target,intent(inout) :: psps
  type(pawrad_type), intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine pspini
end interface

interface
 subroutine radii_ps( vps, rofi, zval, nrval, lmxkb, nrgauss, rgauss, rgauss2)
  use defs_basis
  implicit none
  integer,intent(in) :: lmxkb
  integer,intent(out) :: nrgauss
  integer,intent(in) :: nrval
  real(dp),intent(out) :: rgauss
  real(dp),intent(out) :: rgauss2
  real(dp),intent(in) :: zval
  real(dp),intent(in) :: rofi(:)
  real(dp),intent(in) :: vps(:,0:)
 end subroutine radii_ps
end interface

interface
 subroutine sbf8(nm,xx,sb_out)
  use defs_basis
  implicit none
  integer,intent(in) :: nm
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: sb_out(nm)
 end subroutine sbf8
end interface

interface
 subroutine sincos(iq,irmax,mmax,pspwk,rad,tpiq)
  use defs_basis
  implicit none
  integer,intent(in) :: iq
  integer,intent(in) :: irmax
  integer,intent(in) :: mmax
  real(dp),intent(in) :: tpiq
  real(dp),intent(inout) :: pspwk(mmax,2,2)
  real(dp),intent(in) :: rad(mmax)
 end subroutine sincos
end interface

interface
 subroutine smoothvlocal(  lmax,npts,scale,step,vlocal,vps,zval,rgauss,rgauss2 )
  use defs_basis
  implicit none
  integer, intent(in) :: lmax
  integer, intent(in) :: npts
  real(dp), intent(out) :: rgauss
  real(dp), intent(out) :: rgauss2
  real(dp), intent(in) :: scale
  real(dp), intent(in) :: step
  real(dp), intent(in) :: zval
  real(dp), intent(out) :: vlocal(npts)
  real(dp), intent(in) :: vps(npts,0:lmax)
 end subroutine smoothvlocal
end interface

interface
 subroutine vlocal2( zval, nrval, a, rofi, drdi, s, vps, nrgauss,&  
  &  vlocal,nchloc,chlocal )
  use defs_basis
  implicit none
  integer,  intent(out) :: nchloc
  integer,  intent(inout) :: nrgauss
  integer,  intent(in) :: nrval
  real(dp), intent(in) :: a
  real(dp), intent(in) :: zval
  real(dp), intent(out) :: chlocal(:)
  real(dp), intent(in) :: drdi(:)
  real(dp), intent(in) :: rofi(:)
  real(dp), intent(in) :: s(:)
  real(dp), intent(out) :: vlocal(:)
  real(dp), intent(in) :: vps(:)
 end subroutine vlocal2
end interface

interface
 subroutine vlocal1( zval, nrval, a, rofi, drdi, s, rgauss, vlocal,&  
  &  nchloc, chlocal)
  use defs_basis
  implicit none
  integer,  intent(out) :: nchloc
  integer,  intent(in) :: nrval
  real(dp), intent(in) :: a
  real(dp), intent(inout) :: rgauss
  real(dp), intent(in) :: zval
  real(dp), intent(out) :: chlocal(:)
  real(dp), intent(in) :: drdi(:)
  real(dp), intent(in) :: rofi(:)
  real(dp), intent(in) :: s(:)
  real(dp), intent(out) :: vlocal(:)
 end subroutine vlocal1
end interface

interface
 subroutine upf2abinit (filpsp, znucl, zion, pspxc, lmax_, lloc, mmax,&  
  &  psps, epsatm, xcccrc, indlmn, ekb, ffspl, nproj_l, vlspl, xccc1d)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(out) :: lloc
  integer, intent(out) :: lmax_
  integer, intent(out) :: mmax
  integer, intent(out) :: pspxc
  real(dp), intent(out) :: epsatm
  character(len=fnlen), intent(in) :: filpsp
  type(pseudopotential_type),intent(in) :: psps
  real(dp), intent(out) :: xcccrc
  real(dp), intent(out) :: zion
  real(dp), intent(out) :: znucl
  real(dp), intent(inout) :: ekb(psps%dimekb)
  real(dp), intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
  integer, intent(out) :: indlmn(6,psps%lmnmax)
  integer, intent(out) :: nproj_l(psps%mpssoang)
  real(dp), intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp), intent(inout) :: xccc1d(psps%n1xccc,6)
 end subroutine upf2abinit
end interface

interface
 function vander(a,x) result(f)
  use defs_basis
  implicit none
  real(dp),intent(in) :: a
  real(dp) :: f
  real(dp),intent(in) :: x
 end function vander
end interface

interface
 subroutine vhtnzc(nc,rc,vh_tnzc,mesh,rad,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: mesh
  real(dp),intent(in) :: rc
  real(dp),intent(in) :: znucl
  real(dp),intent(in) :: nc(mesh)
  real(dp),intent(in) :: rad(mesh)
  real(dp),intent(out) :: vh_tnzc(mesh)
 end subroutine vhtnzc
end interface

interface
 subroutine wvlpaw_sin2gauss(basis_size,mparam,nparam,param,wvl)
  use defs_basis
  use m_pawtab
  implicit none
  integer,intent(in) :: basis_size
  integer,intent(in) :: mparam
  type(wvlpaw_type),intent(inout) :: wvl
  integer,intent(in) :: nparam(basis_size)
  real(dp),intent(in) :: param(mparam,basis_size)
 end subroutine wvlpaw_sin2gauss
end interface

end module interfaces_65_psp
!!***
