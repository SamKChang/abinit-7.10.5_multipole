!{\src2tex{textfont=tt}}
!!****p* ABINIT/macroave
!! NAME
!! macroave
!!
!! FUNCTION
!! **********************************************************************
!! The MACROAVE program implements the macroscopic average technique,
!! introduced by A. Baldereschi and coworkers
!! (A. Baldereschi, S. Baroni, and R. Resta, Phys. Rev. Lett. 61, 734 (1988)).
!! This is an extremely powerful method that relates
!! microscopic quantities, typical outputs of first-principles codes,
!! with macroscopic magnitudes, needed to perform electrostatic analysis.
!! Within this methodology, we will be able to wash out all the
!! wiggles of the rapidly-varying functions of position (resembling
!! the underlying atomic structure) of the microscopic quantities,
!! blowing up only the macroscopic features.
!! It can be used to compute band offsets, work functions, effective
!! charges, and high frequency dielectric constants, among others
!! interesting physical properties.
!! Ref: L. Colombo, R. Resta and S. Baroni  Phys Rev B  44, 5572 (1991).
!! Coded by P. Ordejon and J. Junquera, April 1999.
!! Modified by J. Junquera, November 2001.
!! **********************************************************************
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 (P. Ordejon, J. Junquera, J. Soler, A. Garcia)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,four1,hdr_free,hdr_io,io_assign,io_close,iorho
!!      macroav_spline,macroav_splint,thetaft
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program macroave

 use defs_basis
 use defs_abitypes

 use m_header,          only : hdr_free, hdr_io
 use m_xmpi,            only : xmpi_world
 use m_io ! in 21_psiesta_noabirule/m_pspsiesta_io.F90 ... implements as a module the 01_macroave/io.F90. It's a copy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'macroave'
 use interfaces_01_macroavnew_ext
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
! --------- PARAMETERS -------------------------------------------------
! INTEGER   NP    : Parameter needed to define maximum number of points
!              for FFT grid.
!              Number of points = (2**NP)
! INTEGER   N     : Maximum number of complex data point for FFT.
!              MUST be a power of 2
! REAL*8   HARTREE: Conversion factor from Hartrees to Rydbergs
!              1 hartree = 2 Ry
! REAL*8   RYDBERG: Conversion factor from Rydbergs to eV
!              1 Ry = 13.6058 eV
! ----------------------------------------------------------------------
 integer, parameter ::  np=12
 integer, parameter ::  n=2**np
 real(dp), parameter :: hartree=two
 real(dp), parameter :: rydberg=13.6058d0
! --------- VARIABLES --------------------------------------------------
 integer :: i,ii,ij,ip,is,j
 integer :: nconv,npoints,npt,nsm,nspin,nz
 integer :: unit1,unit2,unit3,unit4
 integer :: mesh(3)
 character(len=10) :: code,interp
 character(len=15) :: SNAME,inpdata
 character(len=26) :: fnamerho
 character(len=26) :: fnamedelv
 character(len=26) :: fnameplave
 logical :: siesta,abinit,potential,charge,totalcharge
 logical :: found,linear,splin
 real,allocatable :: rhos(:,:)
 real(dp),allocatable :: rho(:,:)
 real(dp) :: cell(3,3),dcell(3,3)
 real(dp) :: l,sur,ds,length,convfac,qren,qtot,lav1,lav2,vol
 real(dp),allocatable :: z(:),rhoz(:),d2rhoz(:),drhodz(:)
 real(dp) :: data(2*n),th(2*n),re(n),im(n)!,v(2*n)
 real(dp) :: x,delta,yp1,ypn,phi
 complex(dp) :: a,b,c
! ABINIT variables
 integer :: fform,rdwr!,unitfi
 type(hdr_type) :: hdr
! end ABINIT variables

!************************************************************************
!CHARACTER CODE       : First principles-code used to generate the
!electrostatic potential at the points of a grid
!in real space. It is read from file 'macroave.in'
!CHARACTER SNAME      : System Label
!(If code = ABINIT, then SNAME = FNAMERHO)
!CHARACTER INPDATA    : Calculate the band offset from the charge
!density or from the electrostatic potential?
!CHARACTER FNAMERHO   : Name of the file where the electrostatic
!potential at the mesh is stored
!CHARACTER FNAMEPLAVE : Name of the file where the planar average of the
!electronic potential will be stored
!CHARACTER FNAMEDELV  : Name of the file where the profile
!of the electrostatic potential will be stored
!LOGICAL   SIESTA     : Have you used SIESTA to get the electrostatic potential
!or the charge density?
!LOGICAL   ABINIT     : Have you used ABINIT to get the electrostatic potential
!or the charge density?
!LOGICAL   LINEAR     : Linear interpolation to get the charge
!density/potential in the FFT grid.
!LOGICAL   SPLIN      : Cubic spline interpolation to get the charge
!density/potential in the FFT grid.
!LOGICAL   POTENTIAL  : We are going to compute the band offset from
!the electrostatic potential
!LOGICAL   CHARGE     : We are going to compute the band offset from
!the charge density
!LOGICAL   FOUND      : Were data found? (only when task in iorho ='read')
!INTEGER   NATOMS     : Number of atoms in unit cell
!INTEGER   NSPIN      : Number of spin polarizations (1 or 2)
!INTEGER   MESH(3)    : Number of mesh divisions of each lattice vectors,
!INCLUDING subgrid
!INTEGER   NSM        : Number of sub-mesh points per mesh point
!(not used in this version)
!INTEGER   NPOINTS    : Number of mesh subdivisions in the normal direction
!to the interface
!INTEGER   NCONV      : Number of convolutions required to calculate the
!macroscopic average
!INTEGER   NPT        : Total number of mesh points (included subpoints)
!REAL*8    CELL(3,3)  : Unit cell lattice vectors (a.u.) CELL(IXYZ,IVECT)
!REAL*8    DS         : Differential area per point of the mesh
!REAL*8    SUR        : Area of a plane parallel to the interface
!REAL*8    LENGTH     : Distance between two planes parallel to the interface
!REAL*8    L          : Length of the cell in the direction nomal
!to the interface (a.u.)
!REAL*8    LAV1       : Linear period of the electrostatic potential
!in the bulklike region for the first material
!REAL*8    LAV2       : Linear period of the electrostatic potential
!in the bulklike region for the second material
!REAL*8    CONVFAC    : Conversion factor for the output units
!REAL*8    QTOT       : Total electronic charge in the unit cell
!REAL*8    QREN       : Total electronic charge calculated from
!the input data file
!REAL*4    Z(MESH(3)) : Z-coordinate of the planes where the elec. density
!is averaged
!REAL*4    RHO        : Electron density
!Notice single precision in this version
!REAL*8    RHOZ       : Planar average of the electron density
!REAL*8    DATA(2*N)  : Fourier coefficients of the planar average density
!REAL*8    TH(2*N)    : Fourier coefficients of the step functions
!REAL*8    V(2*N)     : Fourier coefficients of the potential
!REAL*4    VREEC(N)   : Real part of the electronic potential in real space
!REAL*4    VIMEC(N)   : Imag. part of the electronic potential in real space
!REAL*4    VRENC(N)   : Real part of the nuclear potential in real space
!REAL*4    VIMNC(N)   : Imag. part of the nuclear potential in real space
!*********************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!Reading input data from a file ---------------------------------------
 CALL IO_ASSIGN(UNIT1)
 OPEN(UNIT=UNIT1, FILE='macroave.in', STATUS='OLD')
 READ(UNIT1,'(A)')CODE
 READ(UNIT1,'(A)')INPDATA
 READ(UNIT1,'(A)')SNAME
 READ(UNIT1,*)NCONV
 READ(UNIT1,*)LAV1
 READ(UNIT1,*)LAV2
 READ(UNIT1,*)QTOT
 READ(UNIT1,*)INTERP
 CALL IO_CLOSE(UNIT1)

!Which code has been used to get the electrostatic potential? ---------
 if ( CODE == 'siesta' .OR. CODE == 'SIESTA' .OR.&
& CODE == 'Siesta' ) then
   SIESTA = .TRUE.
   ABINIT = .FALSE.
 else if ( CODE == 'abinit' .OR. CODE == 'ab-init' .OR.&
&   CODE == 'ABINIT' .OR. CODE == 'AB-INIT' .OR.&
&   CODE == 'Abinit' .OR. CODE == 'Ab-init') then
   SIESTA = .FALSE.
   ABINIT = .TRUE.
 else
   write(std_out,*) 'macroave: Unknown code ', CODE
   STOP
 end if

!Are we going to compute the band offset from the charge density or
!from the electrostatic potential? ------------------------------------
 if ( INPDATA == 'potential' .OR. INPDATA == 'POTENTIAL' .OR.&
& INPDATA == 'Potential' ) then
   POTENTIAL   = .TRUE.
   CHARGE      = .FALSE.
   TOTALCHARGE = .FALSE.
 else if ( INPDATA == 'charge' .OR. INPDATA == 'CHARGE' .OR.&
&   INPDATA == 'Charge' ) then
   POTENTIAL   = .FALSE.
   CHARGE      = .TRUE.
   TOTALCHARGE = .FALSE.
 else if ( INPDATA == 'totalcharge' .OR. &
&   INPDATA == 'Totalcharge' .OR.&
&   INPDATA == 'TotalCharge' .OR.&
&   INPDATA == 'TOTALCHARGE' ) then
   POTENTIAL   = .FALSE.
   CHARGE      = .FALSE.
   TOTALCHARGE = .TRUE.
 else
   write(std_out,*) 'macroave: Unknown input data  ', INPDATA
   STOP
 end if

!What kind of interpolation will we use to get the charge density/
!potential in a FFT grid? ---------------------------------------------
 if ( INTERP == 'linear' .OR. INTERP == 'Linear' .OR.&
& INTERP == 'LINEAR' ) then
   LINEAR = .TRUE.
   SPLIN  = .FALSE.
 else if ( INTERP == 'spline' .OR. INTERP == 'Spline' .OR.&
&   INTERP == 'SPLINE' ) then
   LINEAR = .FALSE.
   SPLIN  = .TRUE.
 end if

!Reading charge density from a file -----------------------------------
 if ( SIESTA ) then
   if (POTENTIAL) then
     FNAMERHO = PASTE(trim(SNAME),'.VH')
   elseif (CHARGE) then
     FNAMERHO = PASTE(trim(SNAME),'.RHO')
   elseif (TOTALCHARGE) then
     FNAMERHO = PASTE(trim(SNAME),'.TOCH')
   end if
 else if ( ABINIT ) then
   FNAMERHO = trim(SNAME)
 end if

 if (SIESTA) then
   NSM   = 1
   NPT   = 0
   NSPIN = 0
   CALL IORHO( 'READ', trim(FNAMERHO), DCELL, MESH, NSM, NPT, NSPIN,&
&   RHOS, FOUND )
   if (FOUND) then
     ABI_ALLOCATE( RHOS,(NPT,NSPIN))
     ABI_ALLOCATE( RHO,(NPT,NSPIN))
     CALL IORHO( 'READ', trim(FNAMERHO), DCELL, MESH, NSM, NPT, NSPIN,&
&     RHOS, FOUND )
     do I = 1, 3
       do J = 1, 3
         CELL(J,I) = DCELL(J,I)
       end do
     end do
!    Transform the density or the potential read from SIESTA
!    from a single precision variable to a double precision variable
     do IS = 1, NSPIN
       do IP = 1, NPT
         RHO(IP,IS) = RHOS(IP,IS) * 1.0D0
       end do
     end do

   else
     write(std_out,*)'macroave: ERROR: file not found: ', FNAMERHO
     STOP
   end if

 else if (ABINIT) then
   CALL IO_ASSIGN(UNIT2)
   OPEN(UNIT=UNIT2, FILE=trim(FNAMERHO), FORM='UNFORMATTED',STATUS='OLD')
   RDWR = 1
   call hdr_io(FFORM,HDR,RDWR,UNIT2)
!  For debugging
!  rdwr=4 ; unitfi=6
!  call hdr_io(fform,hdr,rdwr,unitfi)

   do I = 1, 3
     MESH(I) = HDR%NGFFT(I)
     do J = 1, 3
       CELL(J,I) = HDR%RPRIMD(J,I)
     end do
   end do
   NSPIN = HDR%NSPPOL
   NPT = MESH(1) * MESH(2) * MESH(3)
   ABI_ALLOCATE( RHO,(NPT,NSPIN))

   do IS = 1, NSPIN
     READ(UNIT2) (RHO(IP,IS),IP=1,NPT)
!    Units for the potential in Ab-init are in Hartrees,
!    so we transform them into Ry. No transformation is
!    needed for the charge density
!    (it is directly read in electrons/bohr**3).
     if (POTENTIAL) then
       do IP = 1, NPT
         RHO(IP,IS) = RHO(IP,IS) * HARTREE
       end do
     end if
   end do
   CALL IO_CLOSE(UNIT2)

   call hdr_free(hdr)
 end if

!Initialize some variables (we suppose orthorombic cells) -------------

 L  = CELL(3,3)
 SUR = SURPLA( CELL )
 VOL = VOLCEL( CELL )
 DS = SUR/( MESH(1) * MESH(2) )
 LENGTH = L/DBLE(N)
 NPOINTS = MESH(3)

 ABI_ALLOCATE(Z,(NPOINTS+1))
 ABI_ALLOCATE(RHOZ,(NPOINTS+1))
 ABI_ALLOCATE(D2RHOZ,(NPOINTS+1))
 ABI_ALLOCATE(DRHODZ,(N))

 RHOZ(1:NPOINTS+1)   = 0.D0
 D2RHOZ(1:NPOINTS+1) = 0.D0
 DRHODZ(1:N)         = 0.D0

 if (POTENTIAL) then
   CONVFAC = RYDBERG
 else if (CHARGE) then
   CONVFAC = 1.0D0
 else if (TOTALCHARGE) then
   CONVFAC = 1.0D0
 end if


!Loop over all points and calculate the planar average ----------------
!Warning: The planar average is only done for the first component of
!RHO. Spin polarization is not implemented yet ---------------
 do IP = 1, NPT
   NZ = (IP-1) / (MESH(1)*MESH(2)) + 1
   RHOZ(NZ) =  RHOZ(NZ) + RHO(IP,1)*DS
 end do

 do IP = 1, NPOINTS
   RHOZ(IP) = RHOZ(IP) / SUR
 end do

 do IP  = 1, NPOINTS
   Z(IP) = (IP-1)*CELL(3,3)/DBLE(NPOINTS)
 end do

!Calculate electrostatic potential or electronic charge density -------
!in fft grid, interpolating the planar average calculated before ------
 if (SPLIN) then
   Z(NPOINTS+1)    = L
   RHOZ(NPOINTS+1) = RHOZ(1)
   DELTA = L/DBLE(NPOINTS)
   YP1 = ( RHOZ(2) - RHOZ(NPOINTS) )   / (2.0D0*DELTA)
   YPN = YP1
   CALL MACROAV_SPLINE(DELTA, RHOZ, NPOINTS+1, YP1, YPN, D2RHOZ)
   I = 0
   do II = 1, 2*N-1, 2
     I = I + 1
     X = (I-1)*L/DBLE(N)
     CALL MACROAV_SPLINT( DELTA, RHOZ, D2RHOZ, NPOINTS+1, X, DATA(II), &
&     DRHODZ(I) )
     DATA(II+1) = 0.D0
   end do
 else if (LINEAR) then
   I = 0
   do II = 1,2*N-1,2
     I = I + 1
     X = (I-1)*L/DBLE(N)
     do IJ = 1, NPOINTS
       if (X == Z(IJ)) then
         DATA(II) = RHOZ(IJ)
         DATA(II+1) = 0.D0
         GOTO 20
       end if
       if (Z(IJ) > X) then
         DATA(II) = RHOZ(IJ-1) +&
&         (X-Z(IJ-1))*(RHOZ(IJ)-RHOZ(IJ-1))/&
&         (Z(IJ)-Z(IJ-1))
         DATA(II+1) = 0.D0
         GOTO 20
       end if
     end do
     DATA(II)=RHOZ(NPOINTS) +&
&     (X-Z(NPOINTS))*(RHOZ(1)-RHOZ(NPOINTS))/&
&     (Z(NPOINTS)-Z(NPOINTS-1))
     DATA(II+1) = 0.D0
     20       CONTINUE
   end do
 end if

!Renormalize the charge density ---------------------------------------
 if (CHARGE .OR. TOTALCHARGE) then
   QREN = 0.D0
   do IP = 1, 2*N-1, 2
     QREN = QREN + DATA(IP)*LENGTH*SUR
   end do
   do IP = 1, 2*N-1, 2
     if (CHARGE) then
       DATA(IP) = DATA(IP) * QTOT/QREN
     elseif(TOTALCHARGE) then
       DATA(IP) = DATA(IP) - QREN/VOL
     end if
   end do
   QREN = 0.D0
   do IP = 1, 2*N-1, 2
     QREN = QREN + DATA(IP)*LENGTH*SUR
   end do
!  For debugging
!  write(std_out,*)' QREN = ', QREN
 end if
!...

!Print planar average of the electrostatic potential or ---------------
!the electronic charge density ----------------------------------------
 FNAMEPLAVE = PASTE(SNAME,'.PAV')
 CALL IO_ASSIGN(UNIT3)
 OPEN(UNIT=UNIT3, FILE=FNAMEPLAVE,STATUS='UNKNOWN')
 I = 0
 do II = 1, 2*N-1, 2
   I = I+1
   X=(I-1)*L/DBLE(N)
!  WRITE(UNIT3,'(3F20.12)')X,
!  .           DATA(II)*CONVFAC,DATA(II+1)*CONVFAC
   WRITE(UNIT3,'(2F20.12)')X,&
&   DATA(II)*CONVFAC
 end do
 CALL IO_CLOSE(UNIT3)
!...


!Calculate Fourier transform of the electrostatic potential or
!the electronic density -----------------------------------------------
 CALL FOUR1(DATA,N,1)
!...

!Calculate macroscopic average of the electrostatic potential or the
!electronic charge density taking the convolution with two step functions.
!In Fourier space, it is a product of the Fourier transform components -
!The decompositions in the sum over II is due to the spetial way in which
!the data are stored in subroutine four1( see Fig. 12.2.2, in
!'Numerical Recipes, The Art of Scientific Computing'
!by W.H. Press, S.A. Teukolsky, W.T. Veterling and B.P. Flannery,
!Cambridge U.P. 1987 and 1992.


 CALL THETAFT(N,L,LAV1,TH)

 do II = 1, N+1, 2
   A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
   B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
   C = A*B
   DATA(II) = REAL(C)*L/DBLE(N)
   DATA(II+1) = AIMAG(C)*L/DBLE(N)
 end do

 do II = N+3, 2*N-1, 2
   A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
   B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
   C = A*B
   DATA(II) = REAL(C)*L/DBLE(N)
   DATA(II+1) = AIMAG(C)*L/DBLE(N)
 end do


 if (NCONV == 2) then
   CALL THETAFT(N,L,LAV2,TH)

   do II = 1, N+1, 2
     A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
     B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
     C = A*B
     DATA(II) = REAL(C)*L/DBLE(N)
     DATA(II+1) = AIMAG(C)*L/DBLE(N)
!    if ( POISON ) then
!    IG = (II-1) / 2
!    GSQ= (2.D0*PI*IG/L)**2
!    if(GSQ > 0.D0) then
!    V(II) = DATA(II) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
!    V(II+1) = DATA(II+1) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
!    else
!    V(II) = 0.D0
!    V(II+1) = 0.D0
!    endif
!    endif
   end do

   do II = N+3, 2*N-1, 2
     A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
     B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
     C = A*B
     DATA(II) = REAL(C)*L/DBLE(N)
     DATA(II+1) = AIMAG(C)*L/DBLE(N)
!    if ( POISON ) then
!    IG = (-2*N+II-1) / 2
!    GSQ= (2.D0*PI*IG/L)**2
!    if(GSQ > 0.D0) then
!    V(II) = DATA(II) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
!    V(II+1) = DATA(II+1) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
!    else
!    V(II) = 0.D0
!    V(II+1) = 0.D0
!    endif
!    endif
   end do

 end if
!...

!Transform average electronic density and potential to real space -----
!The decompositions in the sum over J is due to the spetial way in which
!the data are stored in subroutine four1( see Fig. 12.2.2, in
!'Numerical Recipes, The Art of Scientific Computing'
!by W.H. Press, S.A. Teukolsky, W.T. Veterling and B.P. Flannery,
!Cambridge U.P. 1987 and 1992.

 do II = 1, N
   RE(II) = 0.D0
   IM(II) = 0.D0
!  if ( POISON ) then
!  VREEC(II) = 0.D0
!  VIMEC(II) = 0.D0
!  endif
   do J = 1, N+1, 2
     PHI = -2.D0 * PI * (II-1) * ( (J-1)/2 ) / DBLE(N)
     RE(II)=RE(II)+(1.D0/DBLE(N))*(DATA(J)*COS(PHI)&
&     -DATA(J+1)*SIN(PHI))
     IM(II)=IM(II)+(1.D0/DBLE(N))*(DATA(J)*SIN(PHI)&
&     +DATA(J+1)*COS(PHI))
!    if ( POISON ) then
!    VREEC(II)=VREEC(II)+(1.D0/DBLE(N))*(V(J)*COS(PHI)
!    .                           -V(J+1)*SIN(PHI))
!    VIMEC(II)=VIMEC(II)+(1.D0/DBLE(N))*(V(J)*SIN(PHI)
!    .                           +V(J+1)*COS(PHI))
!    endif
   end do

   do J = N+3, 2*N-1, 2
     PHI = -2.0D0 * PI * (II-1) * ((-2*N+J-1)/2) / DBLE(N)
     RE(II)=RE(II)+(1.D0/DBLE(N))*(DATA(J)*COS(PHI)&
&     -DATA(J+1)*SIN(PHI))
     IM(II)=IM(II)+(1.D0/DBLE(N))*(DATA(J)*SIN(PHI)&
&     +DATA(J+1)*COS(PHI))
!    if ( POISON ) then
!    VREEC(II)=VREEC(II)+(1.D0/DBLE(N))*(V(J)*COS(PHI)
!    .                           -V(J+1)*SIN(PHI))
!    VIMEC(II)=VIMEC(II)+(1.D0/DBLE(N))*(V(J)*SIN(PHI)
!    .                           +V(J+1)*COS(PHI))
!    endif
   end do
 end do
!...

!Print averaged electronic charge density and potential ---------------
 FNAMEDELV = PASTE( trim(SNAME),'.MAV')
 CALL IO_ASSIGN(UNIT4)
 OPEN(UNIT=UNIT4, FILE=trim(FNAMEDELV), STATUS='UNKNOWN')
 do I = 1, N
   X=(I-1)*L/DBLE(N)
!  WRITE(UNIT4,'(3F20.5)')X,
!  .           RE(I)*CONVFAC ,IM(I)*CONVFAC
   WRITE(UNIT4,'(2F20.12)')X,&
&   RE(I)*CONVFAC
 end do
 CALL IO_CLOSE(UNIT4)
!...

!Print electrostatic potential ----------------------------------------
!if (POISON) then
!FNAMEVEC = PASTE( SNAME,'.VEC')
!CALL IO_ASSIGN(UNIT5)
!OPEN(UNIT=UNIT5, FILE=FNAMEVEC, STATUS='UNKNOWN')
!do I = 1, N
!X=(I-1)*L/DBLE(N)
!WRITE(UNIT5,'(3F20.12)')X, VREEC(I), VIMEC(I)
!enddo
!CALL IO_CLOSE(UNIT5)
!endif
!...

 ABI_DEALLOCATE(Z)
 ABI_DEALLOCATE(RHOZ)
 ABI_DEALLOCATE(DRHODZ)
 ABI_DEALLOCATE(D2RHOZ)

 end program macroave
!!***
