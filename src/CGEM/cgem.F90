module cgem

!CGEM STATE VARIABLES
use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use grid, only: km
use schism_glbl, only : nea,rkind

implicit none

save

!Write or not
logical sinkwcgem,adjust_ws
real*8 :: adjust_fac

!number of species
integer :: nospA
integer :: nospZ

!misc
integer :: skipcgem,checkwindrad,writecsv,debug
real :: eps,h_massconsv

!constants
real, parameter :: SDay = 86400.

!Sinking
real*8, dimension(:), allocatable :: ws
real, dimension(:), allocatable :: fmin
real :: x !for tiny

!Fixed Stoichiometry:
real :: sx1R,sy1R,sx2R,sy2R
real :: sx1BC,sy1BC,sx2BC,sy2BC

!Module Which_Flux
! =========================================================
! Define which fluxes shall be used
! =========================================================
integer, parameter :: iO2surf  = 1 !O2 surface flux
integer, parameter :: iDICsurf = 2 !DIC surface flux
integer, parameter :: iSOC     = 3 !Sediment Oxygen Consumption
integer, parameter :: iMPB     = 4 !Microphytobethos
integer, parameter :: iNutEx   = 5 !Sediment Nutrient Fluxes
integer, parameter :: iCMAQ    = 6 !CMAQ surface deposition of NH4 and NO3
integer, parameter :: iInRemin = 7 !Instant Remineralization in bottom layer
integer, parameter :: iSDM     = 8 !Sediment Diagenesis Model
integer, parameter :: i_Si      = 9 !Silica (SA, SRP) Fluxes

!Module CGEM_Flux
! =========================================================
! Terms for Flux Calculations
! =========================================================
       REAL :: dT_sed
       REAL,ALLOCATABLE :: pH(:)

!---------------------------------------------------------      
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------  
      integer, dimension(:), allocatable :: iA(:)
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
      integer, dimension(:), allocatable :: iQn(:)
!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
      integer, dimension(:), allocatable :: iQp(:)
!--------------------------------------------------------------------
!-Z: Zooplankton number density (individuals/m3);
!--------------------------------------------------------------------
      integer, dimension(:), allocatable :: iZ(:)
!-------------------------------
!-NO3; Nitrate (mmol-N/m3)
!-------------------------------
      integer :: iNO3
!--------------------------------      
!-NH4; Ammonium (mmol-N/m3)
!--------------------------------
      integer :: iNH4
!-------------------------------------------        
!-PO4: Phosphate (mmol-P/m3)
!--------------------------------------
      integer :: iPO4 
!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!---------------------------------------------------------
      integer :: iDIC 
!-------------------------------------------        
!-O2: Molecular Oxygen (mmol-O2/m3)
!------------------------------
      integer :: iO2 
!-------------------------------------------------------------
!-OM1_A: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton
!-------------------------------------------------------------
      integer :: iOM1_A
!-----------------------------------------------------------------
!-OM2_A: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton 
!------------------------------------------------------------------
      integer :: iOM2_A
!-------------------------------------------------------------
!-OM1_Z:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets.
!-------------------------------------------------------------
      integer :: iOM1_Z
!-------------------------------------------------        
!-OM2_Z:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets.
!-----------------------------------------------
      integer :: iOM2_Z
!--------------------------------------------------------------------
!-OM1_R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!--------------------------------------------------------------------
      integer :: iOM1_R
!-------------------------------------------------      
!-OM2_R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!--------------------------------------------------------------------
      integer :: iOM2_R
!-------------------------------------------
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-------------------------------------------
      integer :: iCDOM
!---------------------------------------------
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-------------------------------------------
      integer :: iSi
!--------------------------------------------------------------------
!-OM1_BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!--------------------------------------------------------------------
      integer :: iOM1_BC
!-------------------------------------------------
!-OM2_BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!--------------------------------------------------------------------
      integer :: iOM2_BC
!-------------------------------------------
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!-------------------------------------------
      integer :: iALK
!Tracer
      integer :: iTr
!Stoichiometry
      integer :: isx1A
      integer :: isy1A
      integer :: isx2A
      integer :: isy2A
      integer :: isx1Z
      integer :: isy1Z
      integer :: isx2Z
      integer :: isy2Z

!Total number of state variables
      integer :: nf

!State Variable Array
      real(rkind),allocatable :: ff(:,:) !state variable array
      real(rkind), allocatable :: ff_new(:,:) !sources array

!----INPUT_VARS_CGEM
!--Switches in GEM---------
integer Which_fluxes(8)
integer Which_uptake
integer Which_quota
integer Which_irradiance
integer Which_photosynthesis
integer Which_growth
integer Which_temperature
integer Which_wind
integer Which_rad
!--Temperature
real, allocatable :: KTg1(:)
real, allocatable :: KTg2(:)
real, allocatable :: Tref(:)
real, allocatable :: Ea(:)
!real, allocatable :: N(:)
!--Optics-----------------------
real Kw
real Kcdom
real Kspm
real Kchla
!--in module LIGHT_VARS
real astar490
real aw490
real astarOMA
real astarOMZ
real astarOMR
real astarOMBC
real PARfac
real sinkCDOM
!---Phytoplankton 
real, allocatable :: ediblevector(:,:)
real, allocatable :: umax(:)
real, allocatable :: CChla(:)
real, allocatable :: alpha(:)
real, allocatable :: beta(:)
real, allocatable :: respg(:)
real, allocatable :: respb(:)
real, allocatable :: QminN(:)
real, allocatable :: QminP(:)
real, allocatable :: QmaxN(:)
real, allocatable :: QmaxP(:)
real, allocatable :: Kn(:)
real, allocatable :: Kp(:)
real, allocatable :: Ksi(:)
real, allocatable :: KQn(:)
real, allocatable :: KQp(:)
real, allocatable :: nfQs(:)
real, allocatable :: vmaxN(:)
real, allocatable :: vmaxP(:)
real, allocatable :: vmaxSi(:)
real, allocatable :: aN(:)
real, allocatable :: volcell(:)
real, allocatable :: Qc(:)
real, allocatable :: Athresh(:)
real, allocatable :: mA(:)
real, allocatable :: A_wt(:)
!---Zooplankton
real, allocatable :: Zeffic(:)
real, allocatable :: Zslop(:)
real, allocatable :: Zvolcell(:)
real, allocatable :: ZQc(:)
real, allocatable :: ZQn(:)
real, allocatable :: ZQp(:)
real, allocatable :: ZKa(:)
real, allocatable :: Zrespg(:)
real, allocatable :: Zrespb(:)
real, allocatable :: Zumax(:)
real, allocatable :: Zm(:)
!---Organic Matter              
real KG1
real KG2
real KG1_R
real KG2_R
real KG1_BC
real KG2_BC
real KNH4
real nitmax
real KO2
real KstarO2
real KNO3
real pCO2
real KGcdom
real CF_SPM
!----Other including Boundary Conditions------------
real KH_coeff
integer Which_Outer_BC
real wt_l, wt_o
real wt_pl, wt_po
real m_OM_init,m_OM_BC,m_OM_sh
real KG_bot

!Light curve parameters
real, allocatable :: alphad(:) ! Initial slope of photosynthesis-irradiance curve / Vmax
real, allocatable :: betad(:)  ! Photoinhibition constant / Vmax

!Initialize variables in cgem_init from cgem_read
  real, dimension(:), allocatable, private :: sinkA
  real, dimension(:), allocatable, private :: A_init,Qn_init,Qp_init
  real, dimension(:), allocatable, private :: Z_init
  real, private :: sinkOM1_A,sinkOM2_A,sinkOM1_Z,sinkOM2_Z,sinkOM1_R,sinkOM2_R,sinkOM1_BC,sinkOM2_BC
  real, private :: NO3_init,NH4_init,PO4_init,DIC_init,O2_init
  real, private :: OM1_A_init,OM2_A_init,OM1_Z_init,OM2_Z_init,OM1_R_init,OM2_R_init,CDOM_init
  real, private :: Si_init,OM1_BC_init,OM2_BC_init,ALK_init,Tr_init

contains

subroutine cgem_dim

  integer           :: istat,iunit
  character(len=1000) :: line
  !http://degenerateconic.com/namelist-error-checking.html
  namelist /nosp/ nospA,nospZ,skipcgem,checkwindrad,writecsv,debug,sinkwcgem,adjust_ws,adjust_fac,h_massconsv

#ifdef DEBUG
write(6,*) "Begin cgem_dim"
#endif

open(action='read',file='cgem.nml',iostat=istat,newunit=iunit)

!namelist /switches/
read(nml=nosp,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist nosp: '//trim(line)
 stop
endif

close(iunit)

if(nospZ.ne.1.and.nospZ.ne.2) then
!write(6,*) "Z's Please, are 1 or 2 for now.  You chose:",nospZ
!write(6,*) "We're going to use nosp2=2 and keep going."
nospZ=2
endif

return
end subroutine cgem_dim



subroutine cgem_read

  integer           :: istat,iunit
  character(len=1000) :: line
  !http://degenerateconic.com/namelist-error-checking.html
  namelist /switches/ Which_fluxes,Which_temperature,Which_uptake,Which_quota,Which_irradiance,&
    Which_photosynthesis,Which_growth,Which_wind,Which_rad
  namelist /sdm/ dT_sed
  namelist /optics/ Kw,Kcdom,Kspm,Kchla,astar490,aw490,astarOMA,astarOMZ,astarOMR,astarOMBC,PARfac,sinkCDOM
  namelist /temperature/ Tref,KTg1,KTg2,Ea
  namelist /phytoplankton/ umax,CChla,alpha,beta,respg,respb,QminN,QminP,QmaxN,QmaxP,Kn,Kp,Ksi,KQn,&
     KQp,nfQs,vmaxN,vmaxP,vmaxSi,aN,volcell,Qc,Athresh,sinkA,mA,A_wt
  namelist /zooplankton/ ediblevector,Zeffic,Zslop,Zvolcell,ZQc,ZQn,ZQp,ZKa,Zrespg,Zrespb,Zumax,Zm
  namelist /OM/ KG1,KG2,KG1_R,KG2_R,KG1_BC,KG2_BC,KNH4,nitmax,KO2,KstarO2,KNO3,pCO2,&
     sx1R,sy1R,sx2R,sy2R,&
     sx1BC,sy1BC,sx2BC,sy2BC,sinkOM1_A,sinkOM2_A,sinkOM1_Z,sinkOM2_Z,sinkOM1_R,&
     sinkOM2_R,sinkOM1_BC,sinkOM2_BC,KGcdom,CF_SPM,KG_bot

#ifdef DEBUG
write(6,*) "Begin cgem_read"
#endif

open(action='read',file='cgem.nml',iostat=istat,newunit=iunit)

!namelist /switches/
read(nml=switches,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist: '//trim(line)
 stop
endif

!namelist /switches/
read(nml=sdm,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist: '//trim(line)
 stop
endif

!namelist /optics/
read(nml=optics,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist optics: '//trim(line)
 stop
endif

!namelist /temperature/
read(nml=temperature,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist temperature: '//trim(line)
 stop
endif

!namelist /phytoplankton/
read(nml=phytoplankton,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist phytoplankton: '//trim(line)
 stop
endif

!namelist /zooplankton/
read(nml=zooplankton,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist zooplankton: '//trim(line)
 stop
endif

!namelist /OM/
read(nml=OM,iostat=istat,unit=iunit)
if (istat /= 0) then
 backspace(iunit)
 read(iunit,fmt='(A)') line
 write(6,'(A)') &
        'Invalid line in namelist OM: '//trim(line)
 stop
endif

!namelist /init/
!read(nml=init,iostat=istat,unit=iunit)
!if (istat /= 0) then
! backspace(iunit)
! read(iunit,fmt='(A)') line
! write(6,'(A)') &
!        'Invalid line in namelist init: '//trim(line)
! stop
!endif


close(iunit)


return
end subroutine cgem_read

subroutine cgem_allocate()

integer i,ierr
integer :: counter = 0

#ifdef DEBUG
write(6,*) "Begin cgem_allocate" 
#endif

!---------------------------------------------------------      
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------  
       allocate (iA(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iA"
       do i=1,nospA
          counter = counter+1
          iA(i) = counter
       enddo
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
       allocate (iQn(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iQn"
       do i=1,nospA
          counter = counter+1
          iQn(i) = counter
       enddo
!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
      allocate (iQp(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iQp"
       do i=1,nospA
          counter = counter+1
          iQp(i) = counter
       enddo
!--------------------------------------------------------------------
!-Z: Zooplankton number density (individuals/m3);
!--------------------------------------------------------------------
      allocate (iZ(nospZ),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iZ"
       do i=1,nospZ
          counter = counter+1
          iZ(i) = counter
       enddo
!-------------------------------
!-NO3; Nitrate (mmol-N/m3)
!-------------------------------
      iNO3 = counter+1
!--------------------------------      
!-NH4; Ammonium (mmol-N/m3)
!--------------------------------
      iNH4 = counter+2
!-------------------------------------------        
!-PO4: Phosphate (mmol-P/m3)
!--------------------------------------
      iPO4 = counter+3
!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!---------------------------------------------------------
      iDIC = counter+4
!-------------------------------------------        
!-O2: Molecular Oxygen (mmol-O2/m3)
!------------------------------
      iO2 = counter+5
!-------------------------------------------------------------
!-OM1_A: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton
!-------------------------------------------------------------
      iOM1_A = counter+6
!-----------------------------------------------------------------
!-OM2_A: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton 
!------------------------------------------------------------------
      iOM2_A = counter+7
!-------------------------------------------------------------
!-OM1_Z:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets.
!-------------------------------------------------------------
      iOM1_Z = counter+8
!-------------------------------------------------        
!-OM2_Z:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets.
!-----------------------------------------------
      iOM2_Z = counter+9
!--------------------------------------------------------------------
!-OM1_R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!--------------------------------------------------------------------
      iOM1_R = counter+10
!-------------------------------------------------      
!-OM2_R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!--------------------------------------------------------------------
      iOM2_R = counter+11
!-------------------------------------------
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-------------------------------------------
      iCDOM = counter+12
!---------------------------------------------
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-------------------------------------------
      iSi = counter+13
!--------------------------------------------------------------------
!-OM1_BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!--------------------------------------------------------------------
      iOM1_BC = counter+14
!-------------------------------------------------
!-OM2_BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!--------------------------------------------------------------------
      iOM2_BC = counter+15
!-------------------------------------------
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!-------------------------------------------
      iALK = counter+16
!Tracer
      iTR = counter+17
!Stoichiometry
      isx1A = counter+18
      isy1A = counter+19
      isx2A = counter+20
      isy2A = counter+21
      isx1Z = counter+22
      isy1Z = counter+23
      isx2Z = counter+24
      isy2Z = counter+25

!How many state variables
      nf = isy2Z

      allocate(ff(km,nf),stat=ierr)
      if(ierr.ne.0) write(6,*) "error in allocating:ff"

      allocate(ff_new(km,nf),stat=ierr)
      if(ierr.ne.0) write(6,*) "error in allocating:ff_new"

      ff = -9999. 
      ff_new = -9999.

!----allocate INPUT_VARS_CGEM

!---Phytoplankton 
allocate( ediblevector(nospZ,nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ediblevector"
allocate( umax(nospA),stat=ierr  )
if(ierr.ne.0) write(6,*) "error in allocating:umax"
allocate( CChla(nospA),stat=ierr  )
if(ierr.ne.0) write(6,*) "error in allocating:CChla"
allocate( alpha(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:alpha"
allocate( beta(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:beta"
allocate( respg(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:respg"
allocate( respb(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:respb"
allocate( QminN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QminN"
allocate( QminP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QminP"
allocate( QmaxN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QmaxN"
allocate( QmaxP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QmaxP"
allocate( Kn(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Kn"
allocate( Kp(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Kp"
allocate( Ksi(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Ksi"
allocate( KQn(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:KQn"
allocate( KQp(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:KQp"
allocate( nfQs(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:nfQs"
allocate( vmaxN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxN"
allocate( vmaxP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxP"
allocate( vmaxSi(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxSi"
allocate( aN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:aN"
allocate( volcell(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:volcell"
allocate( Qc(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Qc"
allocate( Athresh(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Athresh"
allocate( mA(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:mA"
allocate( A_wt(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:A_wt"

!---Zooplankton
allocate( Zeffic(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zeffic"
allocate( Zslop(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zslop"
allocate( Zvolcell(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zvolcell"
allocate( ZQc(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQc"
allocate( ZQn(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQn"
allocate( ZQp(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQp"
allocate( ZKa(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZKa"
allocate( Zrespg(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zrespg"
allocate( Zrespb(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zrespb"
allocate( Zumax(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zumax"
allocate( Zm(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zm"

!Light curve parameters
allocate( alphad(nospA),stat=ierr ) ! Initial slope of photosynthesis-irradiance curve / Vmax
if(ierr.ne.0) write(6,*) "error in allocating:alphad"
allocate( betad(nospA),stat=ierr )  ! Photoinhibition constant / Vmax
if(ierr.ne.0) write(6,*) "error in allocating:betad"


!Temperature parameters for growth rates
allocate(Tref(nospA+nospZ),stat=ierr)                   !Tref(nospA+nospZ): Optimum temperature for growth(C)
if(ierr.ne.0) write(6,*) "error in allocating:Tref"
allocate(KTg1(nospA+nospZ),stat=ierr)                   !KTg1(nospA+nospZ): Effect of T below Topt(C^2)
if(ierr.ne.0) write(6,*) "error in allocating:KTg1"
allocate(KTg2(nospA+nospZ),stat=ierr)                   !KTg2(nospA+nospZ): Effect of T above Topt(C^2)
if(ierr.ne.0) write(6,*) "error in allocating:KTg2"
allocate(Ea(nospA+nospZ),stat=ierr)                     !Ea(nospA+nospZ): Slope of Arrhenius plot(eV)
if(ierr.ne.0) write(6,*) "error in allocating:Ea"

!Sinking
allocate(ws(nf),stat=ierr)
if(ierr.ne.0) write(6,*) "error in allocating:ws"
allocate(fmin(nf),stat=ierr)
if(ierr.ne.0) write(6,*) "error in allocating:fmin"

!pH
  allocate(pH(km),stat=ierr)
  if(ierr.ne.0) write(6,*) "error in allocating:pH"

!sinking
  allocate(sinkA(nospA),stat=ierr)
  if(ierr.ne.0) write(6,*) "error in allocating:sinkA"


#ifdef DEBUG
write(6,*) "End cgem_allocate"
#endif

return
end subroutine cgem_allocate

subroutine cgem_init

integer :: isp
real tot

#ifdef DEBUG
write(6,*) "Begin cgem_init"
#endif


Athresh = Athresh*volcell   ! Threshold for grazing, um^3/m3
eps=0
do isp=1,nospA
   eps=0
   if(umax(isp).eq.0) eps=1.e-18
   alphad(isp) = alpha(isp)/(umax(isp)+eps) ! Initial slope of photosynthesis-irradiance curve / Vmax
   betad(isp)  = beta(isp)/(umax(isp)+eps)  ! Photoinhibition constant / Vmax
enddo

!Convert relative proportions of phytoplankton to percentage of total chlA
tot = SUM(A_wt)
if(tot.le.0) then
 write(6,*) "Error in A_wt, A_wt.le.0"
 stop
endif

do isp=1,nospA
   A_wt(isp) = A_wt(isp)/tot
enddo


!ff is initialized by SCHISM

!namelist /init/ A_init,Qn_init,Qp_init,Z_init,NO3_init,NH4_init,PO4_init,DIC_init,O2_init,&
! OM1_A_init,OM2_A_init,OM1_Z_init,OM2_Z_init,OM1_R_init,OM2_R_init,CDOM_init,Si_init,OM1_BC_init,OM2_BC_init,ALK_init,Tr_init
!Initialize ff for testing
!do isp=1,nospA
! ff(:,iA(isp))  = A_init(isp)
! ff(:,iQn(isp)) = Qn_init(isp)
! ff(:,iQp(isp)) = Qp_init(isp)
!!write(6,*) "Initializing A,isp=",isp,A_init(isp)
!enddo
!do isp=1,nospZ
! ff(:,iZ(isp)) = Z_init(isp)
!! write(6,*) "Initializing Z,isp=",isp,Z_init(isp)  
!enddo
!ff(:,iNO3) = NO3_init
!ff(:,iNH4) = NH4_init 
!ff(:,iPO4) = PO4_init 
!ff(:,iDIC) = DIC_init 
!ff(:,iO2) = O2_init 
!ff(:,iOM1_A) = OM1_A_init
!ff(:,iOM2_A) = OM2_A_init
!ff(:,iOM1_Z) = OM1_Z_init !0. !78.162582                   !OM1_fp 
!ff(:,iOM2_Z) = OM2_Z_init !0. !225.37767                   !OM2_fp 
!ff(:,iOM1_R) = OM1_R_init !0.0000000               !OM1_rp 
!ff(:,iOM2_R) = OM2_R_init !0.0000000               !OM2_rp 
!ff(:,iCDOM) = CDOM_init !2.              !CDOM 
!ff(:,iSi) = Si_init !15.             !Si 
!ff(:,iOM1_BC) = OM1_BC_init !0. !157.09488                   !OM1_bc 
!ff(:,iOM2_BC) = OM2_BC_init !0. !333.65701                   !OM2_bc
!ff(:,iALK) = ALK_init !2134               !ALK 
!ff(:,iTr) = Tr_init !1                  !Tr

pH = -9999.

fmin = tiny(x)
fmin(iA(:)) = 1.
fmin(iZ(:)) = 1.
fmin(iQn(:)) = QminN(:)
fmin(iQp(:)) = QminP(:)


ws = 0.
ws(iA(:))=sinkA(:)
!We didn't do this before, but we should have...
ws(iQn(:)) = sinkA(:)
ws(iQp(:)) = sinkA(:)

#ifdef DEBUG
write(6,*) "ws(ia)",ws(iA(:))
#endif

ws(iCDOM) =  sinkCDOM

ws(iOM1_A) = sinkOM1_A
ws(isx1A) = sinkOM1_A
ws(isy1A) =  sinkOM1_A

ws(iOM2_A) = sinkOM2_A
ws(isx2A) = sinkOM2_A
ws(isy2A) = sinkOM2_A

ws(iOM1_Z) = sinkOM1_Z
ws(isx1Z) = sinkOM1_Z
ws(isy1Z) = sinkOM1_Z

ws(iOM2_Z) = sinkOM2_Z
ws(isx2Z) = sinkOM2_Z
ws(isy2Z) = sinkOM2_Z

ws(iOM1_R) = sinkOM1_R
ws(iOM2_R) = sinkOM2_R

ws(iOM1_BC) = sinkOM1_BC
ws(iOM2_BC) = sinkOM2_BC
#ifdef DEBUG
write(6,*) "ws(iCDOM) =  sinkCDOM",ws(iCDOM)
write(6,*) "ws(iOM1_A) = sinkOM1_A",ws(iOM1_A)
write(6,*) "ws(iOM2_A) = sinkOM2_A",ws(iOM2_A)
write(6,*) "ws(iOM1_Z) = sinkOM1_Z",ws(iOM1_Z)
write(6,*) "ws(iOM2_Z) = sinkOM2_Z",ws(iOM2_Z)
write(6,*) "ws(iOM1_R) = sinkOM1_R",ws(iOM1_R)
write(6,*) "ws(iOM2_R) = sinkOM2_R",ws(iOM2_R)
write(6,*) "ws(iOM1_BC) = sinkOM1_BC",ws(iOM1_BC)
write(6,*) "ws(iOM2_BC) = sinkOM2_BC",ws(iOM2_BC)
write(6,*) "ws",ws
#endif


!Convert per m/d to m/s 
!ws = ws / SDay


#ifdef DEBUG
write(6,*) "ff(1)",ff(:,1)
write(6,*) "End cgem_init"
#endif

return
end subroutine cgem_init


end module cgem
