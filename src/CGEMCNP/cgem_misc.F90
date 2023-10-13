module cgem_misc

use cgem
use grid

implicit none

contains

      Subroutine Command_Line_Args(ivar,ivark)

      integer, intent(out) :: ivar,ivark
      character(10) :: print_var
      character(1) :: print_k
      integer c_count

! ./CGEM Which_code InputFile InitializationFile OutputFileBase 

       !default
       ivar = 0

! --- Command Line Arguments for file names ---
         c_count = command_argument_count()
       !set default no printing
       if (c_count.gt.0) then
        call get_command_argument(1,print_var)  !User selects what to print
        ivark = 1 
        if (c_count.gt.1) then
          call get_command_argument(2,print_k)  !User selects what to print 
          READ (print_k, *) ivark
          !write(6,*) "km,print_k,ivark",km,print_k,ivark
          if(ivark.gt.km) ivark = km 
        endif 

       select case (print_var)
       case ("A")
         ivar = iA(1)
       case ("A1")
         ivar = iA(1)
       case ("A2")
         if(nospA.ge.2) ivar = iA(2)
       case ("A3")
         if(nospA.ge.3) ivar = iA(3)
       case ("A4")
         if(nospA.ge.4) ivar = iA(4)
       case ("A5")
         if(nospA.ge.5) ivar = iA(5)
       case ("A6")
         if(nospA.ge.6) ivar = iA(6)
       case ("Qn")
         ivar = iQn(1)
       case ("Qn1")
         ivar = iQn(1)
       case ("Qn2")
         if(nospA.ge.2) ivar = iQn(2)
       case ("Qn3")
         if(nospA.ge.3) ivar = iQn(3)
       case ("Qn4")
         if(nospA.ge.4) ivar = iQn(4)
       case ("Qn5")
         if(nospA.ge.5) ivar = iQn(5)
       case ("Qn6")
         if(nospA.ge.6) ivar = iQn(6)
       case ("Qp")
         ivar = iQp(1)
       case ("Qp1")
         ivar = iQp(1)
       case ("Qp2")
         if(nospA.ge.2) ivar = iQp(2)
       case ("Qp3")
         if(nospA.ge.3) ivar = iQp(3)
       case ("Qp4")
         if(nospA.ge.4) ivar = iQp(4)
       case ("Qp5")
         if(nospA.ge.5) ivar = iQp(5)
       case ("Qp6")
         if(nospA.ge.6) ivar = iQp(6)
       case ("Z")
         ivar = iZ(1)
       case ("Z1")
         ivar = iZ(1)
       case ("Z2")
         if(nospZ.ge.2) ivar = iZ(2)
       case ("NO3")
         ivar = iNO3
       case ("NH4")
         ivar = iNH4
       case ("PO4")
         ivar = iPO4
       case ("DIC")
         ivar = iDIC
       case ("O2")
         ivar = iO2
       case ("OM1_A")
         ivar = iOM1_A
       case ("OM2_A")
         ivar = iOM2_A
       case ("OM1_Z")
         ivar = iOM1_Z
       case ("OM2_Z")
         ivar = iOM2_Z
       case ("OM1_R")
         ivar = iOM1_R
       case ("OM2_R")
         ivar = iOM2_R
       case ("OM1_BC")
         ivar = iOM1_BC
       case ("OM2_BC")
         ivar = iOM2_BC
       case ("CDOM")
         ivar = iCDOM
       case ("Si")
         ivar = iSi
       case ("Alk")
         ivar = iAlk
       case ("Tr")
         ivar = iTr

       case default
        write(6,*) "These are the only options, case sensitive:"
        write(6,*) "A,A1,Qn,Qp,Z,Z1,Z2,NO3,NH4,PO4,DIC,O2"
        write(6,*) "OM1_A,OM2_A,OM1_Z,OM2_Z,OM1_R,OM2_R"
        write(6,*) "CDOM,Si,Alk,Tr"
        write(6,*) "or you can choose nothing."
        write(6,*) "You chose ", print_var
        STOP
      end select

      else
        ivar=0
      endif

!-A; Phytoplankton number density (cells/m3);
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!-Z: Zooplankton number density (individuals/m3);
!-NO3; Nitrate (mmol-N/m3)
!-NH4; Ammonium (mmol-N/m3)
!-PO4: Phosphate (mmol-P/m3)
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!-O2: Molecular Oxygen (mmol-O2/m3)
!-OM1_A: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton
!-OM2_A: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton 
!-OM1_Z:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets.
!-OM2_Z:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets.
!-OM1_R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!-OM2_R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-OM1_BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!-OM2_BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!- Tr: Tracer
 

       !if (c_count.gt.1) then
       !!  call get_command_argument(2,input_filename)  !User selects input file name
       !endif

       !if (c_count.gt.2) then
       !!  call get_command_argument(3,init_filename) !User selects initial conditions file name
       !endif

       return

       END Subroutine Command_Line_Args


Subroutine Check_InputFile()

integer isp
real x

!read(999,*) Which_fluxes

!read(999,*) Which_temperature
if(Which_temperature.ne.1.and.Which_temperature.ne.2.and.Which_temperature.ne.3) then
  write(6,*) "Which_temperature is outside of range 1-3"
  stop 
endif

!read(999,*) Which_uptake
if(Which_uptake.ne.1.and.Which_uptake.ne.2.and.Which_uptake.ne.3) then
  write(6,*) "Which_uptake is outside of range 1-3"
  stop
endif

!read(999,*) Which_quota
if(Which_quota.ne.1.and.Which_quota.ne.2.and.Which_quota.ne.3) then
  write(6,*) "Which_quota is outside of range 1-3"
  stop
endif


!read(999,*) Which_irradiance 
if(Which_irradiance.ne.1.and.Which_irradiance.ne.2) then
  write(6,*) "Which_irradiance is outside of range 1-2"
  stop
endif

!read(999,*) Which_photosynthesis 
if(Which_photosynthesis.ne.1.and.Which_photosynthesis.ne.2.and.Which_photosynthesis.ne.3) then
  write(6,*) "Which_photosynthesis is outside of range 1-3"
  stop
endif

if(Which_photosynthesis.eq.3.and.Which_growth.ne.3) then
  write(6,*) "If Which_photosynthesis=3 then Which_growth should = 3"
  stop
endif

!read(999,*) Which_growth
if(Which_growth.ne.1.and.Which_growth.ne.2.and.Which_growth.ne.3) then
  write(6,*) "Which_growth is outside of range 1-3"
  stop
endif

if(Which_growth.eq.3.and.Which_photosynthesis.ne.3) then
  write(6,*) "If Which_growth=3 then Which_photosynthesis should =3"
  stop
endif

!read(999,*) SolarRad
!if(SolarRad.ne.0) then
!  write(6,*) "SolarRad will be calculated, currently the only option"
!endif


!Do a test to make sure Qmax is greater than Qmin:
do isp=1,nospA
   if(QmaxN(isp).le.QminN(isp)) then
     write(6,*) "Please set QmaxN greater than QminN"
     stop
   endif
   if(QmaxP(isp).le.QminP(isp)) then
     write(6,*) "Please set QmaxP greater than QminP"
     stop
   endif
enddo

! An hour (3600 secs) must be divisible by dT:
if(mod(3600,dT).ne.0) then
  write(6,*) "Please pick a dT that is a divisor of 3600 seconds (an hour)."
endif

! Diatom/Non Diatom check:
do isp=1,nospA
   if(KSi(isp).le.tiny(x).and.vmaxSi(isp).gt.0.) then
     write(6,*) "If KSi=0, then vmaxSi must =0 (designates non-diatoms)"
     stop
   endif
   if(vmaxSi(isp).le.tiny(x).and.KSi(isp).gt.0.) then
     write(6,*) "If vmaxSi=0, then KSi must =0 (designates non-diatoms)"
     stop
   endif
enddo


#ifdef DEBUG
write(6,*) "=============Check_InputFile============================"
write(6,*) "!Grid parameters"
write(6,*) "=========== grid ============================"
write(6,*) "integer :: km",km
write(6,*) "integer :: nospA",nospA
write(6,*) "integer :: nospZ",nospZ
write(6,*) ""
write(6,*) "!--Run Specifics---------------"
write(6,*) "integer iYrS,iMonS,iDayS,iHrS,iMinS,iSecS",iYrS,iMonS,iDayS,iHrS,iMinS,iSecS
write(6,*) "integer iYrE,iMonE,iDayE,iHrE,iMinE,iSecE",iYrE,iMonE,iDayE,iHrE,iMinE,iSecE
write(6,*) "integer dT,nstep",dT,nstep
write(6,*) ""
write(6,*) "!Constants"
write(6,*) "integer :: StepsPerDay",StepsPerDay
write(6,*) "real :: dTd",dTd 
write(6,*) "integer :: dT_out !Output timestep",dT_out
write(6,*) "integer iout",iout
write(6,*) "!-output"
write(6,*) "integer istep_out",istep_out
write(6,*) ""
write(6,*) "!This should all be from SCHISM"
write(6,*) "!Solar radiation"
write(6,*) "real :: Rad",Rad
write(6,*) "!Lat/lon of elements"
write(6,*) "real :: lat,lon",lat,lon
write(6,*) "!Depth stuff"
write(6,*) "real, allocatable :: d(:) ! depth cell bottom to surface",d
write(6,*) "real, allocatable :: d_sfc(:) ! depth cell center to surface",d_sfc
write(6,*) "real, allocatable :: dz(:)   !Thickness of cell",dz
write(6,*) "!'tracers'"
write(6,*) "real, allocatable :: S(:)",S
write(6,*) "real, allocatable :: T(:)",T
write(6,*) ""
write(6,*) ""
write(6,*) "============= cgem ============================"
!misc
write(6,*) "real :: eps",eps
write(6,*) ""
write(6,*) "!Sinking"
write(6,*) "real, dimension(:), allocatable :: ws",ws
write(6,*) ""
write(6,*) ""
write(6,*) "!Stoichiometry"
write(6,*) "real, dimension(:), allocatable :: s_x1A,s_x2A,s_y1A,s_y2A"
write(6,*) s_x1A,s_x2A,s_y1A,s_y2A
write(6,*) "real, dimension(:), allocatable :: s_z1A,s_z2A"
write(6,*) s_z1A,s_z2A
write(6,*) "real, dimension(:), allocatable :: s_x1Z,s_x2Z,s_y1Z,s_y2Z"
write(6,*) s_x1Z,s_x2Z,s_y1Z,s_y2Z
write(6,*) "real, dimension(:), allocatable :: s_z1Z,s_z2Z"
write(6,*) s_z1Z,s_z2Z
write(6,*) "       REAL :: Esed",Esed
write(6,*) "       REAL :: CBODW",CBODW
write(6,*) "       REAL,ALLOCATABLE :: pH(:)",pH
write(6,*) "!---------------------------------------------------------      "
write(6,*) "!-A; Phytoplankton number density (cells/m3);"
write(6,*) "!---------------------------------------------------------  "
write(6,*) "      integer, dimension(:), allocatable :: iA(:)",iA
write(6,*) "!----------------------------------------------------------------------"
write(6,*) "!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)"
write(6,*) "!----------------------------------------------------------------------"
write(6,*) "      integer, dimension(:), allocatable :: iQn(:)",iQn
write(6,*) "!----------------------------------------------------------------------"
write(6,*) "!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)"
write(6,*) "!----------------------------------------------------------------------"
write(6,*) "      integer, dimension(:), allocatable :: iQp(:)",iQp
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "!-Z: Zooplankton number density (individuals/m3);"
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "      integer, dimension(:), allocatable :: iZ(:)",iZ
write(6,*) "!-------------------------------"
write(6,*) "!-NO3; Nitrate (mmol-N/m3)"
write(6,*) "!-------------------------------"
write(6,*) "      integer :: iNO3",iNO3
write(6,*) "!--------------------------------      "
write(6,*) "!-NH4; Ammonium (mmol-N/m3)"
write(6,*) "!--------------------------------"
write(6,*) "      integer :: iNH4",iNH4
write(6,*) "!-------------------------------------------        "
write(6,*) "!-PO4: Phosphate (mmol-P/m3)"
write(6,*) "!--------------------------------------"
write(6,*) "      integer :: iPO4 ",iPO4
write(6,*) "!---------------------------------------------------------"
write(6,*) "!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) "
write(6,*) "!---------------------------------------------------------"
write(6,*) "      integer :: iDIC ",iDIC
write(6,*) "!-------------------------------------------        "
write(6,*) "!-O2: Molecular Oxygen (mmol-O2/m3)"
write(6,*) "!------------------------------"
write(6,*) "      integer :: iO2 ",iO2
write(6,*) "!-------------------------------------------------------------"
write(6,*) "!-OM1_A: (mmol-C/m3--particulate)"
write(6,*) "!        -- Particulate Organic Matter arising from "
write(6,*) "!           dead Phytoplankton"
write(6,*) "!-------------------------------------------------------------"
write(6,*) "      integer :: iOM1_A",iOM1_A
write(6,*) "!-----------------------------------------------------------------"
write(6,*) "!-OM2_A: (mmol-C/m3--dissolved)"
write(6,*) "!        -- Dissolved Organic Matter arising from "
write(6,*) "!           dead Phytoplankton "
write(6,*) "!------------------------------------------------------------------"
write(6,*) "      integer :: iOM2_A",iOM2_A
write(6,*) "!-------------------------------------------------------------"
write(6,*) "!-OM1_Z:(mmol-C/m3--particulate)"
write(6,*) "!        -- Particulate Organic Matter arising from "
write(6,*) "!           Zooplankton fecal pellets."
write(6,*) "!-------------------------------------------------------------"
write(6,*) "      integer :: iOM1_Z",iOM1_Z
write(6,*) "!-------------------------------------------------        "
write(6,*) "!-OM2_Z:(mmol-C/m3--dissolved)"
write(6,*) "!        -- Dissolved Organic Matter arising from "
write(6,*) "!          Zooplankton fecal pellets."
write(6,*) "!-----------------------------------------------"
write(6,*) "      integer :: iOM2_Z",iOM2_Z
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "!-OM1_R: (mmol-C/m3--particulate)"
write(6,*) "!         -- Particulate Organic Matter arising from river outflow"
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "      integer :: iOM1_R",iOM1_R
write(6,*) "!---------------------------------------------      "
write(6,*) "!-OM2_R: (mmol-C/m3--dissolved)"
write(6,*) "!         -- Dissolved Organic Matter arising from river outflow"
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "      integer :: iOM2_R",iOM2_R
write(6,*) "!-------------------------------------------"
write(6,*) "!-CDOM: (ppb) "
write(6,*) "!        -- Colored Dissolved Organic Matter"
write(6,*) "!-------------------------------------------"
write(6,*) "      integer :: iCDOM",iCDOM
write(6,*) "!---------------------------------------------"
write(6,*) "!-Silica: (mmol-Si/m3) "
write(6,*) "!        -- Silica"
write(6,*) "!-------------------------------------------"
write(6,*) "      integer :: iSi",iSi
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "!-OM1_BC: (mmol-C/m3--particulate)"
write(6,*) "!         -- Particulate Organic Matter in initial and boundary "
write(6,*) "!            conditions "
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "      integer :: iOM1_BC",iOM1_BC
write(6,*) "!-------------------------------------------------"
write(6,*) "!-OM2_BC: (mmol-C/m3--dissolved)"
write(6,*) "!         -- Dissolved Organic Matter in initial and boundary"
write(6,*) "!            conditions"
write(6,*) "!--------------------------------------------------------------------"
write(6,*) "      integer :: iOM2_BC",iOM2_BC
write(6,*) "!-------------------------------------------"
write(6,*) "!-ALK:  (mmol-HCO3/m3)?"
write(6,*) "!        -- Alkalinity"
write(6,*) "!-------------------------------------------"
write(6,*) "      integer :: iALK",iALK
write(6,*) "!Tracer"
write(6,*) "      integer :: iTR",iTR
write(6,*) "!Total number of state variables"
write(6,*) "      integer :: nf",nf
write(6,*) ""
write(6,*) "!State Variable Array"
write(6,*) "      real,allocatable :: ff(:,:) !state variable array"
write(6,*) ""
write(6,*) "!----INPUT_VARS_CGEM"
write(6,*) "!--Switches in GEM---------"
write(6,*) "integer Which_fluxes(8)",Which_fluxes
write(6,*) "integer Which_uptake",Which_uptake
write(6,*) "integer Which_quota",Which_quota
write(6,*) "integer Which_irradiance",Which_irradiance
write(6,*) "integer Which_chlaC",Which_chlaC
write(6,*) "integer Which_photosynthesis",Which_photosynthesis
write(6,*) "integer Which_growth",Which_growth
write(6,*) "integer Which_temperature",Which_temperature
write(6,*) "!--Temperature"
write(6,*) "real, allocatable :: KTg1(:)",KTg1
write(6,*) "real, allocatable :: KTg2(:)",KTg2
write(6,*) "real, allocatable :: Tref(:)",Tref
write(6,*) "real, allocatable :: Ea(:)",Ea
write(6,*) "!--Optics-----------------------"
write(6,*) "real Kw",Kw
write(6,*) "real Kcdom",Kcdom
write(6,*) "real Kspm",Kspm
write(6,*) "real Kchla",Kchla
write(6,*) "!--in module LIGHT_VARS"
write(6,*) "real astar490",astar490
write(6,*) "real aw490",aw490
write(6,*) "real astarOMA",astarOMA
write(6,*) "real astarOMZ",astarOMZ
write(6,*) "real astarOMR",astarOMR
write(6,*) "real astarOMBC",astarOMBC
write(6,*) "real PARfac",PARfac
write(6,*) "!---Phytoplankton "
write(6,*) "real, allocatable :: ediblevector(:,:)",ediblevector
write(6,*) "real, allocatable :: umax(:)",umax
write(6,*) "real, allocatable :: CChla(:)",CChla
write(6,*) "real, allocatable :: alpha(:)",alpha
write(6,*) "real, allocatable :: beta(:)",beta
write(6,*) "real, allocatable :: respg(:)",respg
write(6,*) "real, allocatable :: respb(:)",respb
write(6,*) "real, allocatable :: QminN(:)",QminN
write(6,*) "real, allocatable :: QminP(:)",QminP
write(6,*) "real, allocatable :: QmaxN(:)",QmaxN
write(6,*) "real, allocatable :: QmaxP(:)",QmaxP
write(6,*) "real, allocatable :: Kn(:)",Kn
write(6,*) "real, allocatable :: Kp(:)",Kp
write(6,*) "real, allocatable :: Ksi(:)",Ksi
write(6,*) "real, allocatable :: KQn(:)",KQn
write(6,*) "real, allocatable :: KQp(:)",KQp
write(6,*) "real, allocatable :: nfQs(:)",nfQs
write(6,*) "real, allocatable :: vmaxN(:)",vmaxN
write(6,*) "real, allocatable :: vmaxP(:)",vmaxP
write(6,*) "real, allocatable :: vmaxSi(:)",vmaxSi
write(6,*) "real, allocatable :: aN(:)",aN
write(6,*) "real, allocatable :: volcell(:)",volcell
write(6,*) "real, allocatable :: Qc(:)",Qc
write(6,*) "real, allocatable :: Athresh(:)",Athresh
write(6,*) "real, allocatable :: mA(:)",mA
write(6,*) "real, allocatable :: A_wt(:)",A_wt
write(6,*) "!---Zooplankton"
write(6,*) "real, allocatable :: Zeffic(:)",Zeffic
write(6,*) "real, allocatable :: Zslop(:)",Zslop
write(6,*) "real, allocatable :: Zvolcell(:)",Zvolcell
write(6,*) "real, allocatable :: ZQc(:)",ZQc
write(6,*) "real, allocatable :: ZQn(:)",ZQn
write(6,*) "real, allocatable :: ZQp(:)",ZQp
write(6,*) "real, allocatable :: ZKa(:)",ZKa
write(6,*) "real, allocatable :: Zrespg(:)",Zrespg
write(6,*) "real, allocatable :: Zrespb(:)",Zrespb
write(6,*) "real, allocatable :: Zumax(:)",Zumax
write(6,*) "real, allocatable :: Zm(:)",Zm
write(6,*) "!---Organic Matter        "
write(6,*) "real KG1",KG1
write(6,*) "real KG2",KG2
write(6,*) "real KG1_R",KG1_R
write(6,*) "real KG2_R",KG2_R
write(6,*) "real KG1_BC",KG1_BC
write(6,*) "real KG2_BC",KG2_BC
write(6,*) "real KNH4",KNH4
write(6,*) "real nitmax",nitmax
write(6,*) "real KO2",KO2
write(6,*) "real KstarO2",KstarO2
write(6,*) "real KNO3",KNO3
write(6,*) "real pCO2",pCO2
write(6,*) "real stoich_x1R",stoich_x1R
write(6,*) "real stoich_y1R",stoich_y1R
write(6,*) "real stoich_z1R",stoich_z1R
write(6,*) "real stoich_x2R",stoich_x2R
write(6,*) "real stoich_y2R",stoich_y2R
write(6,*) "real stoich_z2R",stoich_z2R
write(6,*) "real stoich_x1BC",stoich_x1BC
write(6,*) "real stoich_y1BC",stoich_y1BC
write(6,*) "real stoich_z1BC",stoich_z1BC
write(6,*) "real stoich_x2BC",stoich_x2BC
write(6,*) "real stoich_y2BC",stoich_y2BC
write(6,*) "real stoich_z2BC",stoich_z2BC
write(6,*) "real KGcdom",KGcdom
write(6,*) "real CF_SPM",CF_SPM
write(6,*) "!----Other including Boundary Conditions------------"
write(6,*) "real KH_coeff",KH_coeff
write(6,*) "integer Which_Outer_BC",Which_Outer_BC
write(6,*) "real wt_l, wt_o",wt_l,wt_o
write(6,*) "real wt_pl, wt_po",wt_pl,wt_po
write(6,*) "real m_OM_init,m_OM_BC,m_OM_sh",m_OM_init,m_OM_BC,m_OM_sh
write(6,*) "real KG_bot",KG_bot
write(6,*) ""
write(6,*) "!Light curve parameters"
write(6,*) "real, allocatable :: alphad(:) ! Initial slope of photosynthesis-irradiance"
write(6,*) "curve / Vmax      ",alphad
write(6,*) "real, allocatable :: betad(:)  ! Photoinhibition constant / Vmax",betad
write(6,*) ""
write(6,*) "========================================="
write(6,*) ""

!stop
#endif

return
END SUBROUTINE Check_InputFile


end module cgem_misc
