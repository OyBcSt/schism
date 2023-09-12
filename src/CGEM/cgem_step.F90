!======================================================================     
    Subroutine cgem_step( TC_8, dT, istep, inea, myrank )

!======================================================================
    use grid !, only: lat,T,S,Rad
    use cgem
    use cgem_light
    use cgem_growth
    use cgem_utils
    use date_time
    use mvars

    IMPLICIT NONE

!---------------------------------------------
! Interface variables
!---------------------------------------------------------------------
    integer, intent(in)  :: TC_8         ! Model time (seconds from beginning of Jan 1, 2002)
    integer, intent(in)  :: istep     ! Current time step
    integer, intent(in)  :: inea  !grid element for this step
    integer, intent(in) :: myrank !process number
    real, intent(in) :: dT
!---------------------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------
    integer        ::  k, isp, isz ! Loop indicies, isp/isz is for phytoplankton/zooplankton species
    integer        ::  Is_Day            ! Switch for day/night for phytoplankton nutrient uptake only, Is_Day=0 means night
!------------------------------------ 
! Phytoplankton parameters
! Phytoplankton uptake and growth
    real, dimension(km,nospA) :: A     ! Zooplankton number density (indv./m3)
    real, dimension(km,nospA) :: Agrow ! Phytoplankton growth (cells/m3/d)
    real, dimension(km,nospA) :: Aresp ! Total respiration from a phytoplankton group (cells/m3/d)
    real, dimension(km,nospA) :: uA    ! Specific growth rate (1/d)
    real, dimension(km,nospA) :: uN    ! Nitrogen Limited growth rate (1/d)
    real, dimension(km,nospA) :: uP    ! Phosphorus limited growth rate (1/d)
    real, dimension(km,nospA) :: uE    ! Light limited growth rate (1/d)
    real, dimension(km,nospA) :: uSi   ! Silica limited growth rate (1/d)
    real, dimension(km,nospA) :: Qn    ! Phytoplankton Nitrogen Quota (mmol-N/cell)
    real, dimension(km,nospA) :: Qp    ! Phytoplankton Phosphorus Quota (mmol-P/cell)
    real, dimension(km,nospA)    :: f_Qn  ! Quota model for N
    real, dimension(km,nospA)    :: f_Qp  ! Quota model for P
    real, dimension(km)    :: vN    ! Phytoplankton uptake rate of Nitrogen (mmol-N/cell/d)
    real, dimension(km)    :: vP    ! Phytoplankton uptake rate of Phosphorus (mmol-P/cell/d)
    real, dimension(km)    :: vSi   ! Phytoplankton uptake rate of Silica (mmol-Si/cell/d)
    real, dimension(km)    :: AupN  ! Total Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
    real, dimension(km)    :: AupP  ! Total Phytoplankton uptake of Phosphorus (mmol-P/m3/d)
    real, dimension(km)    :: AupSi ! Total Phytoplankton uptake of Silica (mmol-Si/m3/d)
    integer, dimension(km) :: RLN   ! Rate Limiting Nutrient of N, P, and Si
 ! Monod equations for phytoplankton
    real, dimension(nospA)    :: monodN  !Monod term in nitrogen uptake
    real, dimension(nospA)    :: monodP  !Monod term in phosphorus uptake
    real, dimension(nospA)    :: monodSi !Monod term in Si uptake
    real, dimension(km)       :: Ntotal   ! Total N (mmol/m3)
 ! Phytoplankton nutrient loss
    real, dimension(nospA)    :: Amort ! Dead phytoplankton (cells/m3/day)
    real, dimension(km)    :: AexudN          ! Sum of Exudation of N from all phytoplankton groups (mmol-N/m3/d)
    real, dimension(km)    :: AexudP          ! Sum of Exudation of P from all phytoplankton groups (mmol-P/m3/d)
    real, dimension(km)    :: ArespC          ! Phytoplankton equivalent carbon loss from respiration (mmol-C/m3/d)
!------------------------------------------------------------------
! Zooplankton parameters
 !Zooplankton uptake and growth
    real, dimension(km,nospZ)   :: Z         ! Zooplankton number density (indv./m3)
    real, dimension(nospZ)       :: optNP    ! Optimal nutrient ratio for zooplankton
    real, dimension(km,nospZ)       :: Zgrow    ! Zooplankton growth (indv./m3/d)
    real, dimension(nospA,nospZ) :: Zgrazvol ! Grazing rate in units of biovolume (um3/m3/d)
    real, dimension(nospA,nospZ) :: ZgrazA   ! Zooplankton grazing of phytoplankton (cells/m3/d)
    real, dimension(nospA)       :: ZgrazA_tot ! Total zooplankton grazing of phytoplankton (cells/m3/d)
    real, dimension(nospZ)       :: ZgrazN   ! Zooplankton grazing uptake of Nitrogen (mmol-N/m3/d)
    real, dimension(nospZ)       :: ZgrazP   ! Zooplankton grazing uptake of Phosphorus (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZgrazC   ! Zooplankton grazing uptake of Carbon (mmol-C/m3/d)
    real, dimension(nospZ)       :: ZinN     ! Zooplankton ingestion of Nitrogen (mmol-N/m3/d)
    real, dimension(nospZ)       :: ZinP     ! Zooplankton ingestion of Phosphorus (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZinC     ! Zooplankton ingestion of Carbon (mmol-C/m3/d)
 !Monod equations for zooplankton ingestion of phytoplankton
    real, dimension(km,nospA,nospZ) :: monodZ   ! Monod term for zooplankton grazing
    real, dimension(km,nospA)     :: Abiovol  ! Algae biovolume vector (um3/m3)
    real, dimension(km,nospA)     :: Abiovol_thresh  ! Algae biovolume vector (um3/m3)
    real, dimension(nospA,nospZ) :: top_A    ! Monod numerator value for phytoplankton group
    real, dimension(nospA,nospZ) :: bottom_A ! Monod Denominator value for phytoplankton group
    real, dimension(nospZ)       :: bottom   ! Sum of Monod Denominator value for all phytoplankton groups
 !Zooplankton nutrient loss
    real, dimension(km,nospZ)       :: Zresp    ! Zooplankton respiration (individuals/m3/d)
    real, dimension(km)          :: ZrespC   ! Carbon loss from zooplankton respiration (mmol-C/m3/day)
    real, dimension(nospZ)       :: ZunC     ! Unassimilated ingested Carbon (mmol-C/m3/d)
    real, dimension(nospZ)       :: ZunN     ! Unassimilated ingested Nitrogen (mmol-N/m3/d)
    real, dimension(nospZ)       :: ZunP     ! Unassimilated ingested Phosphorus (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZunSi    ! Unassimilated ingested Silica (mmol-Si/m3/d)
    real, dimension(km)          :: sumZunSi
    real, dimension(km,nospZ)     :: Zmort    ! Dead zooplankton (individuals/m3/d)
    real :: ZmortC(nospZ), ZmortC_tot        ! Carbon released from dead zooplankton (mmol-C/m3/d)
    real :: ZmortN(nospZ), ZmortN_tot        ! Nitrogen released from dead zooplankton (mmol-N/m3/d)
    real :: ZmortP(nospZ), ZmortP_tot        ! Phosphorus released from dead zooplankton (mmol-P/m3/d)
    real :: ZslopC(nospZ), ZslopC_tot        ! Carbon lost to sloppy feeding (mmol-C/m3/d)
    real :: ZslopN(nospZ), ZslopN_tot        ! Nitrogen lost to sloppy feeding (mmol-N/m3/d)
    real :: ZslopP(nospZ), ZslopP_tot        ! Phosphorus lost to sloppy feeding (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZexN     ! Excretion from zooplankton (mmol-N/m3/d)
    real, dimension(km)          :: sumZexN
    real, dimension(nospZ)       :: ZexP     ! Excretion from zooplankton (mmol-P/m3/d)
    real, dimension(km)          :: sumZexP
    real, dimension(nospZ)       :: ZegC     ! Egestion from zooplankton (mmol-C/m3/d)
    real, dimension(nospZ)       :: ZegN     ! Egestion from zooplankton (mmol-N/m3/d)
    real, dimension(nospZ)       :: ZegP     ! Egestion from zooplankton (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZegSi    ! Egestion from zooplankton (mmol-Si/m3/d)
    real, dimension(km)          :: sumZegSi
    real, dimension(km) :: OM1_Ratio, OM2_Ratio             ! Separates sloppy feeding into OM1 and OM2
!---------------------------------------------------------------------- 
! Time variables  
    real, parameter :: one_d_365  = 1.0/365.0 ! Convert 1/yr to 1/day
!-----------------------------------------------------------------------
! Organic Matter Calculations
   ! Variables to calculate stoichiometry C:N:P ratios
    real,dimension(km)    :: OM1_CA, OM1_NA, OM1_PA    ! OM from dead phytoplankton
    real,dimension(km)    :: OM2_CA, OM2_NA, OM2_PA   
    real,dimension(km)    :: OM1_CZ, OM1_NZ, OM1_PZ    ! OM from zooplankton
    real,dimension(km)    :: OM2_CZ, OM2_NZ, OM2_PZ
    real,parameter    :: sz1 = 1. !for stoichiometry x,y,z, z is always normalized to 1 
    real, dimension(km) :: sx1A,sy1A,sx2A,sy2A,sx1Z,sy1Z,sx2Z,sy2Z 
    real, dimension(km) :: OM2A, OM2Z, OM2R, OM2BC !POC in g/m3
!---------------------------------------------------------------------------
! reaction and Nitrification subroutine variables
    real    :: R_11                                     ! Nitrification term
    real    :: RNO3_A, RNO3_Z, RNO3_R, RNO3_BC ! Remineralization terms for NO3
    real    :: RNH4_A, RNH4_Z, RNH4_R, RNH4_BC ! Remineralization terms for NH4
    real, dimension(km) :: ROM1_A, ROM1_Z, ROM1_R, ROM1_BC       ! Remineralization terms for POC
    real, dimension(km) :: ROM2_A, ROM2_Z, ROM2_R, ROM2_BC       ! Remineralization terms for DOC
    real    :: RO2_A, RO2_Z, RO2_R, RO2_BC      ! Remineralization terms for O2
    real    :: RPO4_A, RPO4_Z, RPO4_R, RPO4_BC ! Remineralization terms for PO4
    real    :: RDIC_A, RDIC_Z, RDIC_R, RDIC_BC ! Remineralization terms for DIC
    real    :: RSi_A, RSi_Z, RSi_R, RSi_BC      ! Remineralization terms for Si
    real    :: RALK_A, RALK_Z, RALK_R, RALK_BC ! Remineralization terms for ALK
    real    :: RN2_A, RN2_Z, RN2_R, RN2_BC ! Remineralization terms for N2 
    real, dimension(10) :: RC   ! Array that returns remineralization terms for OM
    real, dimension(km) :: RNO3,RNH4,RO2,RPO4,RDIC,RSi,RALK,RN2
    real, dimension(km) :: NO3,NH4,DIC,O2,PO4,Si,ALK
    !---------------------------------------------------------
! Variables needed for light routine and calc_Agrow
    real                :: PARbot             ! Irradiance at sea floor (quanta/cm2/s)
    real, dimension(km) :: PARdepth           ! Irradiance at center of layer k (quanta/cm2/s)
    real                :: PARsurf            ! Irradiance just below sea surface
    real, parameter :: RADCONV = 1./6.0221413*1.e-19 ! Convert quanta/cm2/s to mol/m2/s:
                                               !  = quanta/cm2/s * 1 mol/Avogadro# * 10,000cm2/m2
                                               !  = (1/6.022e23) * 1.0e4 = (1./6.022)e-23 * 1.0e4
                                               !  = (1./6.0221413)*1.e-19
    real, dimension(km) :: Chla_tot  ! Total amount of Chl-a in all the
                                        !  phytoplankton species (mg/m3) per cell
    real, dimension(km) :: OM1A, OM1Z, OM1R, OM1BC !POC in g/m3
    real, dimension(km) :: CDOM    ! CDOM, ppb
!-----------------------------------------------------------------------
! Other variables 
    real, dimension(km) :: PrimProd                     ! Primary production (photosynthesis)
    real, dimension(nospA+nospZ) :: Tadj ! Temperature adjustment factor
!------------------------------------------------------------------    
! SAVE KGs for instant remineralization
    real, save :: KG1_save, KG2_save
!timestep in days
    real :: dTd
    integer :: StepsPerDay
!------------------------------------------------------------------
!Output vars for alkalinity subroutine:
    real :: ph_calc(1), pco2_calc(1), fco2(1), co2(1), hco3(1), co3(1), omegaa(1), omegac(1), betad_calc(1) 
    real :: rhosw(1), p(1), tempis(1)
    real :: patm(1) = 1.
    real :: m_alk(1), m_dic(1), m_si(1), m_po4(1)
    real :: m_lat(1)
!For tiny
    real x
!mocsy needs lat to be an array
    m_lat = lat

!convert to timestep in days
  StepsPerDay = 86400.
  dTd = dt/StepsPerDay

! write(6,*) "SPD,dTd,dt",StepsPerDay,dTd,dt

if(Which_rad.eq.0) call getSolar( iYrS, TC_8, lon, lat, Rad)
if(Which_wind.eq.0) Wind = 5. 
!write(6,*) "TC_8,Rad,Wind",TC_8,Rad,Wind

#ifdef DEBUG
write(6,*) "Begin cgem, TC_8,istep",TC_8,istep
#endif

   optNP = ZQn/ZQp    ! Optimal nutrient ratio for zooplankton

   if(istep.eq.1.and.inea.eq.10.and.writecsv.eq.1) then
     open(unit=6001,file="growth.csv")
     open(unit=6101,file="rates.csv")
     open(unit=6301,file="hydro.csv")
     open(unit=6401,file="light.csv")
    write(6001,'(A58)') "Agrow,Aresp,uA,uN,uP,uE,uSi,f_E,f_N,f_P,f_Si,Tadj,A,min_S"
    write(6101,'(A314)') "RO2,RNO3,RNH4,RPO4,RDIC,RSi,RALK,RN2,ROM1_A,ROM2_A,RO2_A,RNO3_A,RPO4_A,RDIC_A,RNH4_A,RSi_A,RALK_A,RN2_A,ROM1_Z,ROM2_Z,RO2_Z,RNO3_Z,RPO4_Z,RDIC_Z,RNH4_Z,RSi_Z,RALK_Z,RN2_Z,ROM1_R,ROM2_R,RO2_R,RNO3_R,RPO4_R,RDIC_R,RNH4_R,RSi_R,RALK_R,RN2_R,ROM1_BC,ROM2_BC,RO2_BC,RNO3_BC,RPO4_BC,RDIC_BC,RNH4_BC,RSi_BC,RALK_BC,RN2_BC"
    write(6301,'(A22)') "TC_8,rad,wind,sal,temp"
    write(6401,'(A45)') "ksurf,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kbot"
    endif


do k=1,km
  do isp = 1,nf
    if(ff(k,isp).le.0.) then
      write(6,*) "ff.le.0! set to 0 for istep,myrank,inea,k,isp=",istep,myrank,inea,k,isp,ff(k,isp)
      ff(k,isp) = 0.0
    endif
  enddo
enddo

       do isp=1,nospA
         A(:,isp) = ff(:,iA(isp))
       enddo

       ! After Advection and VMixing, return to Q's
       do isp=1,nospA
        Qn(:,isp) = ff(:,iQn(isp)) / A(:,isp) 
        Qp(:,isp) = ff(:,iQp(isp)) / A(:,isp) 
       enddo

       do isz=1,nospZ
        Z(:,isz) = ff(:,iZ(isz))
       enddo

       !! After Advection and VMixing, return to stoich's
        sx1A(:) = ff(:,isx1A) / ff(:,iOM1_A)
        sy1A(:) = ff(:,isy1A) / ff(:,iOM1_A)
        sx2A(:) = ff(:,isx2A) / ff(:,iOM2_A)
        sy2A(:) = ff(:,isy2A) / ff(:,iOM2_A)
        sx1Z(:) = ff(:,isx1Z) / ff(:,iOM1_Z)
        sy1Z(:) = ff(:,isy1Z) / ff(:,iOM1_Z)
        sx2Z(:) = ff(:,isx2Z) / ff(:,iOM2_Z)
        sy2Z(:) = ff(:,isy2Z) / ff(:,iOM2_Z)

       do isz=1,nospZ
        Z(:,isz) = ff(:,iZ(isz))
       enddo

       NO3(:)   = ff(:,iNO3)
       NH4(:)   = ff(:,iNH4)
       PO4(:)   = ff(:,iPO4)
       DIC(:)   = ff(:,iDIC)
       O2(:)    = ff(:,iO2)
       OM1A(:)  = ff(:,iOM1_A)
       OM2A(:)  = ff(:,iOM2_A)
       OM1Z(:)  = ff(:,iOM1_Z)
       OM2Z(:)  = ff(:,iOM2_Z)
       OM1R(:)  = ff(:,iOM1_R)
       OM2R(:)  = ff(:,iOM2_R)
       CDOM(:)  = ff(:,iCDOM)
       Si(:)    = ff(:,iSi)
       OM1BC(:) = ff(:,iOM1_BC)
       OM2BC(:) = ff(:,iOM2_BC)
       ALK(:)   = ff(:,iALK)
       Ntotal(:)    = NO3(:) + NH4(:)


!write(6,*) "istep=",istep
!-----------------------------------------------------------------
!   Begin main ij loop for the biogeochemistry 
!   calculations at time-level istep
!-----------------------------------------------------------------

!----------------------------------------------------------------
! Get chlorophyll-a quantity per layer
  Chla_tot = Fixed_CChla(A)

!----------------------------------------------------------------------
! Execute the desired atmospheric light model.  To calculate PARsurf,
! the effect amount of downward spectrally integrated irradiance 
! just below the sea surface.  'Rad' is just above sea surface. 
!----------------------------------------------------------------------
  ! Rad(i) is short wave generated by NRL is used,
  ! and is multiplied by SWtoPAR: ratio of PAR to
  ! shortwave radiation (hardcoded 4/30/14 to 0.43).
  ! Hardcoded to 0.47 on 2/11/16, Re: Tsubo and Walker, 2005
  ! PARfac is a multiplication factor for testing
  PARsurf = (0.47 * Rad) * PARfac
!--End Calculate atmospheric model ---------------------------------------------

!----------------------------------------------------------------------------
! Underwater light model to calculate the 1-D radiation array PARdepth 
! PARbot is the downward irradiance (photons/cm2/sec) at the sea bottom
!-------------------------------------------------------------------------
  if(PARsurf.le.0.) then
   PARbot = 0.0
   PARdepth = 0.0
  else
  !This is a wrapper that calls converts mmol to g, Chla, etc., internally 
  call calc_PARdepth(TC_8,PARSurf,S,Chla_tot,CDOM,OM1A,OM1Z,OM1R,OM1BC,&
       &             PARdepth,PARbot )
  endif 
!---------------------End Underwater Light Model-----------------------------------

!-------------------------------------------------------------------------
! Phytoplankton growth model to calculate the 1-D arrays Agrow and Aresp
! Agrow(k) is the growth-rate  for vertical grid column (i) at cell (k).
!-----------------------------------------------------------------------
 if(km .gt. 0) then
   call calc_Agrow(PARdepth,T,Qn,Qp,Si,A,Agrow,uA,Aresp,uN,uP,uE,uSi,inea)
!------end phytoplankton growth model-------------------------
 endif

!Zooplankton grazing stuff
   do isp = 1, nospA
     Abiovol(:,isp) = A(:,isp)*volcell(isp)
     Abiovol_thresh(:,isp) = Abiovol(:,isp) - Athresh(isp)
   enddo


   do k=1,km
     do isp = 1,nospA
       do isz = 1,nospZ
         top_A(isp,isz)    = AMAX1(Abiovol_thresh(k,isp),0.0)*ediblevector(isz,isp)
         bottom_A(isp,isz) = Abiovol(k,isp)       *ediblevector(isz,isp)
       enddo
      enddo

     do isz = 1, nospZ
       bottom(isz) = SUM(bottom_A(:,isz))   ! sum over isp
     enddo

     do isp = 1, nospA
       do isz = 1,nospZ
         monodZ(k,isp,isz)  = top_A(isp,isz)/(ZKa(isz) + bottom(isz))
       enddo 
     enddo
   enddo

!--------------------------------------
! Call temperature and growth functions
!-----------------------------------------
      call func_T( T,Tadj )
!     Nutrient dependent growth function
      call func_Qs( Qn, Qp, f_Qn, f_Qp) 

      do k=1,km
!Nutrients only taken up during the day:
     Is_Day = 1
     if(Rad.le.tiny(x)) Is_Day = 0

#ifdef DEBUG
write(6,*) "In cgem, called T,Qx functions" 
#endif

!--------------------------------------------------------------------------
! Initialize counters to zero that are used to accumulate variable values
! over the nospA phytoplankton groups and the nospZ zooplankton groups.
!--------------------------------------------------------------------------
      PrimProd(k)  = 0.0
      ArespC(k)    = 0.0
      AexudN(k)    = 0.0
      AexudP(k)    = 0.0
      AupN(k)      = 0.0
      AupP(k)      = 0.0
      AupSi(k)     = 0.0
      ZgrazC(k)    = 0.0
      ZgrazN(k)    = 0.0
      ZgrazP(k)    = 0.0
! ----------------------------------------------------------------------
      do isp = 1, nospA      
! ----------------------------------------------------------------------   

!---------------------------------------------------------------------      
! Note that the expressions for PrimProd, ArespC, AexudN, AexudP are 
! sums over the isp phytoplankton groups. When the isp loop is complete,
! PrimProd, ArespC. AexudN, and AexudP will represent totals for
! all the nospA phytoplankton groups. 
!--------------------------------------------------------------------
     PrimProd = PrimProd + Agrow(k,isp)*Qc(isp)     ! Phytoplankton primary production (mmol-C/m3/d)
     ArespC   = ArespC + Aresp(k,isp)*Qc(isp)       ! Phytoplankton respiration     (mmol-C/m3/d)	
     AexudN   = AexudN + Aresp(k,isp)*Qn(k,isp)     ! Total Phytoplankton exudation (mmol-N/m3/d)
     AexudP   = AexudP + Aresp(k,isp)*Qp(k,isp)     ! Total Phytoplankton exudation (mmol-P/m3/d)
     Amort(isp)  = A(k,isp) * mA(isp)                     ! Dead phytoplankton (cells/m3/day)

!    Monod Equations
     monodN(isp)  = Ntotal(k)/(Ntotal(k)+Kn(isp))
     monodP(isp)  = PO4(k)/(PO4(k)+Kp(isp))
     monodSi(isp) = Si(k)/(Si(k)+Ksi(isp))

#ifdef DEBUG
write(6,*) "In cgem, set Amort, Ntotal monod" 
#endif

!------------------------------------------------------------------------     
! Nutrient limited uptake:
! Find Rate Limiting Nutrient RLN for N, P, and Si:
! N==1, P==2, Si==3
   if(Ntotal(k).le.PO4(k).and.Ntotal(k).le.Si(k)) then
     RLN = 1
   elseif (PO4(k).le.Ntotal(k).and.PO4(k).le.Si(k)) then
     RLN = 2
   else
     RLN = 3
   endif

!------------------------------------------------------------------------
   if(Is_Day.eq.0) then  !Nutrient uptake only takes place during the day
        vN(k) = 0
        vP(k) = 0
        vSi(k) = 0
   else  !Day
      !Rate limiting nutrient is N
      if(RLN(k).eq.1) then
         vN(k) = Q10_T(T(k),vmaxN(isp))*monodN(isp)*f_Qn(k,isp)

         vP(k) = Q10_T(T(k),vmaxP(isp))*monodP(isp)*f_Qp(k,isp)&
     &      *( Ntotal(k)/(Ntotal(k)+aN(isp)*Kn(isp)) )

         vSi(k) = Q10_T(T(k),vmaxSi(isp))*monodSi(isp)        &
     &      *( Ntotal(k)/(Ntotal(k)+aN(isp)*Kn(isp)) )
      !Rate limiting nutrient is P
      elseif(RLN(k).eq.2) then
         vN(k) = Q10_T(T(k),vmaxN(isp))*monodN(isp)*f_Qn(k,isp)&
     &      *( PO4(k)/(PO4(k)+aN(isp)*Kp(isp)) )

         vP(k) = Q10_T(T(k),vmaxP(isp))*monodP(isp)*f_Qp(k,isp)

         vSi(k) = Q10_T(T(k),vmaxSi(isp))*monodSi(isp)       &
     &      *( PO4(k)/(PO4(k)+aN(isp)*Kp(isp)) )
      !Rate limiting nutrient is Si 
      else
         vN(k) = Q10_T(T(k),vmaxN(isp))*monodN(isp)*f_Qn(k,isp)&
     &      *( Si(k)/(Si(k)+aN(isp)*Ksi(isp)) )

         vP(k) = Q10_T(T(k),vmaxP(isp))*monodP(isp)*f_Qp(k,isp)&
     &      *( Si(k)/(Si(k)+aN(isp)*Ksi(isp)) )

         vSi(k) = Q10_T(T(k),vmaxSi(isp))*monodSi(isp)

      endif

   endif !Endif Is_Day.eq.0

     
#ifdef DEBUG
write(6,*) "In cgem, finished Nutrients"
#endif
 
!--------------------------------------------------------------      
! When isp loop is done, AupN and AupP are totals for all nospA 
! phytoplankton groups in cell k          
!---------------------------------------------------------------   
      AupN(k) = AupN(k) + A(k,isp)*vN(k)     ! Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
      AupP(k) = AupP(k) + A(k,isp)*vP(k)     ! Phytoplankton uptake of Phosphorus (mmol-P/m3/d) 
      AupSi(k) = AupSi(k) + A(k,isp)*vSi(k)   ! Phytoplankton uptake of Silica (mmol-Si/m3/d)
 
!-----------------------------------------------------------------------
! Note that Zumax(1)*monodZ(isp,1) is volume of type isp phytoplankton
! eaten per day by type 1 zooplankton. Therefore
!      Zumax(1)*monodZ(isp,1)/volcell(isp)
! is the number of type isp phytoplankton eaten per day by type 1 
! zooplankton. An analogous statement holds for type 2 zooplankton
!-----------------------------------------------------------------------
      Zgrazvol(isp,:)     = Z(k,:)*Zumax(:)*monodZ(k,isp,:)   ! Grazing of phytoplankton by biovolume (um3/m3/d)
      ZgrazA(isp,:)       = Zgrazvol(isp,:)/volcell(isp)  ! Grazing of phytoplankton (cells/m3/d)
      ZgrazA_tot(isp) = SUM( ZgrazA(isp,:) ) 

!----------------------------------------------------------------------
! When the isp loop is finished, ZgrazC, ZgrazN, and ZgrazP will be total
! carbon, nitrogen, and phosphorous uptake of zooplankton from grazing
! all phytoplankton groups
!---------------------------------------------------------------------
      ZgrazC(:) = ZgrazC(:) + ZgrazA(isp,:) * Qc(isp)     ! Carbon uptake from grazing (mmol-C/m3/day)
      ZgrazN(:) = ZgrazN(:) + ZgrazA(isp,:) * Qn(k,isp)   ! Nitrogen uptake from grazing( mmol-N/m3/day)
      ZgrazP(:) = ZgrazP(:) + ZgrazA(isp,:) * Qp(k,isp)   ! Phosphorus uptake from grazing (mmol-P/m3/day)

!---------------------------------------------------------
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------
      ff_new(k,iA(isp)) = AMAX1(A(k,isp)        &
      & + ( Agrow(k,isp) - Aresp(k,isp) - ZgrazA_tot(isp) - Amort(isp) )*dTd, 1.)
    
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
!! Enforce minima for Droop, also enforce maxima if not equal Droop (which_quota=1)
  if(which_quota.eq.1) then
     ff_new(k,iQn(isp)) = AMAX1(Qn(k,isp) + (vN(k) - Qn(k,isp)*uA(k,isp))*dTd,QminN(isp))
!! , also enforce maxima if not equal Droop (which_quota=1)
  else
     ff_new(k,iQn(isp)) = AMIN1(AMAX1(Qn(k,isp) + (vN(k) - Qn(k,isp)*uA(k,isp))*dTd,QminN(isp)),QmaxN(isp))
  endif

!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
!! Enforce minima for Droop, also enforce maxima if not equal Droop (which_quota=1)
  if(which_quota.eq.1) then
     ff_new(k,iQp(isp)) = AMAX1(Qp(k,isp) + (vP(k) - Qp(k,isp)*uA(k,isp))*dTd,QminP(isp))
!! , also enforce maxima if not equal Droop (which_quota=1)
  else
     ff_new(k,iQp(isp)) = AMIN1(AMAX1(Qp(k,isp) + (vP(k) - Qp(k,isp)*uA(k,isp))*dTd,QminP(isp)),QmaxP(isp))
  endif

!----------------------------------------------------------------------- 
  enddo  ! END OF do isp = 1, nospA 

!-------------------------------------------------------------------
! Now calculate the total ingested ZinC, ZinN, and ZinP of C, N, and P
!-------------------------------------------------------------------
      ZslopC(:)  = Zslop(:)*ZgrazC(:)                      ! Sloppy feeding (mmol-C/m3/d)
      ZslopC_tot = SUM(ZslopC)                          ! Total Sloppy feeding (mmol-C/m3/d)
      ZunC(:)    = (1.-Zeffic(:))*(ZgrazC(:)-ZslopC(:)) ! Unassimilated (mmol-C/m3/d)
      ZinC(:)    = ZgrazC(:) - ZslopC(:) - ZunC(:)         ! Ingested (mmol-C/m3/d)

      ZslopN(:)  = Zslop(:)*ZgrazN(:)                      ! Sloppy feeding (mmol-N/m3/d)
      ZslopN_tot = SUM(ZslopN)                             ! Total Sloppy feeding (mmol-N/m3/d) 
      ZunN(:)    = (1.-Zeffic(:))*(ZgrazN(:)-ZslopN(:)) ! Unassimilated (mmol-N/m3/d)
      ZinN(:)    = ZgrazN(:) - ZslopN(:) - ZunN(:)         ! Ingested (mmol-N/m3/d)

      ZslopP(:)  = Zslop(:)*ZgrazP(:)                      ! Sloppy feeding (mmol-P/m3/d)
      ZslopP_tot = SUM(ZslopP)                             ! Total Sloppy feeding (mmol-P/m3/d)
      ZunP(:)    = (1.-Zeffic(:))*(ZgrazP(:)-ZslopP(:)) ! Unassimilated (mmol-P/m3/d)
      ZinP(:)    = ZgrazP(:) - ZslopP(:) - ZunP(:)         ! Ingested (mmol-P/m3/d)
!-------------------------------------------------


!------------------------------------         
! Liebigs Law for zooplankton group isz 
!------------------------------------
  do isz=1,nospZ

     if (ZinN(isz) .gt. optNP(isz)*ZinP(isz)) then  
        Zgrow(k,isz)= ZinP(isz)/ZQp(isz)                   ! P-limited growth (indv./m3/d) 
        ZegN(isz) = ZinN(isz) - ZinP(isz)*optNP(isz)       ! P-limited N excretion (mmol-N/m3/d) 
                                                   ! determined by subtracting N-equivalent of ZinP
        ZegC(isz) = ZinC(isz) - ZinP(isz)/ZQp(isz)*ZQc(isz)  ! P-limited C excretion (mmol-C/m3/d)
        ZegP(isz) = 0.                        
      else
        Zgrow(k,isz)= ZinN(isz)/ZQn(isz)                   ! N-limited growth (indv./m3/d)
        ZegP(isz) = ZinP(isz) - ZinN(isz)/optNP(isz)       ! N-limited P excretion (mmol-P/m3/d)    
                                                   ! determined by subtracting P-equivalent of ZinN
        ZegC(isz) = ZinC(isz) - ZinN(isz)/ZQn(isz)*ZQc(isz)  ! N-limited C excretion (mmol-C/m3/d)
        ZegN(isz) = 0.
      endif

  enddo

!------------------------------------------------

!-----------------------------------------------------
! ZegC should not be negative 
  do isz=1,nospZ
      if(ZegC(isz).lt.0.) then
          ZegC(isz) = 0.
          !ZegN(isz) = 0.
          !ZegP(isz) = 0.
      endif
  enddo

! Egestion and unassimilated for Si set equivalent to that of N
      ZegSi = ZegN
      ZunSi = ZunN
      sumZegSi(k) = SUM(ZegSi)
      sumZunSi(k) = SUM(ZunSi)

! Zooplankton respiration based on growth and basal metabolism, both modified by a temperature adjustment factor 
      Zresp(k,:) = (Zgrow(k,:)*Zrespg(:) + Z(k,:)*Zrespb(:)) !Zooplankton respiration (indv./m3/d)

      ZrespC(k)   = SUM(Zresp(k,:)*ZQc)                                       !Total Carbon loss from respiration (mmol-C/m3/d)

                                                ! Excretion
      ZexN(:)   = Zresp(k,:)*ZQn(:)               ! (mmol-N/m3/d)
      sumZexN(k) = SUM(ZexN)
      ZexP(:)   = Zresp(k,:)*ZQp(:)               ! (mmol-P/m3/d)
      sumZexP(k) = SUM(ZexP)
                                                ! Mortality
     Zmort(k,:)       = Zm(:) * Z(k,:) * Z(k,:)       ! (indv./m3/d)
     ZmortC(:)      = Zmort(k,:)*ZQc(:)           ! (mmol-C/m3/d)
     ZmortC_tot     = SUM(ZmortC)
     ZmortN(:)      = Zmort(k,:)*ZQn(:)           ! (mmol-N/m3/d)
     ZmortN_tot     = SUM(ZmortN)
     ZmortP(:)      = Zmort(k,:)*ZQp(:)           ! (mmol-P/m3/d)
     ZmortP_tot     = SUM(ZmortP)
!-------------------------------------------------------------------------

!---------------------------------------------------------
!-G; Zooplankton number density (individuals/m3);
!---------------------------------------------------------
      ff_new(k,iZ(:))  = AMAX1( Z(k,:)                         &
      &      + (Zgrow(k,:) - Zresp(k,:) - Zmort(k,:))*dTd, 1.)

#ifdef DEBUG
write(6,*) "In cgem, updated iZ, which_Fluxes(iInRemin),KG_bot=",which_Fluxes(iInRemin),KG_bot
#endif

!------------------------------------------------------------------------

!-----------------------------------------------------------
! Remineralization - reactions
!---------------------------------------------------------------
       ! Instant Remineralization, if on bottom of shelf, redefine KG's
       if(k.eq.km.and.which_fluxes(iInRemin).eq.1) then
           KG1 = KG_bot
           KG2 = KG_bot
       endif
!------------------------------------------------------------
! Nitrification
!--------------------------------------------------------------
        call Nitrification( O2(k), NH4(k), KO2, KNH4, nitmax, T(k), R_11 )

#ifdef DEBUG
write(6,*) "In cgem, called Nitrification" 
#endif

!------------------------------------------------------------
! Carbon Chemistry
!--------------------------------------------------------------
!!! MOCSY alkalinity expressions:
        m_alk = ALK(k)/1000.
        m_dic = DIC(k)/1000.
        m_si  = Si(k)/1000.
        m_po4 = PO4(k)/1000.
        call vars(ph_calc, pco2_calc, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD_calc, rhoSW, p, tempis,&
         &    T(k), S(k), m_alk, m_dic, m_si, m_po4, patm, d_sfc(k), m_lat, 1, &
         &    'mol/m3', 'Tinsitu', 'm ', 'u74', 'l  ', 'pf ', 'Pzero  ')
        pH(k) = ph_calc(1)
#ifdef DEBUG
write(6,*) "In cgem, called mocsy"
#endif

!write(6,*) "In cgem, begin reaction A, istep,k,inea,myrank", istep,k,inea,myrank

!------------------------------------------------------------
! Particulate and Dissolved dead phytoplankton, rate of remineralization
!--------------------------------------------------------------
        call reaction( OM1A(k), OM2A(k), O2(k), NO3(k), KG1, KG2, KO2, KstarO2, KNO3,               &
     &  sx1A(k), sy1A(k), sz1, sx2A(k), sy2A(k), sz1, T(k), RC )

        RC        = one_d_365 * RC  !Change units from /year to /day

        ROM1_A(k)     = RC(1)          ! units are /m3/day
        ROM2_A(k)     = RC(2)
        RO2_A      = RC(3)
        RNO3_A     = RC(4)
        RPO4_A     = RC(5)
        RDIC_A     = RC(6)
        RNH4_A     = RC(7)
        RSi_A      = RC(8)
        RALK_A     = RC(9)
        RN2_A      = RC(10)

#ifdef DEBUG
write(6,*) "In cgem, finished reaction A"
#endif


!write(6,*) "In cgem, begin reaction Z"
!------------------------------------------------------------
! Particulate and Dissolved fecal pellets, rate of remineralization
!--------------------------------------------------------------
        call reaction( OM1Z(k), OM2Z(k), O2(k), NO3(k), KG1, KG2, KO2, KstarO2, KNO3,               &
     &  sx1Z(k), sy1Z(k), sz1, sx2Z(k), sy2Z(k), sz1, T(k), RC )
        RC       = one_d_365 * RC   !Change units from /year to /day

        ROM1_Z(k)     = RC(1)         ! units are /m3/day
        ROM2_Z(k)     = RC(2)
        RO2_Z      = RC(3)
        RNO3_Z     = RC(4)
        RPO4_Z     = RC(5)
        RDIC_Z     = RC(6)
        RNH4_Z     = RC(7)
        RSi_Z      = RC(8)
        RALK_Z     = RC(9)
        RN2_Z      = RC(10)

#ifdef DEBUG
write(6,*) "In cgem, finished reaction Z"
#endif

!write(6,*) "In cgem, begin reaction R"
!------------------------------------------------------------
! Particulate and Dissolved riverine OM, rate of remineralization 
!------------------------------------------------------------
        call reaction( OM1R(k), OM2R(k), O2(k), NO3(k), KG1_R, KG2_R, KO2, KstarO2, KNO3,               &
     &  sx1R, sy1R, sz1, sx2R, sy2R, sz1, T(k), RC )
        RC       = one_d_365 * RC   !Change units from /year to /day

        ROM1_R(k)     = RC(1)         ! units are /m3/day
        ROM2_R(k)     = RC(2)
        RO2_R      = RC(3)
        RNO3_R     = RC(4)
        RPO4_R     = RC(5)
        RDIC_R     = RC(6)
        RNH4_R     = RC(7)
        RSi_R      = RC(8)
        RALK_R     = RC(9)
        RN2_R      = RC(10)

#ifdef DEBUG
write(6,*) "In cgem, finished reaction R"
#endif

!write(6,*) "In cgem, begin reaction BC"
!------------------------------------------------------------
! Particulate and Dissolved initial and boundary OM, rate of remineralization
!------------------------------------------------------------
        call reaction( OM1BC(k), OM2BC(k), O2(k), NO3(k), KG1_BC, KG2_BC, KO2, KstarO2, KNO3,               &
     &  sx1BC, sy1BC, sz1, sx2BC, sy2BC, sz1, T(k), RC )

        RC       = one_d_365 * RC   !Change units from /year to /day

        ROM1_BC(k)     = RC(1)         ! units are /m3/day
        ROM2_BC(k)     = RC(2)
        RO2_BC      = RC(3)
        RNO3_BC     = RC(4)
        RPO4_BC     = RC(5)
        RDIC_BC     = RC(6)
        RNH4_BC     = RC(7)
        RSi_BC      = RC(8)
        RALK_BC     = RC(9)
        RN2_BC      = RC(10)

       ! Instant Remineralization, change KG's back to original
       if(k.eq.km.and.which_fluxes(iInRemin).eq.1) then
           KG1 = KG1_save
           KG2 = KG2_save
       endif

#ifdef DEBUG
write(6,*) "In cgem, finished reaction BC"
#endif


!--------------------------------------------------------------------
! Sum remineralization terms from dead phytoplankton, fecal pellets, and riverine particulate
  RO2(k)   = RO2_A  + RO2_Z  + RO2_R  + RO2_BC - 2.*R_11  ! (mmol-O2/m3/d)
  RNO3(k)  = RNO3_A + RNO3_Z + RNO3_R + RNO3_BC + R_11    ! (mmol-NO3/m3/d)
  RNH4(k)  = RNH4_A + RNH4_Z + RNH4_R + RNH4_BC - R_11    ! (mmol-NH4/m3/d)
  RPO4(k)  = RPO4_A + RPO4_Z + RPO4_R + RPO4_BC           ! (mmol-PO4/m3/d)
  RDIC(k)  = RDIC_A + RDIC_Z + RDIC_R + RDIC_BC           ! (mmol-DIC/m3/d)
  RSi(k)   = RSi_A  + RSi_Z  + RSi_R  + RSi_BC            ! (mmol-Si/m3/d)
  RALK(k)  = RALK_A + RALK_Z + RALK_R + RALK_BC - 2.*R_11 ! (mmol-HCO3/m3/d)
  RN2(k)   = RN2_A + RN2_Z + RN2_R + RN2_BC         ! (mmol-N2/m3/d)
!--------------------------------------------------------------------
if(writecsv==1.and.k.eq.1.and.inea.eq.10) then
    write(6101,'(*(g0,:,", "))') &
    & RO2(k), RNO3(k), RNH4(k), RPO4(k), RDIC(k),RSi(k),RALK(k),RN2(k), &
    & ROM1_A(k), ROM2_A(k), RO2_A, RNO3_A, RPO4_A, RDIC_A, RNH4_A, RSi_A, RALK_A, RN2_A, &
    & ROM1_Z(k), ROM2_Z(k), RO2_Z, RNO3_Z, RPO4_Z, RDIC_Z, RNH4_Z, RSi_Z, RALK_Z, RN2_Z, &
    & ROM1_R(k), ROM2_R(k), RO2_R, RNO3_R, RPO4_R, RDIC_R, RNH4_R, RSi_R, RALK_R, RN2_R, &
    & ROM1_BC(k),ROM2_BC(k),RO2_BC,RNO3_BC,RPO4_BC,RDIC_BC,RNH4_BC,RSi_BC,RALK_BC,RN2_BC
endif

!
!---------------------------------------------------------------------
! Stoichiometry - calculate C:N:P ratios for Remineralization equations
!---------------------------------------------------------------------
!-- Organic Matter from dead phytoplankton --------------------------
OM1_CA(k) = 0.
OM2_CA(k) = 0.
OM1_NA(k) = 0.
OM2_NA(k) = 0.
OM1_PA(k) = 0.
OM2_PA(k) = 0.

do isp=1,nospA
 if ( uN(k,isp) .lt. uP(k,isp)  ) then
!Particulate
   OM1_CA(k) = OM1_CA(k) + Amort(isp)*(Qn(k,isp)-QminN(isp))/Qn(k,isp)*Qc(isp)
   OM1_NA(k) = OM1_NA(k) + Amort(isp)*(Qn(k,isp)-QminN(isp))
   OM1_PA(k) = OM1_PA(k) + Amort(isp)*(Qn(k,isp)-QminN(isp))/Qn(k,isp)*Qp(k,isp)
!!Dissolved
   OM2_CA(k) = OM2_CA(k) + Amort(isp)*QminN(isp)/Qn(k,isp)*Qc(isp)
   OM2_NA(k) = OM2_NA(k) + Amort(isp)*QminN(isp)
   OM2_PA(k) = OM2_PA(k) + Amort(isp)*QminN(isp)/Qn(k,isp)*Qp(k,isp)
 else
!Particulate
   OM1_CA(k) = OM1_CA(k) + Amort(isp)*(Qp(k,isp)-QminP(isp))/Qp(k,isp)*Qc(isp)
   OM1_NA(k) = OM1_NA(k) + Amort(isp)*(Qp(k,isp)-QminP(isp))/Qp(k,isp)*Qn(k,isp)
   OM1_PA(k) = OM1_PA(k) + Amort(isp)*(Qp(k,isp)-QminP(isp))
!!Dissolved
   OM2_CA(k) = OM2_CA(k) + Amort(isp)*QminP(isp)/Qp(k,isp)*Qc(isp)
   OM2_NA(k) = OM2_NA(k) + Amort(isp)*QminP(isp)/Qp(k,isp)*Qn(k,isp)
   OM2_PA(k) = OM2_PA(k) + Amort(isp)*QminP(isp)
 endif
enddo

                                             ! Dissolved
   !This calculates the cumulative stoichiometry ratios for OM1_A
   if(OM1_CA(k).gt.tiny(x)) then
!   if(OM1_CA.ne.0) then
    ff_new(k,isx1A) = (OM1_CA(k)*dTd + OM1A(k)) / (OM1_PA(k)*dTd + (1./sx1A(k))*OM1A(k)) ! C/P
    ff_new(k,isy1A) = (OM1_NA(k)*dTd + (sy1A(k)/sx1A(k))*OM1A(k)) / (OM1_PA(k)*dTd + (1./sx1A(k))*OM1A(k)) !N/P
    !stoich_z1A = 1., using sz1=1
   else
    ff_new(k,isx1A) = sx1A(k)
    ff_new(k,isy1A) = sy1A(k)
   endif

#ifdef DEBUG
write(6,*) "In cgem, finished stoich OM1A"
#endif

   !This calculates the cumulative stoichiometry ratios for OM2_A
   if(OM2_CA(k).gt.tiny(x)) then
!    if(OM2_CA.ne.0) then
    ff_new(k,isx2A) = (OM2_CA(k)*dTd + OM2A(k)) / (OM2_PA(k)*dTd + (1./sx2A(k))*OM2A(k)) ! C/P
    ff_new(k,isy2A) = (OM2_NA(k)*dTd + (sy2A(k)/sx2A(k))*OM2A(k)) / (OM2_PA(k)*dTd + (1./sx2A(k))*OM2A(k)) !N/P
   else
    ff_new(k,isx2A) = sx2A(k)
    ff_new(k,isy2A) = sy2A(k)
   endif

#ifdef DEBUG
write(6,*) "In cgem, finished stoich OM2A"
#endif


!!-- Organic Matter from fecal pellets ---------------------------------
    OM1_Ratio(k) = SUM( (Qn(k,:)-QminN)/Qn(k,:)*A(k,:))/SUM(A(k,:))
    OM2_Ratio(k) = SUM( (QminN/Qn(k,:))*A(k,:))/SUM(A(k,:))

#ifdef DEBUG
write(6,*) "In cgem, OM1,2 Ratio=",OM1_Ratio,OM2_Ratio
#endif

    if(nospZ.eq.1) then 
                                                                  ! Particulate
     OM1_CZ(k)  = .5*(ZegC(1) + ZunC(1) + ZmortC_tot) + OM1_Ratio(k)*ZslopC_tot !  (mmol-C/m3/d)
     OM1_NZ(k)  = .5*(ZegN(1) + ZunN(1) + ZmortN_tot) + OM1_Ratio(k)*ZslopN_tot !  (mmol-N/m3/d)
     OM1_PZ(k)  = .5*(ZegP(1) + ZunP(1) + ZmortP_tot) + OM1_Ratio(k)*ZslopP_tot !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ(k)  = .5*(ZegC(1) + ZunC(1) + ZmortC_tot) + OM2_Ratio(k)*ZslopC_tot              !  (mmol-C/m3/d)
     OM2_NZ(k)  = .5*(ZegN(1) + ZunN(1) + ZmortN_tot) + OM2_Ratio(k)*ZslopN_tot              !  (mmol-N/m3/d)
     OM2_PZ(k)  = .5*(ZegP(1) + ZunP(1) + ZmortP_tot) + OM2_Ratio(k)*ZslopP_tot              !  (mmol-P/m3/d)

    else if(nospZ.eq.2) then
                                                                  ! Particulate
     OM1_CZ(k)  = ZegC(1) + ZunC(1) + ZmortC_tot + OM1_Ratio(k)*ZslopC_tot !  (mmol-C/m3/d)
     OM1_NZ(k)  = ZegN(1) + ZunN(1) + ZmortN_tot + OM1_Ratio(k)*ZslopN_tot !  (mmol-N/m3/d)
     OM1_PZ(k)  = ZegP(1) + ZunP(1) + ZmortP_tot + OM1_Ratio(k)*ZslopP_tot !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ(k)  = ZegC(2) + ZunC(2) + OM2_Ratio(k)*ZslopC_tot              !  (mmol-C/m3/d)
     OM2_NZ(k)  = ZegN(2) + ZunN(2) + OM2_Ratio(k)*ZslopN_tot              !  (mmol-N/m3/d)
     OM2_PZ(k)  = ZegP(2) + ZunP(2) + OM2_Ratio(k)*ZslopP_tot              !  (mmol-P/m3/d)

    else 

                                                                  ! Particulate
     OM1_CZ(k)  = ZegC(1) + ZunC(1) + ZmortC_tot + OM1_Ratio(k)*ZslopC_tot !  (mmol-C/m3/d)
     OM1_NZ(k)  = ZegN(1) + ZunN(1) + ZmortN_tot + OM1_Ratio(k)*ZslopN_tot !  (mmol-N/m3/d)
     OM1_PZ(k)  = ZegP(1) + ZunP(1) + ZmortP_tot + OM1_Ratio(k)*ZslopP_tot !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ(k)  = SUM(ZegC(2:nospZ)) + SUM(ZunC(2:nospZ)) + OM2_Ratio(k)*ZslopC_tot              !  (mmol-C/m3/d)
     OM2_NZ(k)  = SUM(ZegN(2:nospZ)) + SUM(ZunN(2:nospZ)) + OM2_Ratio(k)*ZslopN_tot              !  (mmol-N/m3/d)
     OM2_PZ(k)  = SUM(ZegP(2:nospZ)) + SUM(ZunP(2:nospZ)) + OM2_Ratio(k)*ZslopP_tot              !  (mmol-P/m3/d)

    endif

   !This calculates the cumulative stoichiometry ratios for OM1_Z
   if(OM1_CZ(k).gt.tiny(x)) then
!    if(OM1_CZ.ne.0) then
    ff_new(k,isx1Z) = (OM1_CZ(k)*dTd + OM1Z(k)) / (OM1_PZ(k)*dTd + (1./sx1Z(k))*OM1Z(k)) ! C/P
    ff_new(k,isy1Z) = (OM1_NZ(k)*dTd + (sy1Z(k)/sx1Z(k))*OM1Z(k)) / (OM1_PZ(k)*dTd + (1./sx1Z(k))*OM1Z(k)) !N/P
   else
    ff_new(k,isx1Z) = sx1Z(k)
    ff_new(k,isy1Z) = sy1Z(k)
   endif

   !This calculates the cumulative stoichiometry ratios for OM2_Z
   if(OM2_CZ(k).gt.tiny(x)) then
!    if(OM2_CZ.ne.0) then
    ff_new(k,isx2Z) = (OM2_CZ(k)*dTd + OM2Z(k)) / (OM2_PZ(k)*dTd + (1./sx2Z(k))*OM2Z(k)) !  C/P
    ff_new(k,isy2Z) = (OM2_NZ(k)*dTd + (sy2Z(k)/sx2Z(k))*OM2Z(k)) / (OM2_PZ(k)*dTd + (1./sx2Z(k))*OM2Z(k)) !N/P
   else
    ff_new(k,isx2Z) = sx2Z(k)
    ff_new(k,isy2Z) = sy2Z(k)
   endif

enddo !end k loop

!-------------------------------
!-NO3; (mmol-N/m3)
  ff_new(:,iNO3) = NO3(:)                            &
  &  + ( RNO3(:) - AupN(:)*NO3(:)/Ntotal(:))*dTd
!--------------------------------
!-NH4; Ammonium (mmol-N/m3)
  ff_new(:,iNH4) = ff(:,iNH4)                            &
  & + ( RNH4(:) - AupN(:)*NH4(:)/(Ntotal(:)) + AexudN(:) + sumZexN(:)  )*dTd
!----------------------------
!-Silica: (mmol-Si/m3)
  ff_new(:,iSi) =  Si(:)                             &
  & + ( RSi(:) - AupSi(:) + sumZegSi(:) + sumZunSi(:) )*dTd
!---------------------------------------------
!-PO4: Phosphate (mmol-P/m3)
 ff_new(:,iPO4) = PO4(:)                             &
 & + ( RPO4(:) - AupP(:) + AexudP(:) + sumZexP(:)  )*dTd

!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3)
  ff_new(:,iDIC) = DIC(:)                            &
  &  + ( RDIC(:) - PrimProd(:) + ArespC(:)  + ZrespC(:) )*dTd
!-----------------------------------------------------------------------      
!-O2: Oxygen (mmol O2 m-3) 
  ff_new(:,iO2)  = O2(:)                             &  
  &  + ( PrimProd(:) - ArespC(:) + RO2(:) - ZrespC(:))*dTd
!-----------------------------------------
!-OM1_A: (mmol-C/m3-- Dead Phytoplankton Particulate)
  ff_new(:,iOM1_A) = OM1A(:) + (ROM1_A(:) + OM1_CA(:))*dTd
!-----------------------------------------------
!-OM2_A: (mmol-C/m3-- Dead Phytoplankton Dissolved)
  ff_new(:,iOM2_A) = OM2A(:) + (ROM2_A(:) + OM2_CA(:))*dTd
!-----------------------------------------------
!-OM1_Z:(mmol-C/m3--G particulate)
  ff_new(:,iOM1_Z) = OM1Z(:) + (ROM1_Z(:) + OM1_CZ(:))*dTd
!-----------------------------------------------
!-OM2_Z:(mmol-C/m3--G dissolved)
  ff_new(:,iOM2_Z) = OM2Z(:) + (ROM2_Z(:) + OM2_CZ(:))*dTd
!---------------------------------------------------------------------
!-OM1_R: (mmol-C/m3--SPM particulate)
  ff_new(:,iOM1_R) = OM1R(:) + ROM1_R(:)*dTd
!---------------------------------------------------------------------
!-OM2_R: (mmol-C/m3--SPM dissolved)
  ff_new(:,iOM2_R) = OM2R(:) + ROM2_R(:)*dTd
!---------------------------------------------------------------------
!-OM1_BC: (mmol-C/m3--initial and boundary condition OM particulate)
  ff_new(:,iOM1_BC) = OM1BC(:) + ROM1_BC(:)*dTd
!---------------------------------------------------------------------
!-OM2_BC: (mmol-C/m3--initial and boundary condition OM dissolved)
  ff_new(:,iOM2_BC) = OM2BC(:) + ROM2_BC(:)*dTd
!---------------------------------------------------------------------
!-CDOM: (ppb) 
  ff_new(:,iCDOM) =  CDOM(:)*(1.0 - KGcdom*dTd)
!!---------------------------------------------------------------------
!!-ALK: (mmol-HCO3/m3)
  ff_new(:,iALK) =  ALK(:) +                  &
  & (RALK(:) + AupN(:)*NO3(:)/(Ntotal(:))           &
  &          - AupN(:)*NH4(:)/(Ntotal(:))           &
  &          + AupP(:) + 4.8*AupP(:))*dTd

  ff_new(:,iTr) = ff(:,iTr) 

do k=1,km
  do isp = iNO3,nf
    if(ff_new(k,isp).le.0.) then
      write(6,*) "ff_new.le.0! set to 0 for istep,myrank,inea,k,isp=",istep,myrank,inea,k,isp,ff_new(k,isp)
      ff_new(k,isp) = 0.0
    endif
  enddo
enddo

!-- End Main GEM Calculations ---------------------------------------------------

  !Write csv files
  k=1
  if(writecsv==1.and.k.eq.1.and.inea.eq.10) then
    write(6301,'(*(g0,:,", "))') TC_8,rad,wind,S(k),T(k)
    write(6401,'(*(g0,:,", "))') PARsurf,PARdepth,PARbot
  endif

  !! Before Advection and VMixing, combine A's and Q's
  do isp=1,nospA
    ff_new(:,iQn(isp)) = ff_new(:,iQn(isp)) * ff_new(:,iA(isp))
    ff_new(:,iQp(isp)) = ff_new(:,iQp(isp)) * ff_new(:,iA(isp))
  enddo

  !! Before Advection and VMixing, combine stoich with OMs 
   ff_new(:,isx1A) = ff_new(:,isx1A) * ff_new(:,iOM1_A)
   ff_new(:,isy1A) = ff_new(:,isy1A) * ff_new(:,iOM1_A)
   ff_new(:,isx2A) = ff_new(:,isx2A) * ff_new(:,iOM2_A)
   ff_new(:,isy2A) = ff_new(:,isy2A) * ff_new(:,iOM2_A)
   ff_new(:,isx1Z) = ff_new(:,isx1Z) * ff_new(:,iOM1_Z)
   ff_new(:,isy1Z) = ff_new(:,isy1Z) * ff_new(:,iOM1_Z)
   ff_new(:,isx2Z) = ff_new(:,isx2Z) * ff_new(:,iOM2_Z)
   ff_new(:,isy2Z) = ff_new(:,isy2Z) * ff_new(:,iOM2_Z)


   return
   end subroutine cgem_step
!---------------------------------------------------------------------- 
