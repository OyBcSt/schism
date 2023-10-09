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
    integer, intent(in)  :: myrank !process number
    real, intent(in)     :: dT
!---------------------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------
    integer        ::  k, isp, isz ! Loop indicies, isp/isz is for phytoplankton/zooplankton species
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
    real, dimension(km,nospA) :: f_Qn  ! Quota model for N
    real, dimension(km,nospA) :: f_Qp  ! Quota model for P
    real, dimension(km,nospA) :: vN    ! Phytoplankton uptake rate of Nitrogen (mmol-N/cell/d)
    real, dimension(km,nospA) :: vP    ! Phytoplankton uptake rate of Phosphorus (mmol-P/cell/d)
    real, dimension(km,nospA) :: vSi   ! Phytoplankton uptake rate of Silica (mmol-Si/cell/d)
    real, dimension(km) :: AupN  ! Total Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
    real, dimension(km) :: AupP  ! Total Phytoplankton uptake of Phosphorus (mmol-P/m3/d)
    real, dimension(km) :: AupSi ! Total Phytoplankton uptake of Silica (mmol-Si/m3/d)
 ! Monod equations for phytoplankton
    real, dimension(km,nospA) :: monodN  !Monod term in nitrogen uptake
    real, dimension(km,nospA) :: monodP  !Monod term in phosphorus uptake
    real, dimension(km,nospA) :: monodSi !Monod term in Si uptake
    real, dimension(km)       :: Ntotal   ! Total N (mmol/m3)
 ! Phytoplankton nutrient loss
    real, dimension(km,nospA) :: Amort   ! Dead phytoplankton (cells/m3/day)
    real, dimension(km)       :: AexudN  ! Sum of Exudation of N from all phytoplankton groups (mmol-N/m3/d)
    real, dimension(km)       :: AexudP  ! Sum of Exudation of P from all phytoplankton groups (mmol-P/m3/d)
    real, dimension(km)       :: ArespC  ! Phytoplankton equivalent carbon loss from respiration (mmol-C/m3/d)
!------------------------------------------------------------------
! Zooplankton parameters
 !Zooplankton uptake and growth
    real, dimension(km,nospZ)    :: Z         ! Zooplankton number density (indv./m3)
    real, dimension(nospZ)       :: optNP    ! Optimal nutrient ratio for zooplankton
    real, dimension(km,nospZ)    :: Zgrow    ! Zooplankton growth (indv./m3/d)
    real, dimension(km,nospA,nospZ) :: Zgrazvol ! Grazing rate in units of biovolume (um3/m3/d)
    real, dimension(km,nospA,nospZ) :: ZgrazA   ! Zooplankton grazing of phytoplankton (cells/m3/d)
    real, dimension(km,nospA)    :: ZgrazA_tot ! Total zooplankton grazing of phytoplankton (cells/m3/d)
    real, dimension(km,nospZ)    :: ZgrazN   ! Zooplankton grazing uptake of Nitrogen (mmol-N/m3/d)
    real, dimension(km,nospZ)    :: ZgrazP   ! Zooplankton grazing uptake of Phosphorus (mmol-P/m3/d)
    real, dimension(km,nospZ)    :: ZgrazC   ! Zooplankton grazing uptake of Carbon (mmol-C/m3/d)
    real, dimension(km,nospZ)    :: ZinN     ! Zooplankton ingestion of Nitrogen (mmol-N/m3/d)
    real, dimension(km,nospZ)    :: ZinP     ! Zooplankton ingestion of Phosphorus (mmol-P/m3/d)
    real, dimension(km,nospZ)    :: ZinC     ! Zooplankton ingestion of Carbon (mmol-C/m3/d)
 !Monod equations for zooplankton ingestion of phytoplankton
    real, dimension(km,nospA,nospZ) :: monodZ   ! Monod term for zooplankton grazing
    real                            :: Abiovol  ! Algae biovolume vector (um3/m3)
    real, dimension(nospA,nospZ)    :: top_A    ! Monod numerator value for phytoplankton group
    real, dimension(nospA,nospZ)    :: bottom_A ! Monod Denominator value for phytoplankton group
    real, dimension(nospZ)          :: bottom   ! Sum of Monod Denominator value for all phytoplankton groups
 !Zooplankton nutrient loss
    real, dimension(km,nospZ)       :: Zresp    ! Zooplankton respiration (individuals/m3/d)
    real, dimension(km)             :: ZrespC   ! Carbon loss from zooplankton respiration (mmol-C/m3/day)
    real, dimension(km,nospZ)       :: ZunC     ! Unassimilated ingested Carbon (mmol-C/m3/d)
    real, dimension(km,nospZ)       :: ZunN     ! Unassimilated ingested Nitrogen (mmol-N/m3/d)
    real, dimension(km,nospZ)       :: ZunP     ! Unassimilated ingested Phosphorus (mmol-P/m3/d)
    real, dimension(km,nospZ)       :: ZunSi    ! Unassimilated ingested Silica (mmol-Si/m3/d)
    real, dimension(km)             :: ZunSi_tot
    real, dimension(km,nospZ)       :: Zmort    ! Dead zooplankton (individuals/m3/d)
    real :: ZmortC(km,nospZ), ZmortC_tot(km)  ! Carbon released from dead zooplankton (mmol-C/m3/d)
    real :: ZmortN(km,nospZ), ZmortN_tot(km)  ! Nitrogen released from dead zooplankton (mmol-N/m3/d)
    real :: ZmortP(km,nospZ), ZmortP_tot(km)  ! Phosphorus released from dead zooplankton (mmol-P/m3/d)
    real :: ZslopC(km,nospZ), ZslopC_tot(km)  ! Carbon lost to sloppy feeding (mmol-C/m3/d)
    real :: ZslopN(km,nospZ), ZslopN_tot(km)  ! Nitrogen lost to sloppy feeding (mmol-N/m3/d)
    real :: ZslopP(km,nospZ), ZslopP_tot(km)  ! Phosphorus lost to sloppy feeding (mmol-P/m3/d)
    real, dimension(km,nospZ)       :: ZexN   ! Excretion from zooplankton (mmol-N/m3/d)
    real, dimension(km)             :: ZexN_tot
    real, dimension(km,nospZ)       :: ZexP   ! Excretion from zooplankton (mmol-P/m3/d)
    real, dimension(km)             :: ZexP_tot
    real, dimension(km,nospZ)       :: ZegC   ! Egestion from zooplankton (mmol-C/m3/d)
    real, dimension(km,nospZ)       :: ZegN   ! Egestion from zooplankton (mmol-N/m3/d)
    real, dimension(km,nospZ)       :: ZegP   ! Egestion from zooplankton (mmol-P/m3/d)
    real, dimension(km,nospZ)       :: ZegSi  ! Egestion from zooplankton (mmol-Si/m3/d)
    real, dimension(km)          :: ZegSi_tot
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
    real, dimension(km,nospA+nospZ) :: Tadj ! Temperature adjustment factor
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
     open(unit=6201,file="cgem.csv")
     open(unit=6301,file="hydro.csv")
     open(unit=6401,file="light.csv")
     write(6001,'(A58)') "Agrow,Aresp,uA,uN,uP,uE,uSi,f_E,f_N,f_P,f_Si,Tadj,A,min_S"
    write(6101,'(A319)') "RO2,RNO3,RNH4,RPO4,RDIC,RSi,RALK,RN2,ROM1_A,ROM2_A,RO2_A,RNO3_A,RPO4_A,RDIC_A,RNH4_A,RSi_A,RALK_A,RN2_A,ROM1_Z,ROM2_Z,RO2_Z,RNO3_Z,RPO4_Z,RDIC_Z,RNH4_Z,RSi_Z,RALK_Z,RN2_Z,ROM1_R,ROM2_R,RO2_R,RNO3_R,RPO4_R,RDIC_R,RNH4_R,RSi_R,RALK_R,RN2_R,ROM1_BC,ROM2_BC,RO2_BC,RNO3_BC,RPO4_BC,RDIC_BC,RNH4_BC,RSi_BC,RALK_BC,RN2_BC,R_11"
    write(6201,'(A132)') "A1,Qn1,Qp1,Z1,Z2,NO3,NH4,PO4,DIC,O2,OM1A,OM2A,OM1Z,OM2Z,OM1R,OM2R,CDOM,Si,OM1BC,OM2BC,Alk,sx1A,sy1A,sx2A,sy2A,sx1Z,sy1Z,sx2Z,sy2Z"
    write(6301,'(A22)') "TC_8,rad,wind,sal,temp"
    write(6401,'(A45)') "ksurf,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kbot"
    endif


  do isp = 1,nf
    do k=1,km
      if(ff(k,isp).le.0.) then
        !write(6,*) "ff.le.0! set to 0 for istep,myrank,inea,k,isp=",istep,myrank,inea,k,isp,ff(k,isp)
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

       !! After Advection and VMixing, return to stoich's
        sx1A(:) = ff(:,isx1A) / OM1A(:)
        sy1A(:) = ff(:,isy1A) / OM1A(:)
        sx2A(:) = ff(:,isx2A) / OM2A(:)
        sy2A(:) = ff(:,isy2A) / OM2A(:)
        sx1Z(:) = ff(:,isx1Z) / OM1Z(:)
        sy1Z(:) = ff(:,isy1Z) / OM1Z(:)
        sx2Z(:) = ff(:,isx2Z) / OM2Z(:)
        sy2Z(:) = ff(:,isy2Z) / OM2Z(:)

!-----------------------------------------------------------------
!   Begin biogeochemistry calculations at time-level istep
!-----------------------------------------------------------------

!----------------------------------------------------------------
! Get chlorophyll-a quantity per layer
  Chla_tot = Fixed_CChla(A)

!---- Underwater Light Model ----------------------------------
! PARsurf is downward spectrally integrated irradiance just below the sea surface.  
! Rad is just above sea surface. 
  ! Rad was short wave generated by NRL, multiplied by SWtoPAR: ratio of PAR to
  ! shortwave radiation (hardcoded 4/30/14 to 0.43).
  ! Hardcoded to 0.47 on 2/11/16, Re: Tsubo and Walker, 2005
  ! PARfac is a multiplication factor for testing
  PARsurf = (0.47 * Rad) * PARfac
! PARbot is the downward irradiance (photons/cm2/sec) at the sea bottom
  if(PARsurf.le.0.) then
   PARbot = 0.0
   PARdepth = 0.0
  else
  call calc_PARdepth(TC_8,PARSurf,S,Chla_tot,CDOM,OM1A,OM1Z,OM1R,OM1BC,&
       &             PARdepth,PARbot )
  endif 

!---
!A.1-A1 Phytoplankton
!    dA_i/dt = Agrow_i - Aresp_i - ZgrazA_tot_i - Amort_i
!--

 !- calc_Agrow calculates Agrow_i, Aresp_i
 !    Agrow_i: production (cells/m3/s)
 !    Aresp_i: sum of somatic and basal respirtion (cells/m3/s)
 call calc_Agrow(PARdepth,T,Qn,Qp,Si,A,Agrow,uA,Aresp,uN,uP,uE,uSi,inea)

  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "Agrow,Aresp,uA,uN,uP,uSi"
    write(6,'(*(g0,:,", "))') Agrow(1,1),Aresp(1,1),uA(1,1),uN(1,1),uP(1,1),uSi(1,1)
  endif

 !- ZgrazA_tot_i: total zooplankton grazing on Ai by all zooplankton groups (cells/m3/d)
 do k=1,km

  do isp = 1, nospA
     Abiovol         = A(k,isp)*volcell(isp) 
     top_A(isp,:)    = AMAX1((Abiovol-Athresh(isp))*ediblevector(:,isp),0.0)
     bottom_A(isp,:) = Abiovol * ediblevector(:,isp)
  enddo
 
  do isz = 1, nospZ
     bottom(isz) = SUM(bottom_A(:,isz))   ! sum over isp
  enddo
 
  do isp = 1, nospA
     monodZ(k,isp,:)  = top_A(isp,:)/(ZKa(:) + bottom(:))
  enddo 

 enddo

  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "monodZ1:2,A,volcell,Athresh,top_A,bottom1:2,ediblevector(1:2,1)"
    write(6,'(*(g0,:,", "))') monodZ(1,1,:),A(1,1),volcell(1),Athresh(1),top_A(1,1),bottom(:),ediblevector(:,1)
  endif



!--------------------------------------
!-- Temperature adjustment factor 
   call func_T( T,Tadj )
!-- Nutrient dependent growth function
   call func_Qs( Qn, Qp, f_Qn, f_Qp)
  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "Tadj,f_Qn,f_Qp"
    write(6,'(*(g0,:,", "))') Tadj(1,1),f_Qn(1,1),f_Qp(1,1)
  endif


!! Sum over phytoplankton groups for PrimProd, ArespC, AexudN, AexudP
   do k = 1,km 
     PrimProd(k) = SUM(Agrow(k,:)*Qc(:))    ! Phytoplankton primary production (mmol-C/m3/d)
     ArespC(k)   = SUM(Aresp(k,:)*Qc(:))     ! Phytoplankton respiration     (mmol-C/m3/d)
     AexudN(k)   = SUM(Aresp(k,:)*Qn(k,:))   ! Total Phytoplankton exudation (mmol-N/m3/d)
     AexudP(k)   = SUM(Aresp(k,:)*Qp(k,:))   ! Total Phytoplankton exudation (mmol-P/m3/d)
   enddo

   do isp=1,nospA
   Amort(:,isp)  = A(:,isp) * mA(isp)                         ! Dead phytoplankton (cells/m3/d)
   enddo

  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "PrimProd,ArespC,AexudN,AexudP,Amort"
    write(6,'(*(g0,:,", "))') PrimProd(1),ArespC(1),AexudN(1),AexudP(1),Amort(1,1) 
  endif

!------------------------------------------------------------------------     
! Nutrient limited uptake:
! Find Rate Limiting Nutrient RLN for N, P, and Si:
   if(Rad.le.tiny(x)) then  !Nutrient uptake only takes place during the day
        vN = 0
        vP = 0
        vSi = 0
   else  !Day
      !Rate limiting nutrient is N
      !    Monod Equations
      ! Kx is half saturation coefficient for x
      do isp = 1, nospA
        monodN(:,isp)  = Ntotal(:)/(Ntotal(:)+Kn(isp))
        monodP(:,isp)  = PO4(:)/(PO4(:)+Kp(isp))
        monodSi(:,isp) = Si(:)/(Si(:)+Ksi(isp))
      enddo

      do k=1,km
        !Rate limiting nutrient is N
        if(Ntotal(k).le.PO4(k).and.Ntotal(k).le.Si(k)) then
          do isp = 1, nospA
           vN(k,isp) = Q10_T(T(k),vmaxN(isp))*monodN(k,isp)*f_Qn(k,isp)

           vP(k,isp) = Q10_T(T(k),vmaxP(isp))*monodP(k,isp)*f_Qp(k,isp)&
     &      *( Ntotal(k)/(Ntotal(k)+aN(isp)*Kn(isp)) )

           vSi(k,isp) = Q10_T(T(k),vmaxSi(isp))*monodSi(k,isp)        &
     &      *( Ntotal(k)/(Ntotal(k)+aN(isp)*Kn(isp)) )
          enddo

        !Rate limiting nutrient is P
        elseif(PO4(k).le.Ntotal(k).and.PO4(k).le.Si(k)) then
           do isp = 1, nospA
            vN(k,isp) = Q10_T(T(k),vmaxN(isp))*monodN(k,isp)*f_Qn(k,isp)&
     &         *( PO4(k)/(PO4(k)+aN(isp)*Kp(isp)) )

            vP(k,isp) = Q10_T(T(k),vmaxP(isp))*monodP(k,isp)*f_Qp(k,isp)

            vSi(k,isp) = Q10_T(T(k),vmaxSi(isp))*monodSi(k,isp)       &
     &         *( PO4(k)/(PO4(k)+aN(isp)*Kp(isp)) )
           enddo

        !Rate limiting nutrient is Si 
        else
           do isp = 1, nospA
            vN(k,isp) = Q10_T(T(k),vmaxN(isp))*monodN(k,isp)*f_Qn(k,isp)&
     &         *( Si(k)/(Si(k)+aN(isp)*Ksi(isp)) )

            vP(k,isp) = Q10_T(T(k),vmaxP(isp))*monodP(k,isp)*f_Qp(k,isp)&
     &         *( Si(k)/(Si(k)+aN(isp)*Ksi(isp)) )

            vSi(k,isp) = Q10_T(T(k),vmaxSi(isp))*monodSi(k,isp)
           enddo
        endif
      enddo !k

   endif !End Nutrient Uptake

  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "vN,vP,vSi"
    write(6,'(*(g0,:,", "))') vN(1,1),vP(1,1),vSi(1,1) 
  endif


!! ----------------------------------------------------------------------
!! Sum over phytoplankton groups for PrimProd, ArespC, AexudN, AexudP
   do k = 1, km 
     AupN(k) = SUM(A(k,:)*vN(k,:))     ! Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
     AupP(k) = SUM(A(k,:)*vP(k,:))     ! Phytoplankton uptake of Phosphorus (mmol-P/m3/d)
     AupSi(k) = SUM(A(k,:)*vSi(k,:))  ! Phytoplankton uptake of Silica (mmol-Si/m3/d)
   enddo
  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "AupN,AupP,AupSi"
    write(6,'(*(g0,:,", "))') AupN(1),AupP(1),AupSi(1)
  endif


  do k=1,km
    do isp=1,nospA
      ZgrazA(k,isp,:)   = (Z(k,:)*Zumax(:)*monodZ(k,isp,:))/volcell(isp)   ! Grazing of phytoplankton by plankton (cells/m3/d)
      ZgrazA_tot(k,isp) = SUM( ZgrazA(k,isp,:) )
    enddo
  enddo

    do isz = 1,nospZ
    do k=1,km
      ZgrazC(k,isz) = SUM(ZgrazA(k,:,isz) * Qc(:))     ! Carbon uptake from grazing (mmol-C/m3/day)
      ZgrazN(k,isz) = SUM(ZgrazA(k,:,isz) * Qn(k,:))  ! Nitrogen uptake from grazing( mmol-N/m3/day)
      ZgrazP(k,isz) = SUM(ZgrazA(k,:,isz) * Qp(k,:)) ! Phosphorus uptake from grazing (mmol-P/m3/day)
    enddo
   enddo

  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "ZgrazA,ZgrazA_tot,ZgrazC,ZgrazN,ZgrazP"
    write(6,'(*(g0,:,", "))') ZgrazA(1,1,1),ZgrazA_tot(1,1),ZgrazC(1,1),ZgrazN(1,1),ZgrazP(1,1) 
    write(6,'(*(g0,:,", "))') ZgrazA(1,1,2),ZgrazA_tot(1,1),ZgrazC(1,2),ZgrazN(1,2),ZgrazP(1,2)
  endif


  !---------------------------------------------------------
  !-A; Phytoplankton number density (cells/m3);
  !---------------------------------------------------------
!---
!A.1-A2 Qn: Phytoplankton Nitrogen Quota (mmol-N/cell) 
!    dQn_i/dt = vN_i - Qn_i*uA_i - AexudN_i/A_i 
!     uptake - utilization to support growth - exudation associated with respiration 


  do isp=1,nospA
      ff_new(:,iA(isp)) = A(:,isp)        &
      & + ( Agrow(:,isp) - Aresp(:,isp) - ZgrazA_tot(:,isp) - Amort(:,isp) )*dTd
  enddo

!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
!! Enforce minima for Droop, also enforce maxima if not equal Droop (which_quota=1)
  if(which_quota.eq.1) then
    do isp=1,nospA
      do k=1,km
       ff_new(k,iQn(isp)) = AMAX1(Qn(k,isp) + (vN(k,isp) - Qn(k,isp)*uA(k,isp))*dTd,QminN(isp))
      enddo
    enddo
  !! , also enforce maxima if not equal Droop (which_quota=1)
  else
    do isp=1,nospA
      do k=1,km
        ff_new(k,iQn(isp)) = AMIN1(AMAX1(Qn(k,isp) + (vN(k,isp) - Qn(k,isp)*uA(k,isp))*dTd,QminN(isp)),QmaxN(isp))
      enddo
    enddo
  endif

!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
!! Enforce minima for Droop, also enforce maxima if not equal Droop (which_quota=1)
  if(which_quota.eq.1) then
    do isp=1,nospA
      do k=1,km
        ff_new(k,iQp(isp)) = AMAX1(Qp(k,isp) + (vP(k,isp) - Qp(k,isp)*uA(k,isp))*dTd,QminP(isp))
      enddo
    enddo
!! , also enforce maxima if not equal Droop (which_quota=1)
  else
    do isp=1,nospA
      do k=1,km
          ff_new(k,iQp(isp)) = AMIN1(AMAX1(Qp(k,isp) + (vP(k,isp) - Qp(k,isp)*uA(k,isp))*dTd,QminP(isp)),QmaxP(isp))
      enddo
    enddo
  endif
!----------------------------------------------------------------------- 


!----------------------------------------------------------------------
! ZgrazC, ZgrazN, and ZgrazP are total carbon, nitrogen, and phosphorus
!   uptake of zooplankton from grazing all phytoplankton groups
!---------------------------------------------------------------------
do isz=1,nospZ
!-------------------------------------------------------------------
! Now calculate the total ingested ZinC, ZinN, and ZinP of C, N, and P
!-------------------------------------------------------------------
  ZslopC(:,isz)  = Zslop(isz)*ZgrazC(:,isz)                      ! Sloppy feeding (mmol-C/m3/d)
  ZunC(:,isz)    = (1.-Zeffic(isz))*(ZgrazC(:,isz)-ZslopC(:,isz))    ! Unassimilated (mmol-C/m3/d)
  ZinC(:,isz)    = ZgrazC(:,isz) - ZslopC(:,isz) - ZunC(:,isz)         ! Ingested (mmol-C/m3/d)

  ZslopN(:,isz)  = Zslop(isz)*ZgrazN(:,isz)                      ! Sloppy feeding (mmol-N/m3/d)
  ZunN(:,isz)    = (1.-Zeffic(isz))*(ZgrazN(:,isz)-ZslopN(:,isz)) ! Unassimilated (mmol-N/m3/d)
  ZinN(:,isz)    = ZgrazN(:,isz) - ZslopN(:,isz) - ZunN(:,isz)         ! Ingested (mmol-N/m3/d)

  ZslopP(:,isz)  = Zslop(isz)*ZgrazP(:,isz)                      ! Sloppy feeding (mmol-P/m3/d)
  ZunP(:,isz)    = (1.-Zeffic(isz))*(ZgrazP(:,isz)-ZslopP(:,isz)) ! Unassimilated (mmol-P/m3/d)
  ZinP(:,isz)    = ZgrazP(:,isz) - ZslopP(:,isz) - ZunP(:,isz)         ! Ingested (mmol-P/m3/d)
!-------------------------------------------------
enddo !isz


do k=1,km
!-------------------------------------------------------------------
! Now calculate the total ingested ZinC, ZinN, and ZinP of C, N, and P
!-------------------------------------------------------------------
  ZslopC_tot(k) = SUM(ZslopC(k,:))                             ! Total Sloppy feeding (mmol-C/m3/d)
  ZslopN_tot(k) = SUM(ZslopN(k,:))                             ! Total Sloppy feeding (mmol-N/m3/d) 
  ZslopP_tot(k) = SUM(ZslopP(k,:))                             ! Total Sloppy feeding (mmol-P/m3/d)
enddo !k

!------------------------------------         
! Liebigs Law for zooplankton group isz 
!------------------------------------
     do isz=1,nospZ
     do k=1,km
      if (ZinN(k,isz) .gt. optNP(isz)*ZinP(k,isz)) then  
        Zgrow(k,isz)= ZinP(k,isz)/ZQp(isz)                   ! P-limited growth (indv./m3/d) 
        ZegN(k,isz) = ZinN(k,isz) - ZinP(k,isz)*optNP(isz)       ! P-limited N excretion (mmol-N/m3/d) 
                                                   ! determined by subtracting N-equivalent of ZinP
        ZegC(k,isz) = ZinC(k,isz) - ZinP(k,isz)/ZQp(isz)*ZQc(isz)  ! P-limited C excretion (mmol-C/m3/d)
        ZegP(k,isz) = 0.                        
      else
        Zgrow(k,isz)= ZinN(k,isz)/ZQn(isz)                   ! N-limited growth (indv./m3/d)
        ZegP(k,isz) = ZinP(k,isz) - ZinN(k,isz)/optNP(isz)       ! N-limited P excretion (mmol-P/m3/d)    
                                                   ! determined by subtracting P-equivalent of ZinN
        ZegC(k,isz) = ZinC(k,isz) - ZinN(k,isz)/ZQn(isz)*ZQc(isz)  ! N-limited C excretion (mmol-C/m3/d)
        ZegN(k,isz) = 0.
    endif
    enddo !k
    enddo !isz

!------------------------------------------------

!-----------------------------------------------------
! ZegC should not be negative 
  do isz=1,nospZ
    do k=1,km
      if(ZegC(k,isz).lt.0.) then
          ZegC(k,isz) = 0.
          !ZegN(isz) = 0.
          !ZegP(isz) = 0.
          write(6,*) "ZegC,ZegN,ZegP",ZegC(k,isz),ZegN(k,isz),ZegP(k,isz)
      endif
    enddo 
  enddo

!
  do isz=1,nospZ
    ! Egestion and unassimilated for Si set equivalent to that of N
    ZegSi(:,isz) = ZegN(:,isz)
    ZunSi(:,isz) = ZunN(:,isz)
    Zresp(:,isz) = (Zgrow(:,isz)*Zrespg(isz) + Z(:,isz)*Zrespb(isz)) !Zooplankton respiration (indv./m3/d)
    ZexN(:,isz)    = Zresp(:,isz)*ZQn(isz)               ! (mmol-N/m3/d)
    ZexP(:,isz)    = Zresp(:,isz)*ZQp(isz)               ! (mmol-P/m3/d)
    Zmort(:,isz)     = Zm(isz) * Z(:,isz) * Z(:,isz)     ! (indv./m3/d)
    ZmortC(:,isz)      = Zmort(:,isz)*ZQc(isz)           ! (mmol-C/m3/d)
    ZmortN(:,isz)      = Zmort(:,isz)*ZQn(isz)           ! (mmol-N/m3/d)
    ZmortP(:,isz)      = Zmort(:,isz)*ZQp(isz)           ! (mmol-P/m3/d)
  enddo

  do k=1,km
    ZrespC(k) = SUM(Zresp(k,:)*ZQc(:))
    ZegSi_tot(k) = SUM(ZegSi(k,:))
    ZunSi_tot(k) = SUM(ZunSi(k,:))
    ZexN_tot(k) = SUM(ZexN(k,:))
    ZexP_tot(k) = SUM(ZexP(k,:))
    ZmortC_tot(k) = SUM(ZmortC(k,:))
    ZmortN_tot(k) = SUM(ZmortN(k,:))
    ZmortP_tot(k) = SUM(ZmortP(k,:))
  enddo



!-------------------------------------------------------------------------

!---------------------------------------------------------
!-Z; Zooplankton number density (individuals/m3);
!---------------------------------------------------------

  do isz=1,nospZ
    do k=1,km
      ff_new(k,iZ(:))  = AMAX1( Z(k,:)                         &
     &      + (Zgrow(k,:) - Zresp(k,:) - Zmort(k,:))*dTd, 1.)
    enddo
  enddo


#ifdef DEBUG
write(6,*) "In cgem, updated iZ, which_Fluxes(iInRemin),KG_bot=",which_Fluxes(iInRemin),KG_bot
#endif

!------------------------------------------------------------------------

do k=1,km
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
    & ROM1_BC(k),ROM2_BC(k),RO2_BC,RNO3_BC,RPO4_BC,RDIC_BC,RNH4_BC,RSi_BC,RALK_BC,RN2_BC,R_11
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
   OM1_CA(k) = OM1_CA(k) + Amort(k,isp)*(Qn(k,isp)-QminN(isp))/Qn(k,isp)*Qc(isp)
   OM1_NA(k) = OM1_NA(k) + Amort(k,isp)*(Qn(k,isp)-QminN(isp))
   OM1_PA(k) = OM1_PA(k) + Amort(k,isp)*(Qn(k,isp)-QminN(isp))/Qn(k,isp)*Qp(k,isp)
!!Dissolved
   OM2_CA(k) = OM2_CA(k) + Amort(k,isp)*QminN(isp)/Qn(k,isp)*Qc(isp)
   OM2_NA(k) = OM2_NA(k) + Amort(k,isp)*QminN(isp)
   OM2_PA(k) = OM2_PA(k) + Amort(k,isp)*QminN(isp)/Qn(k,isp)*Qp(k,isp)
 else
!Particulate
   OM1_CA(k) = OM1_CA(k) + Amort(k,isp)*(Qp(k,isp)-QminP(isp))/Qp(k,isp)*Qc(isp)
   OM1_NA(k) = OM1_NA(k) + Amort(k,isp)*(Qp(k,isp)-QminP(isp))/Qp(k,isp)*Qn(k,isp)
   OM1_PA(k) = OM1_PA(k) + Amort(k,isp)*(Qp(k,isp)-QminP(isp))
!!Dissolved
   OM2_CA(k) = OM2_CA(k) + Amort(k,isp)*QminP(isp)/Qp(k,isp)*Qc(isp)
   OM2_NA(k) = OM2_NA(k) + Amort(k,isp)*QminP(isp)/Qp(k,isp)*Qn(k,isp)
   OM2_PA(k) = OM2_PA(k) + Amort(k,isp)*QminP(isp)
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
    OM1_Ratio(k) = SUM( (Qn(k,:)-QminN(:))/Qn(k,:)*A(k,:))/SUM(A(k,:))
    OM2_Ratio(k) = SUM( (QminN(:)/Qn(k,:))*A(k,:))/SUM(A(k,:))

#ifdef DEBUG
write(6,*) "In cgem, OM1,2 Ratio=",OM1_Ratio,OM2_Ratio
#endif

    if(nospZ.eq.1) then 
                                                                  ! Particulate
     OM1_CZ(k)  = .5*(ZegC(k,1) + ZunC(k,1) + ZmortC_tot(k)) + OM1_Ratio(k)*ZslopC_tot(k) !  (mmol-C/m3/d)
     OM1_NZ(k)  = .5*(ZegN(k,1) + ZunN(k,1) + ZmortN_tot(k)) + OM1_Ratio(k)*ZslopN_tot(k) !  (mmol-N/m3/d)
     OM1_PZ(k)  = .5*(ZegP(k,1) + ZunP(k,1) + ZmortP_tot(k)) + OM1_Ratio(k)*ZslopP_tot(k) !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ(k)  = .5*(ZegC(k,1) + ZunC(k,1) + ZmortC_tot(k)) + OM2_Ratio(k)*ZslopC_tot(k)              !  (mmol-C/m3/d)
     OM2_NZ(k)  = .5*(ZegN(k,1) + ZunN(k,1) + ZmortN_tot(k)) + OM2_Ratio(k)*ZslopN_tot(k)              !  (mmol-N/m3/d)
     OM2_PZ(k)  = .5*(ZegP(k,1) + ZunP(k,1) + ZmortP_tot(k)) + OM2_Ratio(k)*ZslopP_tot(k)              !  (mmol-P/m3/d)

    else if(nospZ.eq.2) then
                                                                  ! Particulate
     OM1_CZ(k)  = ZegC(k,1) + ZunC(k,1) + ZmortC_tot(k) + OM1_Ratio(k)*ZslopC_tot(k) !  (mmol-C/m3/d)
     OM1_NZ(k)  = ZegN(k,1) + ZunN(k,1) + ZmortN_tot(k) + OM1_Ratio(k)*ZslopN_tot(k) !  (mmol-N/m3/d)
     OM1_PZ(k)  = ZegP(k,1) + ZunP(k,1) + ZmortP_tot(k) + OM1_Ratio(k)*ZslopP_tot(k) !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ(k)  = ZegC(k,2) + ZunC(k,2) + OM2_Ratio(k)*ZslopC_tot(k)              !  (mmol-C/m3/d)
     OM2_NZ(k)  = ZegN(k,2) + ZunN(k,2) + OM2_Ratio(k)*ZslopN_tot(k)              !  (mmol-N/m3/d)
     OM2_PZ(k)  = ZegP(k,2) + ZunP(k,2) + OM2_Ratio(k)*ZslopP_tot(k)              !  (mmol-P/m3/d)

    else 
                                                                  ! Particulate
     OM1_CZ(k)  = ZegC(k,1) + ZunC(k,1) + ZmortC_tot(k) + OM1_Ratio(k)*ZslopC_tot(k) !  (mmol-C/m3/d)
     OM1_NZ(k)  = ZegN(k,1) + ZunN(k,1) + ZmortN_tot(k) + OM1_Ratio(k)*ZslopN_tot(k) !  (mmol-N/m3/d)
     OM1_PZ(k)  = ZegP(k,1) + ZunP(k,1) + ZmortP_tot(k) + OM1_Ratio(k)*ZslopP_tot(k) !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ(k)  = SUM(ZegC(k,2:nospZ)) + SUM(ZunC(k,2:nospZ)) + OM2_Ratio(k)*ZslopC_tot(k)              !  (mmol-C/m3/d)
     OM2_NZ(k)  = SUM(ZegN(k,2:nospZ)) + SUM(ZunN(k,2:nospZ)) + OM2_Ratio(k)*ZslopN_tot(k)              !  (mmol-N/m3/d)
     OM2_PZ(k)  = SUM(ZegP(k,2:nospZ)) + SUM(ZunP(k,2:nospZ)) + OM2_Ratio(k)*ZslopP_tot(k)              !  (mmol-P/m3/d)

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
  & + ( RNH4(:) - AupN(:)*NH4(:)/(Ntotal(:)) + AexudN(:) + ZexN_tot(:)  )*dTd
!----------------------------
!-Silica: (mmol-Si/m3)
  ff_new(:,iSi) =  Si(:)                             &
  & + ( RSi(:) - AupSi(:) + ZegSi_tot(:) + ZunSi_tot(:) )*dTd
!---------------------------------------------
!-PO4: Phosphate (mmol-P/m3)
 ff_new(:,iPO4) = PO4(:)                             &
 & + ( RPO4(:) - AupP(:) + AexudP(:) + ZexP_tot(:)  )*dTd

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

  if(inea.eq.10.and.debug.eq.1) then
    write(6,*) "NO3,NO3new,RNO3,AupN,Ntotal" 
    write(6,'(*(g0,:,", "))') NO3(1),ff_new(1,iNO3),RNO3(1),AupN(1),Ntotal(1)
    write(6,*) "NH4,NH4new,RNH4,AupN,AexudN,ZexN_tot"
    write(6,'(*(g0,:,", "))') NH4(1),ff_new(1,iNH4),RNH4(1),AupN(1),AexudN(1),ZexN_tot(1)
    write(6,*) "Si,Sinew,RSi,AupSi,ZegSi_tot,ZunSi_tot"
    write(6,'(*(g0,:,", "))') Si(1),ff_new(1,iSi),RSi(1),AupSi(1),ZegSi_tot(1),ZunSi_tot(1)
    write(6,*) "PO4,PO4new,RPO4,AupP,AexudP,ZexP_tot"
    write(6,'(*(g0,:,", "))') PO4(1),ff_new(1,iPO4),RPO4(1),AupP(1),AexudP(1),ZexP_tot(1)
    write(6,*) "DIC,DICnew,RDIC,PrimProd,ArespC,ZrespC"
    write(6,'(*(g0,:,", "))') DIC(1),ff_new(1,iDIC),RDIC(1),PrimProd(1),ArespC(1),ZrespC(1)
    write(6,*) "O2,O2new,RO2"
    write(6,'(*(g0,:,", "))') O2(1),ff_new(1,iO2),RO2(1)
    write(6,*) "OM1A,OM1A_new,ROM1A,OM1_CA"
    write(6,'(*(g0,:,", "))') OM1A(1),ff_new(1,iOM1_A),ROM1_A(1),OM1_CA(1)
    write(6,*) "OM2A,OM2A_new,ROM2A,OM1_CA"
    write(6,'(*(g0,:,", "))') OM2A(1),ff_new(1,iOM2_A),ROM2_A(1),OM2_CA(1)
    write(6,*) "OM1Z,OM1Z_new,ROM1Z,OM1_CZ"
    write(6,'(*(g0,:,", "))') OM1Z(1),ff_new(1,iOM1_Z),ROM1_Z(1),OM1_CZ(1)
    write(6,*) "OM2Z,OM2Z_new,ROM2Z,OM2_CZ"
    write(6,'(*(g0,:,", "))') OM2Z(1),ff_new(1,iOM2_Z),ROM2_Z(1),OM2_CZ(1)
    write(6,*) "OM1R,OM1R_new,ROM1R"
    write(6,'(*(g0,:,", "))') OM1R(1),ff_new(1,iOM1_R),ROM1_R(1)
    write(6,*) "OM2R,OM2R_new,ROM2R"
    write(6,'(*(g0,:,", "))') OM2R(1),ff_new(1,iOM2_R),ROM2_R(1)
    write(6,*) "OM1BC,OM1BC_new,ROM1BC"
    write(6,'(*(g0,:,", "))') OM1BC(1),ff_new(1,iOM1_BC),ROM1_BC(1)
    write(6,*) "OM2BC,OM2BC_new,ROM2BC"
    write(6,'(*(g0,:,", "))') OM2BC(1),ff_new(1,iOM2_BC),ROM2_BC(1)
  endif


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
      !write(6,*) "ff_new.le.0! set to 0 for istep,inea,k,isp=",istep,inea,k,isp,ff_new(k,isp)
      ff_new(k,isp) = 0.0
    endif
  enddo
enddo

!-- End Main GEM Calculations ---------------------------------------------------

  !Write csv files
  k=1
    if(writecsv==1.and.k.eq.1.and.inea.eq.10) then
    write(6201,'(*(g0,:,", "))') ff_new(k,iA(1)),ff_new(k,iQn(1)),ff_new(k,iQp(1)),ff_new(k,iZ(1)),ff_new(k,iZ(2)),ff_new(k,iNO3),ff_new(k,iNH4),ff_new(k,iPO4),ff_new(k,iDIC),ff_new(k,iO2),ff_new(k,iOM1_A),ff_new(k,iOM2_A),ff_new(k,iOM1_Z),ff_new(k,iOM2_Z),ff_new(k,iOM1_R),ff_new(k,iOM2_R),ff_new(k,iCDOM),ff_new(k,iSi),ff_new(k,iOM1_BC),ff_new(k,iOM2_BC),ff_new(k,iAlk),ff_new(k,isx1A),ff_new(k,isy1A),ff_new(k,isx2A),ff_new(k,isy2A),ff_new(k,isx1Z),ff_new(k,isy1Z),ff_new(k,isx2Z),ff_new(k,isy2Z)
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
