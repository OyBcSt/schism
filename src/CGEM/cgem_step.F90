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
    integer        ::  i, inum, k, isp, isz ! Loop indicies, isp/isz is for phytoplankton/zooplankton species
    integer        ::  Is_Day            ! Switch for day/night for phytoplankton nutrient uptake only, Is_Day=0 means night
!------------------------------------ 
! Phytoplankton parameters
! Phytoplankton uptake and growth
    real, dimension(nospA,km) :: A_k      ! Phytoplankton number density (cells/m3)
    real    :: Agrow                       ! Phytoplankton growth (cells/m3/d)
    real, dimension(nospA,km) :: Agrow_k  ! Phytoplankton growth (cells/m3/d)
    real    :: uA                ! Specific growth rate (1/d)
    real    :: uA_k(km,nospA)    ! Specific growth rate (1/d)
    real    :: uN_k(km,nospA)    ! Nitrogen Limited growth rate (1/d)
    real    :: uP_k(km,nospA)    ! Phosphorus limited growth rate (1/d)
    real    :: uE_k(km,nospA)    ! Light limited growth rate (1/d)
    real    :: uSi_k(km,nospA)   ! Silica limited growth rate (1/d)
    real    :: f_Qn(nospA)        ! Quota model for N
    real    :: f_Qp(nospA)        ! Quota model for P
    real    :: Qn,Qn_k(nospA,km) ! Phytoplankton Nitrogen Quota (mmol-N/cell)
    real    :: Qp,Qp_k(nospA,km) ! Phytoplankton Phosphorus Quota (mmol-P/cell)
    real    :: vN    ! Phytoplankton uptake rate of Nitrogen (mmol-N/cell/d)
    real    :: vP    ! Phytoplankton uptake rate of Phosphorus (mmol-P/cell/d)
    real    :: vSi   ! Phytoplankton uptake rate of Silica (mmol-Si/cell/d)
    real    :: AupN  ! Total Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
    real    :: AupP  ! Total Phytoplankton uptake of Phosphorus (mmol-P/m3/d)
    real    :: AupSi ! Total Phytoplankton uptake of Silica (mmol-Si/m3/d)
    integer :: RLN   ! Rate Limiting Nutrient of N, P, and Si
 ! Monod equations for phytoplankton
    real, dimension(nospA)    :: monodN  !Monod term in nitrogen uptake
    real, dimension(nospA)    :: monodP  !Monod term in phosphorus uptake
    real, dimension(nospA)    :: monodSi !Monod term in Si uptake
    real    :: Ntotal   ! Total N (mmol/m3)
 ! Phytoplankton nutrient loss
    real, dimension(nospA)    :: Amort ! Dead phytoplankton (cells/m3/day)
    real, dimension(nospA)    :: AexudN_A    ! Phytoplankton group exudation (mmol-N/cell/d)
    real, dimension(nospA)    :: AexudP_A    ! Phytoplankton group exudation (mmol-P/cell/d) 
    real    :: AexudN          ! Sum of Exudation of N from all phytoplankton groups (mmol-N/m3/d)
    real    :: AexudP          ! Sum of Exudation of P from all phytoplankton groups (mmol-P/m3/d)
    real    :: Aresp           ! Total respiration from a phytoplankton group (cells/m3/d)
    real    :: Aresp_k(nospA,km) ! Total respiration from a phytoplankton group (cells/m3/d)
    real    :: ArespC          ! Phytoplankton equivalent carbon loss from respiration (mmol-C/m3/d)
!------------------------------------------------------------------
! Zooplankton parameters
 !Zooplankton uptake and growth
    real, dimension(nospZ,km)   :: Z_k      ! Zooplankton number density (indv./m3)
    real, dimension(nospZ)       :: optNP    ! Optimal nutrient ratio for zooplankton
    real, dimension(nospZ)       :: Z        ! Zooplankton
    real, dimension(nospZ)       :: Zgrow    ! Zooplankton growth (indv./m3/d)
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
    real, dimension(nospA,nospZ) :: monodZ   ! Monod term for zooplankton grazing
    real                         :: Abiovol  ! Algae biovolume vector (um3/m3)
    real, dimension(nospA,nospZ) :: top_A    ! Monod numerator value for phytoplankton group
    real, dimension(nospA,nospZ) :: bottom_A ! Monod Denominator value for phytoplankton group
    real, dimension(nospZ)       :: bottom   ! Sum of Monod Denominator value for all phytoplankton groups
 !Zooplankton nutrient loss
    real, dimension(nospZ)       :: Zresp    ! Zooplankton respiration (individuals/m3/d)
    real                         :: ZrespC   ! Carbon loss from zooplankton respiration (mmol-C/m3/day)
    real, dimension(nospZ)       :: ZunC     ! Unassimilated ingested Carbon (mmol-C/m3/d)
    real, dimension(nospZ)       :: ZunN     ! Unassimilated ingested Nitrogen (mmol-N/m3/d)
    real, dimension(nospZ)       :: ZunP     ! Unassimilated ingested Phosphorus (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZunSi    ! Unassimilated ingested Silica (mmol-Si/m3/d)
    real, dimension(nospZ)       :: Zmort    ! Dead zooplankton (individuals/m3/d)
    real :: ZmortC(nospZ), ZmortC_tot        ! Carbon released from dead zooplankton (mmol-C/m3/d)
    real :: ZmortN(nospZ), ZmortN_tot        ! Nitrogen released from dead zooplankton (mmol-N/m3/d)
    real :: ZmortP(nospZ), ZmortP_tot        ! Phosphorus released from dead zooplankton (mmol-P/m3/d)
    real :: ZslopC(nospZ), ZslopC_tot        ! Carbon lost to sloppy feeding (mmol-C/m3/d)
    real :: ZslopN(nospZ), ZslopN_tot        ! Nitrogen lost to sloppy feeding (mmol-N/m3/d)
    real :: ZslopP(nospZ), ZslopP_tot        ! Phosphorus lost to sloppy feeding (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZexN     ! Excretion from zooplankton (mmol-N/m3/d)
    real, dimension(nospZ)       :: ZexP     ! Excretion from zooplankton (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZegC     ! Egestion from zooplankton (mmol-C/m3/d)
    real, dimension(nospZ)       :: ZegN     ! Egestion from zooplankton (mmol-N/m3/d)
    real, dimension(nospZ)       :: ZegP     ! Egestion from zooplankton (mmol-P/m3/d)
    real, dimension(nospZ)       :: ZegSi    ! Egestion from zooplankton (mmol-Si/m3/d)
    real :: OM1_Ratio, OM2_Ratio             ! Separates sloppy feeding into OM1 and OM2
!---------------------------------------------------------------------- 
! Time variables  
    real, parameter :: one_d_365  = 1.0/365.0 ! Convert 1/yr to 1/day
!-----------------------------------------------------------------------
! Organic Matter Calculations
   ! Variables to calculate stoichiometry C:N:P ratios
    real    :: OM1_CA, OM1_NA, OM1_PA    ! OM from dead phytoplankton
    real    :: OM2_CA, OM2_NA, OM2_PA   
    real    :: OM1_CZ, OM1_NZ, OM1_PZ    ! OM from zooplankton
    real    :: OM2_CZ, OM2_NZ, OM2_PZ
    real,parameter    :: sz1 = 1. !for stoichiometry x,y,z, z is always normalized to 1 
!---------------------------------------------------------------------------
! reaction and Nitrification subroutine variables
    real    :: NO3, NH4, O2, DIC, Si, PO4               ! Nutrient input to subroutines
    real    :: OM1_A, OM1_Z, OM1_R, OM2_A, OM2_Z, OM2_R ! OM inputs to subroutines
    real    :: OM1_BC, OM2_BC                           ! OM inputs to subroutines
    real    :: R_11                                     ! Nitrification term
    real    :: RNO3, RNO3_A, RNO3_Z, RNO3_R, RNO3_BC ! Remineralization terms for NO3
    real    :: RNH4, RNH4_A, RNH4_Z, RNH4_R, RNH4_BC ! Remineralization terms for NH4
    real    :: ROM1_A, ROM1_Z, ROM1_R, ROM1_BC       ! Remineralization terms for POC
    real    :: ROM2_A, ROM2_Z, ROM2_R, ROM2_BC       ! Remineralization terms for DOC
    real    :: RO2, RO2_A, RO2_Z, RO2_R, RO2_BC      ! Remineralization terms for O2
    real    :: RPO4, RPO4_A, RPO4_Z, RPO4_R, RPO4_BC ! Remineralization terms for PO4
    real    :: RDIC, RDIC_A, RDIC_Z, RDIC_R, RDIC_BC ! Remineralization terms for DIC
    real    :: RSi, RSi_A, RSi_Z, RSi_R, RSi_BC      ! Remineralization terms for Si
    real    :: RALK, RALK_A, RALK_Z, RALK_R, RALK_BC ! Remineralization terms for ALK
    real    :: RN2, RN2_A, RN2_Z, RN2_R, RN2_BC ! Remineralization terms for N2 
    real, dimension(10) :: RC   ! Array that returns remineralization terms for OM
!---------------------------------------------------------
! Variables needed for light routine and calc_Agrow
    real    :: PARbot             ! Irradiance at sea floor (quanta/cm2/s)
    real    :: PARdepth_k(km)    ! Irradiance at center of layer k (quanta/cm2/s)
    real    :: PARsurf            ! Irradiance just below sea surface
    real, parameter :: RADCONV = 1./6.0221413*1.e-19 ! Convert quanta/cm2/s to mol/m2/s:
                                               !  = quanta/cm2/s * 1 mol/Avogadro# * 10,000cm2/m2
                                               !  = (1/6.022e23) * 1.0e4 = (1./6.022)e-23 * 1.0e4
                                               !  = (1./6.0221413)*1.e-19
    real, dimension(km) :: Chla_tot_k  ! Total amount of Chl-a in all the
                                        !  phytoplankton species (mg/m3) per cell
    real, dimension(nospA,km) :: Chl_C_k     ! Chl:C
    real, dimension(km) :: OM1A_k, OM1Z_k, OM1R_k, OM1BC_k !POC in g/m3
    real, dimension(km) :: CDOM_k    ! CDOM, ppb
!-----------------------------------------------------------------------
! Other variables 
    real :: PrimProd                     ! Primary production (photosynthesis)
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
     open(unit=6001,file="cgem_growth.csv")
     open(unit=6101,file="cgem_rates.csv")
     open(unit=6201,file="cgem.csv")
     open(unit=6301,file="hydro.csv")
     open(unit=6401,file="light.csv")
    write(6001,'(A58)') "Agrow_k,Aresp_k,uA_k,uN_k,uP_k,uE_k,uSi_k,f_E,f_N,f_P,f_Si"
    write(6101,'(A277)') "ROM1_A,ROM2_A,RO2_A,RNO3_A,RPO4_A,RDIC_A,RNH4_A,RSi_A,RALK_A,RN2_A,ROM1_Z,ROM2_Z,RO2_Z,RNO3_Z,RPO4_Z,RDIC_Z,RNH4_Z,RSi_Z,RALK_Z,RN2_Z,ROM1_R,ROM2_R,RO2_R,RNO3_R,RPO4_R,RDIC_R,RNH4_R,RSi_R,RALK_R,RN2_R,ROM1_BC,ROM2_BC,RO2_BC,RNO3_BC,RPO4_BC,RDIC_BC,RNH4_BC,RSi_BC,RALK_BC,RN2_BC"
    write(6201,'(A137)') "A,Qn,Qp,Z1,Z2,NO3,NH4,PO4,DIC,O2,OM1_A,OM2_A,OM1_Z,OM2_Z,OM1_R,OM2_R,CDOM,Si,OM1_BC,OM2_BC,Alk,Tr,sx1A,sy1A,sx2A,sy2A,sx1Z,sy1Z,sx2Z,sy2Z"
    write(6301,'(A22)') "TC_8,rad,wind,sal,temp"
    write(6401,'(A45)') "ksurf,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kbot"
    endif

!    if(istep.ne.1) then
       ! After Advection and VMixing, return to Q's
       do k=1,km
        ff(k,iQn(:)) = ff(k,iQn(:)) / ff(k,iA(:))
        ff(k,iQp(:)) = ff(k,iQp(:)) / ff(k,iA(:))
       enddo
!    endif

do k=1,km
  do isp = 1,nf
    if(ff(k,isp).le.0.) then 
      write(6,*) "ff.le.0! set to 0 for istep,myrank,inea,k,isp=",istep,myrank,inea,k,isp,ff(k,isp)
      ff(k,isp) = 0.0
    endif
  enddo
enddo

!write(6,*) "istep=",istep
!-----------------------------------------------------------------
!   Begin main ij loop for the biogeochemistry 
!   calculations at time-level istep
!-----------------------------------------------------------------
 !---------------------------------------------------------
 ! Calculate and convert variables needed for light routine
 !---------------------------------------------------------
      do k = 1, km

#ifdef DEBUG
write(6,*) "init loop, ZQp, nea, km",ZQp,k
#endif


   !Get algae counts and Nitrogen/phosphorus quotas
             do isp = 1, nospA          
               A_k(isp,k) = ff(k,isp) ! Phytoplankton in group isp, cells/m3
               Qn_k(isp,k) = ff(k,iQn(isp))
               Qp_k(isp,k) = ff(k,iQp(isp))
             enddo 
         !Save Zooplanton to k array
             do isp = 1,nospZ
                Z_k(isp,k) = ff(k,iZ(isp)) ! Zooplankton in group isp, ind./m3
             enddo
#ifdef DEBUG
!write(6,*) "after A,QnQp,Z.  k,Z_k(1,k)=",k,Z_k(1,k)
#endif


#ifdef DEBUG
!write(6,*) "T,S=",k,T(k),S(k)
#endif

         ! CDOM is in ppb
           CDOM_k(k)  = ff(k,iCDOM)
         ! Below is mmol/m3 Organic Matter from rivers converted to equivalent g carbon/m3
           OM1R_k(k) = ff(k,iOM1_R)
         ! Below is mmol/m3 Organic Matter from fecal pellets converted to equivalent g carbon/m3
           OM1Z_k(k) = ff(k,iOM1_Z)
         ! Below is mmol/m3 Organic Matter from dead phytoplankton converted to equivalent g carbon/m3
           OM1A_k(k)  = ff(k,iOM1_A)
         ! Below is mmol/m3 Organic Matter from initial and boundary conditions converted to equivalent g carbon/m3
           OM1BC_k(k)  = ff(k,iOM1_BC)
      enddo ! End of the "DO k = 1, km" block DO loop


!----------------------------------------------------------------
!
! Get chlorophyll-a quantity per layer
! 
!----------------------------------------------------------------
!      select case (which_chlaC)
!      case (1) ! Use fixed C:Chla 
          Chla_tot_k = Fixed_CChla(A_k)
!      case default
!         WRITE(6, "(/'The inputfile specified invalid setting: which_chlaC = ', I2/)") which_chlaC
!         WRITE(6, "('which_chlaC determines the method for calculating the quantity of chlorophyll-a.')")
!         WRITE(6, "('Cloern was removed, and fixed is the only option.'/)")
!         WRITE(6, "('Using fixed CChla.'/)")
!      end select ! which_chlaC


!----------------------------------------------------------------------
! Execute the desired atmospheric light model.  To calculate PARsurf,
! the effect amount of downward spectrally integrated irradiance 
! just below the sea surface.  'Rad' is just above sea surface. 
!----------------------------------------------------------------------

 !--Begin Calculate atmospheric model --------------------------------
         ! Rad(i) is short wave generated by NRL is used,
         ! and is multiplied by SWtoPAR: ratio of PAR to
         ! shortwave radiation (hardcoded 4/30/14 to 0.43).
         ! Hardcoded to 0.47 on 2/11/16, Re: Tsubo and Walker, 2005
         ! PARfac is a multiplication factor for testing
                    PARsurf = (0.47 * Rad) * PARfac
 !--End Calculate atmospheric model ---------------------------------------------

!----------------------------------------------------------------------------
! Execute the desired underwater light model to calculate the 1-D radiation
! arrays PARdepth_k, Esed and PAR_percent_k radiation arrays for
! vertical grid column (i).
!
! PARdepth_k(k) is the downward irradiance (photons/cm2/sec) at the middle
!                    of cell(k).
!
! PARbot is the downward irradiance (photons/cm2/sec) at the sea bottom
!
! PAR_percent_k(k)    is the % of incoming irradiance PARsurf that PARdepth_k(k)
!                 represents. PARsurf is the downward irradiance
!                 (photons/cm2/sec) just below the sea surface.
!-------------------------------------------------------------------------
        if(PARsurf.le.0.) then
         PARbot = 0.0
         PARdepth_k = 0.0
        else

        call calc_PARdepth( TC_8,PARSurf,S,Chla_tot_k,ff(k,iCDOM),OM1A_k,OM1Z_k,      &
   &                        OM1R_k,OM1BC_k,PARdepth_k,PARbot )

        endif !Endif for if it has PAR

#ifdef DEBUG_PAR
       if(myrank==0) write(6,*) "PARdepth_k",PARdepth_k
#endif

!---------------------End Underwater Light Model-----------------------------------

!-------------------------------------------------------------------------
! call subroutine calc_Agrow to execute the desired phytoplankton 
! growth model to calculate the one-D array Agrow_k.
!
! Agrow_k(k) is the growth-rate  for
! vertical grid column (i) at cell (k).
!-----------------------------------------------------------------------

       if(km .gt. 0) call calc_Agrow( PARdepth_k, T,              &
                     & Qn_k, Qp_k, ff(:,iSi) , A_k, Agrow_k,   &
                     & uA_k, Aresp_k, uN_k, uP_k, uE_k, uSi_k, inea )

!------end phytoplankton growth model-------------------------

#ifdef DEBUG
write(6,*) "In cgem, finished calc_Agrow"
#endif


!------------------------------------------
! Begin main k loop within the main ij loop
!------------------------------------------
      do k = 1, km      
!---------------------------------------------------------------------	   
  !Initialize variables
!  Zooplankton groups
   Z(:)         = ff(k,iZ(:))

#ifdef DEBUG
write(6,*) "In cgem, ediblevector,volcell=",ediblevector,volcell
#endif

          do isp = 1, nospA
             Abiovol       = A_k(isp,k)*volcell(isp) 
             top_A(isp,:)    = AMAX1((Abiovol-Athresh(isp))*ediblevector(:,isp),0.0)
             bottom_A(isp,:) = Abiovol * ediblevector(:,isp)
          enddo

         do isz = 1, nospZ
          bottom(isz) = SUM(bottom_A(:,isz))   ! sum over isp
         enddo

         do isp = 1, nospA
           monodZ(isp,:)  = top_A(isp,:)/(ZKa(:) + bottom(:))
         enddo 

#ifdef DEBUG
write(6,*) "In cgem, initialized top, bottom, monod" 
#endif


!--------------------------------------------------------------------------
! Initialize counters to zero that are used to accumulate variable values
! over the nospA phytoplankton groups and the nospZ zooplankton groups.
!--------------------------------------------------------------------------
      PrimProd  = 0.0
      ArespC    = 0.0
      AexudN    = 0.0
      AexudP    = 0.0 
      AupN      = 0.0
      AupP      = 0.0
      AupSi     = 0.0
      ZgrazC    = 0.0
      ZgrazN    = 0.0
      ZgrazP    = 0.0

!---------------------------------------------------------------
! Set the scalar variables that are common to all three reaction
! subroutines
!----------------------------------------------------------------
        O2        = ff(k,iO2)
        NO3       = ff(k,iNO3)
        NH4       = ff(k,iNH4)
        Si        = ff(k,iSi)
        PO4       = ff(k,iPO4)
        DIC       = ff(k,iDIC)
        OM1_A     = ff(k,iOM1_A)
        OM2_A     = ff(k,iOM2_A)
        OM1_Z    = ff(k,iOM1_Z)
        OM2_Z    = ff(k,iOM2_Z)
        OM1_R    = ff(k,iOM1_R)
        OM2_R    = ff(k,iOM2_R)
#ifdef DEBUG_R
write(6,*) "k,OM1_R",k,OM1_R
write(6,*) "k,OM2_R",k,OM2_R
#endif
        OM1_BC    = ff(k,iOM1_BC)
        OM2_BC    = ff(k,iOM2_BC)

#ifdef DEBUG
write(6,*) "In cgem, initialized state vars" 
#endif


!--------------------------------------
! Call temperature and growth functions
!-----------------------------------------
      call func_T( T(k),Tadj )
!     Nutrient dependent growth function
      call func_Qs( Qn_k(:,k), Qp_k(:,k), f_Qn, f_Qp) 


!Nutrients only taken up during the day:
     Is_Day = 1
     if(Rad.le.tiny(x)) Is_Day = 0

#ifdef DEBUG
write(6,*) "In cgem, called T,Qx functions" 
#endif

! ----------------------------------------------------------------------
      do isp = 1, nospA      
! ----------------------------------------------------------------------   
       Aresp   = Aresp_k(isp,k)    
       uA      = uA_k(k,isp)    
       Agrow   = Agrow_k(isp,k) 
       Qn      = Qn_k(isp,k) 
       Qp      = Qp_k(isp,k)      

!---------------------------------------------------------------------      
! Note that the expressions for PrimProd, ArespC, AexudN, AexudP are 
! sums over the isp phytoplankton groups. When the isp loop is complete,
! PrimProd, ArespC. AexudN, and AexudP will represent totals for
! all the nospA phytoplankton groups. 
!--------------------------------------------------------------------
      PrimProd       = PrimProd + Agrow*Qc(isp)    ! Phytoplankton 
                                                   ! primary production
                                                   ! (mmol-C/m3/d)
      ArespC  = ArespC + Aresp*Qc(isp)             ! Phytoplankton respiration 
						   ! equivalent carbon loss
						   ! (mmol-C/m3/d)	
      
      AexudN_A(isp) = Aresp*Qn/A_k(isp,k)   ! Phytoplankton group exudation (mmol-N/cell/d)
      AexudP_A(isp) = Aresp*Qp/A_k(isp,k)   ! Phytoplankton group exudation (mmol-P/cell/d)
      AexudN = AexudN + Aresp*Qn            ! Total Phytoplankton exudation (mmol-N/m3/d)
      AexudP = AexudP + Aresp*Qp            ! Total Phytoplankton exudation (mmol-P/m3/d)
      
!-------------------------------------------------------------------------
! Calculate dead phytoplankton, particulate and dissolved  
!-------------------------------------------------------------------------      
     Amort(isp)     = A_k(isp,k) * mA(isp)    ! dead phytoplankton (cells/m3/day)

!    Monod Equations
     Ntotal       = NO3 + NH4 
     monodN(isp)  = Ntotal/(Ntotal+Kn(isp))
     monodP(isp)  = PO4/(PO4+Kp(isp))
     monodSi(isp) = Si/(Si+Ksi(isp))

#ifdef DEBUG
write(6,*) "In cgem, set Amort, Ntotal monod" 
#endif

!------------------------------------------------------------------------     
! Nutrient limited uptake:
! Find Rate Limiting Nutrient RLN for N, P, and Si:
! N==1, P==2, Si==3
   if(Ntotal.le.PO4.and.Ntotal.le.Si) then
     RLN = 1
   elseif (PO4.le.Ntotal.and.PO4.le.Si) then
     RLN = 2
   else
     RLN = 3
   endif

!------------------------------------------------------------------------
   if(Is_Day.eq.0) then  !Nutrient uptake only takes place during the day
        vN = 0
        vP = 0
        vSi = 0

   else

#ifdef DEBUG
write(6,*) "It Is_Day, calculating vNutrients"
#endif


      if(RLN.eq.1) then
         vN = Q10_T(T(k),vmaxN(isp))*monodN(isp)*f_Qn(isp)

         vP = Q10_T(T(k),vmaxP(isp))*monodP(isp)*f_Qp(isp)&
     &      *( Ntotal/(Ntotal+aN(isp)*Kn(isp)) )

         vSi = Q10_T(T(k),vmaxSi(isp))*monodSi(isp)        &
     &      *( Ntotal/(Ntotal+aN(isp)*Kn(isp)) )

      elseif(RLN.eq.2) then
         vN = Q10_T(T(k),vmaxN(isp))*monodN(isp)*f_Qn(isp)&
     &      *( PO4/(PO4+aN(isp)*Kp(isp)) )

         vP = Q10_T(T(k),vmaxP(isp))*monodP(isp)*f_Qp(isp)

         vSi = Q10_T(T(k),vmaxSi(isp))*monodSi(isp)       &
     &      *( PO4/(PO4+aN(isp)*Kp(isp)) )

      else
         vN = Q10_T(T(k),vmaxN(isp))*monodN(isp)*f_Qn(isp)&
     &      *( Si/(Si+aN(isp)*Ksi(isp)) )

         vP = Q10_T(T(k),vmaxP(isp))*monodP(isp)*f_Qp(isp)&
     &      *( Si/(Si+aN(isp)*Ksi(isp)) )

         vSi = Q10_T(T(k),vmaxSi(isp))*monodSi(isp)

      endif

   endif !Endif Is_Day.eq.0

     
#ifdef DEBUG
write(6,*) "In cgem, finished vNutrients"
#endif
 
!--------------------------------------------------------------      
! When isp loop is done, AupN and AupP are totals for all nospA 
! phytoplankton groups in cell k          
!---------------------------------------------------------------   
      AupN = AupN + A_k(isp,k)*vN     ! Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
      AupP = AupP + A_k(isp,k)*vP     ! Phytoplankton uptake of Phosphorus (mmol-P/m3/d) 
      AupSi = AupSi + A_k(isp,k)*vSi  ! Phytoplankton uptake of Silica (mmol-Si/m3/d)
 
!-----------------------------------------------------------------------
! Note that Zumax(1)*monodZ(isp,1) is volume of type isp phytoplankton
! eaten per day by type 1 zooplankton. Therefore
!      Zumax(1)*monodZ(isp,1)/volcell(isp)
! is the number of type isp phytoplankton eaten per day by type 1 
! zooplankton. An analogous statement holds for type 2 zooplankton
!-----------------------------------------------------------------------
      Zgrazvol(isp,:)     = Z(:)*Zumax(:)*monodZ(isp,:)   ! Grazing of phytoplankton by biovolume (um3/m3/d)
      ZgrazA(isp,:)       = Zgrazvol(isp,:)/volcell(isp)  ! Grazing of phytoplankton (cells/m3/d)
      ZgrazA_tot(isp) = SUM( ZgrazA(isp,:) ) 

!----------------------------------------------------------------------
! When the isp loop is finished, ZgrazC, ZgrazN, and ZgrazP will be total
! carbon, nitrogen, and phosphorous uptake of zooplankton from grazing
! all phytoplankton groups
!---------------------------------------------------------------------
      ZgrazC(:) = ZgrazC(:) + ZgrazA(isp,:) * Qc(isp)     ! Carbon uptake from grazing (mmol-C/m3/day)
      ZgrazN(:) = ZgrazN(:) + ZgrazA(isp,:) * Qn          ! Nitrogen uptake from grazing( mmol-N/m3/day)
      ZgrazP(:) = ZgrazP(:) + ZgrazA(isp,:) * Qp          ! Phosphorus uptake from grazing (mmol-P/m3/day)

#ifdef DEBUG
write(6,*) "In cgem, finished uptake and grazing,k,isp=",k,isp
write(6,*) "Agrow, Aresp, ZgrazA_tot(isp),Amort(isp)",Agrow, Aresp, ZgrazA_tot(isp), Amort(isp)
#endif
!---------------------------------------------------------
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------
      ff_new(k,iA(isp)) = AMAX1(ff(k,iA(isp))        &
      & + ( Agrow - Aresp - ZgrazA_tot(isp) - Amort(isp) )*dTd, 1.)
    

#ifdef DEBUG
write(6,*) "In cgem, updated iA"
#endif
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
      Qn = ff(k,iQn(isp)) + (vN - Qn*uA)*dTd
!      
!! Enforce minima, also enforce maxima if not equal Droop (which_quota=1)
      if(which_quota.eq.1) then
           Qn = AMAX1(Qn,QminN(isp))
      else
           Qn = AMIN1(AMAX1(Qn,QminN(isp)),QmaxN(isp))
      endif

      ff_new(k,iQn(isp)) = Qn

#ifdef DEBUG
write(6,*) "In cgem, updated iQn"
#endif

!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
      Qp = ff(k,iQp(isp)) + (vP - Qp*uA)*dTd
!
!! Enforce minima, also enforce maxima if not equal Droop (which_quota=1)
      if(which_quota.eq.1) then
           Qp = AMAX1(Qp,QminP(isp))
      else
           Qp = AMIN1(AMAX1(Qp,QminP(isp)),QmaxP(isp))
      endif

     ff_new(k,iQp(isp)) = Qp      

#ifdef DEBUG
write(6,*) "In cgem, updated iQn"
#endif

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
        Zgrow(isz)= ZinP(isz)/ZQp(isz)                   ! P-limited growth (indv./m3/d) 
        ZegN(isz) = ZinN(isz) - ZinP(isz)*optNP(isz)       ! P-limited N excretion (mmol-N/m3/d) 
                                                   ! determined by subtracting N-equivalent of ZinP
        ZegC(isz) = ZinC(isz) - ZinP(isz)/ZQp(isz)*ZQc(isz)  ! P-limited C excretion (mmol-C/m3/d)
        ZegP(isz) = 0.                        
      else
        Zgrow(isz)= ZinN(isz)/ZQn(isz)                   ! N-limited growth (indv./m3/d)
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


! Zooplankton respiration based on growth and basal metabolism, both modified by a temperature adjustment factor 
      Zresp(:) = (Zgrow(:)*Zrespg(:) + Z(:)*Zrespb(:)) !Zooplankton respiration (indv./m3/d)

      ZrespC   = SUM(Zresp*ZQc)                                       !Total Carbon loss from respiration (mmol-C/m3/d)

                                                ! Excretion
      ZexN(:)   = Zresp(:)*ZQn(:)               ! (mmol-N/m3/d)
      ZexP(:)   = Zresp(:)*ZQp(:)               ! (mmol-P/m3/d)

                                                ! Mortality
     Zmort(:)       = Zm(:) * Z(:) * Z(:)       ! (indv./m3/d)
     ZmortC(:)      = Zmort(:)*ZQc(:)           ! (mmol-C/m3/d)
     ZmortC_tot     = SUM(ZmortC)
     ZmortN(:)      = Zmort(:)*ZQn(:)           ! (mmol-N/m3/d)
     ZmortN_tot     = SUM(ZmortN)
     ZmortP(:)      = Zmort(:)*ZQp(:)           ! (mmol-P/m3/d)
     ZmortP_tot     = SUM(ZmortP)
!-------------------------------------------------------------------------

!---------------------------------------------------------
!-G; Zooplankton number density (individuals/m3);
!---------------------------------------------------------
      ff_new(k,iZ(:))  = AMAX1( ff(k,iZ(:))                         &
      &      + (Zgrow(:) - Zresp(:) - Zmort(:))*dTd, 1.)

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
        call Nitrification( O2, NH4, KO2, KNH4, nitmax, T(k), R_11 )

#ifdef DEBUG
write(6,*) "In cgem, called Nitrification" 
#endif

!------------------------------------------------------------
! Carbon Chemistry
!--------------------------------------------------------------
!!! MOCSY alkalinity expressions:
        m_alk = ff(k,iALK)/1000.
        m_dic = ff(k,iDIC)/1000.
        m_si  = ff(k,iSi)/1000.
        m_po4 = ff(k,iPO4)/1000.
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
        call reaction( OM1_A, OM2_A, O2, NO3, KG1, KG2, KO2, KstarO2, KNO3,               &
     &  ff(k,isx1A), ff(k,isy1A), sz1, ff(k,isx2A), ff(k,isy2A), sz1, T(k), RC )

        RC        = one_d_365 * RC  !Change units from /year to /day

        ROM1_A     = RC(1)          ! units are /m3/day
        ROM2_A     = RC(2)
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
        call reaction( OM1_Z, OM2_Z, O2, NO3, KG1, KG2, KO2, KstarO2, KNO3,               &
     &  ff(k,isx1Z), ff(k,isy1Z), sz1, ff(k,isx2Z), ff(k,isy2Z), sz1, T(k), RC )
        RC       = one_d_365 * RC   !Change units from /year to /day

        ROM1_Z     = RC(1)         ! units are /m3/day
        ROM2_Z     = RC(2)
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
        call reaction( OM1_R, OM2_R, O2, NO3, KG1_R, KG2_R, KO2, KstarO2, KNO3,               &
     &  sx1R, sy1R, sz1, sx2R, sy2R, sz1, T(k), RC )

        RC       = one_d_365 * RC   !Change units from /year to /day

        ROM1_R     = RC(1)         ! units are /m3/day
        ROM2_R     = RC(2)
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
        call reaction( OM1_BC, OM2_BC, O2, NO3, KG1_BC, KG2_BC, KO2, KstarO2, KNO3,               &
     &  sx1BC, sy1BC, sz1, sx2BC, sy2BC, sz1, T(k), RC )

        RC       = one_d_365 * RC   !Change units from /year to /day

        ROM1_BC     = RC(1)         ! units are /m3/day
        ROM2_BC     = RC(2)
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
  RO2   = RO2_A  + RO2_Z  + RO2_R  + RO2_BC - 2.*R_11  ! (mmol-O2/m3/d)
  RNO3  = RNO3_A + RNO3_Z + RNO3_R + RNO3_BC + R_11    ! (mmol-NO3/m3/d)
  RNH4  = RNH4_A + RNH4_Z + RNH4_R + RNH4_BC - R_11    ! (mmol-NH4/m3/d)
  RPO4  = RPO4_A + RPO4_Z + RPO4_R + RPO4_BC           ! (mmol-PO4/m3/d)
  RDIC  = RDIC_A + RDIC_Z + RDIC_R + RDIC_BC           ! (mmol-DIC/m3/d)
  RSi   = RSi_A  + RSi_Z  + RSi_R  + RSi_BC            ! (mmol-Si/m3/d)
  RALK  = RALK_A + RALK_Z + RALK_R + RALK_BC - 2.*R_11 ! (mmol-HCO3/m3/d)
  RN2   = RN2_A + RN2_Z + RN2_R + RN2_BC         ! (mmol-N2/m3/d)
!--------------------------------------------------------------------
if(writecsv==1.and.k.eq.1.and.inea.eq.10) then
    write(6101,'(*(g0,:,", "))') &
    & ROM1_A, ROM2_A, RO2_A, RNO3_A, RPO4_A, RDIC_A, RNH4_A, RSi_A, RALK_A, RN2_A, &
    & ROM1_Z, ROM2_Z, RO2_Z, RNO3_Z, RPO4_Z, RDIC_Z, RNH4_Z, RSi_Z, RALK_Z, RN2_Z, &
    & ROM1_R, ROM2_R, RO2_R, RNO3_R, RPO4_R, RDIC_R, RNH4_R, RSi_R, RALK_R, RN2_R, &
    & ROM1_BC,ROM2_BC,RO2_BC,RNO3_BC,RPO4_BC,RDIC_BC,RNH4_BC,RSi_BC,RALK_BC,RN2_BC
endif


!
!---------------------------------------------------------------------
! Stoichiometry - calculate C:N:P ratios for Remineralization equations
!---------------------------------------------------------------------
!-- Organic Matter from dead phytoplankton --------------------------
OM1_CA = 0.
OM2_CA = 0.
OM1_NA = 0.
OM2_NA = 0.
OM1_PA = 0.
OM2_PA = 0.

do isp=1,nospA
 if ( uN_k(k,isp) .lt. uP_k(k,isp)  ) then
!Particulate
   OM1_CA = OM1_CA + Amort(isp)*(Qn_k(isp,k)-QminN(isp))/Qn_k(isp,k)*Qc(isp)
   OM1_NA = OM1_NA + Amort(isp)*(Qn_k(isp,k)-QminN(isp))
   OM1_PA = OM1_PA + Amort(isp)*(Qn_k(isp,k)-QminN(isp))/Qn_k(isp,k)*Qp_k(isp,k)
!!Dissolved
   OM2_CA = OM2_CA + Amort(isp)*QminN(isp)/Qn_k(isp,k)*Qc(isp)
   OM2_NA = OM2_NA + Amort(isp)*QminN(isp)
   OM2_PA = OM2_PA + Amort(isp)*QminN(isp)/Qn_k(isp,k)*Qp_k(isp,k)
 else
!Particulate
   OM1_CA = OM1_CA + Amort(isp)*(Qp_k(isp,k)-QminP(isp))/Qp_k(isp,k)*Qc(isp)
   OM1_NA = OM1_NA + Amort(isp)*(Qp_k(isp,k)-QminP(isp))/Qp_k(isp,k)*Qn_k(isp,k)
   OM1_PA = OM1_PA + Amort(isp)*(Qp_k(isp,k)-QminP(isp))
!!Dissolved
   OM2_CA = OM2_CA + Amort(isp)*QminP(isp)/Qp_k(isp,k)*Qc(isp)
   OM2_NA = OM2_NA + Amort(isp)*QminP(isp)/Qp_k(isp,k)*Qn_k(isp,k)
   OM2_PA = OM2_PA + Amort(isp)*QminP(isp)
 endif
enddo

!write(6,*) "it,myi,k,OM12CNP,dTd",istep,myrank,k,OM1_CA,OM1_NA,OM1_PA,OM2_CA,OM2_NA,OM2_PA,dTd

#ifdef DEBUG
write(6,*) "In cgem, PD OMs"
#endif

                                             ! Dissolved

 !write(6,*) stoich_x1A
   !This calculates the cumulative stoichiometry ratios for OM1_A
   if(OM1_CA.gt.tiny(x)) then
!   if(OM1_CA.ne.0) then
    ff_new(k,isx1A) = (OM1_CA*dTd + OM1_A) / (OM1_PA*dTd + 1./ff(k,isx1A)*OM1_A) ! C/P
    ff_new(k,isy1A) = (OM1_NA*dTd + ff(k,isy1A)/ff(k,isx1A)*OM1_A) / (OM1_PA*dTd + 1./ff(k,isx1A)*OM1_A) !N/P
    !stoich_z1A = 1., using sz1=1
   else
    ff_new(k,isx1A) = ff(k,isx1A)
    ff_new(k,isy1A) = ff(k,isy1A)
   endif

#ifdef DEBUG
write(6,*) "In cgem, finished stoich OM1A"
#endif

   !This calculates the cumulative stoichiometry ratios for OM2_A
   if(OM2_CA.gt.tiny(x)) then
!    if(OM2_CA.ne.0) then
    ff_new(k,isx2A) = (OM2_CA*dTd + OM2_A) / (OM2_PA*dTd + 1./ff(k,isx2A)*OM2_A) ! C/P
    ff_new(k,isy2A) = (OM2_NA*dTd + ff(k,isy2A)/ff(k,isx2A)*OM2_A) / (OM2_PA*dTd + 1./ff(k,isx2A)*OM2_A) !N/P
   else
    ff_new(k,isx2A) = ff(k,isx2A)
    ff_new(k,isy2A) = ff(k,isy2A)
   endif

#ifdef DEBUG
write(6,*) "In cgem, finished stoich OM2A"
#endif


!!-- Organic Matter from fecal pellets ---------------------------------
    OM1_Ratio = SUM( (Qn_k(:,k)-QminN)/Qn_k(:,k)*A_k(:,k))/SUM(A_k(:,k))
    OM2_Ratio = SUM( (QminN/Qn_k(:,k))*A_k(:,k))/SUM(A_k(:,k))

#ifdef DEBUG
write(6,*) "In cgem, OM1,2 Ratio=",OM1_Ratio,OM2_Ratio
#endif

    if(nospZ.eq.1) then 
                                                                  ! Particulate
     OM1_CZ  = .5*(ZegC(1) + ZunC(1) + ZmortC_tot) + OM1_Ratio*ZslopC_tot !  (mmol-C/m3/d)
     OM1_NZ  = .5*(ZegN(1) + ZunN(1) + ZmortN_tot) + OM1_Ratio*ZslopN_tot !  (mmol-N/m3/d)
     OM1_PZ  = .5*(ZegP(1) + ZunP(1) + ZmortP_tot) + OM1_Ratio*ZslopP_tot !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ  = .5*(ZegC(1) + ZunC(1) + ZmortC_tot) + OM2_Ratio*ZslopC_tot              !  (mmol-C/m3/d)
     OM2_NZ  = .5*(ZegN(1) + ZunN(1) + ZmortN_tot) + OM2_Ratio*ZslopN_tot              !  (mmol-N/m3/d)
     OM2_PZ  = .5*(ZegP(1) + ZunP(1) + ZmortP_tot) + OM2_Ratio*ZslopP_tot              !  (mmol-P/m3/d)

    else if(nospZ.eq.2) then
                                                                  ! Particulate
     OM1_CZ  = ZegC(1) + ZunC(1) + ZmortC_tot + OM1_Ratio*ZslopC_tot !  (mmol-C/m3/d)
     OM1_NZ  = ZegN(1) + ZunN(1) + ZmortN_tot + OM1_Ratio*ZslopN_tot !  (mmol-N/m3/d)
     OM1_PZ  = ZegP(1) + ZunP(1) + ZmortP_tot + OM1_Ratio*ZslopP_tot !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ  = ZegC(2) + ZunC(2) + OM2_Ratio*ZslopC_tot              !  (mmol-C/m3/d)
     OM2_NZ  = ZegN(2) + ZunN(2) + OM2_Ratio*ZslopN_tot              !  (mmol-N/m3/d)
     OM2_PZ  = ZegP(2) + ZunP(2) + OM2_Ratio*ZslopP_tot              !  (mmol-P/m3/d)

    else 

                                                                  ! Particulate
     OM1_CZ  = ZegC(1) + ZunC(1) + ZmortC_tot + OM1_Ratio*ZslopC_tot !  (mmol-C/m3/d)
     OM1_NZ  = ZegN(1) + ZunN(1) + ZmortN_tot + OM1_Ratio*ZslopN_tot !  (mmol-N/m3/d)
     OM1_PZ  = ZegP(1) + ZunP(1) + ZmortP_tot + OM1_Ratio*ZslopP_tot !  (mmol-P/m3/d)
                                                                  ! Dissolved
     OM2_CZ  = SUM(ZegC(2:nospZ)) + SUM(ZunC(2:nospZ)) + OM2_Ratio*ZslopC_tot              !  (mmol-C/m3/d)
     OM2_NZ  = SUM(ZegN(2:nospZ)) + SUM(ZunN(2:nospZ)) + OM2_Ratio*ZslopN_tot              !  (mmol-N/m3/d)
     OM2_PZ  = SUM(ZegP(2:nospZ)) + SUM(ZunP(2:nospZ)) + OM2_Ratio*ZslopP_tot              !  (mmol-P/m3/d)

    endif

   !This calculates the cumulative stoichiometry ratios for OM1_Z
   if(OM1_CZ.gt.tiny(x)) then
!    if(OM1_CZ.ne.0) then
    ff_new(k,isx1Z) = (OM1_CZ*dTd + OM1_Z) / (OM1_PZ*dTd + 1./ff(k,isx1Z)*OM1_Z) ! C/P
    ff_new(k,isy1Z) = (OM1_NZ*dTd + ff(k,isy1Z)/ff(k,isx1Z)*OM1_Z) / (OM1_PZ*dTd + 1./ff(k,isx1Z)*OM1_Z) !N/P
   else
    ff_new(k,isx1Z) = ff(k,isx1Z)
    ff_new(k,isy1Z) = ff(k,isy1Z)
   endif


   !This calculates the cumulative stoichiometry ratios for OM2_Z
   if(OM2_CZ.gt.tiny(x)) then
!    if(OM2_CZ.ne.0) then
    ff_new(k,isx2Z) = (OM2_CZ*dTd + OM2_Z) / (OM2_PZ*dTd + 1./ff(k,isx2Z)*OM2_Z) !  C/P
    ff_new(k,isy2Z) = (OM2_NZ*dTd + ff(k,isy2Z)/ff(k,isx2Z)*OM2_Z) / (OM2_PZ*dTd + 1./ff(k,isx2Z)*OM2_Z) !N/P
   else
    ff_new(k,isx2Z) = ff(k,isx2Z)
    ff_new(k,isy2Z) = ff(k,isy2Z)
   endif

!------------------------------------------------------------------------
#ifdef DEBUG
write(6,*) "In cgem, finished stoich Z"
#endif

 
!-------------------------------
!-NO3; (mmol-N/m3)
!-------------------------------     
       ff_new(k,iNO3) = AMAX1(ff(k,iNO3)                            &
       &  + ( RNO3 - AupN*NO3/(NO3+NH4) )*dTd, 0.0 )              

!--------------------------------
!-NH4; Ammonium (mmol-N/m3)
!--------------------------------
       ff_new(k,iNH4) = AMAX1(ff(k,iNH4)                            &
       & + ( RNH4 - AupN*NH4/(NO3+NH4) + AexudN + SUM(ZexN)  )*dTd, 0.0)          

!----------------------------
!-Silica: (mmol-Si/m3)
!----------------------------
       ff_new(k,iSi) =  AMAX1(ff(k,iSi)                             &
       & + ( RSi - AupSi + SUM(ZegSi) + SUM(ZunSi) )*dTd, 0.0)

!---------------------------------------------
!-PO4: Phosphate (mmol-P/m3)
!--------------------------------------
      ff_new(k,iPO4) = AMAX1(ff(k,iPO4)                             &
      & + ( RPO4 - AupP + AexudP + SUM(ZexP)  )*dTd, 0.0)

!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3)
!---------------------------------------------------------
       ff_new(k,iDIC) = AMAX1(ff(k,iDIC)                            &
       &  + ( RDIC - PrimProd + ArespC  + ZrespC )*dTd, 0.0)  
!       ff_new(k,iDIC) = RDIC - PrimProd + ArespC  + ZrespC
 
!-----------------------------------------------------------------------      
!-O2: Oxygen (mmol O2 m-3) 
!-----------------------------------------------------------------------
       ff_new(k,iO2)  = AMAX1(ff(k,iO2)                             &  
       &  + ( PrimProd - ArespC + RO2 - ZrespC)*dTd, 0.0)

!-----------------------------------------
!-OM1_A: (mmol-C/m3-- Dead Phytoplankton Particulate)
!-----------------------------------------
       ff_new(k,iOM1_A) = AMAX1(ff(k,iOM1_A)                       &
       &   + ( (ROM1_A) + OM1_CA )*dTd, 0.0)       

!----------------------------------------------------------------------
!---------------------------------------------
!-OM2_A: (mmol-C/m3-- Dead Phytoplankton Dissolved)
!---------------------------------------------
       ff_new(k,iOM2_A) = AMAX1(ff(k,iOM2_A)                       &
       &   + (ROM2_A + OM2_CA )*dTd, 0.0)              

!----------------------------------------------------------------------
!------------------------------------------
!-OM1_Z:(mmol-C/m3--G particulate)
!------------------------------------------
       ff_new(k,iOM1_Z) = AMAX1(ff(k,iOM1_Z)                     &
       &   +( ROM1_Z + OM1_CZ)*dTd, 0.0)              

!---------------------------------------------------------------------
!-----------------------------------------------
!-OM2_Z:(mmol-C/m3--G dissolved)
!-----------------------------------------------
       ff_new(k,iOM2_Z) = AMAX1(ff(k,iOM2_Z)                      &
       &   + ( ROM2_Z + OM2_CZ )*dTd, 0.0)              

!---------------------------------------------------------------------
!-------------------------------------------
!-OM1_R: (mmol-C/m3--SPM particulate)
!-------------------------------------------
       ff_new(k,iOM1_R) = AMAX1(ff(k,iOM1_R)                      &
       &   + ( ROM1_R )*dTd, 0.0)              

!---------------------------------------------------------------------
!------------------------------------------------
!-OM2_R: (mmol-C/m3--SPM dissolved)
!------------------------------------------------
       ff_new(k,iOM2_R) =  AMAX1(ff(k,iOM2_R)                     &
       &   + ( ROM2_R )*dTd, 0.0)       

!---------------------------------------------------------------------
!-------------------------------------------
!-OM1_BC: (mmol-C/m3--initial and boundary condition OM particulate)
!-------------------------------------------
       ff_new(k,iOM1_BC) = AMAX1(ff(k,iOM1_BC)                      &
       &   + ( ROM1_BC )*dTd, 0.0)


!---------------------------------------------------------------------
!------------------------------------------------
!-OM2_BC: (mmol-C/m3--initial and boundary condition OM dissolved)
!------------------------------------------------
       ff_new(k,iOM2_BC) =  AMAX1(ff(k,iOM2_BC)                     &
       &   + ( ROM2_BC )*dTd, 0.0)

!---------------------------------------------------------------------
!----------------------------
!-CDOM: (ppb) 
!----------------------------
       ff_new(k,iCDOM) =  AMAX1(ff(k,iCDOM)*(1.0 - KGcdom*dTd), 0.0)  

!!---------------------------------------------------------------------
!!----------------------------
!!-ALK: (mmol-HCO3/m3)
!!----------------------------
       ff_new(k,iALK) =  AMAX1(ff(k,iALK) +                 &
      & (RALK + AupN*NO3/(NO3+NH4) - AupN*NH4/(NO3+NH4) + AupP + 4.8*AupP)*dTd, 0.0) 
 
#ifdef DEBUG
write(6,*) "In cgem, updated ALK"
#endif
 
!Tracer
       ff_new(k,iTr) =  ff(k,iTr)       

if(writecsv==1.and.k.eq.1.and.inea.eq.10) then
  write(6201,'(*(g0,:,", "))')            ff_new(k,iA(1)), ff_new(k,iQn(1)),ff_new(k,iQp(1)), &
    & ff_new(k,iZ(1)),  ff_new(k,iZ(2)),  ff_new(k,iNO3),  ff_new(k,iNH4),  ff_new(k,iPO4),   &
    & ff_new(k,iDIC),   ff_new(k,iO2),    ff_new(k,iOM1_A),ff_new(k,iOM2_A),ff_new(k,iOM1_Z), & 
    & ff_new(k,iOM2_Z), ff_new(k,iOM1_R), ff_new(k,iOM2_R),ff_new(k,iCDOM), ff_new(k,iSi),    &
    & ff_new(k,iOM1_BC),ff_new(k,iOM2_BC),ff_new(k,iAlk),  ff_new(k,iTr),   ff_new(k,isx1A),  &
    & ff_new(k,isy1A),  ff_new(k,isx2A),  ff_new(k,isy2A), ff_new(k,isx1Z), ff_new(k,isy1Z),  &
    & ff_new(k,isx2Z),  ff_new(k,isy2Z)
  write(6301,'(*(g0,:,", "))') TC_8,rad,wind,S(k),T(k)
  write(6401,'(*(g0,:,", "))') PARsurf,PARdepth_k,PARbot
endif

!--------------------------------------------------------------------
        enddo   ! end of  "do k = 1, km" 


       !! Before Advection and VMixing, combine A's and Q's
       do k=1,km
        ff_new(k,iQn(:)) = ff_new(k,iQn(:)) * ff_new(k,iA(:))
        ff_new(k,iQp(:)) = ff_new(k,iQp(:)) * ff_new(k,iA(:))
       enddo

! ----------------------------------------------------------------------


!-- End Main GEM Calculations ---------------------------------------------------

   return
   end subroutine cgem_step
!---------------------------------------------------------------------- 
