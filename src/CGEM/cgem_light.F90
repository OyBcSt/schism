module cgem_light

use grid, only: km,d,dz,d_sfc,iYrS,lat,lon
use cgem, only: which_irradiance,Kw,Kchla,Kspm,Kcdom,       &
 &              aw490,astar490,astarOMA,astarOMZ,astarOMR,astarOMBC

use date_time

implicit none

contains

   subroutine Call_IOP_PAR(                           & 
                 & TC_8, PARsurf, totChl, CDOM_k,     &
                 & OM1A, OM1Z, OM1R, OM1BC,           &
                 & depth, dz, nz, d_sfc,              &
                 & PARdepth, PARbot)

!-----------------------------------------------------------------------
! Code is based on Brad Penta's code
!---------------------------------------------------------------------
  integer, intent(in)    :: TC_8         ! Model time in seconds from iYrS
  integer, intent(in)    :: nz           ! Number of layers
  real   , intent(in)    :: PARsurf      ! Irradiance just below sea surface
  real   , intent(in), dimension(km)    :: totChl(km)   ! total Chl-a (mg/m3)
  real   , intent(in), dimension(km)    :: CDOM_k(km)   ! CDOM (ppb) 
  real   , intent(in), dimension(km)    :: OM1A(km)     ! Concentration of particulate
  real   , intent(in), dimension(km)    :: OM1Z(km)     ! Concentration of particulate
                                                        !  fecal pellets (g/m3)
  real   , intent(in), dimension(km)    :: OM1R(km)     ! Concentration of particulate
                                                        !  river generated SPM (g/m3)
  real   , intent(in), dimension(km)    :: OM1BC(km)    ! Concentration of particulate initial
                                                        !  and boundary condition SPM (g/m3)
  real   , intent(in)                   :: depth        ! depth at bottom of cell k from surface
  real   , intent(in), dimension(km)    :: dz(km)       ! depth of cell
  real   , intent(in), dimension(km)    :: d_sfc(km)    ! depth at center of cell k from surface
  real  , intent(out), dimension(km)    :: PARdepth(km) ! PAR, visible irradiance at the middle 
                                                        !  of layer k (quanta/cm**2/sec)
  real, intent(out)      :: PARbot                      ! PAR at sea bottom (quanta/cm**2/sec)
!----------------------------------------------------------------
! Calculate absorption (490 nm) components: seawater, chl, SPM from rivers, CDOM,
! detritus (dead cells), fecal pellets ...
  real, dimension(km) :: Chla_tot, CDOM_tot, OM1A_tot, OM1Z_tot, OM1R_tot, OM1BC_tot, CDOM
  real, dimension(km) :: Chla_mass, CDOM_mass, OM1A_mass, OM1Z_mass, OM1R_mass, OM1BC_mass
  real :: sun_zenith ! Angle of the sun
  real :: calc_solar_zenith !function
  real a490_mid, aSw_mid, aChl490_mid, aCDOM490_mid, bbChl490_mid, bb490_mid
  real a490_bot, aSw_bot, aChl490_bot, aCDOM490_bot, bbChl490_bot, bb490_bot
  real aOM1A490_mid, aOM1Z490_mid, aOM1R490_mid, aOM1BC490_mid
  real aOM1A490_bot, aOM1Z490_bot, aOM1R490_bot, aOM1BC490_bot
  real cell_depth  !, bd_km1
  integer :: k  
! Time variables  
  real, parameter :: OneD60     = 1.0/60.0  ! Convert 1/min to 1/sec
  real            :: HrTC          ! Decimal hour of day
  integer         :: iYrTC, iMonTC, iDayTC, iHrTC, iMinTC, iSecTC !Time variables
  integer         :: julianDay     ! Holds Julian Day

 ! First calculate the Julian(GMT) model year (iYrTC), month (iMonTC), 
 ! day (iDayTC), hour (iHrTC), minute (iMinTC), and second (iSecTC) 
 ! associated with the midpoint of the present timestep istep TC_8
      CALL DATE_TIMESTAMP( iYrS, TC_8, &
                           iYrTC, iMonTC, iDayTC, iHrTC, iMinTC, iSecTC )
 ! Calc HrTC, the decimal hour of day
       HrTC = real(iHrTC,4) + OneD60*iMinTC + OneD60*iSecTC
 ! Now calculate the Julian Day associated with model time TC_8
      julianDay = JDAY_IN_YEAR(iYrTC, iMonTC, iDayTC)
 ! Now calculate sun_zenith, the solar beam zenith angle in radians
 ! for a given GMT Julian day, hour, longitude and latitude 
     sun_zenith = calc_solar_zenith(lat, lon,  HrTC, julianDay )

#ifdef DEBUG_PAR
write(6,*) "Begin Call_IOP_Par"
write(6,*) "PARsurf,sun_zenith,totChl",PARsurf,sun_zenith,totChl
write(6,*) "nz,km",nz,km
#endif


! First, convert CDOM(ppb) into CDOM, a490 (m-1)
! Once the CDOM (QSE ppb) is in the model domain, we advect and mix using the same 
! routine as for other dissolved constituents. However, to use the CDOM in the light 
! attenuation models, we need to calculate a490 (Penta et al. 2008). 
! 1) convert CDOM(QSE ppb) back to a312: a312 = (CDOM(QSE ppb)-0.538)/2.933 
! 2) convert a312 to a490: a490 = a312*exp(-0.016*(490-312)), where here S = 0.016 
! (mean value from D'Sa and DiMarco (2008)
   do k=1,nz
      CDOM(k) = (CDOM_k(k) - 0.538)/2.933 !ppb to a312
      CDOM(k) = CDOM(k) * exp(-0.016*(490.-312.))
      CDOM(k) = AMAX1(CDOM(k),0.)
   enddo 

#ifdef DEBUG
write(6,*) "Call_IOP_Par: initialized CDOM"
#endif


!Initialize counters for Chla, CDOM, and detritus:
   Chla_tot = 0.
   CDOM_tot = 0.
   OM1A_tot = 0.
   OM1Z_tot = 0.
   OM1R_tot = 0.
   OM1BC_tot = 0.
!   bd_km1 = 0.

!Mass in each cell at layer k (area of volume part cancels out)
!The unit is mg[mmol] / m2
   do k=1,nz
!      cell_depth = bottom_depth(k) - bd_km1
!      bd_km1 = bottom_depth(k)
      cell_depth = dz(k)
      Chla_mass(k) = totChl(k)*cell_depth
      CDOM_mass(k) = CDOM(k)*cell_depth
      OM1A_mass(k) = OM1A(k)*cell_depth
      OM1Z_mass(k) = OM1Z(k)*cell_depth
      OM1R_mass(k) = OM1R(k)*cell_depth
      OM1BC_mass(k) = OM1BC(k)*cell_depth
   enddo


#ifdef DEBUG
write(6,*) "Call_IOP_Par:cell mass"
#endif

!Mass from surface to center of cell at layer k
!Is the sum of the mass of all previous k layers plus 
!half of the current k layer 
!Concentration is that divided by the distance
!from the surface to the center of cell at layer k
!(Division by d_sfc is in the next loop)  
      Chla_tot(1) = 0.5*Chla_mass(1)
      CDOM_tot(1) = 0.5*CDOM_mass(1)
      OM1A_tot(1) = 0.5*OM1A_mass(1)
      OM1Z_tot(1) = 0.5*OM1Z_mass(1)
      OM1R_tot(1) = 0.5*OM1R_mass(1)
      OM1BC_tot(1) = 0.5*OM1BC_mass(1)
   do k=2,nz
      Chla_tot(k)  = 0.5*Chla_mass(k) + SUM(Chla_mass(1:k-1))
      CDOM_tot(k)  = 0.5*CDOM_mass(k) + SUM(CDOM_mass(1:k-1))
      OM1A_tot(k)  = 0.5*OM1A_mass(k) + SUM(OM1A_mass(1:k-1))
      OM1Z_tot(k)  = 0.5*OM1Z_mass(k) + SUM(OM1Z_mass(1:k-1))
      OM1R_tot(k)  = 0.5*OM1R_mass(k) + SUM(OM1R_mass(1:k-1))
      OM1BC_tot(k) = 0.5*OM1BC_mass(k)+ SUM(OM1BC_mass(1:k-1))
   enddo

#ifdef DEBUG
write(6,*) "Call_IOP_Par:total mass"
#endif

#ifdef DEBUG
write(6,*) "Call_IOP_Par:depth,dz,d_sfc:",depth,dz,d_sfc
#endif


   do k=1,nz
!Calculate absorption coefficients:
      aSw_mid = aw490  !Sea water absorption at mid cell
      aChl490_mid = astar490 * Chla_tot(k) / d_sfc(k)        !Chla absorption at mid cell
      aCDOM490_mid = CDOM_tot(k) / d_sfc(k)        !CDOM absorption at mid cell
      aOM1A490_mid = astarOMA * OM1A_tot(k) / d_sfc(k) ! absorption at mid cell
      aOM1Z490_mid = astarOMZ * OM1Z_tot(k) / d_sfc(k) ! absorption at mid cell
      aOM1R490_mid = astarOMR * OM1R_tot(k) / d_sfc(k) ! absorption at mid cell
      aOM1BC490_mid = astarOMBC * OM1BC_tot(k) / d_sfc(k) ! absorption at mid cell
      a490_mid = aSw_mid + aChl490_mid + aCDOM490_mid + aOM1A490_mid + aOM1Z490_mid + aOM1R490_mid + aOM1BC490_mid
!Calculate backscattering coefficients:
      bbChl490_mid = 0.015 * (0.3*((Chla_tot(k) / d_sfc(k))**0.62)*(550./490.)) !Chla backscatter at mid cell
      bb490_mid = bbChl490_mid !Only Chla backscatters for now
! Calculate PAR at depth
      !Why would we check if a490_mid=0??
      !If it is zero, stop.
      !  if(a490_mid.le.0) then
      !     write(6,*) k,CDOM(k),CDOM_k(k)
      !     write(6,*) "a490_mid.le.0, =",a490_mid,aSw_mid,aChl490_mid,aCDOM490_mid
      !     !stop
      !  endif
      call IOP_PARattenuation(a490_mid, bb490_mid, PARsurf, sun_zenith, d_sfc(k), PARdepth(k)) 

#ifdef DEBUG
write(6,*) "Just called IOP_Paratt, PARdepth=:",k,PARdepth(k)
#endif

   enddo

! Calculate PAR at sea bottom
      aSw_bot = aw490  !Sea water absorption at bottom of cell
      aChl490_bot = astar490 * (Chla_tot(nz)+0.5*Chla_mass(nz)) / depth !Chla absorption at bottom
      aCDOM490_bot = CDOM_tot(nz)+(0.5*CDOM_mass(nz)) / depth !CDOM absorption at bottom
      aOM1A490_bot = astarOMA * (OM1A_tot(nz)+0.5*OM1A_mass(nz)) / depth    !A absorption at bottom
      aOM1Z490_bot = astarOMZ * (OM1Z_tot(nz)+0.5*OM1Z_mass(nz)) / depth !FP absorption at bottom
      aOM1R490_bot = astarOMR * (OM1R_tot(nz)+0.5*OM1R_mass(nz)) / depth !SPM absorption at bottom
      aOM1BC490_bot = astarOMbc * (OM1BC_tot(nz)+0.5*OM1BC_mass(nz)) / depth !INIT/BC absorption at bottom
      a490_bot = aSw_bot + aChl490_bot + aCDOM490_bot + aOM1A490_bot + aOM1Z490_bot + aOM1R490_bot + aOM1BC490_bot
      bbChl490_bot = 0.015 * (0.3*(((Chla_tot(nz)+0.5*Chla_mass(nz)) / depth)**0.62)*(550./490.))
!!Chla backscatter at bottom
      bb490_bot = bbChl490_bot !Only Chla backscatters for now
      call IOP_PARattenuation(a490_bot, bb490_bot, PARsurf, sun_zenith, d_sfc(nz), PARbot)

   END subroutine Call_IOP_PAR
!----------------------------------------------------------------------
!*****************************************************************************************
      subroutine IOP_PARattenuation(a490, bb490, PARsurf, sun_zenith, d_sfc, PARdepth)

! Penta et al., 2008 PAR penetration model: based on the IOP (inherent optical properties)
! absorption and backscatter at 490nm and the zenith angle of the sun; modified from Lee 
! et al., 2005 PAR penetration model: based on satellite absorption and backscater at 490nm

!    Lee, Z., K. Du, R. Arnone, S. Liew, and B. Penta, .Penetration of solar radiation
!         in the upper ocean: a numerical model for oceanic and coastal waters,. 
!         J. Geophys. Res. 110, C09019 doi:10.1029/2004JC002780 (2005).
!    Penta, B., Z. Lee, R.M. Kudela, S.L. Palacios, D.J. Gray, J.K. Jolliff, and 
!         I.G. Shulman, .An underwater light attenuation scheme for marine ecosystem 
!         models., Optics Express, 16, 16582-16591 (2008).

! Absorption (a490) is the total absorption at 490 nm calculated from all of the 
! constituents of the model that absorb light (Seawater, chlorophyll, CDOM, detritus, etc.
! Backscatter (bb490) is similar 

! Use the average value from the sea-surface to the depth of the calculation 
! (these calculations moved to a new subroutine(s) and the a490 and bb490 will be passed 
! into this subroutine  

! PAR (photosynthetically active radiation) just below the sea surface is used as a 
! starting value to be multiplied by the attenuation factor computed in this subroutine. 
! The starting value does not affect any of these calculations

! We generally use a factor of 0.4815 to represent the loss of PAR passing across the
! Air/Sea interface - this is included already in the data 

! Define coefficients for light attenuation model       
      real a490, alpha0, alpha1, alpha2, bb490, chi0, chi1, chi2, d_sfc
      real k1, k2, kpar, PARdepth, PARsurf, sun_zenith, zeta0, zeta1
      real zeta2 
     
! sun_zenith => Solar Zenith Angle in radians for time and location
! d_sfc => depth of center of cell from elevated sea surface        

! Set values of coefficients for light attenuation model 
! (Lee et al., 2005; Penta et al., 2008)       
      parameter(alpha0=0.090, alpha1=1.465, alpha2=-0.667, chi0=-0.057, chi1=0.482, & 
     &          chi2=4.221, zeta0=0.183, zeta1=0.702, zeta2=-2.567)
     
! Calculate the attenuation (Equation 2 Penta et al., 2008 errata)
      k1 = ((chi0 + chi1*(a490**0.5) + chi2*bb490) & 
     &     * (1 + alpha0 * sin(sun_zenith)))
     
      k2 = ((zeta0 + zeta1*a490 + zeta2*bb490) &
     &     * (alpha1 + alpha2 * cos(sun_zenith) ))
     
      kpar = k1 + (k2 / (1+(d_sfc))**0.5)
     
      PARdepth = PARsurf * exp(-(kpar*(d_sfc)))

      return
      end
!*****************************************************************************************

subroutine calc_PARdepth( TC_8,PARSurf,S_k,Chla_k,CDOM_k,OM1A_in,OM1Z_in,      &
 &                         OM1R_in,OM1BC_in,PARdepth_k,PARbot )

!---------------------------------------------
! Interface variables
!---------------------------------------------------------------------
    !Inputs
    integer, intent(in)              :: TC_8              ! Model time (seconds from iYrS)
    real, dimension(km), intent(in)  :: S_k               ! Salinity (psu)
    real, dimension(km), intent(in)  :: Chla_k            ! Total amount of Chl-a in all the
                                                          !  phytoplankton species (mg/m3) per cell
    real, dimension(km), intent(in)  :: OM1A_in, OM1Z_in  !  
    real, dimension(km), intent(in)  :: OM1R_in, OM1BC_in ! POC in g/m3
    real, dimension(km), intent(in)  :: CDOM_k            ! CDOM, ppb
    real,                intent(in)  :: PARsurf           ! Irradiance just below the sea surface (quanta/cm2/s) 
    !Outputs
    real, dimension(km), intent(out) :: PARdepth_k(km) ! Irradiance at center of layer k (quanta/cm2/s)
    real,                intent(out) :: PARbot         ! Irradiance at sea floor (quanta/cm2/s)
!---------------------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------
    real, dimension(km)  :: OM1A_k, OM1Z_k    !  
    real, dimension(km)  :: OM1SPM_k, OM1BC_k   ! POC in g/m3

    integer        ::  k    ! Loop index
!------------------------------------ 
! Variables needed for light routine and calc_Agrow
    real    :: Katt               ! Attenuation coefficient for Irradiance model 2 
    real    :: tmpexp             ! Intermediate calculation
    real    :: PARbotkm1          ! Irradiance at bottom of layer k-1 (quanta/cm2/s)
    real    :: PARtopk            ! Irradiance at top of layer k (quanta/cm2/s)
    real, parameter :: C_cf  = 12.0E-3    ! C conversion factor (mmol-C/m3 to g-C/m3) 
!-----------------------------------------------------------------------
#ifdef DEBUG_PAR
write(6,*) "Begin calc_light:calc_PARdepth, TC_8",TC_8
#endif

!-----------------------------------------------------------------
!   Begin main ij loop for the biogeochemistry 
!   calculations at time-level istep
!-----------------------------------------------------------------
 !---------------------------------------------------------
 ! Calculate and convert variables needed for light routine
 !---------------------------------------------------------
      do k = 1, km

         !-------------------------------------------------
         ! -- Convert units for light model 
         !    C_cf == conversion factor (mmol-C/m3 to g-C/m3) 
         !-----
         ! Organic Matter from dead phytoplankton (mmol/m3) 
         !            converted to equivalent (g carbon/m3)
           OM1A_k(k)  = OM1A_in(k) * C_cf        
         !-----
         ! Organic Matter from fecal pellets      (mmol/m3)
         !            converted to equivalent (g carbon/m3)
           OM1Z_k(k) = OM1Z_in(k) * C_cf   
         !-----
         ! Organic Matter from rivers            
         !  Suspended Particulate Matter (SPM)    (mmol/m3) 
         !            converted to equivalent (g carbon/m3)
         ! There is 1.8% Organic Matter in SPM originating from the rivers.
           OM1SPM_k(k) = OM1R_in(k) * C_cf / 0.018
         !-----
         ! Organic Matter from boundary conditions(mmol/m3) 
         !            converted to equivalent (g carbon/m3)
           OM1BC_k(k)  = OM1BC_in(k) * C_cf 

      enddo ! End  "do k = 1, km" loop


!----------------------------------------------------------------------
! Execute the desired atmospheric light model.  To calculate PARsurf,
! the effect amount of downward spectrally integrated irradiance 
! just below the sea surface.  'Rad' is just above sea surface. 
!----------------------------------------------------------------------
#ifdef DEBUG_PAR
write(6,*) "In cgem_light:calc_PARdepth, calculate date_time is next, iYrS=",iYrS
#endif

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
         select case (which_irradiance)

                 !--------------------------------------------
         case (1)! Upgraded form of the underwater light model
                 ! developed by Brad Penta of NRL is used
                 !--------------------------------------------

#ifdef DEBUG_PAR
  write(6,*) "In calc_PARdepth, calling Brad's light model, km,PARsurf=",km,PARsurf
#endif

                call Call_IOP_PAR(                    &
                 & TC_8, PARsurf, Chla_k, CDOM_k,     &
                 & OM1A_k, OM1Z_k, OM1SPM_k, OM1BC_k, &
                 & d(km), dz, km, d_sfc,              &
                 & PARdepth_k, PARbot)

#ifdef DEBUG_PAR
     write(6,*) "In cgem_light:calc_PARdepth."
     write(6,*) "PARdepth_k,PARbot=",PARdepth_k
     write(6,*) "PARbot=",PARbot
#endif

                 !-------------------------------------------------
         case (2)! Upgraded form of the original underwater light
                 ! model of Pete Eldridge is used. Now accounts for
                 ! light attenuation in each k layer rather than for the whole
                 ! mixed layer as in the Eldridge & Roelke(2010) code 
                 !-------------------------------------------------
                 PARbotkm1 = PARsurf             ! initialize at top of
                                                 ! column i., i.e.
                                                 ! at bottom of layer "zero".
                 do k = 1, km
                   !Calculate attenuation coefficient
                     Katt    = Kw                                                   &
                     &       + Kchla * Chla_k(k)                                    &
                     &       + Kspm  * (OM1SPM_k(k)+OM1A_k(k)+OM1Z_k(k)+OM1BC_k(k)) &
                     &       + Kcdom * CDOM_k(k)                                    &
                     &       + (((0.0022*S_k(k))-0.158)*S_k(k)+3.03)
                     PARtopk = PARbotkm1         ! irradiance at top of
                                                 ! layer k is same as
                                                 ! irradiance at bottom
                                                 ! of layer km1
                     tmpexp  = exp(-0.5*Katt*dz(k))
                     PARdepth_k(k) = PARtopk * tmpexp    ! irradiance at middle of layer k
                     PARbot  = PARdepth_k(k) * tmpexp    ! irradiance at bottom of layer k
                     PARbotkm1 = PARbot         ! reinitialize for next top layer
                 enddo 

         case default
             write(6,*) "Error in irradiance switch",which_irradiance
             stop
         end select

!---------------------End Underwater Light Model-----------------------------------
#ifdef DEBUG_PAR
write(6,*) "In cgem_light:calc_PARdepth, Finished Underwater Light Model"
#endif

   return
   end subroutine calc_PARdepth  
!---------------------------------------------------------------------- 

end module
