module cgem_light

use cgem, only: nf,which_irradiance,Kw,Kchla,Kspm,Kcdom,       &
 &              aw490,astar490,astarOMA,astarOMZ,astarOMR,astarOMBC

use date_time

implicit none

contains

subroutine Call_IOP_PAR(ff,TC_8,PARsurf,depth,dz,nz,d_sfc,PARdepth)

  integer, intent(in)    :: TC_8         ! Model time in seconds from iYrS
  integer, intent(in)    :: nz           ! Number of layers
  real   , intent(in)    :: PARsurf      ! Irradiance just below sea surface
  real   , intent(in)                   :: depth   ! depth at bottom of cell k from surface
  real   , intent(in), dimension(nf,nz) :: ff      ! depth of cell
  real   , intent(in), dimension(nz)    :: dz      ! depth of cell
  real  , intent(out), dimension(nz)    :: PARdepth! PAR, visible irradiance at the middle 
                                                   !  of layer k (quanta/cm**2/sec)
  real   , intent(in), dimension(km)    :: totChl  ! total Chl-a (mg/m3)
  real, parameter :: C_cf  = 12.0E-3    ! C conversion factor (mmol-C/m3 to g-C/m3)
!----------------------------------------------------------------
! Calculate absorption (490 nm) components: seawater, chl, SPM from rivers, CDOM,
! detritus (dead cells), fecal pellets ...
  real, dimension(km) :: Chla_tot, CDOM_tot, OM1A_tot, OM1Z_tot, OM1R_tot, OM1BC_tot
  real, dimension(km) :: CDOM !After converting ppb to a490 (m-1)
  real, dimension(km) :: Chla_mass, CDOM_mass, OM1A_mass, OM1Z_mass, OM1R_mass, OM1BC_mass
  real :: sun_zenith ! Angle of the sun
  real :: calc_solar_zenith !function
  real a490_mid, aSw_mid, aChl490_mid, aCDOM490_mid, bbChl490_mid, bb490_mid
  real a490_bot, aSw_bot, aChl490_bot, aCDOM490_bot, bbChl490_bot, bb490_bot
  real aOM1A490_mid, aOM1Z490_mid, aOM1R490_mid, aOM1BC490_mid
  real aOM1A490_bot, aOM1Z490_bot, aOM1R490_bot, aOM1BC490_bot
  real cell_depth  !, bd_km1
! Time variables  
  real, parameter :: OneD60     = 1.0/60.0  ! Convert 1/min to 1/sec
  real            :: HrTC          ! Decimal hour of day
  integer         :: iYrTC, iMonTC, iDayTC, iHrTC, iMinTC, iSecTC !Time variables
  integer         :: julianDay     ! Holds Julian Day

 !---------------------------------------------------------
 ! -- Convert units for light model 
 !    C_cf == conversion factor (mmol-C/m3 to g-C/m3) 
 ! Organic Matter from dead phytoplankton (mmol/m3) 
   OM1A(:) = ff(:,iOM1A) * C_cf
 ! Organic Matter from fecal pellets      (mmol/m3)
   OM1Z(:) = f(:,iOM1CZ) * C_cf
 ! Suspended Particulate Matter (SPM)    (mmol/m3) 
 ! There is 1.8% Organic Matter in SPM originating from the rivers.
   OM1SPM(:) = ff(:,iOM1R) * C_cf / 0.018
 ! Organic Matter from boundary conditions(mmol/m3) 
   OM1BC(:)  = ff(:,iOM1BC) * C_cf
!----------------------------------------------------------------
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

! First, convert CDOM(ppb) into CDOM, a490 (m-1)
! Once the CDOM (QSE ppb) is in the model domain, we advect and mix using the same 
! routine as for other dissolved constituents. However, to use the CDOM in the light 
! attenuation models, we need to calculate a490 (Penta et al. 2008). 
! 1) convert CDOM(QSE ppb) back to a312: a312 = (CDOM(QSE ppb)-0.538)/2.933 
! 2) convert a312 to a490: a490 = a312*exp(-0.016*(490-312)), where here S = 0.016 
! (mean value from D'Sa and DiMarco (2008)
   do k=1,nz
      CDOM(k) = (ff(k,iCDOM) - 0.538)/2.933 !ppb to a312
      CDOM(k) = CDOM(k) * exp(-0.016*(490.-312.))
      CDOM(k) = AMAX1(CDOM(k),0.)
   enddo 

#ifdef DEBUG
write(6,*) "Call_IOP_Par: initialized CDOM"
#endif

do k = 1, nz 
    Chla_mass(k) = 0.0
    do isp = 1, nospA
      Chla_mass(k) =  Chla_mass(k) + ff(iA(isp),k) * ff(iQc(isp),k) * 12. * (1./CChla(isp))
    enddo ! isp = 1, nospA
  enddo ! k = 1, km

!Mass in each cell at layer k (area of volume part cancels out)
!The unit is (g C)
!      cell_depth = bottom_depth(k) - bd_km1
!      bd_km1 = bottom_depth(k)
      Chla_mass(:) = totChl(:)*Vol(:)
      CDOM_mass(:) = CDOM(:)*Vol(:)
      OM1A_mass(:) = OM1A(:)*Vol(:)
      OM1Z_mass(:) = OM1Z(:)*Vol(:)
      OM1R_mass(:) = OM1R(:)*Vol(:)
      OM1BC_mass(:) = OM1BC(:)*Vol(:)

!Mass from surface to center of cell at layer k
!Is the sum of the mass of all previous k layers plus 
!half of the current k layer 
!Concentration is mass/Vol 
      Vol_tot(1)  = 0.5*Vol(1)
      Chla_tot(1) = 0.5*Chla_mass(1)
      CDOM_tot(1) = 0.5*CDOM_mass(1)
      OM1A_tot(1) = 0.5*OM1A_mass(1)
      OM1Z_tot(1) = 0.5*OM1Z_mass(1)
      OM1R_tot(1) = 0.5*OM1R_mass(1)
      OM1BC_tot(1) = 0.5*OM1BC_mass(1)

      Vol_tot(2:nz)  = 0.5*Vol(2:nz) + SUM(Vol(1:k-1))
      Chla_tot(2:nz)  = 0.5*Chla_mass(2:nz) + SUM(Chla_mass(1:k-1))
      CDOM_tot(2:nz)  = 0.5*CDOM_mass(2:nz) + SUM(CDOM_mass(1:k-1))
      OM1A_tot(2:nz)  = 0.5*OM1A_mass(2:nz) + SUM(OM1A_mass(1:k-1))
      OM1Z_tot(2:nz)  = 0.5*OM1Z_mass(2:nz) + SUM(OM1Z_mass(1:k-1))
      OM1R_tot(2:nz)  = 0.5*OM1R_mass(2:nz) + SUM(OM1R_mass(1:k-1))
      OM1BC_tot(2:nz) = 0.5*OM1BC_mass(2:nz)+ SUM(OM1BC_mass(1:k-1))

   do k=1,nz
!Calculate absorption coefficients:
      aSw_mid = aw490  !Sea water absorption at mid cell
      aChl490_mid = astar490 * Chla_tot(k) / Vol_tot(k)        !Chla absorption at mid cell
      aCDOM490_mid = CDOM_tot(k) / Vol_tot(k)        !CDOM absorption at mid cell
      aOM1A490_mid = astarOMA * OM1A_tot(k) / Vol_tot(k) ! absorption at mid cell
      aOM1Z490_mid = astarOMZ * OM1Z_tot(k) / Vol_tot(k) ! absorption at mid cell
      aOM1R490_mid = astarOMR * OM1R_tot(k) / Vol_tot(k) ! absorption at mid cell
      aOM1BC490_mid = astarOMBC * OM1BC_tot(k) / Vol_tot(k) ! absorption at mid cell
      a490_mid = aSw_mid + aChl490_mid + aCDOM490_mid + aOM1A490_mid + aOM1Z490_mid + aOM1R490_mid + aOM1BC490_mid
!Calculate backscattering coefficients:
      bbChl490_mid = 0.015 * (0.3*((Chla_tot(k) / Vol_tot(k))**0.62)*(550./490.)) !Chla backscatter at mid cell
      bb490_mid = bbChl490_mid !Only Chla backscatters for now
      call IOP_PARattenuation(a490_mid, bb490_mid, PARsurf, sun_zenith, d_sfc(k), PARdepth(k)) 
   enddo

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

!---------------------------------------------------------------------- 
end module
