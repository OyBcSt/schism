subroutine cgem_flux(dT)

use grid, only:T,S,dz,Wind
!This is called after cgem_step, which returns ff_new
!This modifies the surface and bottom cells of ff_new
use cgem, only:ff_new,Which_fluxes,iO2surf,iO2
use MOD_UTILITIES

implicit none

real, intent(in) :: dT
real :: T_sfc, Sal_sfc, O2_sfc, Sc, Op_umole, rhow, Op, OsDOp
real :: Vtrans, alpha_O2, O2_atF
real, parameter :: SDay = 86400.0  ! # of sec in 24 hr day


if(Which_fluxes(iO2surf).eq.1) then
!--------------------------------------------------------------
! Calc  O2_atF, the sea surface vertical flux of O2
!--------------------------------------------------------------
   T_sfc    = T(1)       ! Temperature (C)   in sfc layer, k=1
   Sal_sfc  = S(1)       ! Salinity          in sfc layer, k=1
   O2_sfc   = ff_new(1,iO2) ! O2 (mmol-O2/m3) in sfc layer, k=1

   Sc       = SchmidtNumber(Sal_sfc,T_sfc,0)  ! Schmidt number,
                                                          !   0 (zero)
                                                          !   for O2

   Op_umole = o2sat(Sal_sfc,T_sfc)     ! O2 saturation,
                                                   !    (umol-O2/kg)

   rhow     = sw_dens0(Sal_sfc,T_sfc)  ! water density [kg/m3]

   Op       = rhow * Op_umole * 1.0E-3 ! O2 saturation,
                                                   !    (mmol-O2/m3)
   OsDOp    = O2_sfc/Op

!--------------------------------------------------------------
!  Vtrans below is the O2 transfer vel (m/s)
!
!  Vtrans   = (5.9*(kw)*(OsDOp*OsDOp))*(Sc)**X
!    where kw and Sc are dependent on Wind Speed.
!  Values kw and X are from Liss and Merlivat, 1986.
!  Factor of OsDOp**2 is from Justic, et. al 2002 
!   for when saturation levels are above 125%.
!---------------------------------------------------------------
 if(Wind.lt.3.6) then
   Vtrans        = AMAX1((5.9 * (0.17*Wind)         &
   &              * Sc**(-2./3.) / SDay), 0.)
 else if(Wind.le.13.) then
   Vtrans        = AMAX1((5.9 *(2.85*Wind - 9.65 )    &
   &              / SQRT(Sc) / SDay), 0.)
 else
   Vtrans        = AMAX1((5.9 *(5.9*Wind - 49.3 )    &
   &              / SQRT(Sc) / SDay), 0.)
 endif
 if(OsDOp.gt.1.25) Vtrans = Vtrans * (OsDOp*OsDOp)

   alpha_O2       = 1.025

   O2_atF         = Vtrans*(O2_sfc - alpha_O2*Op)
                                           ! flux of O2 thru
                                           ! the
                                           ! sea sfc
                                           ! ((mmol-O2/m2/sec)
                                           ! negative means
                                           ! into
   ff_new(1,iO2) = AMAX1(ff_new(1,iO2) - O2_atF/dz(1)*dT,0.)

endif

end subroutine cgem_flux
