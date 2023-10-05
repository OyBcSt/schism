subroutine cgem_flux(dT)

use grid, only:T,S,dz,Wind
!This is called after cgem_step, which returns ff_new
!This modifies the surface and bottom cells of ff_new
use cgem, only:ff_new,Which_fluxes,iO2surf,iDICsurf,iO2,iDIC,pCO2,iALK,iSi,iPO4,pH
use MOD_UTILITIES
!from mocxy
use gasx

implicit none

real, intent(in) :: dT
real :: T_sfc, Sal_sfc, O2_sfc, Sc, Op_umole, rhow, Op, OsDOp
real :: Vtrans, alpha_O2, O2_atF,zs, DIC_sfc, CO2_atF
real, parameter :: SDay = 86400.0  ! # of sec in 24 hr day
!------------------------------------------------------------------
!Output vars for mocsy subroutine:
      real :: kw660(1), co2flux(1), co2ex(1), dpco2(1)
      real :: ph_calc(1), pco2_calc(1), fco2(1), co2(1), hco3(1), co3(1), omegaa(1) 
      real :: omegac(1), betad_calc(1), rhosw(1), p(1), tempis(1)
      real :: patm(1) = 1., pco2_in(1)
      real :: m_alk(1), m_dic(1), m_si(1), m_po4(1)




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

if(Which_fluxes(iDICsurf).eq.1) then
!--------------------------------------------------------------
! Calc  SFLUX_CO2, the sea surface vertical flux of CO2
!--------------------------------------------------------------
               zs      = dz(1)       ! Thickness (m.) of the water column

               T_sfc   = T(1)        ! Temperature (C) in sfc layer, k=1
               Sal_sfc = S(1)        ! Salinity        in sfc layer, k=1
               DIC_sfc = ff_new(1,iDIC) ! Dissolved Inorganic Carbon
                                         !    (mmol m-3) in sfc layer, k=1

             !----------------------------------------------------------
             ! Units of gas_exchange are mmol CO2 m-2 s-1 
             !----------------------------------------------------------
!use mocsy instead but calculate to compare
             CO2_atF = gas_exchange(T_sfc,Sal_sfc,DIC_sfc,zs,pH(1),pCO2)
                          ff_new(1,iDIC) = AMAX1(ff_new(1,iDIC) - CO2_atF/dz(1)*dT,0.)

elseif(Which_fluxes(iDICsurf).eq.2) then
!---------------using mocsy:------------------------------------
  kw660(1) = 0.01/3600.*0.31*Wind*Wind !*schmidtnumberterm
  pCO2_in(1) = pCO2
            !!! MOCSY alkalinity expressions:
            m_alk = ff_new(1,iALK)/1000.
            m_dic = ff_new(1,iDIC)/1000.
            m_si  = ff_new(1,iSi)/1000.
            m_po4 = ff_new(1,iPO4)/1000.

  call flxco2(co2flux, co2ex, dpco2, &
 &            ph_calc, pco2_calc, fco2, co2, hco3, co3, omegaa, omegac, betad_calc, rhosw, p, tempis,  &
 &            T(1), S(1), m_alk, m_dic, m_si, m_po4, kw660, pCO2_in, patm, zs, 1, &
 &            'mol/m3', 'Tinsitu', 'm ', 'u74', 'l  ', 'pf ', 'Pzero  ')

              ff_new(1,iDIC) = AMAX1(ff_new(1,iDIC) + 1000.*co2flux(1)/dz(1)*dT,0.)

endif

return

end subroutine cgem_flux
