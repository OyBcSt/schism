subroutine cgem_flux(dT,istep,CBODW,Esed)

use grid, only:T,S,dz,Wind,km
!This is called after cgem_step, which returns ff_new
!This modifies the surface and bottom cells of ff_new
use cgem, only:ff_new,Which_fluxes,iO2surf,iDICsurf,iO2,iDIC,pCO2,iALK,iSi,iPO4,pH, &
        & iSOC,iA,iNO3,iNH4,iNutEx,iMPB,nospA,SDay,iSDM,dT_sed,&
        & iOM1_A,iOM2_A,iOM1_Z,iOM2_Z,iOM1_R,iOM2_R,iOM1_BC,iOM2_BC,nf,fmin
use SDM
use MOD_UTILITIES
!from mocxy
use gasx

implicit none

integer, intent(in) :: istep
real, intent(in) :: dT, CBODW, Esed
integer :: nz
real :: T_sfc, Sal_sfc, O2_sfc, Sc, Op_umole, rhow, Op, OsDOp
real :: Vtrans, alpha_O2, O2_atF,zs, DIC_sfc, CO2_atF
!------------------------------------------------------------------
!Output vars for mocsy subroutine:
      real :: kw660(1), co2flux(1), co2ex(1), dpco2(1)
      real :: ph_calc(1), pco2_calc(1), fco2(1), co2(1), hco3(1), co3(1), omegaa(1) 
      real :: omegac(1), betad_calc(1), rhosw(1), p(1), tempis(1)
      real :: patm(1) = 1., pco2_in(1)
      real :: m_alk(1), m_dic(1), m_si(1), m_po4(1)
!Bottom flux
  real :: SOC, DICFlux,tau,O2Flux,NO3Flux,NH4Flux,PO4Flux,SiFlux,ALKFlux
  real :: sedflux_iOM1_A,sedflux_iOM2_A,sedflux_iOM1_Z,sedflux_iOM2_Z
  real :: sedflux_iOM1_BC,sedflux_iOM2_BC,sedflux_iOM1_R,sedflux_iOM2_R
  real, dimension(nf) :: f
  !-- Convert quanta/cm2/s to mol photons/m2/d
      !   N_Av=6.0221413E+23
      !   quanta/cm2/s * 1 mol/N_av quanta * 10,000cm2/m2 * 86400s/d = mol/m2/d
      real, parameter :: convert_I   = 1. / 6.0221413 * 8.64 * 1.e-15
      integer :: nstep_sed

nz = km
nstep_sed = int(dT_sed / dT)  ! number of steps in-between calls to sediment diagenesis

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
   ff_new(1,iO2) = AMAX1(ff_new(1,iO2) - O2_atF/dz(1)*dT,fmin(iO2))

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
                          ff_new(1,iDIC) = AMAX1(ff_new(1,iDIC) - CO2_atF/dz(1)*dT,fmin(iDIC))

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

              ff_new(1,iDIC) = AMAX1(ff_new(1,iDIC) + 1000.*co2flux(1)/dz(1)*dT,fmin(iDIC))

endif

!-- BOTTOM FLUXES -------------------------------------------------------------------------
if(Which_fluxes(iSOC).eq.1) then
!Murrell and Lehrter sediment oxygen consumption
       SOC = - 0.0235*2.**(.1*T(nz))*ff_new(nz,iO2)
               ff_new(nz,iO2) = AMAX1(ff_new(nz,iO2)  + SOC/  &
     & dz(nz)*dT/SDay,fmin(iO2))
       DICFlux = (-3.7*log(AMAX1(ff_new(nz,iO2),1.e-8)) + 19.4)*SOC
               ff_new(nz,iDIC) = AMAX1(ff_new(nz,iDIC) + DICFlux/  &
     & dz(nz)*dT/SDay,fmin(iDIC))
elseif(Which_fluxes(iSOC).eq.2.or.Which_fluxes(iSOC).eq.3) then
        write(6,*) "Error, iSOC is 2 or 3"
        write(6,*) "JW_SOC requires storing CBODS across grid and has not been implemented"
        write(6,*) "Stopping"
        stop
!L3flux !Justic and Wang sediment oxygen consumption
!L3flux      tau=0.
!L3flux      call JW_SOC(O2Flux,NH4Flux,PO4Flux,CBODW,ff_new(nz,iA(1):iA(nospA)),Esed,ff_new(nz,iO2), &
!L3flux        T(nz),tau)
!L3flux !O2
!L3flux                ff_new(nz,iO2) = AMAX1(ff_new(nz,iO2)    + O2Flux/  &
!L3flux      & dz(nz)*dT/SDay,0.)
!L3flux !NH4
!L3flux                ff_new(nz,iNH4) = AMAX1(ff_new(nz,iNH4)  + NH4Flux/  &
!L3flux      & dz(nz)*dT/SDay,0.)
!L3flux !PO4
!L3flux                ff_new(nz,iPO4) = AMAX1(ff_new(nz,iPO4)  + PO4Flux/  &
!L3flux      & dz(nz)*dT/SDay,0.)
elseif(Which_fluxes(iSOC).eq.4) then
!Meta Model
     f(:) = ff_new(nz,:)
     call Meta_SOC(f,T(nz),S(nz),dz(nz),dT)
!Update state variables with flux values
      ff_new(nz,iO2) =  f(iO2)
      ff_new(nz,iNH4) = f(iNH4)
      ff_new(nz,iNO3) = f(iNO3)
      ff_new(nz,iOM1_A) = f(iOM1_A)
      ff_new(nz,iOM1_Z) = f(iOM1_Z)
      ff_new(nz,iOM1_R) = f(iOM1_R)
      ff_new(nz,iOM1_BC) =f(iOM1_BC)
endif


if(Which_fluxes(iNutEx).eq.1) then
!NO3 Exchange
       NO3Flux = 0.0057*ff_new(nz,iO2) - 0.52
               ff_new(nz,iNO3) = AMAX1(ff_new(nz,iNO3) + NO3Flux/ &
     & dz(nz)*dT/SDay,fmin(iNO3))

!NH4 Exchange
       NH4Flux = -1.55*NO3Flux + 0.69
               ff_new(nz,iNH4) = AMAX1(ff_new(nz,iNH4) + NH4Flux/ &
     & dz(nz)*dT/SDay,fmin(iNH4))

!PO4 Exchange
      PO4Flux = 0.094*NH4Flux - 0.0125
               ff_new(nz,iPO4) = AMAX1(ff_new(nz,iPO4) + PO4Flux/ &
     & dz(nz)*dT/SDay,fmin(iPO4))

!Si Exchange
      SiFlux = 1.68 
               ff_new(nz,iSi)  = AMAX1(ff_new(nz,iSi)  + SiFlux/ &
     & dz(nz)*dT/SDay,fmin(iSi))

!ALK Exchange
      ALKFlux = NO3Flux - NH4Flux + PO4Flux
               ff_new(nz,iALK)  = AMAX1(ff_new(nz,iALK)  + ALKFlux/ &
     & dz(nz)*dT/SDay,fmin(iALK))
endif


!MPB O2 Production
if(Which_fluxes(iMPB).eq.1) then
! Gatusso et al. 2006
               ff_new(nz,iO2) = ff_new(nz,iO2) + 120.82*(1.-exp(-convert_I*Esed/2.09))/ &
     & dz(nz)*dT/SDay
elseif(Which_fluxes(iMPB).eq.2) then
! Jahnke et al. 2008
               ff_new(nz,iO2) = ff_new(nz,iO2) + 132./12.*convert_I*Esed**(1.45)/ &
     & dz(nz)*dT/SDay
elseif(Which_fluxes(iMPB).eq.3) then
! Lehrter et al. (2014)
               ff_new(nz,iO2) = ff_new(nz,iO2) + 0.33*convert_I*Esed**(2.93)/ &
     & dz(nz)*dT/SDay
endif


 if (Which_Fluxes(iSDM) .eq. 1) then

         write(6,*) "SDM is in the code and compiles, but it was never checked, and not expected to work"
         write(6,*) "stopping"
         stop

 ! Sediment Diagenesis Model
 !        if(init.eq.1.or.mod(istep,288).eq.0) then  !Call every day, every 288 timesteps, assumes timestep = 5 min
 !           call Sediment_Diagenesis_Flux(dT*288,f(i,j,nz,:),T(i,j,nz),S(i,j,nz),pH(i,j,nz),sedflux(i,j,:),s_x1A(i,j,nz),&
     if (istep == 1 .OR. mod(istep,nstep_sed) == 0) then
         call Sediment_Diagenesis_Flux(dT_sed, ff_new(nz,:), T(nz), S(nz), pH(nz), sedflux(:), &
                                     & YY_init(:), pph_init(:) )
     endif
 
 !DIC Exchange
                ff_new(nz,iDIC) = AMAX1(ff_new(nz,iDIC) - sedflux(sDIC)/dz(nz)*dT/SDay, 0.)
 
 !NH4 Exchange
                ff_new(nz,iNH4) = AMAX1(ff_new(nz,iNH4) - sedflux(sNH4)/dz(nz)*dT/SDay, 0.)
 !NO3 Exchange
                ff_new(nz,iNO3) = AMAX1(ff_new(nz,iNO3) - sedflux(sNO3)/dz(nz)*dT/SDay, 0.)
 !O2 Exchange
                ff_new(nz,iO2)  = AMAX1(ff_new(nz,iO2)  - sedflux(sO2)/dz(nz)*dT/SDay, 0.)
 
 !DOC Exchange
        sedflux_iOM2_A = ff_new(nz,iOM2_A) / (ff_new(nz,iOM2_A) + ff_new(nz,iOM2_Z) + ff_new(nz,iOM2_R) + ff_new(nz,iOM2_bc))&
      &   * sedflux(sDOC)/dz(nz)*dT/SDay
        ff_new(nz,iOM2_A) = AMAX1(ff_new(nz,iOM2_A)-sedflux_iOM2_A, 0.)
 
        sedflux_iOM2_Z = ff_new(nz,iOM2_Z) / (ff_new(nz,iOM2_A) + ff_new(nz,iOM2_Z) + ff_new(nz,iOM2_R) + ff_new(nz,iOM2_bc)) &
      &   * sedflux(sDOC)/dz(nz)*dT/SDay
        ff_new(nz,iOM2_Z) = AMAX1(ff_new(nz,iOM2_Z)-sedflux_iOM2_Z, 0.)
 
        sedflux_iOM2_R = ff_new(nz,iOM2_R) / (ff_new(nz,iOM2_A) + ff_new(nz,iOM2_Z) + ff_new(nz,iOM2_R) + ff_new(nz,iOM2_bc)) &
      &   * sedflux(sDOC)/dz(nz)*dT/SDay
        ff_new(nz,iOM2_R) = AMAX1(ff_new(nz,iOM2_R)-sedflux_iOM2_R, 0.)
 
        sedflux_iOM2_BC = ff_new(nz,iOM2_bc) / (ff_new(nz,iOM2_A) + ff_new(nz,iOM2_Z) + ff_new(nz,iOM2_R) + & 
            ff_new(nz,iOM2_bc)) * sedflux(sDOC)/dz(nz)*dT/SDay
        ff_new(nz,iOM2_bc) = AMAX1(ff_new(nz,iOM2_bc)-sedflux_iOM2_bc, 0.)
 
 !OM1 Exchange for OM1_A and OM1_Z
        sedflux_iOM1_A = ff_new(nz,iOM1_A) / (ff_new(nz,iOM1_A) + ff_new(nz,iOM1_Z)) * sedflux(sOM1)/dz(nz)*dT/SDay
        ff_new(nz,iOM1_A) = AMAX1(ff_new(nz,iOM1_A)-sedflux_iOM1_A, 0.)
 
        sedflux_iOM1_Z = ff_new(nz,iOM1_Z) / (ff_new(nz,iOM1_A) + ff_new(nz,iOM1_Z)) * sedflux(sOM1)/dz(nz)*dT/SDay
        ff_new(nz,iOM1_Z) = AMAX1(ff_new(nz,iOM1_Z)-sedflux_iOM1_Z, 0.)
 
 !OM2 Exchange for OM1_R and OM1_bc 
        sedflux_iOM1_R = ff_new(nz,iOM1_R) / (ff_new(nz,iOM1_R) + ff_new(nz,iOM1_bc)) * sedflux(sOM2)/dz(nz)*dT/SDay
        ff_new(nz,iOM1_R) = AMAX1(ff_new(nz,iOM1_R)-sedflux_iOM1_R, 0.)
 
        sedflux_iOM1_bc = ff_new(nz,iOM1_bc) / (ff_new(nz,iOM1_R) + ff_new(nz,iOM1_bc)) * sedflux(sOM2)/dz(nz)*dT/SDay
        ff_new(nz,iOM1_bc) = AMAX1(ff_new(nz,iOM1_bc)-sedflux_iOM1_bc, 0.)
 
 !ALK Exchange
                ff_new(nz,iALK) = AMAX1(ff_new(nz,iALK) - sedflux(sALK)/dz(nz)*dT/SDay, 0.)
 
 endif



RETURN

end subroutine cgem_flux 
