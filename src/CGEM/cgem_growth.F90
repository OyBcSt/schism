module cgem_growth

use cgem, only: nospA,nospZ,km,umax,QminN,QmaxN,QminP,QmaxP,nfQs,   &
  & respg,respb,alphad,betad,Tref,KTg1,KTg2,Ea,KQn,KQp,KSi,Qc, &
  & which_growth,which_photosynthesis,which_quota,which_temperature,     &
  & which_uptake,CChla

implicit none

contains

! ------------------------------------------------------------------------
subroutine calc_Agrow( E, T_k, Qn, Qp, Si, A_k, Agrow_k, &
  &              uA_k, Aresp_k, uN_k, uP_k, uE_k, uSi_k )       
! ------------------------------------------------------------------------
! Call subroutine calc_Agrow to execute the desired phytoplankton 
! growth model to calculate the 1D array (water column) Agrow_k 
!-----------------------------------------------------------------------
! -- Declare input variables coming thru the interface ---------------------
  real,intent(in)  ::  E(km)           ! Irradiance (quanta/cm2/sec) 
                                       !   at middle of layer k
  real,intent(in)  ::  T_k(km)         ! Water temperature in Celsius
  real,intent(in)  ::  Qn(nospA,km)    ! Phytoplankton Nitrogen Quota (mmol-N/cell)         
  real,intent(in)  ::  Qp(nospA,km)    ! Phytoplankton Phosphorous Quota (mmol-P/cell)     
  real,intent(in)  ::  Si(km)          ! Silica (mmol-Si/m3)
  real,intent(in)  ::  A_k(nospA,km)   ! Number density of phytoplankton group isp 
! -- Declare calculated variables being returned ---------------------
  real,intent(out) ::  Agrow_k(nospA,km)  ! Specific growth rate    
                                          !   of phytoplankton group isp
  real,intent(out) ::  uA_k(km,nospA)     ! Temperature adjusted light factor
                                          !   phytoplankton group isp
  real,intent(out) ::  Aresp_k(nospA,km)  ! Phytoplankton respiration of group       	
                                          !   isp, including dark respiration. 
  real,intent(out) ::  uN_k(km,nospA)     ! Nitrogen limited growth rate (1/d)
  real,intent(out) ::  uP_k(km,nospA)     ! Phosphorus limited growth rate (1/d)
  real,intent(out) ::  uE_k(km,nospA)     ! Light limited growth rate (1/d)
  real,intent(out) ::  uSi_k(km,nospA)    ! Silica limited growth rate (1/d)

! -- Local variables --------------------------------------------------------------   
  integer :: k, isp ! loop indices     
  real,dimension(nospA+nospZ) :: Tadj            ! Temperature adjustment factor, variable and function 
  real,dimension(nospA)       :: uA              ! Specific growth, 1/d      
  real,dimension(nospA)       :: f_E             ! Light growth function 
  real,dimension(nospA)       :: f_N, f_P, f_Si  ! Nutrient growth functions
  real,dimension(nospA)       :: min_S           ! Limiting substrate values
  real,dimension(nospA)       :: respg2          ! Actual respiration coefficient
!------------------------------------------------------------------------

#ifdef DEBUG
write(6,*) "In calc_Agrow: Begin calc_Agrow"
#endif

!-------------------------------
! Begin growth rate calculations
!-------------------------------

  do k = 1, km

    call func_T( T_k(k), Tadj ) ! Temperature adjustment
    call func_S( Qn(:,k), Qp(:,k), Si(k), f_N, f_P, f_Si ) ! Nutrient dependent growth function
    do isp = 1, nospA
      min_S(isp) = AMIN1( f_N(isp), f_P(isp), f_Si(isp) )
    enddo
    call func_E( E(k), min_S, f_E ) ! Light growth function

    !Output variables for netCDF to examine light vs. nutrient limitations 
    uN_k(k,:)   = f_N(:)  * umax(:) * Tadj(1:nospA) 
    uP_k(k,:)   = f_P(:)  * umax(:) * Tadj(1:nospA)
    uE_k(k,:)   = f_E(:)  * umax(:) * Tadj(1:nospA) 
    uSi_k(k,:)  = f_Si(:) * umax(:) * Tadj(1:nospA)

#ifdef DEBUG
  if(k.eq.1) write(6,*) "f_E,umax,uE_k", f_E,umax,uE_k(k,:)
#endif

    if(which_growth.eq.1) then
      do isp=1,nospA
        uA(isp) = umax(isp) * Tadj(isp) * AMIN1(min_S(isp),f_E(isp)) ! Minimum Formulation
      enddo
    else if(which_growth.eq.2) then
      uA(:) = umax(:) * Tadj(1:nospA) * f_E(:) * min_S(:)   ! Product Formulation
    else if(which_growth.eq.3) then
      uA(:) = umax(:) * Tadj(1:nospA) * f_E(:) * min_S(:) ! Nutrient dependence is in f_E
    else !Let default be Minimum Formulation
      do isp=1,nospA
        uA(isp) = umax(isp) * Tadj(isp) * AMIN1(min_S(isp),f_E(isp)) ! Minimum Formulation
      enddo
    endif

    uA_k(k,:)      = uA(:)           ! Save specific growth rate to array for netCDF, 1/d
    Agrow_k(:,k)   = A_k(:,k)*uA(:)  ! Phytoplankton growth, cells/m3/d

    ! If uA < 0.25d-1, set respiration to zero; Laws and Bannister(1980) 
    do isp=1,nospA
      if(uA(isp).lt.0.25) then
        respg2(isp) = 0.
      else
        respg2(isp) = respg(isp)
      endif
    enddo

    !-----------------------------------------      
    ! Calculate the total respiration Aresp
    !-----------------------------------------
    Aresp_k(:,k) =  Agrow_k(:,k) * respg2(:) &  ! Growth dependent respiration (loss of cells), cells/m3/d
      & + Tadj(1:nospA)  * respb(:) * A_k(:,k)  ! Basal respiration (loss of cells) , cells/m3/d

  enddo    
  
return 
end subroutine calc_Agrow
!-------------------------------------------------------------

!-------------------------------------------------------------
subroutine func_E( E, min_S, f_E )   
!-------------------------------------------------------------
! INPUT:  
!   E = Irradiance at cell k (quanta/cm**2/sec)
!   min_S = Minimum substrate value (of f_N, f_P, and Si)
!
! OUTPUT:
!   f_E = Dimensionless factor for light dependent phytoplankton growth rate 
!
! REFERENCES:
!   nospA = Number of phytoplankton groups
!   which_photosynthesis  ==  1 : With photoinhibition, Platt et al. (1980)
!                         ==  2 : Without photoinhibition
!                         ==  3 : Nutrient dependent, Flynn (2003)
!   alphad =
!   betad  =
!------------------------------------------------------------------------
! -- Declare input variables coming thru the interface ---------------------
  real, intent(IN) :: E    ! Irradiance (quanta/cm**2/sec)
  real, intent(IN),  dimension(nospA) :: min_S ! Function of rate limiting nutrient
! -- Declare calculated variables being returned ---------------------
  real, intent(out),dimension(nospA)  :: f_E   ! Growth rate factor (dimensionless) 

  if (which_photosynthesis.eq.1) then         !With photoinhibition 
    f_E(1:nospA) = ( 1.0 - exp(-alphad(1:nospA) * E) ) * exp(-betad(1:nospA)*E)
  else if (which_photosynthesis.eq.2) then    !Without photoinhibition
    f_E(1:nospA) = ( 1.0 - exp(-alphad(1:nospA) * E) )
  else if (which_photosynthesis.eq.3) then    !Nutrient dependent
    f_E(1:nospA) = ( 1.0 - exp(-alphad(1:nospA) * E / min_S) )
  else
    write(6,*) "Error in func_E"
      stop
  endif
 
return
end subroutine func_E  
!------------------------------------------------------------

!------------------------------------------------------------
subroutine func_Qs( Qn, Qp, f_Qn, f_Qp)
!-- func_Qs is for a function of substrate 'S' --------------
!--------------------------------------------------------------------------
! INPUT:  
!   Qn - Phytoplankton Nitrogen Quota (mmol-N/cell) 
!   Qp - Phytoplankton Nitrogen Quota (mmol-P/cell)
!
! OUTPUT:
!   f_Qn  - Nitrogen dependent growth function  
!   f_Qp  - Phosphorus dependent growth function
!
! REFERENCES:
!   nospA = Number of phytoplankton groups
!   which_photosynthesis  ==  1 : With photoinhibition, Platt et al. (1980)
!                         ==  2 : Without photoinhibition
!                         ==  3 : Nutrient dependent, Flynn (2003)
!   nospA,nfQs,QmaxN,QminN,QmaxP,QminP,which_uptake 
!   Qmin/max_X - minimum and maximum nutrient cell quota (mmol/cell)
!------------------------------------------------------------------------
! -- Declare input variables coming thru the interface ---------------------
  real, intent(IN), dimension(nospA)  :: Qn    ! Phytoplankton Nitrogen Quota (mmol-N/cell)
  real, intent(IN), dimension(nospA)  :: Qp    ! Phytoplankton Phosporus Quota (mmol-P/cell) 
! -- Declare calculated variables being returned ---------------------
  real, intent(out), dimension(nospA) :: f_Qn  ! Function based on N
  real, intent(out), dimension(nospA) :: f_Qp  ! Function based on P

  if (which_uptake.eq.1) then !Michaelis-Menten
    f_Qn(:) = 1. 
    f_Qp(:) = 1. 
  else if (which_uptake.eq.2) then !Geider(1998), Lehman(1975) is nfQs=1
    f_Qn(:) = ( (QmaxN(:) - Qn(:))/(QmaxN(:) - QminN(:)) ) ** nfQs(:)
    f_Qp(:) = ( (QmaxP(:) - Qp(:))/(QmaxP(:) - QminP(:)) ) ** nfQs(:)
  else if (which_uptake.eq.3) then !Flynn(2003)
    f_Qn(:) = QmaxN(:)/Qn(:) 
    f_Qp(:) = QmaxP(:)/Qp(:)
  else  
    write(6,*) "Error in func_Qs"
    stop
  endif

return
end subroutine func_Qs  
!------------------------------------------------------------

!------------------------------------------------------------
subroutine func_S( Qn, Qp, Si, f_N, f_P, f_Si)
!-- func_S is for a function of substrate 'S' --------------- 
!  USE cgem_vars, only: nospA,which_quota,QminN,QminP,QmaxN,QmaxP,&
!      KQn,KQp,KSi
      
!--------------------------------------------------------------------------
! INPUT:  
!   Qn - Phytoplankton Nitrogen Quota (mmol-N/cell) 
!   Qp - Phytoplankton Nitrogen Quota (mmol-P/cell)
!   Si - Silica concentration in seawater (mmol-Si/m3)
!   N  - NO3+NH4 concentration in seawater (mmol-N/m3)
!   P  - PO4 concentration in seawater (mmol-P/m3)
!
! OUTPUT:
!   f_N  - Nitrogen dependent growth function  
!   f_P  - Phosphorus dependent growth function
!   f_Si - Silica dependent growth function
!
!  USE cgem_vars, only: nospA,which_quota,QminN,QminP,QmaxN,QmaxP,&
!      is_diatom,KQn,KQp,KSi
!   K_X - half saturation constants for phytoplankton group  (X mmol/m3)
!   Qmin/max_X - minimum and maximum nutrient cell quota (mmol/cell)
! 
!------------------------------------------------------------------------

! -- Declare input variables coming thru the interface ---------------------
  real, intent(in), dimension(nospA)  :: Qn    ! Phytoplankton Nitrogen Quota (mmol-N/cell)
  real, intent(in), dimension(nospA)  :: Qp    ! Phytoplankton Phosporus Quota (mmol-P/cell) 
  real, intent(in)                    :: Si    ! Silica concentration in seawater (mmol-Si/m3)
! -- Declare calculated variables being returned ---------------------
  real, intent(out), dimension(nospA) :: f_N   ! Function based on N
  real, intent(out), dimension(nospA) :: f_P   ! Function based on P
  real, intent(out), dimension(nospA) :: f_Si  ! Function based on Si

  if (which_quota.eq.1) then !Droop(1968)
    f_N(:) = ( Qn(:) - QminN(:) ) / Qn(:)
    f_P(:) = ( Qp(:) - QminP(:) ) / Qp(:)
    f_Si(:) = Si / ( Si + Ksi(:) ) !Monod        
  else if (which_quota.eq.2) then !Nyholm(1978)
    f_N(:) = ( Qn(:) - QminN(:) ) / ( QmaxN(:) - QminN(:) )
    f_P(:) = ( Qp(:) - QminP(:) ) / ( QmaxP(:) - QminP(:) )
    f_Si(:) = Si / ( Si + Ksi(:) ) !Monod 
  else if (which_quota.eq.3) then !Flynn(2003)
    f_N(:) = ( 1. + KQn(:) ) * ( Qn(:) - QminN(:) ) /        &
     &           ( Qn(:) - QminN(:) + KQn(:)*( QmaxN(:) - QminN(:) ) )
    f_P(:) = ( 1. + KQp(:) ) * ( Qp(:) - QminP(:) ) /        &
     &           ( Qp(:) - QminP(:) + KQp(:)*( QmaxP(:) - QminP(:) ) )
    f_Si(:) = Si / ( Si + Ksi(:) ) !Monod 
  else
    write(6,*) "Error in func_S. which_quota=",which_quota
    stop
  endif

return
end subroutine func_S  
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
subroutine func_T( T, Tadj )   
!---------------------------------------------------------------------------

!--------------------------------------------------------------------------
! INPUT:  
!   T = temperature [degree C]
!
! OUTPUT:
!   Temperature Adjustment 
! 
! REFERENCES:
!
! nospA,nospZ,Tref,KTg1,KTg2,which_temperature,Ea
!------------------------------------------------------------------------
! -- Declare input variables coming thru the interface ------------------
  real, intent(in) :: T    ! Temperature (deg C)
! -- Declare calculated variables being returned ---------------------
  real, intent(out), dimension(nospA+nospZ) :: Tadj 
! -- Local variables ------------------------------------------------------   
  real, parameter  :: f0    = 0.1
  real, parameter  :: r     = 0.3
  real, parameter  :: f1    = 1.0/f0 - 1.0  
  real, parameter  :: r1    = r*(46.5/18.0) 
  real, parameter  :: k_b    = 8.6173303e-5 !Boltzmann constant in eV/K
  real             :: denom(nospA+nospZ)
  real             :: T_in_K,Tref_in_K(nospA+nospZ)
  integer          :: isp

#ifdef DEBUG
write(6,*) "func_T, nospA, which_temperature",nospA,which_temperature
#endif

  if (which_temperature.eq.1) then !Sigmoidal 
    denom(:) = 1.0 + f1*exp(-r1*( T - Tref(:))) 
    Tadj(:)  = 0.3 *(1.0/denom(:)) + 0.7            
  else if (which_temperature.eq.2) then !Optimum temperature threshold T (Cerco and Noel, 2004)
    do isp=1,nospA+nospZ
      if(T.le.Tref(isp)) then
        Tadj(isp) = exp( -KTg1(isp) * (T - Tref(isp))**2 )           
      else
        Tadj(isp) = exp( -KTg2(isp) * (Tref(isp) - T)**2 )
      endif
    enddo
  else if (which_temperature.eq.3) then !Decrease in growth rate at threshold T (Arrhenius form, Geider 1997)
    T_in_K  = T + 273.15 !Temp. in Kelvin
    Tref_in_K(:) = Tref(:) + 273.15 !Temp. in Kelvin 
    Tadj(:) = exp ( -(Ea(:)/k_b) * ( 1./T_in_K - 1./Tref_in_K(:) ) ) 
  else  
    write(6,*) "Error in func_T"
    stop
  endif
 
return
end subroutine func_T  


! These subroutines contain code for calculating mg chlorophyll-a from phytoplankton abundance.  Cloern is removed, so there is only:
! 1. Fixed C:Chla ratio 

function Fixed_CChla(A_k) RESULT(Chla_tot)

  ! Input parameters
  real, intent(in) :: A_k(nospA,km)  ! A's number density, cells/m3
  ! Function return value
  real :: Chla_tot(km)
  
  ! Local variables
  integer :: k, isp

  do k = 1, km 
    Chla_tot(k) = 0.0
    do isp = 1, nospA
      Chla_tot(k) =  Chla_tot(k) + A_k(isp,k) * Qc(isp) * 12. * (1./CChla(isp)) 
    enddo ! isp = 1, nospA
  enddo ! k = 1, km 

return
end function Fixed_CChla 


end module cgem_growth
