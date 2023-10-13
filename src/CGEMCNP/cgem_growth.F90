module cgem_growth

use cgem, only: nospA,nospZ,writecsv,km,umax,QminN,QmaxN,QminP,QmaxP,nfQs,   &
  & respg,respb,alphad,betad,Tref,KTg1,KTg2,Ea,KQn,KQp,KSi,Qc, &
  & which_growth,which_photosynthesis,which_quota,which_temperature,     &
  & which_uptake,CChla

implicit none

contains

! ------------------------------------------------------------------------
subroutine calc_Agrow( E, T, Qn, Qp, Si, A, Agrow, &
  &              uA, Aresp, uN, uP, uE, uSi, inea )       
! ------------------------------------------------------------------------
! Call subroutine calc_Agrow to execute the desired phytoplankton 
! growth model to calculate the 1D array (water column) Agrow 
!-----------------------------------------------------------------------
! -- Declare input variables coming thru the interface ---------------------
  integer, intent(in) :: inea !Grid number, for writing
  real,intent(in)  ::  E(km)           ! Irradiance (quanta/cm2/sec) 
                                       !   at middle of layer k
  real,intent(in)  ::  T(km)         ! Water temperature in Celsius
  real,intent(in)  ::  Qn(km,nospA)    ! Phytoplankton Nitrogen Quota (mmol-N/cell)         
  real,intent(in)  ::  Qp(km,nospA)    ! Phytoplankton Phosphorous Quota (mmol-P/cell)     
  real,intent(in)  ::  Si(km)          ! Silica (mmol-Si/m3)
  real,intent(in)  ::  A(km,nospA)   ! Number density of phytoplankton group isp 
! -- Declare calculated variables being returned ---------------------
  real,intent(out) ::  Agrow(km,nospA)  ! Specific growth rate    
                                          !   of phytoplankton group isp
  real,intent(out) ::  uA(km,nospA)     ! Temperature adjusted light factor
                                          !   phytoplankton group isp
  real,intent(out) ::  Aresp(km,nospA)  ! Phytoplankton respiration of group       	
                                          !   isp, including dark respiration. 
  real,intent(out) ::  uN(km,nospA)     ! Nitrogen limited growth rate (1/d)
  real,intent(out) ::  uP(km,nospA)     ! Phosphorus limited growth rate (1/d)
  real,intent(out) ::  uE(km,nospA)     ! Light limited growth rate (1/d)
  real,intent(out) ::  uSi(km,nospA)    ! Silica limited growth rate (1/d)

! -- Local variables --------------------------------------------------------------   
  integer :: k, isp ! loop indices     
  real,dimension(km,nospA+nospZ) :: Tadj            ! Temperature adjustment factor, variable and function 
  real,dimension(km,nospA)       :: f_E             ! Light growth function 
  real,dimension(km,nospA)       :: f_N, f_P, f_Si  ! Nutrient growth functions
  real,dimension(km,nospA)       :: min_S           ! Limiting substrate values
  real,dimension(nospA)       :: respg2          ! Actual respiration coefficient
!------------------------------------------------------------------------

#ifdef DEBUG
write(6,*) "In calc_Agrow: Begin calc_Agrow"
#endif

!-------------------------------
! Begin growth rate calculations
!-------------------------------
    call func_S( Qn, Qp, Si, f_N, f_P, f_Si ) ! Nutrient dependent growth function
    call func_T( T, Tadj ) ! Temperature adjustment
    do isp = 1, nospA
       min_S(:,isp) = AMIN1( f_N(:,isp), f_P(:,isp), f_Si(:,isp) )
    enddo
    call func_E( E, min_S, f_E ) ! Light growth function

    !Output variables for netCDF to examine light vs. nutrient limitations 
    do isp = 1,nospA
      uN(:,isp)   = f_N(:,isp)  * umax(isp) * Tadj(:,isp) 
      uP(:,isp)   = f_P(:,isp)  * umax(isp) * Tadj(:,isp)
      uE(:,isp)   = f_E(:,isp)  * umax(isp) * Tadj(:,isp) 
      uSi(:,isp)  = f_Si(:,isp) * umax(isp) * Tadj(:,isp)
    enddo

    if(which_growth.eq.1) then
      do isp=1,nospA
        uA(:,isp) = umax(isp) * Tadj(:,isp) * AMIN1(min_S(:,isp),f_E(:,isp)) ! Minimum Formulation
      enddo
    else if(which_growth.eq.2) then
      do isp=1,nospA
        uA(:,isp) = umax(isp) * Tadj(:,isp) * f_E(:,isp) * min_S(:,isp)   ! Product Formulation
      enddo
    else if(which_growth.eq.3) then
      do isp=1,nospA
        uA(:,isp) = umax(isp) * Tadj(:,isp) * f_E(:,isp) * min_S(:,isp) ! Nutrient dependence is in f_E
      enddo
    else !Let default be Minimum Formulation
      do isp=1,nospA
        uA(:,isp) = umax(isp) * Tadj(:,isp) * AMIN1(min_S(:,isp),f_E(:,isp)) ! Minimum Formulation
      enddo
    endif

    do isp=1,nospA
      Agrow(:,isp)   = A(:,isp)*uA(:,isp)  ! Phytoplankton growth, cells/m3/d
    enddo

    !-----------------------------------------      
    ! Calculate the total respiration Aresp
    !-----------------------------------------
    do k=1,km
      do isp=1,nospA
      ! If uA < 0.25d-1, set respiration to zero; Laws and Bannister(1980) 
        if(uA(k,isp).lt.0.25) then
          respg2(isp) = 0.
        else
          respg2(isp) = respg(isp)
        endif
        Aresp(k,isp) =  Agrow(k,isp) * respg2(isp) &  ! Growth dependent respiration (loss of cells), cells/m3/d
          & + Tadj(k,isp)  * respb(isp) * A(k,isp)  ! Basal respiration (loss of cells) , cells/m3/d
    enddo
    enddo

k=1
if(writecsv==1.and.k.eq.1.and.inea.eq.10) then
    write(6001,'(*(g0,:,", "))') Agrow(k,1),Aresp(k,1),uA(k,1),uN(k,1),uP(k,1),uE(k,1),uSi(k,1),f_E(k,1),f_N(k,1),f_P(k,1),f_Si(k,1),Tadj(k,1),A(k,1),min_S(k,1)
endif


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
  real, dimension(km), intent(in) :: E    ! Irradiance (quanta/cm**2/sec)
  real, dimension(km,nospA), intent(in) :: min_S ! Function of rate limiting nutrient
! -- Declare calculated variables being returned ---------------------
  real, intent(out),dimension(km,nospA)  :: f_E   ! Growth rate factor (dimensionless) 
  integer isp

  if (which_photosynthesis.eq.1) then         !With photoinhibition 
    do isp=1,nospA
      f_E(:,isp) = ( 1.0 - exp(-alphad(isp) * E(:)) ) * exp(-betad(isp)*E(:))
    enddo
  else if (which_photosynthesis.eq.2) then    !Without photoinhibition
    do isp=1,nospA
      f_E(:,isp) = ( 1.0 - exp(-alphad(isp) * E(:)) )
    enddo
  else if (which_photosynthesis.eq.3) then    !Nutrient dependent
    do isp=1,nospA
      f_E(:,isp) = ( 1.0 - exp(-alphad(isp) * E(:) / min_S(:,isp)) )
    enddo
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
  real, intent(in), dimension(km,nospA)  :: Qn    ! Phytoplankton Nitrogen Quota (mmol-N/cell)
  real, intent(in), dimension(km,nospA)  :: Qp    ! Phytoplankton Phosporus Quota (mmol-P/cell) 
! -- Declare calculated variables being returned ---------------------
  real, intent(out), dimension(km,nospA) :: f_Qn  ! Function based on N
  real, intent(out), dimension(km,nospA) :: f_Qp  ! Function based on P
  integer isp

  if (which_uptake.eq.1) then !Michaelis-Menten
    f_Qn = 1. 
    f_Qp = 1. 
  else if (which_uptake.eq.2) then !Geider(1998), Lehman(1975) is nfQs=1
    do isp=1,nospA
      f_Qn(:,isp) = ( (QmaxN(isp) - Qn(:,isp))/(QmaxN(isp) - QminN(isp)) ) ** nfQs(isp)
      f_Qp(:,isp) = ( (QmaxP(isp) - Qp(:,isp))/(QmaxP(isp) - QminP(isp)) ) ** nfQs(isp)
    enddo
  else if (which_uptake.eq.3) then !Flynn(2003)
    do isp=1,nospA
      f_Qn(:,isp) = QmaxN(isp)/Qn(:,isp) 
      f_Qp(:,isp) = QmaxP(isp)/Qp(:,isp)
    enddo
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
  real, intent(in), dimension(km,nospA)  :: Qn    ! Phytoplankton Nitrogen Quota (mmol-N/cell)
  real, intent(in), dimension(km,nospA)  :: Qp    ! Phytoplankton Phosporus Quota (mmol-P/cell) 
  real, intent(in), dimension(km)        :: Si    ! Silica concentration in seawater (mmol-Si/m3)
! -- Declare calculated variables being returned ---------------------
  real, intent(out), dimension(km,nospA) :: f_N   ! Function based on N
  real, intent(out), dimension(km,nospA) :: f_P   ! Function based on P
  real, intent(out), dimension(km,nospA) :: f_Si  ! Function based on Si
! -- local
  integer :: isp

  if (which_quota.eq.1) then !Droop(1968)
    do isp=1,nospA
      f_N(:,isp)  = ( Qn(:,isp) - QminN(isp) ) / Qn(:,isp)
      f_P(:,isp)  = ( Qp(:,isp) - QminP(isp) ) / Qp(:,isp)
      f_Si(:,isp) = Si(:) / ( Si(:) + Ksi(isp) ) !Monod        
    enddo
  else if (which_quota.eq.2) then !Nyholm(1978)
    do isp=1,nospA
      f_N(:,isp)  = ( Qn(:,isp) - QminN(isp) ) / ( QmaxN(isp) - QminN(isp) )
      f_P(:,isp)  = ( Qp(:,isp) - QminP(isp) ) / ( QmaxP(isp) - QminP(isp) )
      f_Si(:,isp) = Si(:) / ( Si(:) + Ksi(isp) ) !Monod 
    enddo
  else if (which_quota.eq.3) then !Flynn(2003)
    do isp=1,nospA
      f_N(:,isp)  = ( 1. + KQn(isp) ) * ( Qn(:,isp) - QminN(isp) ) /        &
       &           ( Qn(:,isp) - QminN(isp) + KQn(isp)*( QmaxN(isp) - QminN(isp) ) )
      f_P(:,isp)  = ( 1. + KQp(isp) ) * ( Qp(:,isp) - QminP(isp) ) /        &
       &           ( Qp(:,isp) - QminP(isp) + KQp(isp)*( QmaxP(isp) - QminP(isp) ) )
      f_Si(:,isp) = Si(:) / ( Si(:) + Ksi(isp) ) !Monod 
    enddo
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
  real, dimension(km), intent(in) :: T    ! Temperature (deg C)
! -- Declare calculated variables being returned ---------------------
  real, intent(out), dimension(km,nospA+nospZ) :: Tadj
! -- Local variables ------------------------------------------------------
  real, parameter  :: f0    = 0.1
  real, parameter  :: r     = 0.3
  real, parameter  :: f1    = 1.0/f0 - 1.0
  real, parameter  :: r1    = r*(46.5/18.0)
  real, parameter  :: k_b    = 8.6173303e-5 !Boltzmann constant in eV/K
  real, dimension(km,nospA+nospZ)           :: denom
  real, dimension(nospA+nospZ)              :: Tref_in_K
  integer          :: isp,k

#ifdef DEBUG
write(6,*) "func_T, nospA, which_temperature",nospA,which_temperature
#endif

  if (which_temperature.eq.1) then !Sigmoidal
    do isp=1,nospA+nospZ
      denom(:,isp) = 1.0 + f1*exp(-r1*( T(:) - Tref(isp)))
    enddo
    do isp=1,nospA+nospZ
      Tadj(:,isp)  = 0.3 *(1.0/denom(:,isp)) + 0.7
    enddo
  else if (which_temperature.eq.2) then !Optimum temperature threshold T (Cerco and Noel, 2004)
    do k=1,km
      do isp=1,nospA+nospZ
        if(T(k).le.Tref(isp)) then
          Tadj(k,isp) = exp( -KTg1(isp) * (T(k) - Tref(isp))**2 )
        else
          Tadj(k,isp) = exp( -KTg2(isp) * (Tref(isp) - T(k))**2 )
        endif
      enddo
    enddo
  else if (which_temperature.eq.3) then !Decrease in growth rate at threshold T (Arrhenius form, Geider 1997)
    Tref_in_K(:) = Tref(:) + 273.15 !Temp. in Kelvin
    do k=1,km
      Tadj(k,:) = exp ( -(Ea(:)/k_b) * ( 1./(T(k)+273.15) - 1./Tref_in_K(:) ) )
    enddo
  else
    write(6,*) "Error in func_T"
    stop
  endif

return
end subroutine func_T

! These subroutines contain code for calculating mg chlorophyll-a from phytoplankton abundance.  Cloern is removed, so there is only:
! 1. Fixed C:Chla ratio 

function Fixed_CChla(A) RESULT(Chla_tot)

  ! Input parameters
  real, intent(in) :: A(km,nospA)  ! A's number density, cells/m3
  ! Function return value
  real :: Chla_tot(km)
  
  ! Local variables
  integer :: k, isp

  do k = 1, km 
    Chla_tot(k) = 0.0
    do isp = 1, nospA
      Chla_tot(k) =  Chla_tot(k) + A(k,isp) * Qc(isp) * 12. * (1./CChla(isp)) 
    enddo ! isp = 1, nospA
  enddo ! k = 1, km 

return
end function Fixed_CChla 


end module cgem_growth
