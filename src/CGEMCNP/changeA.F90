      Subroutine changeA(A,f,T,S,pH)

      use cgem, only: nf,iOM1_A,iOM2_A,iOM1_Z,iOM2_Z,iOM1_R,iOM2_R,iOM1_BC,iOM2_BC,&
      & isx1A,isy1A,isx2A,isy2A,isx1Z,isx2Z,isy1Z,isy2Z, &
      & sx1R,sy1R,sx2R,sy2R,sx1BC,sx2BC,sy1BC,sy2BC,iO2,iNH4,KO2,KNH4, &
      & KG1,KG1_R,KG1_BC,KG2,KG2_R,KG2_BC,nitmax,KNO3,iNO3,iDIC,ws
      use cgem_utils

      IMPLICIT NONE

      REAL(kind=8), dimension(100), intent(INOUT) :: A
      real, intent(IN) :: f(nf),T,S,pH ! State variables, Temp, Salinity
      real :: R_11


! Calculate nitrification, time units are in years, R_11:
      call Nitrification(f(iO2),f(iNH4),KO2,KNH4,nitmax,T,R_11)

! Change A values according to state variables, T, and S
      A(1) = T   ! Temperature in C
      A(2) = S   ! Salinity 
      A(4) = pH  ! pH

! Weighted average of KG1 for OM1_A and OM1_Z (KG1 is the same)
      A(8) = KG1
! Weighted average of KGX for OM1_R and OM1_
      if((f(iOM1_R) + f(iOM1_BC)).gt.0.) then
      A(9) = (KG1_R * f(iOM1_R) + KG1_BC * f(iOM1_BC)) / &
     &                         (f(iOM1_R) + f(iOM1_BC))
      else
      A(9) = KG1 
      endif

!DOM: Weighted Averages of K_OM for all OM2
      if((f(iOM2_A) + f(iOM2_Z) + f(iOM2_R) + f(iOM2_BC)).gt.0.) then
      A(10) = (KG2*(f(iOM2_A)+f(iOM2_Z)) + KG2_R*f(iOM2_R) + KG2_BC*f(iOM2_BC)) / & 
     &             (f(iOM2_A) + f(iOM2_Z) + f(iOM2_R) + f(iOM2_BC))
      else
      A(10) = KG2 
      endif

      A(16) = R_11 

! OM1=OM1_A + OM1_Z
      if((f(iOM1_A)+f(iOM1_Z)).gt.0.) then
      A(23) = (f(isx1A)/1.*f(iOM1_A) + f(isx1Z)/1.*f(iOM1_Z))/(f(iOM1_A)+f(iOM1_Z)) 
      A(24) = (f(isy1A)/1.*f(iOM1_A) + f(isy1Z)/1.*f(iOM1_Z))/(f(iOM1_A)+f(iOM1_Z)) 
      else
      A(23) = 105.
      A(24) = 25.
      endif

      A(25) = 1. 

! OM2=OM1_R + OM1_BC
      if((f(iOM1_R)+f(iOM1_BC)).gt.0.) then
      A(26) = (sx1R/1.*f(iOM1_R) + sx1BC/1.*f(iOM1_BC))/(f(iOM1_R)+f(iOM1_BC))
      A(27) = (sy1R/1.*f(iOM1_R) + sy1BC/1.*f(iOM1_BC))/(f(iOM1_R)+f(iOM1_BC))
      else
      A(26) = 105.
      A(27) = 25.
      endif

      A(28) = 1. 

! DOM = Sum of OM
      if((f(iOM1_A)+f(iOM1_Z)+f(iOM1_R)+f(iOM1_BC)).gt.0.) then
      A(29) = (f(isx2A)/1.*f(iOM2_A) + f(isx2Z)/1.*f(iOM2_Z) +     &
     &       sx2R/1.*f(iOM2_R) + sx2BC/1.*f(iOM2_BC)) / &
     &        (f(iOM2_A) + f(iOM2_Z) + f(iOM2_R) + f(iOM2_BC))
      A(30) = (f(isy2A)/1.*f(iOM2_A) + f(isy2Z)/1.*f(iOM2_Z) +     &
     &       sy2R/1.*f(iOM2_R) + sy2BC/1.*f(iOM2_BC)) / &
     &        (f(iOM2_A) + f(iOM2_Z) + f(iOM2_R) + f(iOM2_BC))
      else
      A(29) = 105.
      A(30) = 10.
      endif

      A(31) = 1. 

      A(32) = kO2
      A(33) = kNO3
      A(47) = f(iO2)
      A(48) = f(iNO3)
      A(49) = f(iNH4) 
      A(57) = f(iOM2_A) + f(iOM2_Z) + f(iOM2_R) + f(iOM2_BC)
      A(58) = f(iDIC)
      A(59) = A(57)
      A(60) = A(49) 
      A(62) = (f(iOM1_A)*ws(iOM1_A)+f(iOM1_Z)*ws(iOM1_Z))*12.
      A(63) = (f(iOM1_R)*ws(iOM1_R)+f(iOM1_BC)*ws(iOM1_BC))*12.


      RETURN

      END Subroutine changeA 
