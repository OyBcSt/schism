module grid

!CGEM STATE VARIABLES
use date_time
use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use schism_glbl, only: nea,rkind

implicit none

save

!Grid parameters
integer :: km

!Although nospA, nospZ seem like 'cgem' parameters, they are needed to allocate
!schism arrays  (to calculate number of tracers)

!Simulation parameters
integer :: START_SECONDS

!--Run Specifics---------------
integer iYrS,iMonS,iDayS,iHrS,iMinS,iSecS

!This should all be from SCHISM
!Solar radiation and Wind
real :: Rad, Wind
!Lat/lon of elements
real :: lat,lon
!Depth stuff
real(rkind), allocatable :: d(:) ! depth cell bottom to surface
real(rkind), allocatable :: d_sfc(:) ! depth cell center to surface 
real(rkind), allocatable :: dz(:)   !Thickness of cell
!'tracers'
real, allocatable :: S(:),T(:)


contains


subroutine grid_read(lat_in,lon_in)
  real, intent(out) :: lat_in,lon_in
  integer           :: istat,iunit
  character(len=1000) :: line
  !http://degenerateconic.com/namelist-error-checking.html
  namelist /hydro/ lon_in,lat_in 

#ifdef DEBUG
write(6,*) "Begin grid_init"
#endif

  open(action='read',file='grid.nml',iostat=istat,newunit=iunit)

  !namelist /hydro/
  read(nml=hydro,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
          'Invalid line in namelist hydro: '//trim(line)
  endif

  close(iunit)

return
end subroutine grid_read


subroutine grid_allocate()
integer ierr

#ifdef DEBUG
write(6,*) "Begin grid_allocate" 
#endif
!Depth stuff
allocate(d(km),stat=ierr) ! depth cell bottom to surface  
if(ierr.ne.0) write(6,*) "error in allocating:d"
allocate(d_sfc(km),stat=ierr) ! depth cell center to surface  
if(ierr.ne.0) write(6,*) "error in allocating:d_sfc"
allocate(dz(km),stat=ierr) !Thickness of cell
if(ierr.ne.0) write(6,*) "error in allocating:dz"
!SCHISM 'tracers'
allocate(S(km),stat=ierr) ! Salinity (psu) 
if(ierr.ne.0) write(6,*) "error in allocating:S"
allocate(T(km),stat=ierr) ! Temperature (C) 
if(ierr.ne.0) write(6,*) "error in allocating:T"

#ifdef DEBUG
write(6,*) "End grid_allocate"
#endif

return
end subroutine grid_allocate

subroutine grid_init(lat_in,lon_in) 

  integer :: k

  real, intent(in) :: lat_in,lon_in


  ! Compute starting time of run in seconds since Model_dim::iYrS:
  START_SECONDS = &
  TOTAL_SECONDS( iYrS, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )
#ifdef DEBUG
write(6,*) "In ReadInput"
write(6,*) "StartSeconds",START_SECONDS
#endif

  !Set lat/lon
  lat = lat_in
  lon = lon_in
  !Radiation
  !Rad = Rad_in
  call getSolar( iYrS, START_SECONDS, lon, lat, Rad)

#ifdef DEBUG
write(6,*) "lat,lon,Rad",lat,lon,Rad
#endif


end subroutine grid_init



!----------------------------------------------------------------------
      Subroutine getSolar( iYrS, TC_8, lon, lat, Rad) 
!----------------------------------------------------------------------
!     Written by  ::  D.S.Ko/NRL
!
!     Calculate visible solar radiation irradiance
!                        using just the solar zenith and ignore
!                        the effects of clouds and atmospheric
!                        absorption and scattering.
! ---------------------------------------------------------------------
     USE DATE_TIME

!----------------------
! Interface variables:
!----------------------
      real   , intent(in)  :: lon  ! longitude (deg E) at center of cell 
      real   , intent(in)  :: lat  ! latitude (deg N) at center of cell 
      integer, intent(in)  :: iYrS
      integer, intent(in) :: TC_8 ! Current time in seconds since Model_dim::iYrS
      real   , intent(out) :: Rad ! Solar Radiation
      integer  :: iYr      ! Year that Time_8 corresponds to
      integer  :: iMon     ! Month that Time_8 corresponds to
      integer  :: iDay     ! Day that Time_8 corresponds to
      integer  :: iHr      ! Hour that Time_8 corresponds to
      integer  :: iMin     ! Minute that Time_8 corresponds to
      integer  :: iSec     ! Second that Time_8 corresponds to
 
!-----------------------
! Local variables 
!-----------------------
      real, parameter :: cv        = 2.77e14 ! multiplicative factor used
                                                 ! to convert from watts/m2 
                                                 ! to photons/cm2/sec
                                                 ! Morel and Smith (1974)
      real, parameter :: OneD60   =  1./60.
      real, parameter :: OneD3600 =  1./3600.              
      integer  :: jday
      real                :: rhr         ! decimal hr in the Julian Day 
      real                :: Z   ! solar zenith angle 
      real calc_solar_zenith
      real                :: solconst      

         ! Note in next line that 1200 is the average clear-sky solar constant
         ! in watts/m^2
         solconst = 1200.00 * cv  ! in photons/cm2/sec
!-----------------------------------------------------------------------

         !Calculate Yr/Mon/Day from Model time in seconds	 
         call date_timestamp(iYrS,TC_8,iYr,iMon,iDay,iHr,iMin,iSec)

         !Hours in day
         rhr = real(iHr,4) + real(iMin,4)*OneD60 + real(iSec,4)*OneD3600

         !Day of the year
         jday = JDAY_IN_YEAR(iYr, iMon, iDay)

         Z =  calc_solar_zenith(lat,lon,rhr,jday) !in rad
         Rad = solconst * AMAX1( COS(Z), 0.0)    ! COS(suna)<= 0 means night 


        RETURN
           
      END Subroutine getSolar

end module grid
