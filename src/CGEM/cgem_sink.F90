subroutine cgem_sink(dT,area,dowrite)

!This is called after cgem_step and cgem_flux, which returns ff_new
use grid, only:km,Vol
use cgem, only:ws,ff_new,nf,fmin,adjust_ws

implicit none

logical, intent(in) :: dowrite
real, intent(in) :: dT
real, dimension(km), intent(in) :: area
integer :: km1
integer :: isp,k
real :: cmin,x
real, dimension(km) :: C
real :: wsc
real, dimension(km) :: mass_in, mass_out

  km1 = km-1

  !if(adjust_ws) then
     !Courant condition
     ! ws*dT/dz < 1 (use .1)
     ! wsc < dz/dT
  !   cmin = 0.1*MINVAL(dz)/dT
  !endif

  do isp = 1,nf

   if(ws(isp).gt.tiny(x)) then

   ! !If sinking rate violates courant condition, change it
   ! if(adjust_ws) then
   !   wsc = AMIN1(cmin,ws(isp)) 
   ! else  !or not
      wsc = ws(isp)
   ! endif

    C(:) = ff_new(:,isp)

    mass_in(1) = 0.
    !mass_in(2:km) = C(1:km1)*area(1:km1)*wsc*dT
    mass_out(1:km1) = C(2:km)*area(2:km)*wsc*dT
    mass_out(km) = 0.
    !Don't sink out more than you started with
    mass_out(1:km1) = AMIN1(mass_out(1:km1),C(1:km1)*Vol(1:km1))
    !mass in should be mass out
    mass_in(2:km) = mass_out(1:km)
    ff_new(:,isp) = ff_new(:,isp) + (mass_in(:) - mass_out(:))/Vol(:)

    do k=1,km
     if(ff_new(k,isp).lt.0) then
      ff_new(k,isp) = 0. 
      !write(6,'(*(g0,:,", "))') "k,isp,fold,fnew,wsc,mass_in,mass_out,m2/vol=",k,isp,C(k),ff_new(k,isp),wsc,mass_in(k),mass_out(k),(mass_in(k) - mass_out(k))/Vol(k)
     endif
    enddo

   endif !end sink state variable

  enddo

  return

end subroutine cgem_sink
