subroutine cgem_sink(dT,i,istep,myrank)

!This is called after cgem_step and cgem_flux, which returns ff_new
use grid, only:km,dz
use cgem, only:ws,ff_new,nf,adjust_ws,adjust_fac
use schism_glbl, only : rkind,area,ze,nvrt,kbe

implicit none

integer, intent(in) :: i,istep,myrank
real(rkind), intent(in) :: dT
integer :: km1
integer :: isp,k
real(rkind)   :: x
real(rkind) :: cmin
real(rkind) :: wsc
real(rkind), dimension(km) :: C,mass_in,mass_out,d_mass

km1 = km-1

  do isp = 1,nf

   if(ws(isp).gt.tiny(x)) then

   ! !If sinking rate violates courant condition, change it
    if(adjust_ws) then
      cmin = adjust_fac*MINVAL(dz)*86.4
      wsc = DMIN1(cmin,ws(isp))
    else  !or not
      wsc = ws(isp)
    endif

      C(:) = ff_new(:,isp)
      wsc = wsc/86400.d0

     km1 = km-1 
     !Nothing sinking in
     mass_in(1) = 0.d0
     !area is constant, cancels
     mass_out(1:km1) = C(1:km1)*wsc*dz(1:km1)
     !Don't sink out
     mass_out(km) = 0.d0
     !mass in should be mass out
     mass_in(2:km) = mass_out(1:km1)
     d_mass = mass_in - mass_out
     ff_new(:,isp) = C(:) + d_mass(:)/dz(:)*dt

!    if(isp.eq.1.and.i.eq.10.and.myrank.eq.1) then
!            write(6,'(*(g0,:,", "))') istep,i,isp,SUM(C),SUM(ff_new(:,isp))
            !do k=1,km
            !write(6,'(*(g0,:,", "))') C(k),ff_new(k,isp),mass_in(k),mass_out(k),d_mass(k),dz(k)
            !enddo
!    endif

!    if(ABS(SUM(C*dz)-SUM(ff_new(:,isp)*dz(:))).gt.1.e-8) then
!        write(6,'(*(g0,:,", "))') istep,isp,SUM(ff_new(:,isp)),SUM(C),SUM(C*dz),SUM(ff_new(:,isp)*dz(:)),MINVAL(dz)
!    endif

    do k=1,km
     if(ff_new(k,isp).le.0) then
       write(6,'(*(g0,:,", "))') "k,isp,fold,fnew,wsc,dz=",k,isp,C(k),ff_new(k,isp),wsc,dz(k)
        ff_new(k,isp) = 0.
      stop 
     endif
    enddo

   endif !end sink state variable

  enddo

  return

end subroutine cgem_sink
