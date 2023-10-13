!==============================================================================
!   CGEM
!   Info
!   Copyright
!   License
!==============================================================================
subroutine cgem_run(istep,myrank)
  use schism_glbl, only : rkind,nea,idry_e,irange_tr,flx_sf,flx_bt,bdy_frc,&
   & nvrt,kbe,dpe,tr_el,dt,srad,elnode,i34,windx,windy,area,ze,wsett
  use grid, only : T,S,km,dz,Vol,d,d_sfc,START_SECONDS,Wind,Rad
  use cgem, only: ws,ff,ff_new,skipcgem,checkwindrad,sinkwcgem

  implicit none

  integer, intent(in) :: istep,myrank
  integer :: itmp1,itmp2,i,m,im,mm,k,TC_8
  logical :: dowrite 
  real :: cgemdt
  real :: mindz
  real, dimension(km) :: cgemarea
  real, parameter :: cv        = 2.77e14 ! multiplicative factor used
                                             ! to convert from watts/m2 
                                             ! to photons/cm2/sec
                                             ! Morel and Smith (1974)
  external :: cgem_step, cgem_flux, cgem_sink
                                             
!Just say Hi in mirror.out
  if(myrank==0) write(16,*) "In cgem_run: istep,dt=",istep,dt

!dt does not pass through cgem_step function properly by itself
  cgemdt = dt

!Time in seconds since start of run
  TC_8 = START_SECONDS + istep*int(dt)

!Range of tracers, first 2 are S,T, then CGEM are next
  itmp1=irange_tr(1,3)
  itmp2=irange_tr(2,3)

!Loop over elements
  do i=1,nea

  !Skip if element is dry
    if(idry_e(i)==1) cycle

  !If not dry, run CGEM

    !Set surface and bottom flux, and body forces to zero
    flx_sf(itmp1:itmp2,i)=0.d0
    flx_bt(itmp1:itmp2,i)=0.d0
    bdy_frc(itmp1:itmp2,:,i)=0.d0

  !Set cgem state variable array, ff, to tracer variables (tr_el) that schism 
  !  has 'transported' for us. 
    mm = 1 
    do m=itmp1,itmp2
     im = km                    !cgem, k=1 is surface
     do k=kbe(i)+1,nvrt         !schism, k=1 is bottom
        ff(im,mm)=tr_el(m,k,i)  !I will avoid rewriting cgem until
       im = im-1                ! performance profiling forces the issue 
     enddo !k                  
     mm = mm+1
    enddo !m

  !Set temperature and salinity, depth and volume
    im = km
    do k=kbe(i)+1,nvrt
      T(im)= tr_el(1,k,i)
      S(im)= tr_el(2,k,i)
      dz(im) = ze(k,i)-ze(k-1,i)
      Vol(im) = (ze(k,i)-ze(k-1,i))*area(i)
      cgemarea(im) = area(i)
      im = im-1
    enddo

    d(1) = dz(1)          !d is distance from surface to bottom of cell
    d_sfc(1) = dz(1)/2.   !d_sfc is surface to mid cell 
                          !First cell is half of total thickness of first cell
    do k=2,km
      d(k) = d(k-1) + dz(k) 
      d_sfc(k) = d_sfc(k-1) + dz(k)
    enddo

  !Wind/Rad expressions copied directly from cosine.F90
    Wind=sqrt((sum(windx(elnode(1:i34(i),i)))/real(i34(i),rkind))**2.0+(sum(windy(elnode(1:i34(i),i)))/real(i34(i),rkind))**2.0)
  !Convert Rad from W/m2 to photons/cm2/sec
    Rad = max(sum(srad(elnode(1:i34(i),i)))/i34(i),0.d0)*cv

  !This will write Wind/Rad for every single timestep, don't run for long or
  !the resulting text file will be enormous.
  if(myrank.eq.1.and.i.eq.1.and.checkwindrad.eq.1) write(16,*) "Wind,Rad,Minutes",Wind,Rad/cv,istep,istep*int(dt)/60./60.

  dowrite=.FALSE.
  if(i.eq.10.and.myrank.eq.1) dowrite=.TRUE.

  !The option to skip cgem is for verifying initial and boundary conditions,
!  sinking, and loading without cgem complicating the process 
if(skipcgem.eq.1) then
  !Don't call cgem
   ff_new = ff
else
  !Call CGEM for a column
  !Input is ff, output is ff_new
  call cgem_step(TC_8,cgemdt,istep,i,myrank)
endif 
  
  if(sinkwcgem) then
    call cgem_sink(cgemdt,cgemarea,dowrite)
    do m=itmp1,itmp2
      wsett(m,:,i) = 0.
    enddo
  else
    !--Since this is constant, I'd rather define it somewhere else, but
    !- it was not registering properly.  I'll get back to this.
    !wsett is settling velocity
    mm = 1                 !cgem tracers are mm=1:nf
    do m=itmp1,itmp2       !schism's are m=itmp1:itmp2 
      wsett(m,:,i)= ws(mm) !ws(nf) is cgem sinking array
      mm = mm+1
    enddo
  endif


!Update schism tracer variables with newly calculated cgem variables
    mm = 1
    do m=itmp1,itmp2
     im = km
     do k=kbe(i)+1,nvrt
       tr_el(m,k,i)=ff_new(im,mm)
       im = im-1
     enddo !k
     mm = mm+1
    enddo !m


  enddo !i=1,nea


return
end subroutine cgem_run
