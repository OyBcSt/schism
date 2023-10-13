subroutine grid_setup(nvrts)

  use grid
  implicit none

  integer, intent(in) :: nvrts
  real :: lat_in,lon_in

  km=nvrts-1

  call grid_read(lat_in,lon_in)
  call grid_allocate
  call grid_init(lat_in,lon_in)

return
end subroutine

