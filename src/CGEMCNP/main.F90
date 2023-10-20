  PROGRAM main

  implicit none

  call grid_setup(nvrt)
  call cgem_setup(ntrs(3))
  call cgem_run(it,myrank)

  return
  END 
