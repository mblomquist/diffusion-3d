! Diffusion Solver
!
! Written by Matt Blomquist
! Last Update: 2018-01-21 (YYYY-MM-DD)
!
! This program solves one-, two-, and three-dimensional diffusions problems
! using the BiCG algorithm.

subroutine output_setup2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  open(unit=2, file="output/output_setup.txt")

  ! Add write statements here...

  close(2)

  return

end subroutine output_setup2d
