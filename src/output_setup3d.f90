! output_setup3d
!
! Written by Matt Blomquist
! Last Update: 2018-02-06 (YYYY-MM-DD)
!
! This program outputs the problem setup data file for a 3D diffusion
! problem.

subroutine output_setup3d

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  open(unit=2, file="output/output_setup.txt")

  ! Add write statements here...

  close(2)

  return

end subroutine output_setup3d
