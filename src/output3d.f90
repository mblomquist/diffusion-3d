! output3d Subroutine for 3D Diffusion Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-19 (YYYY-MM-DD)
!

subroutine output3d

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! Write field data to file
  open(unit=2, file="output/output3d_results.dat")

  write(2, *), 'TITLE = "Example: Simple 3D Volume Data"'
  write(2, *), 'VARIABLES = "X", "Y", "Z", "Temperature"'
  write(2, *), 'ZONE I=10, J=10, K=10, DATAPACKING=POINT'

  do k = 1,l
    do j = 1,n
      do i = 1,m

        write(2, *), length*i*dx, length*j*dy, length*k*dz, T(i,j,k)

      end do
    end do
  end do

  close(2)

  return

end subroutine
