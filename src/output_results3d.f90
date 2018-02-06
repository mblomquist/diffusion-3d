! output_results subroutine
!
! Written by Matt Blomquist
! Last Update: 2018-01-24 (YYYY-MM-DD)
!
! This program solves one-, two-, and three-dimensional diffusions problems
! using the BiCG algorithm.

subroutine output_results2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  open(unit=2, file="output/output_results.txt")

  do j = 1, n
    do i = 1, m
      write(2,*) "Node: ", i, j, "Temperature: ", T(i,j)
    end do
  end do

  close(2)

  open(unit=3, file="output/output_results_r_norm.txt")

!  do i = 1,maxit
!    write(3,*) "Iteration:", i, "BiCG Residual: ", residual(i,1), "BiCGStab Residual: ", residual(i,2)
!  end do

  close(3)

  return

end subroutine output_results2d
