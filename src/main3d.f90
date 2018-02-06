! 3D Diffusion Solver
!
! Written by Matt Blomquist
! Last Update: 2018-02-06 (YYYY-MM-DD)
!
! This program solves three-dimensional diffusions problems using the
! BiCG algorithm and Line-by-Line TDMA.

program main3d

  ! Include standard variable header
  include "var3d.dec"

  ! Initialize boundary conditions
  print *, "Initialize 3D Diffusion Problem"
  call initialize3d

  ! Output problem setup datafile
  print *, "Output Setup Datafile"
  call output_setup3d

  ! Update source terms
  print *, "Update source terms"
  call diffusion_source3d

  ! Set residual vector
  residual = 0

  ! Presolver Notice
  print *, "Problem size (i,j): ", m, n
  print *, "....................."

  T = 0

  ! Run BiCG solver
  print *, "Run BiCGStab Solver"
  call cpu_time(start_time)
  call solver3d_bicgstab(Ab, As, Aw, Ap, Ae, An, At, b, T, m, n, l, tol, maxit)
  call cpu_time(finish_time)

  print *, "Total runtime:", finish_time-start_time, "seconds."
  print *, "....................."

  T = 0

  ! Run BiCG solver
  print *, "Run Line by Line TDMA Solver"
  call cpu_time(start_time)
  call solver23_tdma(Ab, As, Aw, Ap, Ae, An, At, b, T, m, n, l, tol, maxit)
  call cpu_time(finish_time)

  print *, "Total runtime:", finish_time-start_time, "seconds."
  print *, "....................."

  ! Output results datafile
  print *, "Output Results datafile"
  call output_results3d

end program main3d
