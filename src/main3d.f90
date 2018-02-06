! Diffusion Solver
!
! Written by Matt Blomquist
! Last Update: 2018-02-02 (YYYY-MM-DD)
!
! This program solves two-dimensional diffusions problems using the
! BiCG algorithm.

program main2d

  ! Include standard variable header
  include "var2d.dec"

  ! Initialize boundary conditions
  print *, "Initialize 2D Diffusion Problem"
  call initialize2d

  ! Output problem setup datafile
  print *, "Output Setup Datafile"
  call output_setup2d

  ! Update source terms
  print *, "Update source terms"
  call diffusion_source2d

  ! Set residual vector
  residual = 0

  ! Presolver Notice
  print *, "Problem size (i,j): ", m, n
  print *, "....................."

  T = 0

  ! Run BiCG solver
  print *, "Run BiCGStab Solver"
  call cpu_time(start_time)
  call solver2d_bicgstab(As, Aw, Ap, Ae, An, b, T, m, n, tol, maxit)
  call cpu_time(finish_time)

  print *, "Total runtime:", finish_time-start_time, "seconds."
  print *, "....................."

  T = 0

  ! Run BiCG solver
  print *, "Run Line by Line TDMA Solver"
  call cpu_time(start_time)
  call solver2d_tdma(As, Aw, Ap, Ae, An, b, T, m, n, tol, maxit)
  call cpu_time(finish_time)

  print *, "Total runtime:", finish_time-start_time, "seconds."
  print *, "....................."

  ! Output results datafile
  print *, "Output Results datafile"
  call output_results2d

end program main2d
