! 3D Diffusion Solver
! 
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This program solves three-dimensional diffusions problems.

program main3d

  implicit none

  ! Include standard variable header
  include "var3d.dec"

  ! Initialize Problem
  call initialize3d

  ! Set coefficient values
  call boundary3d
  call source3d

  ! Start Timer
  call cpu_time(start_time)

  ! Solve Diffusion Problem
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab, As, Aw, Ap, Ae, An, At, b, u_star, m, n, l, solver_tol, maxit)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab, As, Aw, Ap, Ae, An, At, b, u_star, m, n, l, solver_tol, maxit)
  elseif (solver .eq. 2) then
    call solver3d_gmres(Ab, As, Aw, Ap, Ae, An, At, b, u_star, m, n, l, solver_tol, maxit)
  elseif (solver .eq. 3) then
    call solver3d_paradiso(Ab, As, Aw, Ap, Ae, An, At, b, u_star, m, n, l, solver_tol, maxit)
  else
    call solver3d_tdma(Ab, As, Aw, Ap, Ae, An, At, b, u_star, m, n, l, solver_tol, maxit)
  end if

  ! Stop Timer
  call cpu_time(end_time)

  ! Output Results
  call output3d

end program main3d