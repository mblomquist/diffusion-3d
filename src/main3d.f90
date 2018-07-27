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

  ! Define Internal Variables
  real(8) :: start_time, end_time
  integer :: i, j, k, fault

  ! Initialize Problem
  call initialize3d

  ! Set coefficient values
  call boundary3d
  call source3d

  ! Manually Input Source (West = 1, East = 0)
  !Aw(1,:,:) = 0.
  !Ae(1,:,:) = 0.
  !As(1,:,:) = 0.
  !An(1,:,:) = 0.
  !Ab(1,:,:) = 0.
  !At(1,:,:) = 0.
  !Ap(1,:,:) = 1.
  !b(1,:,:) = 1.

  !Aw(m,:,:) = 0.
  !Ae(m,:,:) = 0.
  !As(m,:,:) = 0.
  !An(m,:,:) = 0.
  !Ab(m,:,:) = 0.
  !At(m,:,:) = 0.
  !Ap(m,:,:) = 1.
  !b(m,:,:) = 0.1

  ! Print coefficient
  !do k = 1,l
  !  do j = 1,n
  !    do i = 1,m
  !      print *,i,j,k, Ap(i,j,k), b(i,j,k)
  !    end do
  !  end do
  !end do

  res_vec = 0.

  ! Start Timer
  call cpu_time(start_time)

  ! Solve Diffusion Problem
  if (solver .eq. 0) then
    call solver3d_bicgstab(Ab, As, Aw, Ap, Ae, An, At, b, T, m, n, l, solver_tol, maxit, res_vec)
  elseif (solver .eq. 1) then
    call solver3d_bicgstab2(Ab, As, Aw, Ap, Ae, An, At, b, T, m, n, l, solver_tol, maxit, res_vec)
  elseif (solver .eq. 2) then

    fault = 0

    do i = 3,maxit

      if (fault .eq. 0) then

        call solver3d_gmres(Ab, As, Aw, Ap, Ae, An, At, b, T, m, n, l, solver_tol, i, fault, res_vec)

        if (fault .eq. 0) then
          print *, "Restarting GMRES."
        end if

      end if

    end do

  elseif (solver .eq. 3) then
    call solver3d_bicg(Ab, As, Aw, Ap, Ae, An, At, b, T, m, n, l, solver_tol, maxit, res_vec)
  else
    call solver3d_tdma(Ab, As, Aw, Ap, Ae, An, At, b, T, m, n, l, solver_tol, maxit, res_vec)
  end if

  ! Stop Timer
  call cpu_time(end_time)

  ! Output Results
  print *, "Total Time: ", end_time-start_time

  call output3d

end program main3d
