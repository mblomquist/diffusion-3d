! solver3d_bicgstab_mcfb
!
! Written by Matt Blomquist
! Last Update: 2018-01-31 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the bi-conjugate gradients stabilized algorithm.
!
! Definition of input arguments
! Inputs:
!   Ab, As, Aw, Ap, Ae, An, At :: These arrays represent the coefficients for adjacent nodes
!   b :: This array represents the right-hand side of the equation Ax=b
!   phi :: This value represents the appropriate solution array (pressure, velocity, temperature)
!   m, n, l :: These values represent the number of nodes for i, j, and k for the phi value
!   tol :: represents the solution tolerance
!   maxit :: represents the maximum number of iterations of the BiCGStab Algorithm
!
! Outputs:
!   phi :: on exit, this value contains the updated solution
!   maxit :: on exit, this value contains the number of iterations of the BiCGStab algorithm
!   tol :: on exit, this value represents the normalized residual

subroutine solver3d_bicgstab_mcfb(Ab, As_in, Aw, Ap_in, Ae, An, At, b, phi, m, n, l, tol, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n, l
  integer, intent(inout) :: maxit
  real(8), intent(inout) :: tol
  real(8), dimension(m,n,l), intent(in) :: Ab, As_in, Aw, Ap_in, Ae, An, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  real(8), dimension(m*n*l,7) :: A_values
  integer, dimension(7) :: A_distance
  real(8), dimension(m*n*l) :: r, r0, r1, x, p, Ap, As, s, Ax, b_values
  real(8) :: alpha, omega, beta, r_norm


  A_distance = (/-m*n, -m, -1, 0, 1, m, m*n/)

  ! Convert values into CDS Format
  do k = 1,l
    do j = 1,n
      do i = 1,m

        ! Compress stiffness matrix values
        A_values(i+(j-1)*m+(k-1)*m*n,1) = Ab(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,2) = As_in(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,3) = Aw(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,4) = Ap_in(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,5) = Ae(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,6) = An(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,7) = At(i,j,k)

        ! Compress right-hand side values
        b_values(i+(j-1)*m+(k-1)*m*n) = b(i,j,k)

        ! Compress preconditioning values
        x(i+(j-1)*m+(k-1)*m*n) = phi(i,j,k)
      end do
    end do
  end do

  ! ================================================================= !
  ! ==================== Start BiCGStab Algoritm ==================== !
  ! ================================================================= !

  ! Set preconditioning array
  if (sum(x) .eq. 0) then
    x = 1
  end if

  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Ax)
  r = b_values - Ax
  r0 = 1
  p = r

  do itr = 1,maxit


    call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, p, Ap)
    alpha = ddot(m*n*l,r,1,r0,1)/ddot(m*n*l,Ap,1,r0,1)

    s = r - alpha*Ap

    call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, s, As)
    omega = ddot(m*n*l,As,1,s,1) / ddot(m*n*l,As,1,As,1)

    x = x + alpha*p + omega*s

    r1 = s - omega*As

    beta = ddot(m*n*l,r1,1,r0,1) / ddot(m*n*l,r,1,r0,1) * alpha/omega

    p = r1 + beta*(p-omega*Ap)

    r = r1

    ! Check convergence of BiCG
    call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Ax)
	r_norm = abs(dnrm2(m*n*l, b_values-Ax, 1))

    if (r_norm < tol) then
        print *, 'BiCGStab Algorithm successfully converged!'
        print *, 'Number of Iterations: ', itr
        print *, 'Relative residual: ', r_norm
		tol = r_norm
		maxit = i
        exit
    end if

    if (itr .eq. maxit) then
      print *, 'BiCGStab Algorithm did not converge!'
      print *, 'Number of Iterations: ', itr
      print *, 'Relative residual: ', r_norm
	  tol = r_norm
    end if

  end do

  ! ================================================================= !
  ! ====================== End BiCG Algoritm ======================== !
  ! ================================================================= !


  ! Update phi with the solution
  do k = 1,l
    do j = 1,n
      do i = 1,m
        phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
      end do
    end do
  end do

  return

end subroutine solver3d_bicgstab_mcfb
