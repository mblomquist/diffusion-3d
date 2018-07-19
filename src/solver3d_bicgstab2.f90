! solver3d_bicgstab2
!
! Written by Matt Blomquist
! Last Update: 2018-07-17 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the bi-conjugate gradients stabilized (2) algorithm by Sleijpen and Vorst.
! reference: Sleijpen et al. 1995
!
! Definition of input arguments
! Inputs:
!   Ab, As, Aw, Ap, Ae, An, At :: These arrays represent the coefficients for adjacent nodes
!   b :: This array represents the right-hand side of the equation Ax=b
!   phi :: This value represents the appropriate solution array (pressure, velocity, temperature)
!   m, n :: These values represent the number of nodes for i and j for the phi value
!   tol :: represents the solution tolerance
!   maxit :: represents the maximum number of iterations of the BiCGStab Algorithm
!
! Outputs:
!   phi :: on exit, this value contains the updated solution
!   maxit :: on exit, this value contains the number of iterations of the BiCGStab algorithm
!   tol :: on exit, this value represents the normalized residual

subroutine solver3d_bicgstab2(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n, l, maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n,l), intent(in) :: As, Aw, Ap, Ae, An, Ab, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  integer, dimension(7) :: A_distance
  real(8) :: alpha, beta, gamma, mu, nu, rho, rho_1, tau, omega_1, omega_2, r_norm
  real(8), dimension(m*n*l) :: Ax, p, r, r0, r0_hat, s, t, w, v, x, x0, b_values
  real(8), dimension(m*n*l, 7) :: A_values

  !  Set A_distance
  A_distance = (/-m*n, -m, -1, 0, 1, m, m*n/)

  ! Convert values into CDS Format
  do k = 1,l
    do j = 1,n
      do i = 1,m

        ! Compress stiffness matrix values
		A_values(i+(j-1)*m+(k-1)*m*n,1) = -Ab(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,2) = -As(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,3) = -Aw(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,4) = Ap(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,5) = -Ae(i,j,k)
        A_values(i+(j-1)*m+(k-1)*m*n,6) = -An(i,j,k)
		A_values(i+(j-1)*m+(k-1)*m*n,7) = -At(i,j,k)

        ! Compress right-hand side values
        b_values(i+(j-1)*m+(k-1)*m*n) = b(i,j,k)

        ! Compress preconditioning values
        x0(i+(j-1)*m+(k-1)*m*n) = phi(i,j,k)

	  end do
    end do
  end do

  ! ======================================================================== !
  ! ========== Start Bi-conjugate Gradients Stabilized (2) Method ========== !
  ! ======================================================================== !

  ! Check inital guess
  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 5, x0, Ax)
  r0 = b_values - Ax

  r_norm = abs(dnrm2(m*n*l, r0, 1))

  if (r_norm < tol) then
    print *, 'Initial guess is a sufficient solution'
	print *, 'relative residual: ', r_norm
    return
  end if

  ! Set r0_hat
  r0_hat = r0

  ! Check that dot(r0, r0_hat) .ne. 0
  rho = ddot(m*n*l, r0, 1, r0_hat, 1)
  if (rho .eq. 0) then
    r0_hat = r0 + 1.
  end if

  ! Set scalars
  rho = 1.
  alpha = 1.
  omega_1 = 1.
  omega_2 = 1.

  ! Set vectors
  w = 0.
  v = 0.
  p = 0.

  ! Update r and x
  x = x0
  r = r0

  ! Start BiCGSTAB(2) Loop
  do itr = 1, maxit+1, 2

    rho_1 = -omega_2*rho

	! Even Bi-CG step:
	rho = ddot(m*n*l,r,1,r0_hat,1)
	beta = alpha * rho / rho_1
	rho_1 = rho
	p = r - beta * (p - omega_1 * v - omega_2 * w)
	call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, p, v)
	gamma = ddot(m*n*l, v, 1, r0_hat, 1)
	alpha = rho / gamma
	r = r - alpha * v
	call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, r, s)
	x = x + alpha * p

	! Check solution
	call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Ax)
	r_norm = abs(dnrm2(m*n*l, b_values - Ax, 1))

	if (r_norm < tol) then
      !print *, 'BiCGSTAB(2) Algorithm successfully converged!(mid)'
      !print *, 'Number of Iterations: ', itr
      !print *, 'Relative residual: ', r_norm

	  do k = 1,l
        do j = 1,n
          do i = 1,m
            phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
		  end do
        end do
      end do

      return
    end if

	! Odd Bi-CG step:
	rho = ddot(m*n*l, s, 1, r0_hat, 1)

	beta = alpha * rho / rho_1

	rho_1 = rho
	v = s - beta * v

	call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, v, w)

	gamma = ddot(m*n*l, w, 1, r0_hat, 1)

	alpha = rho/gamma


	p = r - beta * p

	r = r - alpha * v

	s = s - alpha * w

	call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, s, t)

	! GMRES(2)-part
	omega_1 = ddot(m*n*l, r, 1, s, 1)

	mu = ddot(m*n*l, s, 1, s, 1)

	nu = ddot(m*n*l, s, 1, t, 1)

	tau = ddot(m*n*l, t, 1, t, 1)

	omega_2 = ddot(m*n*l, r, 1, t, 1)

	tau = tau - nu**2 / mu

	omega_2 = (omega_2 - nu * omega_1 / mu) / tau

	omega_1 = (omega_1 - nu*omega_2) / mu

	x = x + alpha * p + omega_1 * r + omega_2 * s
	r = r - omega_1 * s - omega_2 * t


	! Check solution
	call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Ax)
	r_norm = abs(dnrm2(m*n*l, b_values - Ax, 1))

	if (r_norm < tol) then
      print *, 'BiCGSTAB(2) Algorithm successfully converged! (end)'
      print *, 'Number of Iterations: ', itr+1
      print *, 'Relative residual: ', r_norm

	  do k = 1,l
        do j = 1,n
          do i = 1,m
            phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
		  end do
        end do
      end do

      return
    end if

	if (itr .gt. maxit) then
      print *, '************************************'
      print *, '************************************'
      print *, 'BiCGStab Algorithm did not converge!'
      print *, 'Number of Iterations: ', itr
      print *, 'Relative residual: ', r_norm
      print *, '************************************'
      print *, '************************************'
    end if

  end do

  ! ======================================================================== !
  ! ============ End Bi-conjugate Gradients Stabilized (2) Method ========== !
  ! ======================================================================== !

  ! Update phi with the solution
  do k = 1,l
    do j = 1,n
      do i = 1,m
        phi(i,j,k) = x(i+(j-1)*m+(k-1)*m*n)
	  end do
    end do
  end do

  return

end subroutine solver3d_bicgstab2
