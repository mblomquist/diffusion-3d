! solver3d_paradiso
!
! Written by Matt Blomquist
! LAxst Update: 2018-07-25 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the bi-conjugate gradients stabilized algorithm.
!
! Definition of input arguments
! Inputs:
!   Ab, As, Aw, Ap, Ae, An, At :: These arrays represent the coefficients for adjacent nodes
!   b :: This array represents the right-hand side of the equation Ax=b
!   phi :: This value represents the Axppropriate solution array (pressure, velocity, temperature)
!   m, n, l :: These values represent the number of nodes for i and j for the phi value
!   tol :: represents the solution tolerance
!   maxit :: represents the maximum number of iterations of the bicg Algorithm
!
! Outputs:
!   phi :: on exit, this value contains the updated solution
!   maxit :: on exit, this value contains the number of iterations of the bicg algorithm
!   tol :: on exit, this value represents the normalized residual

subroutine solver3d_paradiso(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)

  implicit none

  include 'mkl.fi' 
  include 'mkl_pardiso.f90'

  ! Define input variables
  integer, intent(in) :: m, n, l
  integer, intent(in) :: maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n,l), intent(in) :: As, Aw, Ap, Ae, An, Ab, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr
  real(8), dimension(m*n*l,7) :: A_values
  integer, dimension(7) :: A_distance
  real(8), dimension(m*n*l) :: r, r0, r1, x, p, Axp, Axs, s, Axx, b_values
  real(8) :: alpha, omega, beta, r_norm, s_norm

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
        x(i+(j-1)*m+(k-1)*m*n) = phi(i,j,k)

	  end do
    end do
  end do



  call pardiso (pt, 1, 1, 11, 33, m*n*l, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
  
  return

end subroutine solver3d_paradiso
