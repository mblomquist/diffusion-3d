! solver3d_gmres
!
! Written by Matt Blomquist
! LAxst Update: 2018-07-23 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the generalized minimum residual (GMRES) algorithm.
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

subroutine solver3d_gmres(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)
  
  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"
  include "blas.f90"

  ! Define input variables
  integer, intent(in) :: m, n, l
  integer, intent(in) :: maxit
  real(8), intent(in) :: tol
  real(8), dimension(m,n,l), intent(in) :: As, Aw, Ap, Ae, An, Ab, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables :: Matrix Conversion
  integer :: i, j, k
  real(8), dimension(m*n*l,7) :: A_values
  integer, dimension(7) :: A_distance
  real(8), dimension(m*n*l) :: b_values

  ! Define internal variables :: GMRES


  ! ================================================================= !
  ! ==================== Start Matrix Conversion ==================== !
  ! ================================================================= !

  ! Define distance between diagonals 
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

  ! ================================================================= !
  ! ==================== End Matrix Conversion ====================== !
  ! ================================================================= !

  ! ================================================================= !
  ! ====================== Start GMRES Algoritm ===================== !
  ! ================================================================= !

  ! Compute residual vector
  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Axx)
  r = b_values - Axx

  ! Compute norm of values vector
  b_norm = dnrm2(m*n*l, b_values, 1)
  
  ! Compute norm of r vector
  r_norm = dnrm2(m*n*l, r, 1)

  ! Compute inital error
  error(1) = r_norm/b_norm

  ! Set standard basis
  e1(1) = 1

  ! Initialize Q matrix
  Q(:,1) = r/r_norm

  ! Initalize beta
  beta = r_norm*e1

  ! Start GMRES Loop
  do k = 1,maxit

    ! Run Arnoldi Loop
	q =  A*Q(:,k)

	do i = 1,k

	  h(i) = q*Q(:,i)
	  q = q - h(i)*Q(:,i)

	end do

	h(k+1) = drnm2(m*n*l, q, 1)
	q = q / h(k+1)

	! Apply Given's Rotation
	call arnoldi(A, Q, k, H, m)
	
	! Calculate Rotation
	call given_rotation(H, cs, sn, k)

	beta(k+1) = -sn(k) * beta(k)
	beta(k) = cs(k)*beta(k)
	error(k) = abs(beta(k+1)) / b_norm

	! Check solution
	if (error(k) .le. tol)
	  
	  do ii = 1,k
	    diag(ii) = H(ii,ii)
	  end do

	  y(1:k) = beta(1:k)

	  do ii = 1,k
	    do jj = ii,k
	      H(ii,jj) = H(ii,jj)/diag(ii)
		  y(ii) = y(ii)/diag(ii)
	    end do
      end do

	  call dtrsv('U', 'N', diag, k, H(1:k,1:k), k, y, 1)
	  call dgemv('N', m*n*l, k, 1, Q(:,1:k), m*n*l, y, 1, 1, x, 1)

	end if

  ! ================================================================= !
  ! ======================== End GMRES Algoritm ===================== !
  ! ================================================================= !

  return

end subroutine solver3d_gmres


subroutine arnoldi(A, Qmat, k, H, m)

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"
  include "blas.f90"

  ! Define input variables
  integer, intent(in) :: k
  integer, dimension(7) :: A_distance
  real(8), dimension(m,7) :: A_values
  real(8), dimension(m,k) :: Q
  

  call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, Qmat(:,k), AxQ)

  do i = 1,k
    
	call dgemv('N', m*n*l, k, 1, Qmat(:,1:k), m*n*l, q, 1, 1, h(i), 1)
	
	q = q - h(i)*Qmat(:,i)

  end do

  h(k+1) = dnrm2(m,q,1)
  q = q / h(k+1)

  return
end subroutine arnoldi

subroutine given_rotation(h, cs, sn, k)

  ! Define implicit
  implicit none
  
  ! Define input variables
  integer, intent(in) :: k
  real(8), dimension(k+1) :: h, cs, sn

  ! Define internal variables
  integer :: i
  real(8) :: temp, cs_k, sn_k

  do i = 1,k-1

    temp = cs(i) * h(i) + sn(i)*h(i+1)
	h(i+1) = -sn(i)*h(i) + cs(i)*h(i+1)
	h(i) = temp

  end do

  call calc_rotation(h(k), h(k+1), cs_k, sn_k, k)

  h(k) = cs_k*h(k) + sn_k*h(k+1)
  h(k+1) = 0.

  return
end subroutine given_rotation


subroutine calc_rotation(v1, v2, cs_k, sn_k, k)

  ! Define implicit
  implicit none

  ! Define Input Variables 
  integer, intent(in) :: k
  real(8), intent(in) :: v1, v2
  real(8), intent(inout) :: cs_k, sn_k

  ! Define internal variables
  real(8) :: t

  if (v1 .eq. 0.) then

    cs_k = 0.
	sn_k = 1.

  else

    t = (v1**2.0+v2**2.0)**(0.5)
	cs_k = abs(v1) / t
	sn_k = cs_k * v2 / v1

  end if

  return
end subroutine calc_rotation
