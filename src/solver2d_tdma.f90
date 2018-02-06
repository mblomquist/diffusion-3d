! solver3d_tdma
!
! Written by Matt Blomquist
! Last Update: 2018-02-06 (YYYY-MM-DD)
!
! This program solves a three-dimensional finite volume discretization problem
! using the line-by-line tri-diagonal matrix algorithm.
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

subroutine solver3d_tdma

  ! Define implicit
  implicit none

  ! Include mkl functions
  include "mkl.fi"

  ! Define input variables
  integer, intent(in) :: m, n, l
  integer, intent(inout) :: maxit
  real(8), intent(inout) :: tol
  real(8), dimension(m,n,l), intent(in) :: Ab, As, Aw, Ap, Ae, An, At, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  ! Define internal variables
  integer :: i, j, k, itr_main, itr_inner
  integer, dimension(7) :: A_distance

  real(8) :: r_norm
  real(8), dimension(m) :: bx_line
  real(8), dimension(n) :: by_line
  real(8), dimension(l) :: bz_line
  real(8), dimension(m*n*l,7) :: A_values
  real(8), dimension(m*n*l) :: x_compressed, b_values

  ! ----------------------------------------------------------------------------------------- !
  ! -------------------------------------- Start TDMA --------------------------------------- !
  ! ----------------------------------------------------------------------------------------- !
  do itr_main = 1,maxit

    ! Start west to east loop
    do itr_inner = 1,itr_main
	  
	  k = 1

	  do j = 1,n
	    
		do i = 1,m
		
		  if (j .eq. 1) then
		    bx_line(i) =
		  elseif (j .eq. n) then
		    bx_line(i) =
		  else
		    bx_line(i) =
		  end if
		
		end do

	    call solver1d_tdma()

	  end do

	  do k = 2:l-1
	    do j = 1,n
	    
		  do i = 1,m
		
		    if (j .eq. 1) then
			  bx_line(i) =
		    elseif (j .eq. n) then
			  bx_line(i) =
		    else
			  bx_line(i) =
		    end if
		
		  end do

	      call solver1d_tdma()

	    end do
	  end do
	  
	  k = l

	  do j = 1,n
	    
		do i = 1,m
		
		  if (j .eq. 1) then
		    bx_line(i) =
		  elseif (j .eq. n) then
		    bx_line(i) =
		  else
		    bx_line(i) =
		  end if
		
		end do

	    call solver1d_tdma()

	  end do

	end do

    ! Start south to north loop
    do itr_inner = 1,itr_main

	  i = 1

	  do k = 1:l

	    do j = 1:n

		  if (k .eq. 1) then
		    by_line(j) = 
		  elseif (k .eq. l) then
		    by_line(j) = 
		  else
		    by_line(j) = 
		  end if

		end do

		call solver1d_tdma()

	  end do

	  do i = 2:m-1

	  	do k = 1:l

	      do j = 1:n

		    if (k .eq. 1) then
		      by_line(j) = 
		    elseif (k .eq. l) then
		      by_line(j) = 
		    else
		      by_line(j) = 
		    end if

  		  end do

		  call solver1d_tdma()

	    end do

	  end do

	  i = m

	  do k = 1:l

	    do j = 1:n

		  if (k .eq. 1) then
		    by_line(j) = 
		  elseif (k .eq. l) then
		    by_line(j) = 
		  else
		    by_line(j) = 
		  end if

		end do

		call solver1d_tdma()

	  end do

	end do

    ! Start bottom to top loop
    do itr_inner = 1,itr_main

	  j = 1

	  do i = 1,m

	    do k = 1,l

		  if (i .eq. 1) then
		    bz_line(k) = 
		  elseif (i .eq. m) then
		    bz_line(k) = 
		  else
		    bz_line(k) = 
		  end if

		end do

		call solver1d_tdma()

	  end do

	  do j = 2:n-1

	    do i = 1,m

	      do k = 1,l

		    if (i .eq. 1) then
		      bz_line(k) = 
		    elseif (i .eq. m) then
		      bz_line(k) = 
		    else
		      bz_line(k) = 
		    end if

		  end do

		  call solver1d_tdma()

	    end do

	  end do

	  j = n

	  do i = 1,m

	    do k = 1,l

		  if (i .eq. 1) then
		    bz_line(k) = 
		  elseif (i .eq. m) then
		    bz_line(k) = 
		  else
		    bz_line(k) = 
		  end if

		end do

		call solver1d_tdma()

	  end do

	end do

    ! ----------------------------------------------------------------------------------------- !
    ! ----------------------------------- Check Convergence ----------------------------------- !
    ! ----------------------------------------------------------------------------------------- !

    ! Convert values into CDS Format
    do k = 1,l
      do j = 1,n
        do i = 1,m

          ! Compress stiffness matrix values
          A_values(i+(j-1)*m+(k-1)*m*n,1) = Ab(i,j,k)
          A_values(i+(j-1)*m+(k-1)*m*n,2) = As(i,j,k)
          A_values(i+(j-1)*m+(k-1)*m*n,3) = Aw(i,j,k)
          A_values(i+(j-1)*m+(k-1)*m*n,4) = Ap(i,j,k)
          A_values(i+(j-1)*m+(k-1)*m*n,5) = Ae(i,j,k)
          A_values(i+(j-1)*m+(k-1)*m*n,6) = An(i,j,k)
          A_values(i+(j-1)*m+(k-1)*m*n,7) = At(i,j,k)

          ! Compress right-hand side values
          b_values(i+(j-1)*m+(k-1)*m*n) = b(i,j,k)

          ! Compress preconditioning values
          x_compressed(i+(j-1)*m+(k-1)*m*n) = x(i,j,k)
        end do
      end do
    end do

    ! Check convergence of BiCG
    call mkl_ddiagemv('N', m*n*l, A_values, m*n*l, A_distance, 7, x, Axx)
    r_norm = abs(dnrm2(m*n*l, b_values-Axx, 1))

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


    ! ----------------------------------------------------------------------------------------- !
    ! ------------------------------- End Check Convergence ----------------------------------- !
    ! ----------------------------------------------------------------------------------------- !

  end do

  ! ----------------------------------------------------------------------------------------- !
  ! --------------------------------------- End TDMA ---------------------------------------- !
  ! ----------------------------------------------------------------------------------------- !

  return

end subroutine solver3d_tdma
