! solver3d_tdma
!
! Written by Matt Blomquist
! Last Update: 2018-07-18 (YYYY-MM-DD)
!
! This program solves a three-dimensional discretization problem utilizing a line-by-line
! TDMA (tri-diagonal matrix algorithm).
!
subroutine solver3d_tdma(Ab, As, Aw, Ap, Ae, An, At, b, phi, m, n, l, tol, maxit)

  integer, intent(in) :: m, n, l, maxit
  real(8), dimension(m,n,l), intent(in) :: Aw, Ae, As, An, At, Ab, Ap, b
  real(8), dimension(m,n,l), intent(inout) :: phi

  integer :: i, j, k, itr
  real(8), dimension(m) :: awe, bwe, cwe, dwe, phiwe
  real(8), dimension(n) :: asn, bsn, csn, dsn, phisn
  real(8), dimension(l) :: abt, bbt, cbt, dbt, phibt
  real(8), dimension(m,n,l) :: r
  real(8) :: r_sum, tol

  do itr = 1,maxit

    ! ==================== West - East ==================== !
	do k = 1,l
	  do j = 1,n
	    do i = 1,m

		  awe(i) = Ap(i,j,k)
		  bwe(i) = -Ae(i,j,k)
		  cwe(i) = -Aw(i,j,k)

		  if (j .eq. 1) then

		    if (k .eq. 1) then
			  dwe(i) = b(i,j,k)+An(i,j,k)*phi(i,j+1,k)+At(i,j,k)*phi(i,j,k+1)
			elseif (k .eq. l) then
			  dwe(i) = b(i,j,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)
			else
			  dwe(i) = b(i,j,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
			end if

		  elseif (j .eq. n) then

		    if (k .eq. 1) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+At(i,j,k)*phi(i,j,k+1)
			elseif (k .eq. l) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+Ab(i,j,k)*phi(i,j,k-1)
			else
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
			end if

		  else

		    if (k .eq. 1) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)+At(i,j,k)*phi(i,j,k+1)
			elseif (k .eq. l) then
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)
			else
			  dwe(i) = b(i,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
			end if

		  end if

		end do

    call solver1d_tdma(awe, bwe, cwe, dwe, phiwe, m)

    phi(:,j,k) = phiwe(:)

	  end do
	end do

	! =================== South - North =================== !
	do i = 1,m
	  do k = 1,l
	    do j = 1,n

		  asn(j) = Ap(i,j,k)
		  bsn(j) = -An(i,j,k)
		  csn(j) = -As(i,j,k)

		  if (k .eq. 1) then

		    if (i .eq. 1) then
		      dsn(j) = b(i,j,k)+Ae(i,j,k)*phi(i-1,j,k)+At(i,j,k)*phi(i,j,k+1)
		    elseif (i .eq. m) then
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+At(i,j,k)*phi(i,j,k+1)
		    else
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i-1,j,k)+At(i,j,k)*phi(i,j,k+1)
		    end if

		  elseif (k .eq. l) then

		    if (i .eq. 1) then
		      dsn(j) = b(i,j,k)+Ae(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)
		    elseif (i .eq. m) then
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)
		    else
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)
		    end if

		  else

		    if (i .eq. 1) then
		      dsn(j) = b(i,j,k)+Ae(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
		    elseif (i .eq. m) then
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
		    else
		      dsn(j) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i-1,j,k)+Ab(i,j,k)*phi(i,j,k-1)+At(i,j,k)*phi(i,j,k+1)
		    end if

		  end if

		end do

    call solver1d_tdma(asn, bsn, csn, dsn, phisn, n)

    phi(i,:,k) = phisn(:)

	  end do
	end do


	! ===================== Bottom - Top ================== !
	do j = 1,n
	  do i = 1,m
	    do k = 1,l

		  abt = Ap(i,j,k)
		  bbt = -At(i,j,k)
		  cbt = -Ab(i,j,k)

		  if (j .eq. 1) then

			if (i .eq. 1) then
			  dbt(k) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+An(i,j,k)*phi(i,j+1,k)
			elseif (i .eq. m) then
			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+An(i,j,k)*phi(i,j+1,k)
			else
			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+An(i,j,k)*phi(i,j+1,k)
			end if

		  elseif (i .eq. n) then

			if (i .eq. 1) then
			  dbt(k) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)
			elseif (i .eq. m) then
			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+As(i,j,k)*phi(i,j-1,k)
			else
			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)
			end if

		  else

			if (i .eq. 1) then
			  dbt(k) = b(i,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)
			elseif (i .eq. m) then
			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)
			else
			  dbt(k) = b(i,j,k)+Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)+As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)
			end if

		  end if

		end do

    call solver1d_tdma(abt, bbt, cbt, dbt, phibt, l)

    phi(i,j,:) = phibt(:)

	  end do
	end do

	! =================== Check Solution ================== !
	do i = 1,m
	  do j = 1,n
	    do k = 1,l

		  if (i .eq. 1) then
			if (j .eq. 1) then
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			elseif (j .eq. n) then
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			else
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			end if
		  elseif (i .eq. m) then
			if (j .eq. 1) then
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			elseif (j .eq. n) then
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			else
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			end if
		  else
			if (j .eq. 1) then
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			elseif (j .eq. n) then
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			else
			  if (k .eq. 1) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  elseif (k .eq. l) then
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 b(i,j,k))
			  else
			    r(i,j,k) = Ap(i,j,k)*phi(i,j,k)-(Aw(i,j,k)*phi(i-1,j,k)+ &
				                                 Ae(i,j,k)*phi(i+1,j,k)+ &
				                                 As(i,j,k)*phi(i,j-1,k)+ &
				                                 An(i,j,k)*phi(i,j+1,k)+ &
				                                 Ab(i,j,k)*phi(i,j,k-1)+ &
				                                 At(i,j,k)*phi(i,j,k+1)+ &
				                                 b(i,j,k))
			  end if
			end if
		  end if

        end do
      end do
    end do

	r_sum = 0.

	do i = 1,m
	  do j = 1,n
	    do k = 1,l
		  r_sum = r_sum + abs(r(i,j,k))
		end do
	  end do
	end do

	if (r_sum .le. tol) then
	  print *, "TDMA Compelete."
    print *, "r_sum:", r_sum
    print *, "itrs:", itr
      return
  else
    print *, "r_sum:", r_sum
    print *, "itrs:", itr
	end if

  end do

  return

end subroutine solver3d_tdma
