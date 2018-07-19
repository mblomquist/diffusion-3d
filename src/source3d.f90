! source3d Subroutine for 3D Diffusion Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-18 (YYYY-MM-DD)
!
! This subourtine calculates the coefficents of the stiffness matrix.

subroutine source3d

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k  

  ! Update Coefficients
  do i = 2, m-1
    do j = 2, n-1
	  do k = 2, l-1

	    if (i .eq. 1) then
	      Aw(i,j,k) = 0.
		else
		  Aw(i,j,k) = dy*dz/dx
		end if

		if (i .eq. m) then
		  Ae(i,j,k) = 0.
		else
		  Ae(i,j,k) = dy*dz/dx
		end if 

		if (j .eq. 1) then
		  As(i,j,k) = 0.
		else
		  As(i,j,k) = dz*dx/dy
		end if

		if (j .eq. m) then
		  An(i,j,k) = 0.
		else
		  An(i,j,k) = dz*dx/dy
		end if

		if (k .eq. 1) then
		  Ab(i,j,k) = 0.
		else
		  Ab(i,j,k) = dx*dy/dz
		end if

		if (k .eq. l) then
		  At(i,j,k) = 0.
		else
		  At(i,j,k) = dx*dy/dz
		end if

	    Ap(i,j,k) = Aw(i,j,k)+Ae(i,j,k)+As(i,j,k)+An(i,j,k)+Ab(i,j,k)+At(i,j,k)+Sp(i,j,k)*dx*dy*dz
		b(i,j,k) = Su(i,j,k)*dx*dy*dz

	  end do
	end do
  end do

  return

end subroutine source3d