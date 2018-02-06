! diffusion_source3d
!
! Written by Matt Blomquist
! Last Update: 2018-02-06 (YYYY-MM-DD)
!
! This program updates the sources terms for a 3D diffusion problem.

subroutine diffusion_source3d

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k

  ! Update source terms
  do k = 1,l
    do j = 1,n
	  do i = 1,m

	    ! Set default parameters
	    Ab(i,j,k) = k*A_z/dz
		As(i,j,k) = k*A_y/dy
		Aw(i,j,k) = k*A_x/dx
		Ae(i,j,k) = k*A_x/dx
		An(i,j,k) = k*A_y/dy
		At(i,j,k) = k*A_z/dz
		
		! Check bottom
		if (k .eq. 1) then
		  Ab(i,j,k) = 0
		end if

		! Check top
		if (k .eq. l) then
		  At(i,j,k) = 0
		end if

		! Check south
		if (j .eq. 1) then
		  As(i,j,k) = 0
		end if

		! Check north
		if (j .eq. n) then
		  An(i,j,k) = 0
		end if
		
		! Check west
		if (i .eq. 1) then
		  Aw(i,j,k) = 0
		end if

		! Check east
		if (i .eq. m) then
		  Ae(i,j,k) = 0
		end if

		! Update Ap
		Ap(i,j,k) = Ab(i,j,k)+As(i,j,k)+Aw(i,j,k)+Ae(i,j,k)+An(i,j,k)+At(i,j,k)-Sp_t(i,j,k)

		! Update b
		b(i,j,k) = Sc_t(i,j,k)

	  end do
	end do
  end do

  ! Update A values for solver
  Ab = -Ab
  As = -As
  Aw = -Aw
  Ae = -Ae
  An = -An
  At = -At

  return

end subroutine diffusion_source3d
