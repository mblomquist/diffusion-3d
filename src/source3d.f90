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
  do i = 1, m
    do j = 1, n
      do k = 1, l

        Aw(i,j,k) = dy*dz/dx
        Ae(i,j,k) = dy*dz/dx
        As(i,j,k) = dz*dx/dy
        An(i,j,k) = dz*dx/dy
        Ab(i,j,k) = dx*dy/dz
        At(i,j,k) = dx*dy/dz

        if (i .eq. 1) then

          Aw(i,j,k) = 0.

        end if
        if (i .eq. m) then

          Ae(i,j,k) = 0.

        end if
        if (j .eq. 1) then

          As(i,j,k) = 0.

        end if
        if (j .eq. n) then

          An(i,j,k) = 0.

        end if
        if (k .eq. 1) then

          Ab(i,j,k) = 0.

        end if
        if (k .eq. l) then

          At(i,j,k) = 0.

        end if

        Ap(i,j,k) = Aw(i,j,k)+Ae(i,j,k)+As(i,j,k)+An(i,j,k)+Ab(i,j,k)+At(i,j,k)+Sp(i,j,k)*dx*dy*dz
        b(i,j,k) = Su(i,j,k)*dx*dy*dz

      end do
    end do
  end do

  return

end subroutine source3d
