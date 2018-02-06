! diffusion_source2d
!
! Written by Matt Blomquist
! Last Update: 2018-02-02 (YYYY-MM-DD)
!
! This program updates the sources term for a 2D diffusion problem.

subroutine diffusion_source2d

  ! Include standard variable header
  include "var2d.dec"

  ! Define internal variables
  integer :: i, j

  ! Run loop for internal nodes
  do j = 2, n-1
    do i = 2, m-1

      As(i,j) = k*A_y/dy
      Aw(i,j) = k*A_x/dx
      Ae(i,j) = k*A_x/dx
      An(i,j) = k*A_y/dy
      Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

      b(i,j) = Sc_t(i,j)

      x_initial(i+(j-1)*m) = T(i,j)

    end do
  end do

  ! Run loop for south and north edges
  do j = 2,n-1

    i = 1

    As(i,j) = k*A_x/dx
    Aw(i,j) = 0
    Ae(i,j) = k*A_x/dx
    An(i,j) = k*A_y/dy
    Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

    b(i,j) = Sc_t(i,j)

    x_initial(i+(j-1)*m) = T(i,j)

    i = m

    As(i,j) = k*A_y/dy
    Aw(i,j) = k*A_x/dx
    Ae(i,j) = 0
    An(i,j) = k*A_x/dx
    Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

    b(i,j) = Sc_t(i,j)

    x_initial(i+(j-1)*m) = T(i,j)

  end do

  ! Run loop for west and east edges
  do i = 2,m-1

    j = 1

    As(i,j) = 0
    Aw(i,j) = k*A_x/dx
    Ae(i,j) = k*A_x/dx
    An(i,j) = k*A_y/dy
    Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

    b(i,j) = Sc_t(i,j)

    x_initial(i+(j-1)*m) = T(i,j)

    j = n

    As(i,j) = k*A_y/dy
    Aw(i,j) = k*A_x/dx
    Ae(i,j) = k*A_x/dx
    An(i,j) = 0
    Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

    b(i,j) = Sc_t(i,j)

    x_initial(i+(j-1)*m) = T(i,j)

  end do

  ! Define corners
  ! West - South
  i = 1
  j = 1

  As(i,j) = 0
  Aw(i,j) = 0
  Ae(i,j) = k*A_x/dx
  An(i,j) = k*A_y/dy
  Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

  b(i,j) = Sc_t(i,j)

  x_initial(i+(j-1)*m) = T(i,j)

  ! West-North
  i = 1
  j = n

  As(i,j) = k*A_y/dy
  Aw(i,j) = 0
  Ae(i,j) = k*A_x/dx
  An(i,j) = 0
  Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

  b(i,j) = Sc_t(i,j)

  x_initial(i+(j-1)*m) = T(i,j)

  ! East-South
  i = m
  j = 1

  As(i,j) = 0
  Aw(i,j) = k*A_y/dy
  Ae(i,j) = 0
  An(i,j) = k*A_y/dy
  Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

  b(i,j) = Sc_t(i,j)

  x_initial(i+(j-1)*m) = T(i,j)

  ! East-North
  i = m
  j = n

  As(i,j) = k*A_y/dy
  Aw(i,j) = k*A_y/dy
  Ae(i,j) = 0
  An(i,j) = 0
  Ap(i,j) = As(i,j)+Aw(i,j)+Ae(i,j)+An(i,j)-Sp_t(i,j)

  b(i,j) = Sc_t(i,j)

  x_initial(i+(j-1)*m) = T(i,j)

  ! Update A values for solver
  As = -As
  Aw = -Aw
  Ae = -Ae
  An = -An

  return

end subroutine diffusion_source2d
