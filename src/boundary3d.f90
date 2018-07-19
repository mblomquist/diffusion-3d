! boundary3d Subroutine for 3D Diffusion Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-18 (YYYY-MM-DD)
!

subroutine boundary3d

  ! Include standard variable header
  include "var3d.dec"

  ! Define internal variables
  integer :: i, j, k 

 ! Set Boundary Condition :: West Plane
  Aw(1,2:n-1,2:l-1) = 0.
  Ae(1,2:n-1,2:l-1) = 0.
  As(1,2:n-1,2:l-1) = 0.
  An(1,2:n-1,2:l-1) = 0.
  Ab(1,2:n-1,2:l-1) = 0.
  At(1,2:n-1,2:l-1) = 0.
  Ap(1,2:n-1,2:l-1) = 1.
  b(1,2:n-1,2:l-1) = 0.

  ! Check Boundary Condition Type
  if (T_bc_wc .eq. 1) then
    ! Given Input Condition
    b(1,2:n-1,2:l-1) = (T_bc_wv-T_c)/(T_h-T_c)
  elseif (T_bc_wc .eq. 2) then
    ! Symmetry Condition
    Ae(1,2:n-1,2:l-1) = 1.
    b(1,2:n-1,2:l-1) = 0.
  else
    ! Wall Condition
    b(1,2:n-1,2:l-1) = 0.
  end if

  ! Set Boundary Condition :: East Plane
  Aw(m,2:n-1,2:l-1) = 0.
  Ae(m,2:n-1,2:l-1) = 0.
  As(m,2:n-1,2:l-1) = 0.
  An(m,2:n-1,2:l-1) = 0.
  Ab(m,2:n-1,2:l-1) = 0.
  At(m,2:n-1,2:l-1) = 0.
  Ap(m,2:n-1,2:l-1) = 1.
  b(m,2:n-1,2:l-1) = 0.

  ! Check Boundary Condition Type
  if (T_bc_ec .eq. 1) then
    ! Given Input Condition
    b(m,2:n-1,2:l-1) = (T_bc_ev-T_c)/(T_h-T_c)
  elseif (T_bc_ec .eq. 2) then
    ! Symmetry Condition
    Aw(m,2:n-1,2:l-1) = 1.
    b(m,2:n-1,2:l-1) = 0.
  else
    ! Wall Condition
    b(m,2:n-1,2:l-1) = 0.
  end if

  ! Set Boundary Condition :: South Plane
  Aw(2:m-1,1,2:l-1) = 0.
  Ae(2:m-1,1,2:l-1) = 0.
  As(2:m-1,1,2:l-1) = 0.
  An(2:m-1,1,2:l-1) = 0.
  Ab(2:m-1,1,2:l-1) = 0.
  At(2:m-1,1,2:l-1) = 0.
  Ap(2:m-1,1,2:l-1) = 1.
  b(2:m-1,1,2:l-1) = 0.

  ! Check Boundary Condition Type
  if (T_bc_sc .eq. 1) then
    ! Given Input Condition
    b(2:m-1,1,2:l-1) = (T_bc_sv-T_c)/(T_h-T_c)
  elseif (T_bc_sc .eq. 2) then
    ! Symmetry Condition
    An(2:m-1,1,2:l-1) = 1.
    b(2:m-1,1,2:l-1) = 0.
  else
    ! Wall Condition
    b(2:m-1,1,2:l-1) = 0.
  end if

  ! Set Boundary Condition :: North Plane
  Aw(2:m-1,n,2:l-1) = 0.
  Ae(2:m-1,n,2:l-1) = 0.
  As(2:m-1,n,2:l-1) = 0.
  An(2:m-1,n,2:l-1) = 0.
  Ab(2:m-1,n,2:l-1) = 0.
  At(2:m-1,n,2:l-1) = 0.
  Ap(2:m-1,n,2:l-1) = 1.
  b(2:m-1,n,2:l-1) = 0.

  ! Check Boundary Condition Type
  if (T_bc_nc .eq. 1) then
    ! Given Input Condition
    b(2:m-1,n,2:l-1) = (T_bc_nv-T_c)/(T_h-T_c)
  elseif (T_bc_nc .eq. 2) then
    ! Symmetry Condition
    As(2:m-1,n,2:l-1) = 1.
    b(2:m-1,n,2:l-1) = 0.
  else
    ! Wall Condition
    b(2:m-1,n,2:l-1) = 0.
  end if

  ! Set Boundary Condition :: Bottom Plane
  Aw(2:m-1,2:n-1,1) = 0.
  Ae(2:m-1,2:n-1,1) = 0.
  As(2:m-1,2:n-1,1) = 0.
  An(2:m-1,2:n-1,1) = 0.
  Ab(2:m-1,2:n-1,1) = 0.
  At(2:m-1,2:n-1,1) = 0.
  Ap(2:m-1,2:n-1,1) = 1.
  b(2:m-1,2:n-1,1) = 0.

  ! Check Boundary Condition Type
  if (T_bc_bc .eq. 1) then
    ! Given Input Condition
    b(2:m-1,2:n-1,1) = (T_bc_bv-T_c)/(T_h-T_c)
  elseif (T_bc_bc .eq. 2) then
    ! Symmetry Condition
    At(2:m-1,2:n-1,1) = 1.
    b(2:m-1,2:n-1,1) = 0.
  else
    ! Wall Condition
    b(2:m-1,2:n-1,1) = 0.
  end if

  ! Set Boundary Condition :: Top Plane
  Aw(2:m-1,2:n-1,l) = 0.
  Ae(2:m-1,2:n-1,l) = 0.
  As(2:m-1,2:n-1,l) = 0.
  An(2:m-1,2:n-1,l) = 0.
  Ab(2:m-1,2:n-1,l) = 0.
  At(2:m-1,2:n-1,l) = 0.
  Ap(2:m-1,2:n-1,l) = 1.
  b(2:m-1,2:n-1,l) = 0.

  ! Check Boundary Condition Type
  if (T_bc_tc .eq. 1) then
    ! Given Input Condition
    b(2:m-1,2:n-1,l) = (T_bc_tv-T_c)/(T_h-T_c)
  elseif (T_bc_tc .eq. 2) then
    ! Symmetry Condition
    Ab(2:m-1,2:n-1,l) = 1.
    b(2:m-1,2:n-1,l) = 0.
  else
    ! Wall Condition
    b(2:m-1,2:n-1,l) = 0.
  end if

  ! Set Boundary Condition :: West-South Line
  Aw(1,1,2:l-1) = 0.
  Ae(1,1,2:l-1) = 1.
  As(1,1,2:l-1) = 0.
  An(1,1,2:l-1) = 1.
  Ab(1,1,2:l-1) = 0.
  At(1,1,2:l-1) = 0.
  Ap(1,1,2:l-1) = 2.
  b(1,1,2:l-1) = 0.

  ! Set Boundary Condition :: West-North Line
  Aw(1,n,2:l-1) = 0.
  Ae(1,n,2:l-1) = 1.
  As(1,n,2:l-1) = 1.
  An(1,n,2:l-1) = 0.
  Ab(1,n,2:l-1) = 0.
  At(1,n,2:l-1) = 0.
  Ap(1,n,2:l-1) = 2.
  b(1,n,2:l-1) = 0.

  ! Set Boundary Condition :: West-Bottom Line
  Aw(1,2:n-1,1) = 0.
  Ae(1,2:n-1,1) = 1.
  As(1,2:n-1,1) = 0.
  An(1,2:n-1,1) = 0.
  Ab(1,2:n-1,1) = 0.
  At(1,2:n-1,1) = 1.
  Ap(1,2:n-1,1) = 2.
  b(1,2:n-1,1) = 0.

  ! Set Boundary Condition :: West-Top Line
  Aw(1,2:n-1,l) = 0.
  Ae(1,2:n-1,l) = 1.
  As(1,2:n-1,l) = 0.
  An(1,2:n-1,l) = 0.
  Ab(1,2:n-1,l) = 1.
  At(1,2:n-1,l) = 0.
  Ap(1,2:n-1,l) = 2.
  b(1,2:n-1,l) = 0.

  ! Set Boundary Condition :: East-South Line
  Aw(m,1,2:l-1) = 1.
  Ae(m,1,2:l-1) = 0.
  As(m,1,2:l-1) = 0.
  An(m,1,2:l-1) = 1.
  Ab(m,1,2:l-1) = 0.
  At(m,1,2:l-1) = 0.
  Ap(m,1,2:l-1) = 2.
  b(m,1,2:l-1) = 0.

  ! Set Boundary Condition :: East-North Line
  Aw(m,n,2:l-1) = 1.
  Ae(m,n,2:l-1) = 0.
  As(m,n,2:l-1) = 1.
  An(m,n,2:l-1) = 0.
  Ab(m,n,2:l-1) = 0.
  At(m,n,2:l-1) = 0.
  Ap(m,n,2:l-1) = 2.
  b(m,n,2:l-1) = 0.

  ! Set Boundary Condition :: East-Bottom Line
  Aw(m,2:n-1,1) = 1.
  Ae(m,2:n-1,1) = 0.
  As(m,2:n-1,1) = 0.
  An(m,2:n-1,1) = 0.
  Ab(m,2:n-1,1) = 0.
  At(m,2:n-1,1) = 1.
  Ap(m,2:n-1,1) = 2.
  b(m,2:n-1,1) = 0.

  ! Set Boundary Condition :: East-Top Line
  Aw(m,2:n-1,l) = 1.
  Ae(m,2:n-1,l) = 0.
  As(m,2:n-1,l) = 0.
  An(m,2:n-1,l) = 0.
  Ab(m,2:n-1,l) = 1.
  At(m,2:n-1,l) = 0.
  Ap(m,2:n-1,l) = 2.
  b(m,2:n-1,l) = 0.

  ! Set Boundary Condition :: South-Bottom Line
  Aw(2:m-1,1,1) = 0.
  Ae(2:m-1,1,1) = 0.
  As(2:m-1,1,1) = 0.
  An(2:m-1,1,1) = 1.
  Ab(2:m-1,1,1) = 0.
  At(2:m-1,1,1) = 1.
  Ap(2:m-1,1,1) = 2.
  b(2:m-1,1,1) = 0.

  ! Set Boundary Condition :: South-Top Line
  Aw(2:m-1,1,l) = 0.
  Ae(2:m-1,1,l) = 0.
  As(2:m-1,1,l) = 0.
  An(2:m-1,1,l) = 1.
  Ab(2:m-1,1,l) = 1.
  At(2:m-1,1,l) = 0.
  Ap(2:m-1,1,l) = 2.
  b(2:m-1,1,l) = 0.

  ! Set Boundary Condition :: North-Bottom Line
  Aw(2:m-1,n,1) = 0.
  Ae(2:m-1,n,1) = 0.
  As(2:m-1,n,1) = 1.
  An(2:m-1,n,1) = 0.
  Ab(2:m-1,n,1) = 0.
  At(2:m-1,n,1) = 1.
  Ap(2:m-1,n,1) = 2.
  b(2:m-1,n,1) = 0.

  ! Set Boundary Condition :: North-Top Line
  Aw(2:m-1,n,l) = 0.
  Ae(2:m-1,n,l) = 0.
  As(2:m-1,n,l) = 1.
  An(2:m-1,n,l) = 0.
  Ab(2:m-1,n,l) = 1.
  At(2:m-1,n,l) = 0.
  Ap(2:m-1,n,l) = 2.
  b(2:m-1,n,l) = 0.

  ! Set Boundary Condition :: West-South-Bottom Corner
  Aw(1,1,1) = 0.
  Ae(1,1,1) = 1.
  As(1,1,1) = 0.
  An(1,1,1) = 1.
  Ab(1,1,1) = 0.
  At(1,1,1) = 1.
  Ap(1,1,1) = 3.
  b(1,1,1) = 0.

  ! Set Boundary Condition :: West-North-Bottom Corner
  Aw(1,n,1) = 0.
  Ae(1,n,1) = 1.
  As(1,n,1) = 1.
  An(1,n,1) = 0.
  Ab(1,n,1) = 0.
  At(1,n,1) = 1.
  Ap(1,n,1) = 3.
  b(1,n,1) = 0.

  ! Set Boundary Condition :: East-South-Bottom Corner
  Aw(m,1,1) = 1.
  Ae(m,1,1) = 0.
  As(m,1,1) = 0.
  An(m,1,1) = 1.
  Ab(m,1,1) = 0.
  At(m,1,1) = 1.
  Ap(m,1,1) = 3.
  b(m,1,1) = 0.

  ! Set Boundary Condition :: East-North-Bottom Corner
  Aw(m,n,1) = 1.
  Ae(m,n,1) = 0.
  As(m,n,1) = 1.
  An(m,n,1) = 0.
  Ab(m,n,1) = 0.
  At(m,n,1) = 1.
  Ap(m,n,1) = 3.
  b(m,n,1) = 0.

  ! Set Boundary Condition :: West-South-Top Corner
  Aw(1,1,l) = 0.
  Ae(1,1,l) = 1.
  As(1,1,l) = 0.
  An(1,1,l) = 1.
  Ab(1,1,l) = 1.
  At(1,1,l) = 0.
  Ap(1,1,l) = 3.
  b(1,1,l) = 0.

  ! Set Boundary Condition :: West-North-Top Corner
  Aw(1,n,l) = 0.
  Ae(1,n,l) = 1.
  As(1,n,l) = 1.
  An(1,n,l) = 0.
  Ab(1,n,l) = 1.
  At(1,n,l) = 0.
  Ap(1,n,l) = 3.
  b(1,n,l) = 0.

  ! Set Boundary Condition :: East-South-Top Corner
  Aw(m,1,l) = 1.
  Ae(m,1,l) = 0.
  As(m,1,l) = 0.
  An(m,1,l) = 1.
  Ab(m,1,l) = 1.
  At(m,1,l) = 0.
  Ap(m,1,l) = 3.
  b(m,1,l) = 0.

  ! Set Boundary Condition :: East-North-Top Corner
  Aw(m,n,l) = 1.
  Ae(m,n,l) = 0.
  As(m,n,l) = 1.
  An(m,n,l) = 0.
  Ab(m,n,l) = 1.
  At(m,n,l) = 0.
  Ap(m,n,l) = 3.
  b(m,n,l) = 0.

  return

end subroutine boundary3d