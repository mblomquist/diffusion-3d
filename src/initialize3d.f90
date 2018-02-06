! initialize2d Subroutine for Diffusion Problems
!
! Written by Matt Blomquist
! Last Update: 2018-02-02 (YYYY-MM-DD)
!
! This subroutine runs the static calculations for geometry properties,
! pressure properties, velocity, properties, temperature properties, and
! sets the boundary conditions for pressure, velocity, and temperature.

subroutine initialize2d

  ! Pull in standard variable header
  include "var2d.dec"

  real(8) :: T_west, q_south

  ! Calculate geometry properties
  dx = length/m
  dy = width/n
  dz = depth/l

  A_x = dy*dz
  A_y = dx*dz
  A_z = dx*dy

  ! Assign temperature properties
  k = 1000.0
  T_west = 250.0
  T_east = 100.0
  q_south = 500.0

  ! Set initial values to 0
  T = 0
  Sc_t = 0
  Sp_t = 0

  ! Boundary Conditions
  ! West - constant temperature - 100 C
  ! South - Heat Flux - 100 W.m2
  ! North, East - Adiabatic

  ! Assign boundary conditions
  ! West Nodes
  Sc_t(1,:) = Sc_t(1,:) + 2*k*A_x/dx*T_west
  Sp_t(1,:) = Sp_t(1,:) - 2*k*A_x/dx

  ! East Nodes
  Sc_t(m,:) = Sc_t(m,:) + 2*k*A_x/dx*T_east
  Sp_t(m,:) = Sp_t(m,:) - 2*k*A_x/dx

  ! South Nodes
  Sc_t(:,1) = Sc_t(:,1) + 0
  Sp_t(:,1) = Sp_t(:,1) + 0

  ! North Nodes
  Sc_t(:,n) = Sc_t(:,n) + 0
  Sp_t(:,n) = Sp_t(:,n) + 0

  ! Inerior Nodes
  Sc_t(2:m-1,2:n-1) = 0
  Sp_t(2:m-1,2:n-1) = 0

  ! Set T values
  T(1,:) = T_west
  T(m,:) = T_east

  return

end subroutine initialize2d
