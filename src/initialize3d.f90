! initialize3d Subroutine for 3D Diffusion Problems
!
! Written by Matt Blomquist
! Last Update: 2018-07-18 (YYYY-MM-DD)
!

subroutine initialize3d

  ! Read Input File ..........
  open(unit = 2, file = "input3d.txt")
  read(2,*)
  read(2,*)
  read(2,*) length, width, depth
  read(2,*)
  read(2,*) T_bc_wv, T_bc_ev, T_bc_nv, T_bc_sv, T_bc_bv, T_bc_tv
  read(2,*)
  read(2,*) T_bc_wc, T_bc_ec, T_bc_nc, T_bc_sc, T_bc_bc, T_bc_tc
  read(2,*)
  read(2,*) maxit, solver_tol, solver
  close(2)

  ! Determine T_h and T_c
  T_h = maxval((/T_bc_wv, T_bc_ev, T_bc_sv, T_bc_nv, T_bc_bv, T_bc_tv/))
  T_c = minval((/T_bc_wv, T_bc_ev, T_bc_sv, T_bc_nv, T_bc_bv, T_bc_tv/))

  return

end subroutine initialize3d
