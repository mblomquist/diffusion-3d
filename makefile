# makefile for 3D Diffusion Repository
#
# Written by Matt Blomquist
# Last Update: 2018-07-18 (YYYY-MM-DD)
#
# This file compiles and links the diffusion-3d repository source code into the
# executable file main3d.out
#
main3d.out : boundary3d.o initialize3d.o main3d.o output3d.o solver1d_tdma.o solver3d_bicgstab.o solver3d_bicgstab2.o solver3d_gmres.o solver3d_paradiso.o solver3d_tdma.o source3d.out
	ifort -o build/main3d.out -mkl build/boundary3d.o build/initialize3d.o build/main3d.o build/output3d.o build/solver1d_tdma.o build/solver3d_bicgstab.o build/solver3d_bicgstab2.o build/solver3d_gmres.o build/solver3d_paradiso.o build/solver3d_tdma.o build/source3d.out

boundary3d.o : boundary3d.f90
	ifort -o build/boundary3d.o -c src/boundary3d.f90

initialize3d.o : initialize3d.f90
	ifort -o build/initialize3d.o -c src/initialize3d.f90

main3d.o : main3d.f90
	ifort -o build/main3d.o -c src/main3d.f90

output3d.o : output3d.f90
	ifort -o build/output3d.o -c src/output3d.f90

solver1d_tdma.o : solver1d_tdma.f90
	ifort -o build/solver1d_tdma.o -c src/solver1d_tdma.f90

solver3d_bicgstab.o : solver3d_bicgstab.f90
	ifort -o build/solver3d_bicgstab.o -c src/solver3d_bicgstab.f90

solver3d_bicgstab2.o : solver3d_bicgstab2.f90
	ifort -o build/solver3d_bicgstab2.o -c src/solver3d_bicgstab2.f90

solver3d_gmres.o : solver3d_gmres.f90
	ifort -o build/solver3d_gmres.o -c src/solver3d_gmres.f90

solver3d_paradiso.o : solver3d_paradiso.f90
	ifort -o build/solver3d_paradiso.o -c src/solver3d_paradiso.f90

solver3d_tdma.o : solver3d_tdma.f90
	ifort -o build/solver3d_tdma.o -c src/solver3d_tdma.f90

source3d.o : source3d.f90
	ifort -o build/source3d.o -c src/source3d.f90