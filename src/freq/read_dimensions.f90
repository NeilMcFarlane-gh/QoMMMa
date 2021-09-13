SUBROUTINE read_dimensions(ns)
implicit none

INTEGER, INTENT(out) :: ns
character(len=80) :: dummy

OPEN(UNIT=8,FILE="hessian_input")

read (UNIT=8, FMT=*) dummy
read (UNIT=8, FMT=*) ns
close(8)
END SUBROUTINE read_dimensions


