SUBROUTINE read_converg()
use nrtype ; use coordinates 
implicit none

! this routine is used to read convergence criteria set by user
! file converg.data is created by qommma.py
! TO-DO : Add GSM convergence criteria here.

character(80) :: dummy
open(unit=8,file="converg.data",position="rewind")

read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) tolde_org
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) tolgmax_org
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) tolgrms_org
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) toldxmax_org
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) toldxrms_org
read(unit=8,fmt=*) dummy
read(unit=8,fmt=*) tolper
close(8)
END SUBROUTINE read_converg
