SUBROUTINE ReadInput_prj(ns,nx,q,h,amx,at,gd)
implicit none

INTEGER, INTENT(IN) :: ns, nx
DOUBLE PRECISION, INTENT(OUT) :: Q(nx), amx(nx), h(nx,nx),gd(nx)
character(len=2), intent(out) :: at(ns)

character(len=80) :: dummy,atfind
character(len=2) :: labels(54)
integer :: I, ii, J, jj, kk, K
double precision :: masses(54)
double precision, parameter:: amu_to_au= 1.822887d3

labels=(/" H","He","Li","Be"," B"," C"," N"," O"," F","Ne","Na","Mg","Al","Si"," P"," S","Cl","Ar", &
	& " K","Ca","Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br", &
	& "Kr","Rb","Sr"," Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te"," I","Xe" /)
masses = (/ 1.007825D+00,4.0026D+00,7.01600D+00,9.01218D+00,11.00931D+00, &
     &   12.0D+00,14.00307D+00,15.99491D+00,18.99840D+00,19.99244D+00, &
     &   22.9898D+00,23.98504D+00,26.98153D+00,27.97693D+00, &
     &   30.97376D+00,31.97207D+00,34.96885D+00,39.948D+00, &
     &   38.96371D+00,39.96259D+00,44.95592D+00,47.90D+00,50.9440D+00, &
     &   51.9405119D+00,54.9349421D+00,55.9349421D+00,58.9332002D+00, &
     &   57.9353479D+00, &
     &   0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
     &   0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
     &   0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
     &   0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
     &   0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
     &   131.9041545d0 /)

OPEN(UNIT=8,FILE="hessian_input",POSITION="rewind")

read (UNIT=8, FMT=*) dummy
read (UNIT=8, FMT=*) dummy
read (UNIT=8, FMT=*) dummy
do i =1,ns
	j=3*(i-1)
	read(8,*) at(i), (q(k),k=j+1,j+3)
        atfind='FALSE'  
	do j=1,54
		if (ADJUSTL(at(i)).eq.ADJUSTL(labels(j))) then
			amx(3*i-2)=masses(j)*amu_to_au
			amx(3*i-1)=amx(3*i-2)
			amx(3*i)=amx(3*i-2)
                        atfind='TRUE'
			exit
		end if
	end do
        IF (atfind.eq.'FALSE') THEN
           write(*,*) 'Atomic mass is not assigned in frequency code for atom  ', i,'   ',at(i)
           stop
        END IF
end do
read (UNIT=8, FMT=*) dummy
do i=1, nx
	j=i/4
	if (j.gt.0) then
	do k=1,j
		read (8,*) kk,(h(i,ii),ii=4*k-3,4*k)
	end do
	end if
	if (i.ne.4*j) then
		read (8,*) kk,(h(i,ii),ii=4*j+1,i)
	end if
	h(1:i-1,i)=h(i,1:i-1)
end do
read (UNIT=8, FMT=*) dummy
do i=1,ns
    ii=(i-1)*3+1
    read (8,*) dummy,(gd(k),k=ii,ii+2)
enddo
CLOSE(8)
END SUBROUTINE ReadInput_prj


