SUBROUTINE Jacobi(a,n,d,v,nrot)
 implicit none

INTEGER, INTENT(OUT) :: nrot
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(INOUT) :: a(n,n)
DOUBLE PRECISION, INTENT(OUT) :: v(n,n), d(n)

INTEGER :: i,ip,iq,j
DOUBLE PRECISION :: c, g, h, s, sm, t, tau, theta, tresh, b(n), z(n)
DOUBLE PRECISION, PARAMETER :: ZERO = 0.00D+00

v = ZERO
DO ip=1,n
	v(ip,ip) = 1.D+00
	b(ip) = a(ip,ip)
END DO
d = b
z = ZERO

nrot = 0
do i = 1, 50
	sm = ZERO
	DO ip=1,n-1
		sm = sm + SUM(ABS(a(ip,ip+1:n)))
	END DO
	IF (sm .eq. ZERO) RETURN
	IF (i .lt. 4) THEN
		tresh = 0.2 * sm / n**2
	ELSE
		tresh=ZERO
	END IF
	DO ip = 1, n-1
		DO iq = ip+1, n
			g = 100.*ABS(a(ip,iq))
			IF ((i.gt.4).and.((abs(d(ip))+g.eq.abs(d(ip)))).and.((abs(d(iq))+g.eq.abs(d(iq))))) THEN
				a(ip,iq) = ZERO
			ELSE IF(abs(a(ip,iq)) .gt. tresh) THEN
				h = d(iq) - d(ip)
				IF (abs(h)+g.eq.abs(h)) THEN
					t = a(ip,iq) / h
				ELSE
					theta = 0.5d0 * h/a(ip,iq)
					t = 1.d0/(abs(theta)+sqrt(1.d0 + theta**2))
					IF (theta .lt. ZERO) t = -t
				END IF
				c = 1.d0 / sqrt(1.d0 + t**2)
				s = t * c
				tau = s / (1.d0 + c)
				h = t * a(ip,iq)
				z(ip) = z(ip) - h
				z(iq) = z(iq) + h
				d(ip) = d(ip) - h
				d(iq) = d(iq) + h
				a(ip,iq) = ZERO
				DO j = 1, ip-1
					g = a(j,ip) ; h = a(j,iq)
					a(j,ip) = g - s * (h + g * tau)
					a(j,iq) = h + s * (g - h * tau)
				END DO
				DO J = ip + 1, iq -1
					g = a(ip,j) ; h = a(j,iq)
					a(ip,j) = g - s * (h + g * tau)
					a(j,iq) = h + s * (g - h * tau)
				END DO
				DO j = iq + 1, n
					g = a(ip,j) ; h = a(iq,j)
					a(ip,j) = g - s * (h + g * tau)
					a(iq,j) = h + s * (g - h * tau)
				END DO
				DO j = 1, n
					g = v(j,ip) ; h = v(j,iq)
					v(j,ip) = g - s * (h + g * tau)
					v(j,iq) = h + s * (g - h * tau)
				END DO
				nrot = nrot + 1
			END IF
		END DO
	END DO
	b(:) = b(:) + z(:)
	d(:) = b(:)
	z(:) = ZERO
END DO

END SUBROUTINE Jacobi



