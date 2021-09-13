c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine extra  --  user defined extra potentials  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "extra" calculates any additional user defined potential
c     energy contribution
c
c
      subroutine extra
      implicit none
      integer freeunit
      integer i,j,igrd
      real*8 edifqm
      include 'sizes.i'
      include 'atoms.i'
      include 'energi.i'
c
c
c     zero out the energy due to extra potential terms
c
      ex = 0.0d0
      
      edifqm=0.0d0
      do i = 1, n
         edifqm = edifqm +(x(i)-xini(i))*ddeqm(1,i)+(y(i)-yini(i))*
     &      ddeqm(2,i) + (z(i)-zini(i))*ddeqm(3,i)
      end do
c
c     add any user-defined extra potentials below here
c
c     e = ......
c     ex = ex + e
c
      ex = edifqm
      return
      end

