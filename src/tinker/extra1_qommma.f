c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra1  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra1" calculates any additional user defined potential
c     energy contribution and its first derivatives
c
c
      subroutine extra1
      implicit none
      integer freeunit
      integer i,j,igrd
      real*8 edifqm
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'energi.i'
c
c
c     this is simple - just make a copy!
c
      edifqm=0.0d0
      do i = 1, n
         edifqm =edifqm +(x(i)-xini(i))*ddeqm(1,i)+(y(i)-yini(i))*
     &      ddeqm(2,i) + (z(i)-zini(i))*ddeqm(3,i)
         dex(1,i) = ddeqm(1,i)
         dex(2,i) = ddeqm(2,i)
         dex(3,i) = ddeqm(3,i)
      end do

      ex = edifqm
c
c     add any user-defined extra potentials and derivatives;
c     also increment intermolecular energy and virial as needed
c
c     e = ......
c     ex = ex + e
c     do i = 1, n
c        dex(1,i) = ......
c        dex(2,i) = ......
c        dex(3,i) = ......
c     end do
c
      return
      end
