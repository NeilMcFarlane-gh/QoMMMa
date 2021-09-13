c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readdiffgrad - qm/mm   gradient correction   ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readdiffgrad" gets a vector of gradient corrections from file 'gradcorrection'
c
c
      subroutine readdiffgrad (silly)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,igrd,freeunit
      logical silly
c
c     read the ref. geometry & gradient
c
      igrd = freeunit ()
      open (unit=igrd,file="gradcorrection")
      read (unit=igrd,fmt=*) i
      do i = 1, n
         read (unit=igrd,fmt='(6F16.8)') xini(i),yini(i),zini(i),
     &              ddeqm(1,i),ddeqm(2,i),ddeqm(3,i)
      end do
      return
      end
