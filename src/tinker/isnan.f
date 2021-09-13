c
c    ##################################
c    subroutine isnan
c    ##################################
c
c
      subroutine is 
      real b 
      integer ib 
      logical isnan 
      equivalence(b,ib) 
c 
      logical function isnan(a) 
      real a 
      if (a.ne.a) then 
        isnan = .true. 
      else 
        isnan = .false. 
      end if 
      return 
      end 
