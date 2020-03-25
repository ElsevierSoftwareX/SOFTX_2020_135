!>Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %\n
!>Junior University.  All rights reserved.                             %\n
!>
!>The programs in GSLIB are distributed in the hope that they will be  %\n
!>useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %\n
!>responsibility to anyone for the consequences of using them or for   %\n
!>whether they serve any particular purpose or work at all, unless he  %\n
!>says so in writing.  Everyone is granted permission to copy, modify  %\n
!>and redistribute the programs in GSLIB, but only under the condition %\n
!>that this notice and the above copyright notice remain intact.       %\n
      double precision function powint(xlow,xhigh,ylow,yhigh,xval,pow)
c-----------------------------------------------------------------------
c
c Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
c                 for a value of x and a power pow.
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(EPSLON=1.0e-20)

      if((xhigh-xlow).lt.EPSLON) then
            powint = (yhigh+ylow)/2.0d0
      else
            powint = ylow + (yhigh-ylow)* 
     +               (((xval-xlow)/(xhigh-xlow))**pow)
      end if

      return
      end
