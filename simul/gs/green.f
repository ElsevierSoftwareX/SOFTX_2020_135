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
      subroutine green(value,hexrep,gfrac)
c-----------------------------------------------------------------------
c
c Provided with a real value ``value'' this subroutine returns the green
c portion of the color specification.
c
c Note common block "color" and call to "hexa"
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include 'gslib.inc'
      double precision value
      character (len=2) :: hexrep
      character (len=2) :: hexa
      hexrep = '00'
      if(value.lt.cint(1))then
c
c Scale it between (0,0):
c
            integ  = 0
      else if((value.ge.cint(1)).and.(value.lt.cint(2)))then
c
c Scale it between (0,255):
c
            integ = int((value-cint(1))/(cint(2)-cint(1))*255.)
            if(integ.gt.255) integ = 255
            if(integ.lt.0)   integ = 0
      else if((value.ge.cint(2)).and.(value.lt.cint(4)))then
c
c Scale it between (255,255):
c
            integ  = 255
      else if(value.ge.cint(4))then
c
c Scale it between (255,0):
c
            integ=int((ccmax-value)/(ccmax-cint(4))*255.)
            if(integ.gt.255) integ = 255
            if(integ.lt.0)   integ = 0
      end if
c
c Establish coding and return:
c
      gfrac  = dble(integ) / 255.
      hexrep = hexa(integ)
      return
      end
