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
      subroutine strlen(str,MAXLEN,lostr)
c-----------------------------------------------------------------------
c
c      Determine the length of the string minus trailing blanks
c
c
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character (len=132) :: str
      lostr = MAXLEN
      do i=1,MAXLEN
            j = MAXLEN - i + 1
            if(str(j:j).ne.' ') return
            lostr = lostr - 1
      end do
      return
      end
