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
      subroutine chknam(str,len)
c-----------------------------------------------------------------------
c
c                   Check for a Valid File Name
c                   ***************************
c
c This subroutine takes the character string "str" of length "len" and
c removes all leading blanks and blanks out all characters after the
c first blank found in the string (leading blanks are removed first).
c
c
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (MAXLEN=132)
      character (len=1), dimension (MAXLEN) ::str
c
c Remove leading blanks:
c
      do i=1,len-1
            if(str(i).ne.' ') then
                  if(i.eq.1) go to 1
                  do j=1,len-i+1
                        k = j + i - 1
                        str(j) = str(k)
                  end do
                  do j=len,len-i+2,-1
                        str(j) = ' '
                  end do
                  go to 1
            end if
      end do
 1    continue
c
c Find first blank and blank out the remaining characters:
c
      do i=1,len-1
            if(str(i).eq.' ') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 2
            end if
      end do
 2    continue
c
c Return with modified file name:
c
      return
      end
