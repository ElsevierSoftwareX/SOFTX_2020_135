! MIT License
!
! Copyright (c) 2020 SHEMAT-Suite
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!>    @brief adaptive relaxation parameter
!>    @param[in] iter iteration counter
!>    @param[in] r <description>
!>    @param[in] rold <description>
!>    @param[in] emax <description>
!>    @param[out] relax current relaxation value
!>    @param[in,out] relaxold old relaxation value
!>    @param[in] ismpl local sample index
!>    @details
!> calculate adaptive relaxation parameter\n
!> see Cooley (1983), WRR 19(5), 1271-1285\n
      SUBROUTINE nl_relax(iter,r,rold,emax,relax,relaxold,ismpl)
        use mod_genrl
        use mod_genrlc
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        INTEGER iter
        DOUBLE PRECISION r, rold, emax, relax, relaxold, s

!  step 1
        IF (iter==1) THEN
          s = 1.0D0
        ELSE
          s = r/(relaxold*rold)
        END IF
        WRITE(*,*) 'r, rold =', r, rold, '    s= ', s
! step 2
        IF (s>=-1.D0) THEN
          relax = (3.D0+s)/(3.D0+abs(s))
        ELSE
          relax = 0.5D0/abs(s)
        END IF
! step 3
        IF (relax*abs(r)>emax) relax = emax/abs(r)
        relax = min(1.D0,max(0.05D0,relax))
        WRITE(*,*) '**** relax= ', relax, ' (', relaxold, ')'
!
        relaxold = relax
!
        IF (linfos(3)>=1) THEN
          WRITE(*,'(1a,1e10.3)') ' relaxation factor: ', relax
        END IF
!
        RETURN
      END
