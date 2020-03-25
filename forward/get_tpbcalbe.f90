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

!>    @brief time depended boundary condition modificators
!>    @param[out] malfa alfa modificator
!>    @param[out] mbeta beta modificator
!>    @param[in] tpbcu time period BC table index
!>    @param[in] ismpl local sample index
!>    @details
!> "GET Time Periods Boundary Condition ALfa & BEta"\n
!> get the alfa and beta modificators for time dependend bc-values\n
      SUBROUTINE get_tpbcalbe(malfa,mbeta,tpbcu,ismpl)
        use arrays
        use mod_genrl
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: k

        INTEGER tpbcu, imt
        DOUBLE PRECISION malfa, mbeta, mtime
        INTRINSIC abs

!       default - when not time depended
        malfa = 0.0D0
        mbeta = 1.0D0
!
        IF (tpbcu>0) THEN
          k = 0
          imt = 0
!         next bc-tp entry
100       imt = imt + 1
          mtime = bcperiod(imt,1,tpbcu,ismpl)
          IF (mtime<=simtime(ismpl)) THEN
            malfa = bcperiod(imt,2,tpbcu,ismpl)
            mbeta = bcperiod(imt,3,tpbcu,ismpl)
!           save index of the valid start time
            k = imt
          END IF
!
!         !!! this IF statement (with the GOTO) needs to be outside of the
!             other "mtime<=simtime(ismpl)" scopes, to generate reverse-mode code !!!
          IF (mtime<=simtime(ismpl).AND.imt<ibcperiod(tpbcu)) GO TO 100
          IF (k>0) THEN
            IF ( .NOT. lbcperiod(k,tpbcu)) THEN
!              disable this BC (switched off)
              tpbcu = -abs(tpbcu)
            END IF
          END IF
        END IF
!
        RETURN
      END
