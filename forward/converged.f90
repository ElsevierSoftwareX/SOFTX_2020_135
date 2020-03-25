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

!>    @brief compute average for break condition
!>    @param[in] enough test for break
!>    @param[in] depsilon break condition
!>    @param[in] Htype history (number)
!>    @param[in] maxhistlen max history length
!>    @param[in] history HISTORY buffer
!>    @param[in] hlen History Length
!>    @param[in] ipos Position Index
!>    @param[in] Hmax history length
!>    @return "true" when converged
!>    @details
!> proof "enough" < "depsilon" or the average of the last "hlen" steps\n
      LOGICAL FUNCTION converged(enough,depsilon,htype,maxhistlen,hmax,history,hlen,ipos)

        USE mod_genrl
        use mod_linfos

        IMPLICIT NONE

        ! step diff. or precision
        DOUBLE PRECISION enough

        ! break condition
        DOUBLE PRECISION depsilon

        ! max history length
        INTEGER maxhistlen

        ! number of probes for each history length
        INTEGER maxprob
        PARAMETER (maxprob=4)

        ! History (selection, counter, number)
        INTEGER htype, h, hmax

        ! HISTORY buffer
        DOUBLE PRECISION history(maxhistlen,hmax)

        ! History Length
        INTEGER hlen(hmax)

        ! Position Index
        INTEGER ipos(hmax)

        ! locale SUM for average computation
        DOUBLE PRECISION lsum, smin, smax

        ! loop Index
        INTEGER i

        ! Wrong htype values
        IF (htype==0) THEN
          WRITE(*,*) &
            'error: "history type"=0 in "converged" not allowed !'
          STOP
        END IF
        IF (htype>hmax) THEN
          WRITE(*,*) 'error: "history number" (Hmax<', htype, &
            ') in "converged" to low !'
          STOP
        END IF

        ! nonlinear convergence test is disabled by input file
        IF (.NOT. nlconverge .eq. 0) THEN
           converged = .FALSE.
           RETURN
        END IF

        ! quick test: maximum difference smaller than tolerance
        IF (abs(enough)<depsilon) THEN
          converged = .TRUE.
          RETURN
        END IF

        ! set defaults
        converged = .FALSE.
        h = htype
        IF (htype<0) THEN
          ! init iterations, then end
          h = -htype
          hlen(h) = 1
          ipos(h) = -1
          RETURN
        END IF

        ! break, if history buffer too small
        IF (ipos(h)>=maxprob*maxhistlen) THEN
          WRITE(*,*) 'warning: history buffer in "converged.f"', &
            ' to small ! <break outer iterations>'
          converged = .TRUE.
        END IF

        ! fill history with new element
        IF (ipos(h)>=maxprob*hlen(h)) THEN

          ! increase history length
          hlen(h) = hlen(h) + 1
          ipos(h) = -1
          history(hlen(h),h) = enough

        ELSE

          ! Add position
          ipos(h) = ipos(h) + 1
          history(1+mod(ipos(h),hlen(h)),h) = enough

        END IF

        ! compute sum for average
        lsum = 0.0D0
        smin = history(1,h)
        smax = history(1,h)
        DO i = 1, hlen(h)
          lsum = lsum + history(i,h)
          smin = min(smin,history(i,h))
          smax = max(smax,history(i,h))
        END DO

!     can we break iterations (enough precision)
        IF ((abs(lsum)/dble(hlen(h)))<depsilon) converged = .TRUE.

        IF ((abs(enough)>depsilon) .AND. (converged)) THEN
          IF (h==1) WRITE(*,'(2(A,1e15.8),A,I3)') &
            'warning: oscillating convergence in HEAD, [min:max] = [', &
            smin, ':', smax, '], period length=', hlen(h)
          IF (h==2) WRITE(*,'(2(A,1e15.8),A,I3)') &
            'warning: oscillating convergence in TEMP, [min:max] = [', &
            smin, ':', smax, '], period length=', hlen(h)
          IF (h>=4) WRITE(*,'(2(A,1e15.8),A,I3)') &
            'warning: oscillating convergence in CONC, [min:max] = [', &
            smin, ':', smax, '], period length=', hlen(h)
        END IF

        RETURN
      END
