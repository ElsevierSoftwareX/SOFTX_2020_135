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

!>    @brief probe for break criteria, (OpenMP version)
!>    @param[in] enough value for ||[res]|| or max(abs(res)), see 'criteria'
!>    @param[in] iter count iterations
!>    @param[in] max_It maximum of iterations, counted with 'iter'
!>    @param[in] depsilon precision criteria to break iterations
!>    @param[in] need_Ax switch to compute an extra MVP:([A]x[x]) in [ax], see 'criteria'
!>    @param[in] criteria switch to set when should break\n
!>    - 0 : relative stopping crit. : ||[res]|| < depsilon*||[res0]||\n
!>    - 1 : absolute stopping crit. : ||[res]|| < depsilon\n
!>    - 2 : maximum  stopping crit. : max(abs([res])) < depsilon\n
!>    - 3 : abs. and rel. stopping crit. : ( ||[res]|| < depsilon ) and ( ||[res]|| < 0.5d0*||[res0]|| )\n
!>    0.5d0 is a constant for testing only and is named 'minRel'\n
!>    first [res] ^= [r_0] and later (if precise enough) : [res] ^= ([A]x[x]-[b])
!>    @param[in] res0 res0 ^= ||[res0]||, should be start residue
!>    @param[in] divide_zero (>0) division by zero ?
!>    @param[in,out] e_count counts the number of reaching the criteria
!>    @param[in,out] e_old last good `enough` value
!>    @return false : continue (not enough precision), true : break - enough precision or max_It reached
      LOGICAL FUNCTION omp_abbruch(enough,iter,max_it,depsilon, &
          need_ax,criteria,res0,divide_zero,e_count,e_old)
        use mod_OMP_TOOLS
        use mod_linfos
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     break with enough precision
        DOUBLE PRECISION depsilon, enough
        DOUBLE PRECISION res0, e_old
        INTEGER e_count
!     see above
        DOUBLE PRECISION minrel
        PARAMETER (minrel=0.99D0)
!     need_Ax  : switch to compute an extra MVP:([A]x[x]) in [ax]
!     criteria : switch to break
        LOGICAL need_ax, flag_abbruch
        INTEGER criteria
        character (len=3) :: rstr
!     iter   : count iterations,
!     max_It : max iterations,
!     divide_zero <> 0 : break should be done too
        INTEGER iter, max_it, divide_zero
        INTRINSIC dsqrt


        rstr = '[r]'
        IF (need_ax) rstr = '(R)'

!     default -> break iterations
        flag_abbruch = .TRUE.

!     max-iterations reached ?
        IF (iter>max_it) GO TO 100

!     same result the last three times ?
        IF ((e_count>3) .AND. (e_old==enough)) THEN
!$OMP master
          IF (linfos(4)>=1) WRITE(*,*) &
            ' Warning : stopping before precision reached'
!$OMP end master
          GO TO 100
        END IF

!     enough precision ?

        IF (criteria==1) THEN
!        absolute
!        write(*,*) "Testing absolute ",enough,"<=",depsilon*depsilon,"?"
          IF (enough<=(depsilon*depsilon)) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        ELSE IF (criteria==2) THEN
!        maximum
!        write(*,*) "Testing maximum ",enough,"<=",depsilon,"?"
        IF (enough<=depsilon) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        ELSE IF (criteria==3) THEN
!        absolute and relative
!        write(*,*) "Testing absolut and relative ",enough,"<=",res0*res0*minrel*minrel," and",depsilon*depsilon,"?"
        IF ((enough<=(res0*res0*minrel*minrel)) .AND. (enough<=( &
              depsilon*depsilon))) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        ELSE IF (criteria==0) THEN
!        relative
!        write(*,*) "Testing relative ",enough,"<=",depsilon*res0*depsilon*res0
          IF (enough<=(depsilon*res0*depsilon*res0)) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        END IF

!     continue, no break
        IF (divide_zero==0) flag_abbruch = .FALSE.

        IF (e_old==enough) THEN
          e_count = e_count + 1
        ELSE
          e_old = enough
          e_count = 0
        END IF

!     stop label
100     CONTINUE

333     FORMAT (3A,1D20.13,1A)
!$OMP master
!     print latest
        IF ((linfos(4)>=1) .AND. (flag_abbruch)) THEN
          IF (criteria==2) THEN
            WRITE(*,333,advance='NO') '  damax', rstr, ' =', enough, &
              ' '
          ELSE
            WRITE(*,333,advance='NO') '  nrm2', rstr, ' =', &
              dsqrt(enough), ' '
          END IF
        ELSE IF (linfos(4)>=3 .AND. linfos(4)/=100) THEN
          IF (criteria==2) THEN
            WRITE(*,333,advance='NO') '  damax', rstr, ' =', enough, &
              ','
          ELSE
            WRITE(*,333,advance='NO') '  nrm2', rstr, ' =', &
              dsqrt(enough), ','
          END IF
        END IF
!$OMP end master

        omp_abbruch = flag_abbruch
        RETURN
      END

!>    @brief probe for break criteria, serial (no OpenMP) implementation of "omp_abbruch", see above
!>    @param[in] enough value for ||[res]|| or max(abs(res)), see 'criteria'
!>    @param[in] iter count iterations
!>    @param[in] max_It maximum of iterations, counted with 'iter'
!>    @param[in] depsilon precision criteria to break iterations
!>    @param[in] need_Ax switch to compute an extra MVP:([A]x[x]) in [ax], see 'criteria'
!>    @param[in] criteria switch to set when should break\n
!>    - 0 : relative stopping crit. : ||[res]|| < depsilon*||[res0]||\n
!>    - 1 : absolute stopping crit. : ||[res]|| < depsilon\n
!>    - 2 : maximum  stopping crit. : max(abs([res])) < depsilon\n
!>    - 3 : abs. and rel. stopping crit. : ( ||[res]|| < depsilon ) and ( ||[res]|| < 0.5d0*||[res0]|| )\n
!>    0.5d0 is a constant for testing only and is named 'minRel'\n
!>    first [res] ^= [r_0] and later (if precise enough) : [res] ^= ([A]x[x]-[b])
!>    @param[in] res0 res0 ^= ||[res0]||, should be start residue
!>    @param[in] divide_zero (>0) division by zero ?
!>    @param[in,out] e_count counts the number of reaching the criteria
!>    @param[in,out] e_old last good `enough` value
!>    @return false : continue (not enough precision), true : break - enough precision or max_It reached
      LOGICAL FUNCTION s_abbruch(enough,iter,max_it,depsilon,need_ax, &
          criteria,res0,divide_zero,e_count,e_old)
        use mod_linfos
        IMPLICIT NONE
!     break with enough precision
        DOUBLE PRECISION depsilon, enough
        DOUBLE PRECISION res0, e_old
        INTEGER e_count
!     see above
        DOUBLE PRECISION minrel
        PARAMETER (minrel=0.99D0)
!     need_Ax  : switch to compute an extra MVP:([A]x[x]) in [ax]
!     criteria : switch to break
        LOGICAL need_ax, flag_abbruch
        INTEGER criteria
        character (len=3) :: rstr
!     iter   : count iterations,
!     max_It : max iterations,
!     divide_zero <> 0 : break should be done too
        INTEGER iter, max_it, divide_zero
        INTRINSIC dsqrt


        rstr = '[r]'
        IF (need_ax) rstr = '(R)'

!     default -> break iterations
        flag_abbruch = .TRUE.

!     max-iterations reached ?
        IF (iter>max_it) GO TO 100

!     same result the last three times ?
        IF ((e_count>3) .AND. (e_old==enough)) THEN
          IF (linfos(4)>=1) WRITE(*,*) &
            ' Warning : stopping before precision reached'
          GO TO 100
        END IF

!     enough precision ?

        IF (criteria==1) THEN
!        absolute
          IF (enough<=(depsilon*depsilon)) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        ELSE IF (criteria==2) THEN
!        maximum
          IF (enough<=depsilon) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        ELSE IF (criteria==3) THEN
!        absolute and relative
          IF ((enough<=(res0*res0*minrel*minrel)) .AND. (enough<=( &
              depsilon*depsilon))) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        ELSE IF (criteria==0) THEN
!        relative
          IF (enough<=(depsilon*res0*depsilon*res0)) THEN
            IF (need_ax) GO TO 100
            need_ax = .TRUE.
          END IF
        END IF

!     continue, no break
        IF (divide_zero==0) flag_abbruch = .FALSE.

        IF (e_old==enough) THEN
          e_count = e_count + 1
        ELSE
          e_old = enough
          e_count = 0
        END IF

!     stop label
100     CONTINUE

333     FORMAT (3A,1D20.13,1A)
!     print latest
        IF ((linfos(4)>=1) .AND. (flag_abbruch)) THEN
          IF (criteria==2) THEN
            WRITE(*,333,advance='NO') '  damax', rstr, ' =', enough, &
              ' '
          ELSE
            WRITE(*,333,advance='NO') '  nrm2', rstr, ' =', &
              dsqrt(enough), ' '
          END IF
        ELSE IF (linfos(4)>=3 .AND. linfos(4)/=100) THEN
          IF (criteria==2) THEN
            WRITE(*,333,advance='NO') '  damax', rstr, ' =', enough, &
              ','
          ELSE
            WRITE(*,333,advance='NO') '  nrm2', rstr, ' =', &
              dsqrt(enough), ','
          END IF
        END IF

        s_abbruch = flag_abbruch
        RETURN
      END
