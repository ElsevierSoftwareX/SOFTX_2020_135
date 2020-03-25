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

!>    @brief init the seed matrix for divided difference method
!>    @param[in] seedi seeding component/index
!>    @param[in,out] hh step size h
!>    @param[out] hh_bak backup value of the step size h
!>    @param[in] ismpl local sample index
!>    @details
!>      kx-,ky-,kz-,lx-,ly-,lz-,q-units and more are seeding parameters\n
      SUBROUTINE dd_seeding(seedi,hh,hh_bak,ismpl)
        use arrays
#ifndef AD_RM
        use g_arrays
#else
        use arrays_ad
#endif
        use mod_genrl
        use mod_inverse
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        INTEGER seedi, s_k, s_u, p_offs, p_end
        DOUBLE PRECISION hh, hh_bak, dnrm2
        DOUBLE PRECISION get_optip
        EXTERNAL get_optip
        LOGICAL test_null
        EXTERNAL dnrm2, test_null

        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error: in "dd_seeding", no valid "seeding=', &
            seedi, '"'
          STOP
        END IF
!
        hh_bak = get_optip(seedi,ismpl)
        IF ( .NOT. test_null(hh_bak)) THEN
          CALL set_optip(seedi,hh_bak*(1.0D0+hh),ismpl)
        ELSE
          IF (hh==1.0D-2) hh = 1.0D-20
          CALL set_optip(seedi,hh,ismpl)
          WRITE(*,'(A,1e16.8,A)') 'warning: step size h=', hh, &
            ' may be to hard, please choose a nonzero apriori value !'
        END IF
!
        RETURN
      END
