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

!>    @brief compute reference value of "rhof"
!>    @param[in] ismpl local sample index
!>    @details
!> recompute the reference value [rref] and\n
!> makes a sanity proof with the given (default) one\n
      SUBROUTINE check_props(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION rmin, rmax, ravrg, rlocal, rhof
        EXTERNAL rhof

#ifdef head_base
!     check rref
        rmin = rhof(1,1,1,ismpl)
        rmax = rhof(1,1,1,ismpl)
        ravrg = 0.D0
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              rlocal = rhof(i,j,k,ismpl)
              ravrg = ravrg + rlocal
              rmin = min(rlocal,rmin)
              rmax = max(rlocal,rmax)
            END DO
          END DO
        END DO
        ravrg = ravrg/dble(i0*j0*k0)
!
        IF (rref>rmax) THEN
          WRITE(*,'(1A,1e12.4,1A)') 'error: rref=', rref, &
            ' not consistent'
          WRITE(*,'(3A,4(1e12.4,1A))') &
            'You should explicitely define a ',&
            'section "'//key_char//' rref" with ', &
            'a value like ', rmax, &
            ', calulated range [', rmin, ', ', ravrg, ', ', rmax, &
            '] !!!'
          rref=rmax
!vr         rref=ravrg
          WRITE(*,'(1A,1e12.4,1A)') 'rref set to:', rref, &
            ' (= rmax)'
        END IF
!
#endif

        RETURN
      END
