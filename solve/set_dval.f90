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

!>    @brief initialise a double precision vector
!>    @param[in] N length of [x]-vectors
!>    @param[in] alpha value for initialisation
!>    @param[out] x initialised vector [x]
      SUBROUTINE set_dval(n,alpha,x)
        IMPLICIT NONE
!     N : length of [x]-vector
!     i : loop variable
        integer :: n, i
!     vector [x]
        double precision, dimension (n) :: x
        double precision :: alpha

        DO i = 1, n
          x(i) = alpha
        END DO
        RETURN
      END

!>    @brief initialise a double precision vector, (OpenMP version)
!>    @param[in] N length of [x]-vectors
!>    @param[in] alpha value for initialisation
!>    @param[out] x initialised vector [x]
      SUBROUTINE omp_set_dval(n,alpha,x)
        IMPLICIT NONE
!     N : length of [x]-vector
!     i : loop variable
        integer :: n, i, tpos, tanz
!     vector [x]
        double precision, dimension (n) :: x
        double precision :: alpha

        CALL omp_part(n,tpos,tanz)
        DO i = tpos, tpos + tanz - 1
          x(i) = alpha
        END DO
        RETURN
      END
