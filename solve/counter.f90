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

!>    @brief initialisation of the floating point performance counter (only for benchmarking)
      SUBROUTINE dp_init()
        IMPLICIT NONE
        INTEGER dp
        COMMON /floating_points/dp
        dp = 0
        RETURN
      END

!>    @brief increase the floating point performance counter (only for benchmarking)
!>    @param[in] n additional number of floating point operations
      SUBROUTINE dp_count(n)
        IMPLICIT NONE
        INTEGER n, dp
        COMMON /floating_points/dp
!$OMP atomic
        dp = dp + n
        RETURN
      END

!>    @brief write out the floating point performance (only for benchmarking)
!>    @param[in] time run time during benchmark
      SUBROUTINE dp_comp(time)
        IMPLICIT NONE
        DOUBLE PRECISION time
        INTEGER dp
        COMMON /floating_points/dp
        INTRINSIC dble
        WRITE(*,'(1A,1F10.2)') ' [I] : MBytes/sec = ', &
          dble(dp*8)/(time*1.0D6)
        RETURN
      END
