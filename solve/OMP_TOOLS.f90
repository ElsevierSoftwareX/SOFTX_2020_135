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

!>    @brief OpenMP wrapper for "omp_get_thread_num"
!>    @return thread index number
      INTEGER FUNCTION omp_get_his_thread_num()
        IMPLICIT NONE
        INTEGER n
!$      integer (kind=4) :: OMP_GET_THREAD_NUM
!$      external OMP_GET_THREAD_NUM

        n = 0
!$      N=OMP_GET_THREAD_NUM()
        omp_get_his_thread_num = n

        RETURN
      END

!>    @brief OpenMP wrapper for "omp_get_num_threads"
!>    @return number of threads
      INTEGER FUNCTION omp_get_num_of_threads()
        IMPLICIT NONE
        INTEGER n
!$      integer (kind=4) :: OMP_GET_NUM_THREADS
!$      external OMP_GET_NUM_THREADS

        n = 1
!$      N=OMP_GET_NUM_THREADS()
!AW-TESTC$      N=OMP_GET_MAX_THREADS()
        omp_get_num_of_threads = n

        RETURN
      END
