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

!>    @brief gives start position for proc. 'my_rank'
!>    @param[in] N number of elements
!>    @param[in] my_rank my rank to part the vector (lengths given by N)
!>    @param[in] omp_P maximum rank to part the vector
!>    @return start position of elements for this rank
      INTEGER FUNCTION ppos(n,my_rank,omp_p)
        IMPLICIT NONE

        INTEGER n
!     OMP stuff
        INTEGER my_rank, rank, omp_p
!     tmp variables
        INTEGER delta, schwelle, iall

!     checking border
        rank = my_rank
        IF (rank>(omp_p-1)) rank = omp_p - 1
        IF (rank<0) rank = 0

        delta = n/omp_p
        schwelle = n - delta*omp_p
        iall = delta*rank + rank + 1
        IF ((rank+1)>schwelle) iall = iall + schwelle - rank

        ppos = iall

        RETURN
      END

!>    @brief gives number of elements for proc. 'my_rank'
!>    @param[in] N number of elements
!>    @param[in] my_rank my rank to part the vector (lengths given by N)
!>    @param[in] omp_P maximum rank to part the vector
!>    @return number of elements for this rank
      INTEGER FUNCTION panz(n,my_rank,omp_p)
        IMPLICIT NONE

        INTEGER n
!     OMP stuff
        INTEGER rank, my_rank, omp_p
!     tmp variables
        INTEGER delta, schwelle

!     checking border
        rank = my_rank
        IF (rank>(omp_p-1)) rank = omp_p - 1
        IF (rank<0) rank = 0

        delta = n/omp_p
        schwelle = n - delta*omp_p
        IF ((rank+1)<(schwelle+1)) delta = delta + 1

        panz = delta

        RETURN
      END

!>    @brief compute position an number of local elements (OpenMP based)
!>    @param[in] N number of elements to parting
!>    @param[out] tpos start position of local elements
!>    @param[out] tanz local number of elements
      SUBROUTINE omp_part(n,tpos,tanz)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

        INTEGER n, panz, ppos
        EXTERNAL panz, ppos
!     thread stuff
        INTEGER tpos, tanz
        INTEGER my_thd, thd_p

        my_thd = omp_get_his_thread_num()
        thd_p = omp_get_num_of_threads()
        tanz = panz(n,my_thd,thd_p)
        tpos = ppos(n,my_thd,thd_p)

        RETURN
      END
