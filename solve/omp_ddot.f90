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

!>    @brief compute ([x]^T*[y]) distributed and build a global sum, (OpenMP version)
!>    @param[in] N length of all vectors
!>    @param[in] X vector [x]
!>    @param[in] Y vector [y]
!>    @param[in] sh_help openmp-shared vector [# threads *1]
!>    @param[out] S solution ( [x]^T*[y] )
      SUBROUTINE omp_ddot(n,x,y,s,sh_help)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N    : length of all vector x,y
        INTEGER n
!     lsum : local-summary
!     lsum : source summary
!     vectors [x], [y]
        DOUBLE PRECISION lsum(1), slsum(1), x(n), y(n), s, sh_help(*)
        DOUBLE PRECISION qddot
        EXTERNAL qddot
        DOUBLE PRECISION ddot
        EXTERNAL ddot

!     compute local-sum
#ifdef USE_QDDOT
        lsum(1) = qddot(n,x,1,y,1)
#else
        lsum(1) = ddot(n,x,1,y,1)
#endif
!     compute global sum (OpenMP)
        CALL xsum_0(1,lsum,slsum,sh_help)
        s = slsum(1)
        RETURN
      END

!>    @brief compute ([x]^T*[y]) and build the sum, serial (no OpenMP) implementation
!>    @param[in] N length of all vectors
!>    @param[in] X vector [x]
!>    @param[in] Y vector [y]
!>    @param[out] S solution ( [x]^T*[y] )
      SUBROUTINE s_ddot(n,x,y,s)
        IMPLICIT NONE
!     N    : length of all vector x,y
        INTEGER n
!     vectors [x], [y]
        DOUBLE PRECISION x(n), y(n), s
        DOUBLE PRECISION qddot
        EXTERNAL qddot
        DOUBLE PRECISION ddot
        EXTERNAL ddot

!     compute local-sum
#ifdef USE_QDDOT
        s = qddot(n,x,1,y,1)
#else
        s = ddot(n,x,1,y,1)
#endif
        RETURN
      END

!>    @brief 2x compute ([x]^T*[y]) distributed and build a global sum, (OpenMP version)
!>    @param[in] N length of all vectors
!>    @param[in] X vector [x]
!>    @param[in] X2 vector [x]
!>    @param[in] Y vector [y]
!>    @param[in] Y2 vector [y]
!>    @param[in] sh_help openmp-shared vector [# threads *2]
!>    @param[out] S solution ( [x]^T*[y] )
!>    @param[out] S2 solution ( [x]^T*[y] )
!>    @details
!>    Function : compute ([x]^T*[y]) distributed and build a global sum,\n
!>               2 products at once\n
      SUBROUTINE omp_2ddot(n,x,y,s,x2,y2,s2,sh_help)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N    : length of all vector x,y
        INTEGER n
!     lsum : local-summary
!     lsum : source summary
!     vectors [x], [y]
        DOUBLE PRECISION lsum(2), slsum(2), sh_help(*)
        DOUBLE PRECISION x(max(n,1)), y(max(n,1)), s
        DOUBLE PRECISION x2(max(n,1)), y2(max(n,1)), s2
        INTEGER von, bis
        DOUBLE PRECISION qddot
        EXTERNAL qddot
        DOUBLE PRECISION ddot
        EXTERNAL ddot

        lsum(1) = 0.0D0
        lsum(2) = 0.0D0
!     blocking
        DO von = 1, n, int(bl_size/bldiv_dot(1))
          bis = min(n,von+int(bl_size/bldiv_dot(1))-1)
!        compute local-sum
#ifdef USE_QDDOT
          lsum(1) = lsum(1) + qddot(bis-von+1,x(von),1,y(von),1)
          lsum(2) = lsum(2) + qddot(bis-von+1,x2(von),1,y2(von),1)
#else
          lsum(1) = lsum(1) + ddot(bis-von+1,x(von),1,y(von),1)
          lsum(2) = lsum(2) + ddot(bis-von+1,x2(von),1,y2(von),1)
#endif
        END DO
!     compute global sum (OpenMP)
        CALL xsum_0(2,lsum,slsum,sh_help)
        s = slsum(1)
        s2 = slsum(2)
        RETURN
      END

!>    @brief 3x compute ([x]^T*[y]) distributed and build a global sum, (OpenMP version)
!>    @param[in] N length of all vectors
!>    @param[in] X vector [x]
!>    @param[in] X2 vector [x]
!>    @param[in] X3 vector [x]
!>    @param[in] Y vector [y]
!>    @param[in] Y2 vector [y]
!>    @param[in] Y3 vector [y]
!>    @param[in] sh_help openmp-shared vector [# threads *3]
!>    @param[out] S solution ( [x]^T*[y] )
!>    @param[out] S2 solution ( [x]^T*[y] )
!>    @param[out] S3 solution ( [x]^T*[y] )
!>    @details
!>    Function : compute ([x]^T*[y]) distributed and build a global sum,\n
!>               3 products at once\n
      SUBROUTINE omp_3ddot(n,x,y,s,x2,y2,s2,x3,y3,s3,sh_help)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N    : length of all vector x,y
        INTEGER n
!     lsum : local-summary
!     lsum : source summary
!     vectors [x], [y]
        DOUBLE PRECISION lsum(3), slsum(3), sh_help(*)
        DOUBLE PRECISION x(max(n,1)), y(max(n,1)), s
        DOUBLE PRECISION x2(max(n,1)), y2(max(n,1)), s2
        DOUBLE PRECISION x3(max(n,1)), y3(max(n,1)), s3
        INTEGER von, bis
        DOUBLE PRECISION qddot
        EXTERNAL qddot
        DOUBLE PRECISION ddot
        EXTERNAL ddot

        lsum(1) = 0.0D0
        lsum(2) = 0.0D0
        lsum(3) = 0.0D0
!     blocking
        DO von = 1, n, int(bl_size/bldiv_dot(2))
          bis = min(n,von+int(bl_size/bldiv_dot(2))-1)
!        compute local-sum
#ifdef USE_QDDOT
          lsum(1) = lsum(1) + qddot(bis-von+1,x(von),1,y(von),1)
          lsum(2) = lsum(2) + qddot(bis-von+1,x2(von),1,y2(von),1)
          lsum(3) = lsum(3) + qddot(bis-von+1,x3(von),1,y3(von),1)
#else
          lsum(1) = lsum(1) + ddot(bis-von+1,x(von),1,y(von),1)
          lsum(2) = lsum(2) + ddot(bis-von+1,x2(von),1,y2(von),1)
          lsum(3) = lsum(3) + ddot(bis-von+1,x3(von),1,y3(von),1)
#endif
        END DO
!     compute global sum (OpenMP)
        CALL xsum_0(3,lsum,slsum,sh_help)
        s = slsum(1)
        s2 = slsum(2)
        s3 = slsum(3)
        RETURN
      END

!>    @brief 4x compute ([x]^T*[y]) distributed and build a global sum, (OpenMP version)
!>    @param[in] N length of all vectors
!>    @param[in] X vector [x]
!>    @param[in] X2 vector [x]
!>    @param[in] X3 vector [x]
!>    @param[in] X4 vector [x]
!>    @param[in] Y vector [y]
!>    @param[in] Y2 vector [y]
!>    @param[in] Y3 vector [y]
!>    @param[in] Y4 vector [y]
!>    @param[in] sh_help openmp-shared vector [# threads *4]
!>    @param[out] S solution ( [x]^T*[y] )
!>    @param[out] S2 solution ( [x]^T*[y] )
!>    @param[out] S3 solution ( [x]^T*[y] )
!>    @param[out] S4 solution ( [x]^T*[y] )
!>    @details
!>    Function : compute ([x]^T*[y]) distributed and build a global sum,\n
!>               4 products at once\n
      SUBROUTINE omp_4ddot(n,x,y,s,x2,y2,s2,x3,y3,s3,x4,y4,s4,sh_help)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N    : length of all vector x,y
        INTEGER n
!     lsum : local-summary
!     lsum : source summary, don't do that
!     vectors [x], [y]
        DOUBLE PRECISION lsum(4), slsum(4), sh_help(*)
        DOUBLE PRECISION x(max(n,1)), y(max(n,1)), s
        DOUBLE PRECISION x2(max(n,1)), y2(max(n,1)), s2
        DOUBLE PRECISION x3(max(n,1)), y3(max(n,1)), s3
        DOUBLE PRECISION x4(max(n,1)), y4(max(n,1)), s4
        INTEGER von, bis
        DOUBLE PRECISION qddot
        EXTERNAL qddot
        DOUBLE PRECISION ddot
        EXTERNAL ddot

        lsum(1) = 0.0D0
        lsum(2) = 0.0D0
        lsum(3) = 0.0D0
        lsum(4) = 0.0D0
!     blocking
        DO von = 1, n, int(bl_size/bldiv_dot(3))
          bis = min(n,von+int(bl_size/bldiv_dot(3))-1)
!        compute local-sum
#ifdef USE_QDDOT
          lsum(1) = lsum(1) + qddot(bis-von+1,x(von),1,y(von),1)
          lsum(2) = lsum(2) + qddot(bis-von+1,x2(von),1,y2(von),1)
          lsum(3) = lsum(3) + qddot(bis-von+1,x3(von),1,y3(von),1)
          lsum(4) = lsum(4) + qddot(bis-von+1,x4(von),1,y4(von),1)
#else
          lsum(1) = lsum(1) + ddot(bis-von+1,x(von),1,y(von),1)
          lsum(2) = lsum(2) + ddot(bis-von+1,x2(von),1,y2(von),1)
          lsum(3) = lsum(3) + ddot(bis-von+1,x3(von),1,y3(von),1)
          lsum(4) = lsum(4) + ddot(bis-von+1,x4(von),1,y4(von),1)
#endif
        END DO
!     compute global sum (OpenMP)
        CALL xsum_0(4,lsum,slsum,sh_help)
        s = slsum(1)
        s2 = slsum(2)
        s3 = slsum(3)
        s4 = slsum(4)
        RETURN
      END
