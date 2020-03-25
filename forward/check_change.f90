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

!>    @brief check changes between vectors [new] and [old]
!>    @param[in] mode switch absolute/relative
!>    @param[out] rms return value
!>    @param[out] difmax maximal difference
!>    @param[in] ni I-dimension for vectors [new], [old]
!>    @param[in] nj J-dimension for vectors [new], [old]
!>    @param[in] nk K-dimension for vectors [new], [old]
!>    @param[in] new vector with new values
!>    @param[in] old vector with old values
!>    @param[in] pv_idx index number (physical value), only needed for AD code generation/modification
!>    @param[out] loc_nltol tolerance criteria, only needed for AD code generation/modification
!>    @param[in] ismpl local sample index
!>    @details
!> Computes two difference metrics between the physical variable array
!> of this iteration (new) and the one from the previous iteration
!> (old): \n\n
!>
!> 1. difmax  : max. difference of fields new and old\n
!>    mode= 0 : difmax is absolute maximum difference\n
!>    else    : difmax is relative maximum difference\n\n
!>
!> 2. rms     : root mean square difference of fields new and old
!>    mode= 0 : rms is root mean square of absolute differences\n
!>    else    : rms is root mean square of relative differences\n\n
      SUBROUTINE check_change(mode,pv_idx,loc_nltol,rms,difmax,ni,nj,nk,new,old,ismpl)

        use arrays
        use mod_linfos

        IMPLICIT NONE

        INTEGER i, j, k, ni, nj, nk, ijk, ipt, jpt, kpt, mode, pv_idx, ismpl
        DOUBLE PRECISION new(ni,nj,nk), old(ni,nj,nk)
        ! DOUBLE PRECISION dif
        DOUBLE PRECISION rms, difmax, loc_nltol
        INTEGER idamax
        EXTERNAL idamax
        INTRINSIC dabs, dble, sqrt


        IF (linfos(3)>=2) WRITE(*,*) ' ... check_change'

        ! Initial values for output
        rms = 0.D0
        difmax = 0.D0
        ipt = 0
        jpt = 0
        kpt = 0

        IF (mode==0) THEN
          ! Absolute difference

          ! Number of cells
          ijk = ni*nj*nk

          ! Copy new values to x
          CALL dcopy(ijk,new,1,x(1,1,1,ismpl),1)

          ! Absolute difference array: new - old
          CALL daxpy(ijk,-1.0d0,old,1,x(1,1,1,ismpl),1)

          ! Indices of maximum absolute difference
          i = idamax(ijk,x(1,1,1,ismpl),1)
          CALL ijk_m(i,ipt,jpt,kpt)

          ! Maximum absolute difference
          difmax = x(ipt,jpt,kpt,ismpl)

          ! Sum of squares
          CALL s_ddot(ijk,x(1,1,1,ismpl),x(1,1,1,ismpl),rms)

        ELSE
          ! Relative difference (currently not used)

          ! Number of cells
          ijk = ni*nj*nk

          ! Copy new values to x
          CALL dcopy(ijk,new,1,x(1,1,1,ismpl),1)

          ! Relative difference array (new - old)/old
          CALL daxpy(ijk,-1.0d0,old,1,x(1,1,1,ismpl),1)
          DO k = 1, nk
            DO j = 1, nj
              DO i = 1, ni
                IF (dabs(old(i,j,k))>1.D-200) x(i,j,k,ismpl) = x(i,j,k,ismpl)/old(i,j,k)
              END DO
            END DO
          END DO

          ! Indices of maximum relative difference
          i = idamax(ijk,x(1,1,1,ismpl),1)
          CALL ijk_m(i,ipt,jpt,kpt)

          ! Maximum relative difference
          difmax = x(ipt,jpt,kpt,ismpl)

          ! Sum of squares
          CALL s_ddot(ijk,x(1,1,1,ismpl),x(1,1,1,ismpl),rms)

        END IF

        ! Root-Mean part of the RMSE
        rms = sqrt(rms/dble(ijk))

        ! Standard output
        IF (linfos(3)>=2) WRITE(*,'(A,1e12.5,A,1e12.5,3(A,i5),A)') &
          '   nl iteration rms =', rms, ',  difmax =', difmax, ', [', &
          ipt, ',', jpt, ',', kpt, ']'

        RETURN
      END
