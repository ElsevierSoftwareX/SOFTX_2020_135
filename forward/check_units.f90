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

!>    @brief prove the thickness for each unit layer
!>    @param[in] ismpl local sample index
!>    @details
!> count for each rock layer the number of cell neighbours\n
!> and warns about small layers (currently disabled !!!)\n
      SUBROUTINE check_units(ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
!     counts the directions (thickness large enough)
        INTEGER icount
!     current unit, min. thickness criteria, ()
        INTEGER un, xdcrit, lines
!     delta offset for neighbours
        INTEGER d_i, d_j, d_k
!     interval 'i':[iv,ib]
        INTEGER iv_, ib_, id_
!     interval 'j':[jv_,jb_]
        INTEGER jv_, jb_, jd_
!     interval 'k':[kv_,kb_]
        INTEGER kv_, kb_, kd_
!     counts all to small areas for each unit
        INTEGER, ALLOCATABLE :: warn(:,:)
!     : in percentage
        DOUBLE PRECISION warn_p
        integer :: i
        integer :: j
        integer :: k
        integer :: ismpl

! --- routine disabled !!! ---
        RETURN
! --- routine disabled !!! ---

!     init.
        lines = 0
        ALLOCATE(warn(maxunits,2))
        DO i = 1, maxunits
          warn(i,1) = 0
          warn(i,2) = 0
        END DO

!     setup all ranges
        iv_ = 3
        ib_ = i0 - 2
        id_ = 1
        IF (i0<5) THEN
!       disable this dimension
          iv_ = 1
          ib_ = i0
          id_ = 0
        END IF
        jv_ = 3
        jb_ = j0 - 2
        jd_ = 1
        IF (j0<5) THEN
!       disable this dimension
          jv_ = 1
          jb_ = j0
          jd_ = 0
        END IF
        kv_ = 3
        kb_ = k0 - 2
        kd_ = 1
        IF (k0<5) THEN
!       disable this dimension
          kv_ = 1
          kb_ = k0
          kd_ = 0
        END IF

!     setup right criteria
        IF (id_+jd_+kd_==1) THEN
!       1D criteria
          xdcrit = 1
        ELSE IF (id_+jd_+kd_==2) THEN
!       2D criteria
          xdcrit = 3
        ELSE
!       3D criteria
          xdcrit = 7
        END IF

!     analyse all inner elements
        DO k = kv_, kb_
          DO j = jv_, jb_
            DO i = iv_, ib_
!       current unit number
              un = uindex(i,j,k)
              icount = 0
!       prove for all directions
              DO d_k = -kd_, + kd_
                DO d_j = -jd_, + jd_
                  DO d_i = -id_, + id_
                    IF ((d_i/=0) .OR. (d_j/=0) .OR. (d_k/=0)) THEN
                      IF (((uindex(i+d_i,j+d_j,k+d_k)== &
                        un) .AND. (uindex(i-d_i,j-d_j,k-d_k)== &
                        un)) .OR. ((uindex(i+d_i,j+d_j,k+d_k)== &
                        un) .AND. (uindex(i+2*d_i,j+2*d_j,k+2*d_k)== &
                        un)) .OR. ((uindex(i-d_i,j-d_j,k-d_k)== &
                        un) .AND. (uindex(i-2*d_i,j-2*d_j,k-2*d_k)== &
                        un))) icount = icount + 1
                    END IF
                  END DO
                END DO
              END DO

              IF ((icount/2)<xdcrit) THEN
!         thickness to small !!!
                warn(un,1) = warn(un,1) + 1
              ELSE
!         enough thickness !!!
                warn(un,2) = warn(un,2) + 1
              END IF
            END DO
          END DO
        END DO

!     make an output for all units layers which are to small
        DO i = 1, maxunits
          IF (warn(i,1)+warn(i,2)>27) THEN
            warn_p = dble(warn(i,1))/dble(warn(i,1)+warn(i,2))
            IF (warn_p>0.01D0) THEN
              IF (lines==0) THEN
                WRITE(*,*)
                lines = 1
              END IF
              WRITE(*,'(A,I7,A,F5.1,A)') '  *** warning: unit layer ' &
                , i, ' is too fine (about ', warn_p*100.0D0, '%) ***'
            END IF
          END IF
        END DO
        IF (lines==1) WRITE(*,*)

!     free memory
        DEALLOCATE(warn)

        RETURN
      END
