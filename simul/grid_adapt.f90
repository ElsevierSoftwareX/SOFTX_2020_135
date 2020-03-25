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

!> @brief check grid adapter usage, sanity check
!> @param[in] fname SIM file name
!> @details
!> compares the physical domain of Shemat-Suite definition
!> (`<projectname>`) with the SIMUL parameter file (`<fname>`)
      SUBROUTINE grid_adapt_check(fname)
        use arrays
        USE simul_arrays
        use mod_genrl
        use mod_simul
        use mod_linfos
        IMPLICIT NONE
        integer :: i
        character (len=*) :: fname
        DOUBLE PRECISION dtmp, dtmin
        INTRINSIC dble, min, dabs, trim

        IF (linfos(2)>=1) WRITE(*,*) ' [I] : check grid adapting'
!     sanity checks
        dtmp = delxa(i0) + 0.5D0*delx(i0)
        dtmin = delx(1)
        DO i = 2, i0
          dtmin = min(delx(i),dtmin)
        END DO
        IF (dabs(sm_delx*dble(sm_i0)-dtmp)>dtmin*0.01D0) THEN
          WRITE(*,'(3A)') &
            'error: grid mismatch beween "delx" and SIM "', &
            trim(fname), '" configuration!'
          WRITE(*,*) 'SM-delx:', sm_delx, ', SM-#:', sm_i0, &
            ', set-sum:', dtmp
          STOP
        END IF
        IF (sm_i0>i0) THEN
          WRITE(*,'(3A)') 'error: SIM "', trim(fname), &
            '" I0-dimension must be lesser or equal than model-I0!'
          STOP
        END IF
!
        dtmp = delya(j0) + 0.5D0*dely(j0)
        dtmin = dely(1)
        DO i = 2, j0
          dtmin = min(dely(i),dtmin)
        END DO
        IF (dabs(sm_dely*dble(sm_j0)-dtmp)>dtmin*0.01D0) THEN
          WRITE(*,'(3A)') &
            'error: grid mismatch beween "dely" and SIM "', &
            trim(fname), '" configuration!'
          WRITE(*,*) 'SM-dely:', sm_dely, ', SM-#:', sm_j0, &
            ', set-sum:', dtmp
          STOP
        END IF
        IF (sm_j0>j0) THEN
          WRITE(*,'(3A)') 'error: SIM "', trim(fname), &
            '" J0-dimension must be lesser or equal than model-J0!'
          STOP
        END IF
!
        dtmp = delza(k0) + 0.5D0*delz(k0)
        dtmin = delz(1)
        DO i = 2, k0
          dtmin = min(delz(i),dtmin)
        END DO
        IF (dabs(sm_delz*dble(sm_k0)-dtmp)>dtmin*0.01D0) THEN
          WRITE(*,'(3A)') &
            'error: grid mismatch beween "delz" and SIM "', &
            trim(fname), '" configuration!'
          STOP
        END IF
        IF (sm_k0>k0) THEN
          WRITE(*,'(3A)') 'error: SIM "', trim(fname), &
            '" K0-dimension must be lesser or equal than model-K0!'
          WRITE(*,*) 'SM-delz:', sm_delz, ', SM-#:', sm_k0, &
            ', set-sum:', dtmp
          STOP
        END IF
!
        RETURN
      END

!> @brief converts the [I0,J0,K0]-system coords into a *SIM-system index
!> @param[in] i I0-dimension index
!> @param[in] j J0-dimension index
!> @param[in] k K0-dimension index
!> @param[out] sm_ijk *SIM-index position for [i,j,k]
      SUBROUTINE grid_adapt_ijk(i,j,k,sm_ijk)
        use arrays
        USE simul_arrays
        use mod_genrl
        use mod_simul
        use mod_linfos
        IMPLICIT NONE
        integer :: i, j, k
        integer :: sm_ijk
        INTRINSIC int

!     the absolute position (center of a cell) divided by the
!     equidistant *SIM cell-length (plus one) is the target index
        sm_ijk = int(delxa(i)/sm_delx) + 1 + int(delya(j)/sm_dely)* &
          sm_i0 + int(delza(k)/sm_delz)*sm_i0*sm_j0
!
        IF (linfos(2)>=3) WRITE(*,'(4(1A,1I6,1A,1I6))') &
          ' [debug] : ', int(delxa(i)/sm_delx) + 1, '==', i, ', ', &
          int(delya(j)/sm_dely) + 1, '==', j, ', ', &
          int(delza(k)/sm_delz) + 1, '==', k, ', ', sm_ijk, '/', &
          sm_i0*sm_j0*sm_k0
!
        RETURN
      END
