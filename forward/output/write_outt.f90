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

!>    @brief determine output times, interpolate, and call the write routine
!>    @param[in] deltt time step length
!>    @param[in] ismpl local sample index
      SUBROUTINE write_outt(deltt,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_time
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION deltt, talfa, malfa, numdiff
        DOUBLE PRECISION, ALLOCATABLE :: tmp_new(:,:)
        INTEGER ijk0, id, ii, i1, i2, i3, i4
        character (len=80) :: sfx_old
        INTRINSIC trim


        IF (write_disable) RETURN

!     allowed numerical difference
        numdiff = 1.D-14*tunit

        id = 0
        DO i = 1, noutt + 1
          IF (outt(i)>simtime(ismpl)+numdiff) THEN
            id = i - 1
            GO TO 100
          END IF
        END DO

100     IF (id==0) THEN
!        no output
          RETURN
        ELSE IF (outt(id)>simtime(ismpl)-deltt+numdiff .AND. &
            outt(id)<=simtime(ismpl)+numdiff) THEN
!        time interval match
          ijk0 = i0*j0*k0
!        save current main values
          ALLOCATE(tmp_new(ijk0,5))
          CALL dcopy(ijk0,head(1,1,1,ismpl),1,tmp_new(1,1),1)
          CALL dcopy(ijk0,temp(1,1,1,ismpl),1,tmp_new(1,2),1)
          CALL dcopy(ijk0,pres(1,1,1,ismpl),1,tmp_new(1,3),1)
          CALL dcopy(ijk0,conc(1,1,1,1,ismpl),1,tmp_new(1,5),1)

          CALL old_restore(cgen_time,ismpl)
          ii = 0
!        interpolate
          talfa = (simtime(ismpl)-outt(id))/deltt
          malfa = 1.0D0 - talfa
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                ii = ii + 1
                head(i,j,k,ismpl) = talfa*head(i,j,k,ismpl) + malfa*tmp_new(ii,1)
                temp(i,j,k,ismpl) = talfa*temp(i,j,k,ismpl) + malfa*tmp_new(ii,2)
                pres(i,j,k,ismpl) = talfa*pres(i,j,k,ismpl) + malfa*tmp_new(ii,3)
                conc(i,j,k,1,ismpl) = talfa*conc(i,j,k,1,ismpl) + malfa*tmp_new(ii,5)
              END DO
            END DO
          END DO

          sfx_old = project_sfx(ismpl)
          project_sfx(ismpl) = trim(sfx_old) // '_time' // &
              '_out'
          CALL forward_write(id,ismpl)
          project_sfx(ismpl) = sfx_old


!        restore current main values
          CALL dcopy(ijk0,tmp_new(1,1),1,head(1,1,1,ismpl),1)
          CALL dcopy(ijk0,tmp_new(1,2),1,temp(1,1,1,ismpl),1)
          CALL dcopy(ijk0,tmp_new(1,3),1,pres(1,1,1,ismpl),1)
          CALL dcopy(ijk0,tmp_new(1,5),1,conc(1,1,1,1,ismpl),1)
          DEALLOCATE(tmp_new)
        END IF

        RETURN
      END
