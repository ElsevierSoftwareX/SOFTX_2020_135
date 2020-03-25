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

!>    @brief determine output times and call the write routine
!>    @param[in] iseed seeding component/index
!>    @param[in] ismpl local sample index
      SUBROUTINE write_joutt(iseed,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_time
        IMPLICIT NONE
        integer :: ismpl
        integer :: i

        DOUBLE PRECISION deltt, talfa, malfa, numdiff
        INTEGER id, ii, i1, i2, i3, i4, iseed
        CHARACTER snumber*30, snumber2*30
        DOUBLE PRECISION deltat
        EXTERNAL deltat

!     allowed numerical difference
        numdiff = 1.D-14*tunit
        deltt = deltat(simtime(ismpl),ismpl)
!
        id = noutt + 1
        DO i = 1, noutt + 1
          IF (outt(i)>simtime(ismpl)+numdiff) THEN
            id = i 
            GO TO 100
          END IF
        END DO
!
100     IF (id == noutt + 1) THEN
          IF ( .NOT. transient) THEN
            project_sfx(ismpl) = '_out'
            CALL write_joutt_hdf(iseed,ismpl)
            project_sfx(ismpl) = ' '
          END IF
          RETURN
        ELSE IF (outt(id)>simtime(ismpl)+numdiff .AND. &
            outt(id)<=simtime(ismpl)+deltt+numdiff) THEN
!
          WRITE(snumber,'(1e16.8)') (simtime(ismpl)+deltt)/tunit
          CALL chln(snumber,i1,i2)
          WRITE(snumber2,'(1I7)') id
          CALL chln(snumber2,i3,i4)
          project_sfx(ismpl) = '_time_' // snumber(i1:i2) // &
            '_out_' // snumber2(i3:i4)
          CALL write_joutt_hdf(iseed,ismpl)
          project_sfx(ismpl) = ' '
        END IF

        RETURN
      END
