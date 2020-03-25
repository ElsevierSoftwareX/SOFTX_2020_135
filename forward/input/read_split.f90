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

!>    @brief read unit spliting
!>    @param[in] filename file name
!>    @param[in] ismpl local sample index
!>    @details
!> read configuration for splitting some units,\n
!> split each former unit (layer) into different units for each cell\n
      SUBROUTINE read_split(filename,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        character (len=80) :: filename
        character (len=80) :: line
        INTEGER nsplit, msplit
        ! INTEGER munits
        INTEGER nn

        INTEGER, ALLOCATABLE :: isplit(:), sindex(:)
!     copy fields
        ! INTEGER, ALLOCATABLE :: itmp(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: dtmp(:,:)

        INTEGER lblank
        LOGICAL found
        EXTERNAL lblank, found


!     open project Config file
        OPEN(79,file=filename,status='old')

        IF (found(79,key_char//' split units',line,.FALSE.)) THEN
          WRITE(*,*) ' '
          WRITE(*,*) '  reading splitting parameter:'
          WRITE(*,*) '    from file "', filename(:lblank(filename)), &
            '"'
          WRITE(*,*) ' '

          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) nsplit
          ELSE
            READ(line(i:j),*) nsplit
          END IF
          ALLOCATE(isplit(nsplit))
          READ(79,*) (isplit(i),i=1,nsplit)
          WRITE(*,'(1A,1I4,1A)') '  [I] : splitting for ', nsplit, &
            ' unit(s)'

          ALLOCATE(sindex(i0*j0*k0+nunits))
          DO l = 1, nunits
            sindex(l) = l
          END DO
!        count the number of unit-elements
          msplit = nunits
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                DO l = 1, nsplit
                  IF (uindex(i,j,k)==isplit(l)) THEN
                    msplit = msplit + 1
                    sindex(msplit) = uindex(i,j,k)
                    uindex(i,j,k) = msplit
!                 break the search loop
                    GO TO 100
                  END IF
                END DO
100             CONTINUE
              END DO
            END DO
          END DO
          WRITE(*,'(1A,1I9)') '  [I] : (splitted) units = ', msplit

!        # propunit
          nn = nprop*nsmpl
          ALLOCATE(dtmp(nunits,nn))
!        save values
          CALL dcopy(nunits*nn,propunit,1,dtmp,1)
!        reallocate space
          DEALLOCATE(propunit)
          ALLOCATE(propunit(msplit,nprop,nsmpl))
          memory = memory - nunits*nn
          memory = memory + msplit*nn
!        split elements, execption for "propunit", because of the bc-parts
!        over all elements, normal units
          DO k = 1, nsmpl
            DO j = firstidx, lastidx
              l = j + nprop*(k-1)
              DO i = 1, msplit
!              copy old values for new unit
                propunit(i,j,k) = dtmp(sindex(i),l)
              END DO
            END DO
          END DO
!        over all elements, bc units
          DO k = 1, nsmpl
            DO j = bc_firstidx, bc_lastidx
              l = j + nprop*(k-1)
              DO i = 1, bc_maxunits
!              copy old values for new unit
                propunit(i,j,k) = dtmp(i,l)
              END DO
            END DO
          END DO
          DEALLOCATE(dtmp)

          DEALLOCATE(sindex)
          DEALLOCATE(isplit)
          maxunits = msplit
          nunits = msplit
        END IF

!     close project config file
        CLOSE(79)

        RETURN
      END
