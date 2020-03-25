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

!>    @brief write parameter depending on "simul"
!>    @param[in] index_i index number
!>    @param[in] ismpl local sample index
!>    @details
!> Multiple runs at the same time (different parameters for different samples) -> OpenMP "ordered" needed -> different sorted lines!\n
      SUBROUTINE write_hdfparameter2(index_i,ismpl)
        use arrays
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_simul
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i

        INTEGER i1, i2, index_i, lblank, anzahl, anzi
        EXTERNAL lblank
        character (len=256) :: filename
        DOUBLE PRECISION get_optip
        EXTERNAL get_optip

        character (len=16) :: pname
        character (len=16), dimension (:), allocatable :: ctmp
        DOUBLE PRECISION, ALLOCATABLE :: ptmp(:)

#ifndef noHDF

        IF ((mpara<=0) .OR. ( .NOT. write_param)) RETURN
#ifdef NOOUT
        RETURN
#endif

        IF (index_i<=0) THEN
          WRITE(*,'(A)') &
            'error: wrong index number in "write_hdfparameter.f" !!!'
          STOP
        END IF

! !!! to avoid compiler bugs
!--- C$OMP ordered
#ifdef fOMP
        CALL omp_ordered_begin(index_i)
#endif

        anzahl = mpara
        ALLOCATE(ptmp(anzahl))
        ALLOCATE(ctmp(anzahl))

        CALL chln(project,i1,i2)
        filename = project(i1:i2) // '_parameters.h5'

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : Parameter information to "', &
            filename(1:lblank(filename)), '"'
        END IF

        anzi = 0

!     collect values
        DO i = 1, mpara
          ptmp(i) = get_optip(i,ismpl)
        END DO

!     collect names
        DO i = 1, mpara
          CALL param_name(i,pname,ismpl)
          anzi = anzi + 1
          WRITE(ctmp(anzi),'(1A16)') pname
        END DO

        IF ((anzi)/=anzahl) THEN
          WRITE(*,*) &
            'error, "anzi"<>"anzahl" in "write_hdfparameter.f" !!!'
          STOP
        END IF

#ifdef fOMP
!$OMP critical
#endif

        CALL h5open_f(error)
        IF (index_i<=1) THEN
          CALL create_hdf5(filename)
          CALL write2_hdf5_char(16,anzahl,ctmp,'title',filename)
        END IF
        CALL add_line(anzahl,index_i,ptmp,'parameter',filename)
        CALL h5close_f(error)


#ifdef fOMP
!$OMP end critical
#endif

        DEALLOCATE(ctmp)
        DEALLOCATE(ptmp)


! !!! to avoid compiler bugs
!--- C$OMP end ordered
#ifdef fOMP
        CALL omp_ordered_end(index_i)
#endif

#endif
        RETURN
      END
