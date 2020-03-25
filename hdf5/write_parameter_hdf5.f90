!>    @brief write parameter depending on "inverse"
!>    @param[in] index_i index number
!>    @param[in] ismpl local sample index
!>    @details
!> Only one single run at one time (one global parameter result) -> no OpenMP "ordered" needed!\n
      SUBROUTINE write_hdfparameter(index_i,ismpl)
        use arrays
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_inverse
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i

        INTEGER i1, i2, index_i, lblank, anzahl, anzi
        EXTERNAL lblank
        CHARACTER filename*80
        DOUBLE PRECISION get_optip
        EXTERNAL get_optip

        CHARACTER*16 pname
        CHARACTER*16, ALLOCATABLE :: ctmp(:)
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

        anzahl = mpara + 3
        ALLOCATE(ptmp(anzahl))
        ALLOCATE(ctmp(anzahl))

        CALL chln(project,i1,i2)
        filename = project(i1:i2) // '_parameter.h5'

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

!     rms_data
        anzi = anzi + 1
        ptmp(anzi) = 0.0D0
        IF (ndata>0) ptmp(anzi) = dsqrt(rms_data)/dble(ndata)
        WRITE(ctmp(anzi),'(A16)') '        rms_data'
!     rms_para
        anzi = anzi + 1
        ptmp(anzi) = dsqrt(rms_para)/dble(mpara)
        WRITE(ctmp(anzi),'(A16)') '        rms_para'
!     quality function
        anzi = anzi + 1
        ptmp(anzi) = 0.5D0*(rms_data+rms_para)
        WRITE(ctmp(anzi),'(A16)') '    quality_func'

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

#endif
        RETURN
      END
