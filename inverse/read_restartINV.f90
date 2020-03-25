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

!>    @brief read data & parameter for restart
!>    @param[in] filename file name
!>    @param[out] r_iter inverse iteration counter
!>    @param[in] ismpl local sample index
      SUBROUTINE read_restartinv(filename,r_iter,ismpl)
        use arrays
#ifdef AD
        use g_arrays
#endif
#ifdef AD_RM
        use arrays_ad
#endif
        use mod_genrl
        use mod_data
        use mod_inverse
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, l
        CHARACTER filename*80, line*80
        DOUBLE PRECISION, ALLOCATABLE :: dptmp(:,:)
        INTEGER r_iter, lout
        INTEGER lblank
        LOGICAL found, no_ext_link
        EXTERNAL found, no_ext_link, lblank


        WRITE(*,*)
        WRITE(*,*) &
          '  reading inverse model input parameter for restart:'
        WRITE(*,*) '    from file "', filename(:lblank(filename)), &
          '"'
        WRITE(*,*)

        CALL omp_new_file_handler(lout,1)
!     read restart file
        OPEN(lout,file=filename,status='old')

        IF (found(lout,key_char//' iter_inv',line,.TRUE.)) THEN
          READ(lout,*) r_iter
          WRITE(*,*) ' [R] : iteration number - for restart'
        END IF

        IF (mpara>0) THEN
!     init HDF5 support, when available
          CALL open_hdf5(' ')

#ifndef JACOBI_FREE
          IF (found(lout,key_char//' jacoby',line,.TRUE.)) THEN
            IF (no_ext_link(ndata,mpara,1,jac,'jacoby',line)) &
              READ(lout,*) ((jac(i,l),i=1,ndata),l=1,mpara)
            WRITE(*,*) ' [R] : Jacoby matrix - for restart'
          END IF
#endif

!     finish HDF5 support, when available
          CALL close_hdf5()
        END IF

        IF (found(lout,key_char//' units',line,.TRUE.)) THEN
          ALLOCATE(dptmp(nprop_load,maxunits))
          IF (no_ext_link(nprop_load,maxunits,1,dptmp,'units',line)) &
              THEN
            DO i = 1, maxunits
              READ(lout,*) (dptmp(j,i),j=1,nprop_load)
            END DO
          END IF
          IF (nprop_load==lastidx-firstidx+1) THEN
!          load all, dense entries (no specific index)
            DO i = 1, maxunits
              DO j = 1, nprop_load
                propunit(i,firstidx-1+j,ismpl) = dptmp(j,i)
              END DO
            END DO
          ELSE
!          needs to handle manual reading with specific index
            WRITE(*,'(1A)') &
              'error in "read_restartINV.f": reading units: specific index handling not implemented, use dense structure instead!'
            STOP
          END IF
          WRITE(*,*) ' [R] : units - for restart'
          DEALLOCATE(dptmp)
        END IF

        IF (found(lout,key_char//' bcunits',line,.TRUE.)) THEN
          ALLOCATE(dptmp(nbc,bc_maxunits))
          IF (no_ext_link(nbc,bc_maxunits,1,dptmp,'bcunits',line)) &
              THEN
            DO i = 1, bc_maxunits
              READ(lout,*) (dptmp(j,i),j=1,nbc)
            END DO
          END IF
          DO i = 1, bc_maxunits
            DO j = 1, nbc
              propunit(i,j,ismpl) = dptmp(j,i)
            END DO
          END DO
          WRITE(*,*) ' [R] : bcunits - for restart'
          DEALLOCATE(dptmp)
        END IF

!     close restart file
        CLOSE(lout)
        CALL omp_del_file_handler(lout)

        RETURN
      END

!>    @brief write data & parameter for restart
!>    @param[in] filename file name
!>    @param[in] r_iter inverse iteration counter
!>    @param[in] ismpl local sample index
      SUBROUTINE write_restartinv(filename,r_iter,ismpl)
        use arrays
#ifdef AD
        use g_arrays
#endif
#ifdef AD_RM
        use arrays_ad
#endif
        use mod_genrl
        use mod_data
        use mod_inverse
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        CHARACTER filename*80, hfilename*80, pform*10
        INTEGER r_iter, lout
        INTEGER lblank
        LOGICAL test_hdf5
        EXTERNAL lblank, test_hdf5
#ifdef JACOBI_FREE
        write(*,'(3A)') ' [E] : Jacobian can not be written for restart'
        RETURN
#endif
!$OMP master
        CALL omp_new_file_handler(lout,1)
!     write restart file
        OPEN(lout,file=filename,status='old',position='append')

        WRITE(lout,'(A)') key_char//' iter_inv'
        WRITE(lout,'(I8)') r_iter
!      write(*,*) ' [W] : iteration number - for restart'

        IF (mpara>0 .AND. ndata>0) THEN
          IF (test_hdf5()) THEN
!         init HDF5 support, when available
            CALL open_hdf5(' ')
            hfilename = filename(:lblank(filename)) // '.h5'

            WRITE(lout,'(A)') key_char//' jacoby, HDF5 = ' // hfilename
            CALL write2_hdf5(ndata,mpara,jac,'jacoby',hfilename)

            WRITE(*,'(3A)') '  [W] : Jacoby matrix to "', &
              filename(:lblank(filename)), '.h5" - for restart'

!         finish HDF5 support, when available
            CALL close_hdf5()
          ELSE
            WRITE(lout,'(A)') key_char//' jacoby'
            CALL write_dense3d(jac,ndata,mpara,1,lout,ismpl)

            WRITE(*,'(3A)') '  [W] : Jacoby matrix to "', &
              filename(:lblank(filename)), '" - for restart'
          END IF
        END IF

        WRITE(lout,'(A)') key_char//' units'
!     store all, dense entries (no specific index)
        DO i = 1, maxunits
          WRITE(lout,'('//c_npropunit//'e24.16)') (propunit(i,j,ismpl),j=firstidx, &
            lastidx)
        END DO

        WRITE(lout,'(A)') key_char//' bcunits'
!     store all, dense entries (no specific index)
        DO i = 1, bc_maxunits
          WRITE(lout,'('//c_nbcunit//'e24.16)') (propunit(i,j,ismpl),j=bc_firstidx, &
            bc_lastidx)
        END DO

        WRITE(*,'(3A)') '  [W] : inverse parameters to "', &
          filename(:lblank(filename)), '" - for restart'

!     close restart file
        CLOSE(lout)
        CALL omp_del_file_handler(lout)
!$OMP end master

        RETURN
      END
