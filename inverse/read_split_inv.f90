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

!>    @brief read unit spliting (support of inverse mode)
!>    @param[in] filename file name
!>    @param[in] ismpl local sample index
!>    @details
!> active units (layer) not supported,\n
!> only units with not active parameters can be splitted\n
      SUBROUTINE read_split_inv(filename,ismpl)
        use arrays
#ifdef AD
        use g_arrays
#endif 
#ifdef AD_RM 
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_inverse
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l


        CHARACTER filename*80, line*80
        INTEGER nsplit, msplit, munits, nn, npara
        LOGICAL ltake

        INTEGER, ALLOCATABLE :: isplit(:), sindex(:)
!     copy fields
        INTEGER, ALLOCATABLE :: itmp(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: dtmp(:,:)

        INTEGER lblank
        LOGICAL found
        EXTERNAL lblank, found


!     open projeCt Config file
        OPEN(79,file=filename,status='old')

        IF (found(79, key_char//' split units',line,.FALSE.)) THEN
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
          memory = memory - nunits*nn + msplit*nn
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

!        # a_propunit
          nn = nprop
          ALLOCATE(dtmp(nunits,nn))
!        save values
          CALL dcopy(nunits*nn,a_propunit,1,dtmp,1)
!        reallocate space
          DEALLOCATE(a_propunit)
          ALLOCATE(a_propunit(msplit,nprop))
          memory = memory - nunits*nn + msplit*nn
!        split elements, execption for "a_propunit", because of the bc-parts
!        over all elements, normal units
          DO j = firstidx, lastidx
            DO i = 1, msplit
!              copy old values for new unit
              a_propunit(i,j) = dtmp(sindex(i),j)
            END DO
          END DO
!        over all elements, bc units
          DO j = bc_firstidx, bc_lastidx
            DO i = 1, bc_maxunits
!              copy old values for new unit
              a_propunit(i,j) = dtmp(i,j)
            END DO
          END DO
          DEALLOCATE(dtmp)

!        # e_propunit
          nn = nprop
          ALLOCATE(dtmp(nunits,nn))
!        save values
          CALL dcopy(nunits*nn,e_propunit,1,dtmp,1)
!        reallocate space
          DEALLOCATE(e_propunit)
          ALLOCATE(e_propunit(msplit,nprop))
          memory = memory - nunits*nn + msplit*nn
!        split elements, execption for "e_propunit", because of the bc-parts
!        over all elements, normal units
          DO j = firstidx, lastidx
            DO i = 1, msplit
!              copy old values for new unit
              e_propunit(i,j) = dtmp(sindex(i),j)
            END DO
          END DO
!        over all elements, bc units
          DO j = bc_firstidx, bc_lastidx
            DO i = 1, bc_maxunits
!              copy old values for new unit
              e_propunit(i,j) = dtmp(i,j)
            END DO
          END DO
          DEALLOCATE(dtmp)

!        # d_propunit
          nn = nprop
          ALLOCATE(dtmp(nunits,nn))
!        save values
          CALL dcopy(nunits*nn,d_propunit,1,dtmp,1)
!        reallocate space
          DEALLOCATE(d_propunit)
          ALLOCATE(d_propunit(msplit,nprop))
          memory = memory - nunits*nn + msplit*nn
!        split elements, execption for "d_propunit", because of the bc-parts
!        over all elements, normal units
          DO j = firstidx, lastidx
            DO i = 1, msplit
!              copy old values for new unit
              d_propunit(i,j) = dtmp(sindex(i),j)
            END DO
          END DO
!        over all elements, bc units
          DO j = bc_firstidx, bc_lastidx
            DO i = 1, bc_maxunits
!              copy old values for new unit
              d_propunit(i,j) = dtmp(i,j)
            END DO
          END DO
          DEALLOCATE(dtmp)

#ifndef AD_RM 
!        # g_propunit
          nn = nprop*nsmpl
          ALLOCATE(dtmp(nunits,nn))
!        save values
          CALL dcopy(nunits*nn,g_propunit,1,dtmp,1)
!        reallocate space
          DEALLOCATE(g_propunit)
          ALLOCATE(g_propunit(msplit,nprop,nsmpl))
          memory = memory - nunits*nn + msplit*nn
!        split elements, execption for "g_propunit", because of the bc-parts
!        over all elements, normal units
          DO k = 1, nsmpl
            DO j = firstidx, lastidx
              l = j + nprop*(k-1)
              DO i = 1, msplit
!              copy old values for new unit
                g_propunit(i,j,k) = dtmp(sindex(i),l)
              END DO
            END DO
          END DO
!        over all elements, bc units
          DO k = 1, nsmpl
            DO j = bc_firstidx, bc_lastidx
              l = j + nprop*(k-1)
              DO i = 1, bc_maxunits
!              copy old values for new unit
                g_propunit(i,j,k) = dtmp(i,l)
              END DO
            END DO
          END DO
          DEALLOCATE(dtmp)
#else
!        # propunit_ad
          nn = nprop*nsmpl
          ALLOCATE(dtmp(nunits,nn))
!        save values
          CALL dcopy(nunits*nn,propunit_ad,1,dtmp,1)
!        reallocate space
          DEALLOCATE(propunit_ad)
          ALLOCATE(propunit_ad(msplit,nprop,nsmpl))
          memory = memory - nunits*nn + msplit*nn
!        split elements, execption for "g_propunit", because of the bc-parts
!        over all elements, normal units
          DO k = 1, nsmpl
            DO j = firstidx, lastidx
              l = j + nprop*(k-1)
              DO i = 1, msplit
!              copy old values for new unit
                propunit_ad(i,j,k) = dtmp(sindex(i),l)
              END DO
            END DO
          END DO
!        over all elements, bc units
          DO k = 1, nsmpl
            DO j = bc_firstidx, bc_lastidx
              l = j + nprop*(k-1)
              DO i = 1, bc_maxunits
!              copy old values for new unit
                propunit_ad(i,j,k) = dtmp(i,l)
              END DO
            END DO
          END DO
          DEALLOCATE(dtmp)
#endif
!        # propunitold
!        no save/restore step, because of later initialising
          DEALLOCATE(propunitold)
          ALLOCATE(propunitold(msplit,nprop))
          memory = memory - nunits*nprop + msplit*nprop
!        # opti_bc
          ALLOCATE(itmp(nbc,bc_maxunits))
!        save values
          DO j = 1, bc_maxunits
            DO i = 1, nbc
              itmp(i,j) = opti_bc(i,j)
            END DO
          END DO
!        reallocate space
          DEALLOCATE(opti_bc)
          ALLOCATE(opti_bc(nbc,msplit))
          memory = memory - nbc*nunits + nbc*msplit
!        over all elements, normal units
          DO j = 1, bc_maxunits
            DO i = 1, nbc
              opti_bc(i,j) = itmp(i,j)
            END DO
          END DO
          DEALLOCATE(itmp)

!        # opti_props
          ALLOCATE(itmp(nprop,nunits))
!        save values
          DO j = 1, maxunits
            DO i = firstidx, lastidx
              itmp(i,j) = opti_props(i,j)
            END DO
          END DO
!        reallocate space
          DEALLOCATE(opti_props)
          ALLOCATE(opti_props(nprop,msplit))
          memory = memory - nunits*nprop + msplit*nprop
!        split elements
!        over all elements, normal units
          DO i = 1, msplit
            DO j = firstidx, lastidx
!              copy old values for new unit
              opti_props(j,i) = itmp(j,sindex(i))
            END DO
          END DO
!        delete all old units-numbers that are now splitted
          DO i = 1, nsplit
            DO j = firstidx, lastidx
              opti_props(j,isplit(i)) = 0
            END DO
          END DO
          DEALLOCATE(itmp)

!        first count new "mpara"
          npara = 0
          DO j = 1, mpara
            DO i = 1, msplit
              IF (seed_para(2,j)==sindex(i)) THEN
                ltake = .TRUE.
                IF (i<=maxunits) THEN
!                    excepting out-dated unit
                  DO l = 1, nsplit
                    IF (i==isplit(l)) ltake = .FALSE.
                  END DO
                END IF
                IF (ltake) npara = npara + 1
              END IF
            END DO
          END DO
          IF (npara/=mpara) THEN
            WRITE(*,'(1A)') &
              'error: unit splitting not supported for inversion!'
            STOP
! ### for inversion following splittings needed ! ###
! seed_para(2,mpara)
! apri_input(mpara))
! linlog_input(mpara))
! grad(mpara))
! grad_sec(mpara))
! p_delta(mpara))
! main_input(mpara,nsmpl))
! main_input_master(mpara))
#ifndef JACOBI_FREE
! jac(ndata,mpara))
! covar_p(mpara,mpara))
! resmat_p(mpara,mpara))
! ginv(mpara,ndata))
! tmp_mat(mpara,mpara))
#endif
! tmp_vec(mpara,max(ndata,mpara),5)) <-> tmp_vec(max(ndata,mpara),1,4))
! covar_prior_p(mpara,mpara))
! tmp_p_k(mpara))
! ###################################################
          END IF

          DEALLOCATE(sindex)
          DEALLOCATE(isplit)
          maxunits = msplit
          nunits = msplit
        END IF

!     close project config file
        CLOSE(79)

        RETURN
      END
