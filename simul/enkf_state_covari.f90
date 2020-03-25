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

!> @brief Compute and write state covariance matrix
!> @param[in] a_before_after switch before or after assimilation ('bef'/'aft')
!> @details
!> The state and sample covariance matrix is computed from mem
!> and the an output is written.
      SUBROUTINE enkf_state_covari(a_before_after)

        use arrays, only: &
             project_sfx

        use mod_simul, only: &
             senkf_outdir

        use mod_enkf, only: &
             mem,&
             ave,&
             nstate,&
             nrens

        implicit none

        character (len=3), intent(in) :: a_before_after

        double precision, ALLOCATABLE :: memdot(:,:)
        double precision, ALLOCATABLE :: p(:,:)

        double precision :: fac
        integer :: i, j, iens
        character (len=5) :: teci
        INTRINSIC trim

        ALLOCATE(memdot(nstate,nrens))
        ALLOCATE(p(nstate,nstate))

        fac=1.d0/dble(float(nrens-1))

        DO iens=1,nrens
          memdot(:,iens)= mem(:,iens) - ave(:)
        END DO
        call dgemm('n','t',nstate,nstate,nrens,fac,memdot,nstate,memdot,nstate,  &
                   0.d0,p,nstate)

        if( .not.(a_before_after == 'bef' .or. a_before_after == 'aft') ) then
           write(unit = *, fmt = *) '[E1] Error in enkf_state_covari()'
        end if

        OPEN(unit=32,file=senkf_outdir//'state_covariance'//'_'//&
             trim(a_before_after)//trim(project_sfx(1))//'.plt')
        WRITE(32,*)'TITLE = '//'"state covariance"'
        WRITE(32,*)'VARIABLES = '
        WRITE(32,*)'"I"'//','//'"J"'//','//'"P"'
        WRITE(teci,'(i5)')nstate
        WRITE(32,*)'ZONE  I= '//teci//',J= '//teci//', f=POINT'
        DO i=1,nstate
          DO j=1,nstate
            write(32,'(i5,2x,i5,2x,e12.4)')i,j,p(i,j)
          END DO
        END DO
        CLOSE(32)

        DEALLOCATE(memdot)
        DEALLOCATE(p)

      END SUBROUTINE enkf_state_covari
