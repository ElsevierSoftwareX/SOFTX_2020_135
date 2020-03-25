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

!> @brief Save overall standard deviation for the state vector variables
!> @param[in] a_before_after bef/aft assimilation
!> @details
!> __USE__: The vector containing the variances for each single state
!> vector variables, the usual numbers specifying the number of
!> observations and length of the state vector, as well as number
!> of active species in the state vector.
!>
!> __SET__: The OVERALL ENSEMBLE Standard deviation for the whole grid,
!> but for each variable seperately in stdvar. The square-root of
!> the arithmetic mean of the variances across the grid is
!> calculated.
subroutine enkf_stddev(a_before_after)

  use mod_genrl, only: &
       i0,&
       j0,&
       k0

  use mod_enkf, only: &
       act_s,&
       lstate,&
       var,&
       irobs,&
       num_enkf_vars,&
       stdvar

  implicit none

  character (len=3), intent(in) :: a_before_after

  integer :: i, j, k, imem
  integer :: istd, ivar

  integer, external :: index_loc_to_mem

  if(a_before_after == 'bef') then
     istd = 2*irobs -1
  else if(a_before_after == 'aft') then
     istd = 2*irobs
  else
     write(unit = *, fmt = *) '[E1] Error in enkf_stddev()'
     stop 1
  end if

  ! Square root of mean variance
  do ivar = 1, num_enkf_vars 
     if(act_s(ivar)==1) then

        do k = 1, k0
           do j = 1, j0
              do i = 1, i0
                 imem = index_loc_to_mem(i,j,k,ivar)
                 if (imem > 0) then

                    stdvar(istd,ivar) = stdvar(istd,ivar) &
                         + var(imem)

                 end if
              end do
           end do
        end do

        stdvar(istd,ivar) = stdvar(istd,ivar)/dble(float(lstate))
        stdvar(istd,ivar) = dsqrt(stdvar(istd,ivar))

     end if
  end do

end subroutine enkf_stddev
