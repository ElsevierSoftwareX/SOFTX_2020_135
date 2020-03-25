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

!> @brief Generate state vector for Pilot Point EnKF
!> @details
!> Generate smaller pilot point version of `mem`. For the other
!> variables in the state vector, values are copied. For the
!> Pilot-Point variable, only those values at Pilot-Point locations
!> are copied.
subroutine enkf_state_vector_pp()

  use mod_enkf, only:&
       mem,&
       num_pp,&
       mem_pp,&
       mem_rp,&
       pp_state_inds,&
       nrens,&
       nstate,&
       nstate_pp_temp,&
       act_s,&
       lstate,&
       pp_ivar,&
       ipp_start,&
       ipp_end

  implicit none

  integer :: ipp, irens, imem, jpp, irp

  allocate(mem_pp(nstate,nrens))
  allocate(mem_rp(lstate-num_pp,nrens))

  if(.not. act_s(pp_ivar) == 1) then
     write (unit = *, fmt = *) "[E] Error in enkf_state_vector_pp"
     stop 1
  end if

  do irens = 1, nrens
     irp = 1
     jpp = 1
     ipp = 1
     do imem = 1, nstate_pp_temp
        if(imem<ipp_start) then
           mem_pp(ipp,irens) = mem(imem,irens)
           ipp = ipp + 1
        else if(imem>ipp_end) then
           mem_pp(ipp,irens) = mem(imem,irens)
           ipp = ipp + 1
        else
           if ((jpp <= num_pp) .and. (jpp >= 1)) then
              if (imem == pp_state_inds(jpp)) then
                 mem_pp(ipp,irens) = mem(imem,irens)
                 ipp = ipp + 1
                 jpp = jpp + 1
              else
                 mem_rp(irp,irens) = mem(imem,irens)
                 irp = irp + 1
              end if
           else
              mem_rp(irp,irens) = mem(imem,irens)
              irp = irp + 1
           end if
        end if
     end do
  end do

  deallocate(mem)
  allocate(mem(nstate,nrens))

  do irens = 1, nrens
     do ipp = 1, nstate
        mem(ipp,irens) = mem_pp(ipp,irens)
     end do
  end do

end subroutine enkf_state_vector_pp
