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

!> @brief forward model initialisation
!> @param[in] ismpl local sample index
!> @details
!> - For the head-based model computation, the pressure may need to be
!>   computed from head input. For the pres-based model computation,
!>   the head may need to be computed from pressure input. \n\n
!>
!> - Call user and props initializing subroutines \n\n
!>
!> - doing additional stuff for the multi-sample computation:\n
!>   - prepare monitor output files (headers)\n
!>   - prepare main input/output vector variables\n\n
!>
!> - output memory information\n
      subroutine forward_init(ismpl)

        use arrays
        use mod_genrl
        use mod_linfos
        use mod_data

        implicit none

        ! Local sample index
        integer :: ismpl

        ! Local variable for array memory
        double precision :: memloc

        integer :: i

        double precision, external :: get_optip


        if (linfos(3)>=2) write(*,*) ' ... forward_init'

        ! If needed: Transformation to non-computed flow variable.
#ifdef head_base
        if (is_init_flow_trafo_needed) call head2pres(0,ismpl)
#endif
#ifdef pres_base
        if (is_init_flow_trafo_needed) call pres2head(0,ismpl)
#endif

        ! USER model init
        call user_init(ismpl)

        ! PROPS model init
        call props_init(ismpl)
!
!       setup main input/output variables for "forward_compute"
        ALLOCATE(main_input(max(mpara,1),nsmpl))
        memory = memory +max(mpara,1)*nsmpl
        ALLOCATE(main_output(max(ndata,1),nsmpl))
        memory = memory +max(ndata,1)*nsmpl
!       init master input - parameter vector
        DO i = 1, mpara
          main_input(i,idx_master) = get_optip(i,idx_master)
        END DO
!
        ! init master copy, save state from file init (necessary for transient models, INVERSE and SIMUL/ENKF)
        CALL old_save(cgen_opti,ismpl)

        ! Memory information output
        if (linfos(1)>=0) then
          ! Assume double precision data type
          memloc = memory*8.0D0
          ! ... in KByte
          memloc = memloc/1024.0D0
          ! ... in MByte
          memloc = memloc/1024.0D0
          write(*,*) ''
          write(*,'(a,f11.3,a)') '  [I] : memory: ', memloc, ' MByte data'
          write(*,*) ''
        end if

        return

      end subroutine forward_init
