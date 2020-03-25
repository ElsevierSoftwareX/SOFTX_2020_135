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

!> @brief Output of variables/parameters at single cells
!> @param[in] a_before_after switch before or after assimilation ('bef'/'aft')
!> @details
!> Output to directory
!>
!>     model_output/single_cell_output
!>
!> to the file
!>
!>     single_cell_E1_i_j_k_var_bef.plt
subroutine enkf_output_single_cell(a_before_after)

  use arrays, only: &
       project_sfx,&
       head,&
       temp,&
       conc

  use mod_simul, only: &
       ssingle_cell_outdir

  use mod_enkf, only: &
       irobs,&
       enkf_variable_names,&
       nrens,&
       mat_single_out,&
       num_single_out,&
       iassout_single_start
  
  implicit none

  character (len=4) :: ax, ay, az

  character (len=4) :: avariable

  character (len=3), intent(in) :: a_before_after

  integer :: i_so, j_so, k_so, ivar_so
  character (len=3) :: single_out_lin_log

  integer :: irens

  double precision :: kz, lz, por
  external kz, lz, por
  double precision, dimension(nrens) :: varholder
  
  logical :: single_cell_ex, single_cell_file_ex

  integer :: i

  do i = 1, num_single_out

     i_so = mat_single_out(i,1)
     j_so = mat_single_out(i,2)
     k_so = mat_single_out(i,3)
     ivar_so = mat_single_out(i,4)
     select case (mat_single_out(i,5))
     case(0)
        single_out_lin_log = 'lin'
     case(1)
        single_out_lin_log = 'log'
     case default
        write(unit = *, fmt = *) '[E4] Error in enkf_output_single_cell()'
     end select

     !Error message if more than four digits
     if( i_so > 9999 &
          .or. j_so > 9999 &
          .or. k_so > 9999 &
          .or. irobs > 9999) then
        write(unit = *, fmt = *) "[E1] Error in enkf_output_single_cell()"
        stop
     end if

     !Error message if bef/aft/ini wrong
     if( a_before_after /= 'bef' &
          .and. a_before_after /= 'aft' &
          .and. a_before_after /= 'ini') then
        write(unit = *, fmt = *) "[E2] Error in enkf_output_single_cell()"
        stop
     end if

     !Write location to character variables
     write(ax, fmt = '(i4.4)') i_so
     write(ay, fmt = '(i4.4)') j_so
     write(az, fmt = '(i4.4)') k_so
     write(avariable, fmt = '(i4.4)') ivar_so

     !Make directory (in case it does not exist)
     single_cell_ex = .false.
     inquire(file="single_cell_output/exist.txt", exist = single_cell_ex)

     if(.not. single_cell_ex) then
        call sys_mkdir(ssingle_cell_outdir)
        open(unit = 79, file = "single_cell_output/exist.txt", status = "new")
        write(unit = 79, fmt = *) "SHEMAT knows: single_cell_output exists"
        close(unit = 79)
     end if

     ! Open/Reopen file
     inquire(file=ssingle_cell_outdir//'single_cell'//&
          trim(project_sfx(1))//'_'//trim(adjustl(ax))//&
          '_'//trim(adjustl(ay))//'_'//trim(adjustl(az))//&
          '_'//trim(adjustl(avariable))//'_'//a_before_after//&
          trim('.txt'),&
          exist = single_cell_file_ex)
     if (single_cell_file_ex) then
        open(unit=19,file=ssingle_cell_outdir//'single_cell'//&
             trim(project_sfx(1))//'_'//trim(adjustl(ax))//&
             '_'//trim(adjustl(ay))//'_'//trim(adjustl(az))//&
             '_'//trim(adjustl(avariable))//'_'//a_before_after//&
             trim('.txt'),&
             status="old",position="append",action="write")
     else 
        open(unit=19,file=ssingle_cell_outdir//'single_cell'//&
             trim(project_sfx(1))//'_'//trim(adjustl(ax))//&
             '_'//trim(adjustl(ay))//'_'//trim(adjustl(az))//&
             '_'//trim(adjustl(avariable))//'_'//a_before_after//&
             trim('.txt'),&
             status="new", action="write")
     end if

     ! Write header
     if(a_before_after == 'ini' .or. irobs == iassout_single_start) then
        write(unit = 19, fmt = *) "Variable:    ", &
             enkf_variable_names(mod(ivar_so,10))
        write(unit = 19, fmt = *) "N:    ", nrens
        write(unit = 19, fmt = *) "i:    ", i_so
        write(unit = 19, fmt = *) "j:    ", j_so
        write(unit = 19, fmt = *) "k:    ", k_so
     end if

     ! Write values
     select case (ivar_so)
     case(1,11)
        write(unit = 19, fmt = '(100es18.10)') &
             (head(i_so,j_so,k_so,irens), irens=1,nrens)
        
     case(2,12)
        write(unit = 19, fmt = '(100es18.10)') &
             (temp(i_so,j_so,k_so,irens), irens=1,nrens)
        
     case(3,13)
        select case (single_out_lin_log)
        case('lin')
           write(unit = 19, fmt = '(100es18.10)') &
                (conc(i_so,j_so,k_so,1,irens) , irens=1,nrens)
        case('log')
           do irens=1,nrens
              varholder(irens) = log(conc(i_so,j_so,k_so,1,irens))/log(10.0d0)
           end do
           write(unit = 19, fmt = '(100es18.10)') (varholder(irens), irens=1,nrens)
        case default
           write(unit = *, fmt = *) "[E5] Error in enkf_output_single_cell"
        end select
        
     case(4,14)
        select case (single_out_lin_log)
        case('lin')
           write(unit = 19, fmt = '(100es18.10)') &
                (kz(i_so,j_so,k_so,irens),irens=1,nrens)
        case('log')
           do irens=1,nrens
              varholder(irens) = log(kz(i_so,j_so,k_so,irens))/log(10.0d0)
           end do
           write(unit = 19, fmt = '(100es18.10)') (varholder(irens), irens=1,nrens)
        case default
           write(unit = *, fmt = *) "[E4] Error in enkf_output_single_cell"
        end select
        
     case(5,15)
        write(unit = 19, fmt = '(100es18.10)') &
             (lz(i_so,j_so,k_so,irens) , irens=1,nrens)

     case(6,16)
        write(unit = 19, fmt = '(100es18.10)') &
             (por(i_so,j_so,k_so,irens) , irens=1,nrens)
     case default
        write(unit = *, fmt = *) "[E3] Error in enkf_output_single_cell()"
        stop
     end select

     close(unit = 19)

  end do

end subroutine enkf_output_single_cell
