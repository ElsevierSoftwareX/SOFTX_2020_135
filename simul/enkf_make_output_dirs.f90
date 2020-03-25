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

!> @brief Make output directores enkf_output and 05_output
!> @param[in] is_sample If it is samples dir or not.
!> @details
!> Check for the existence using the files `exist.txt` and if the
!> directories do not exist, make them.
subroutine enkf_make_output_dirs(is_sample)

  use mod_genrl, only: &
       runmode

  use mod_simul, only: &
       s05_outdir,&
       senkf_outdir,&
       ssample_outdir

  implicit none

  logical, intent(in) :: is_sample

  logical :: enkf_ex
  logical :: dir_05_ex
  logical :: samples_ex

  if(is_sample) then

     ! Creates the samples output directory
     samples_ex = .false.
     inquire(file="samples_output/exist.txt", exist = samples_ex)

     if(.not. samples_ex) then
        CALL sys_mkdir(ssample_outdir)
        open(unit = 79, file = "samples_output/exist.txt", status = 'new')
        write(unit = 79, fmt = *) "SHEMAT knows: samples_output exists"
        close(unit = 79)
     end if

  else

     ! Check existence of directories
     enkf_ex = .false.
     dir_05_ex = .false.
     inquire(file="enkf_output/exist.txt", exist = enkf_ex)
     inquire(file="05_output/exist.txt", exist = dir_05_ex)

     ! Creates ENKF output directories
     IF (runmode>=2) THEN
        if(.not. enkf_ex) then
           CALL sys_mkdir(senkf_outdir)
           open(unit = 79, file = "enkf_output/exist.txt", status = 'new')
           write(unit = 79, fmt = *) "SHEMAT knows: enkf_output exists"
           close(unit = 79)
        end if
        if(.not. dir_05_ex) then
           CALL sys_mkdir(s05_outdir)
           open(unit = 79, file = "05_output/exist.txt", status = 'new')
           write(unit = 79, fmt = *) "SHEMAT knows: 05_output exists"
           close(unit = 79)
        end if
     END IF

  end if

end subroutine enkf_make_output_dirs
