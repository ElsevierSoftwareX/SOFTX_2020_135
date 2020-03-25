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

!> @brief Output of overall standard deviations.
!> @details
!> Output of overall standard deviations.
subroutine enkf_output_stddev()

  use arrays, only: &
       project_sfx

  use mod_simul, only: &
       senkf_outdir

  use mod_genrl, only: &
       key_char,&
       tec_out

  use mod_enkf, only: &
       nrobs_int,&
       stdvar,&
       obst,&
       dual_enkf_switch,&
       vtk_out_stddev
       
  implicit none 

  integer :: ii, i, nint
  character (len=12), dimension(12), parameter :: variable_names_vtk=&
       (/ "std_head_bef", "std_temp_bef",&
       "std_conc_bef", "std_kz_bef  ", "std_lz_bef  ", "std_por_bef ",&
       "std_head_aft", "std_temp_aft", "std_conc_aft", "std_kz_aft  ",&
       "std_lz_aft  ", "std_por_aft " /)

  !-------------------------------------------------------------------
  !-------------------------- TECPLOT --------------------------------
  !-------------------------------------------------------------------
  if(tec_out) then
     
     OPEN(unit=18,file=senkf_outdir//'std-dev'//trim(project_sfx(1))//'.plt')

     !----------------------------------------------------
     !-------------- WRITING THE OUTPUT ------------------
     !----------------------------------------------------

     WRITE(18,216)
216  FORMAT ('TITLE = ','"Standard deviation"')
     WRITE(18,'(3A)') 'VARIABLES = "obsnr", "obstime", ', &
          ' "std-head_bef", "std_temp_bef", "std_conc_bef", "std_kz_bef", "std_lz_bef", "std_por_bef"', &
          ' "std-head_aft", "std_temp_aft", "std_conc_aft", "std_kz_aft", "std_lz_aft", "std_por_aft"'

110  FORMAT (2X,I5,5X,13(D11.5,2X))
     DO ii = 1, nrobs_int
        WRITE(18,110) ii, obst(ii), &
             stdvar(2*ii-1,1), stdvar(2*ii-1,2), stdvar(2*ii-1,3), &
             stdvar(2*ii-1,4), stdvar(2*ii-1,5), stdvar(2*ii-1,6), &
             stdvar(2*ii,1), stdvar(2*ii,2), stdvar(2*ii,3), &
             stdvar(2*ii,4), stdvar(2*ii,5), stdvar(2*ii,6)
        write(unit = 18, fmt = *) 

     END DO
     CLOSE(18)


  end if

     !-------------------------------------------------------------------
     !-------------------------- PARAVIEW -------------------------------
     !-------------------------------------------------------------------

  if(vtk_out_stddev) then

     !Variable names
     
     !Open the file
     open(unit = 18, file = senkf_outdir//'stddev'//&
          trim(project_sfx(1))//'.vtk')

     if(dual_enkf_switch) then
        nint = nrobs_int/2
     else
        nint = nrobs_int
     end if
     
     !----------------------------------------------------
     !-------------- WRITING THE OUTPUT ------------------
     !----------------------------------------------------

     !HEADER
     WRITE(unit = 18, fmt = '(a/a/a/a/a,3I7)') key_char//&
          ' vtk DataFile Version 2.0', &
          'stddevs',&
          'ASCII',&
          'DATASET RECTILINEAR_GRID', &
          'DIMENSIONS  ', nint+1, 1, 1

     !X_COORDINATES
     WRITE(18,'(a,I7,a)') 'X_COORDINATES  ', nint+1, ' float'
     write(unit = 18, fmt = '(10e16.6)') (real(i,8),i=0,nint)

     !Y_COORDINATES
     WRITE(18,'(a,I7,a)') 'Y_COORDINATES  ', 1, ' float'
     write(unit = 18, fmt = '(10e16.6)') 1.0d0


     !Z_COORDINATES
     WRITE(18,'(a,I7,a)') 'Z_COORDINATES  ', 1, ' float'
     WRITE(18,'(10e16.6)') 1.0d0

     !VARIABLES
     WRITE(18,'(/a,i8)') 'CELL_DATA', nrobs_int

     WRITE(18,'(a)') 'SCALARS '//'obsnr       '//' float 1'
     WRITE(18,'(a)') 'LOOKUP_TABLE default:mean'
     write(unit = 18, fmt = '(10i5)') (ii,ii=1,nint)

     WRITE(18,'(a)') 'SCALARS '//'obstime     '//' float 1'
     WRITE(18,'(a)') 'LOOKUP_TABLE default:mean'
     if(.not. dual_enkf_switch) then
        write(unit = 18, fmt = '(10e16.6)') (obst(ii),ii=1,nint)
     else
        write(unit = 18, fmt = '(10e16.6)') (obst(2*ii-1),ii=1,nint)
     end if

     do i=1,12

        ! Input of variables
        if(.not. dual_enkf_switch) then
           if(i<7) then
              WRITE(18,'(a)') 'SCALARS '//variable_names_vtk(i)//' float 1'
              WRITE(18,'(a)') 'LOOKUP_TABLE default:mean'
              write(unit = 18, fmt = '(10e16.6)') (stdvar(2*ii-1,i),ii=1,nint)
           else if(i<13) then
              WRITE(18,'(a)') 'SCALARS '//variable_names_vtk(i)//' float 1'
              WRITE(18,'(a)') 'LOOKUP_TABLE default:mean'
              write(unit = 18, fmt = '(10e16.6)') (stdvar(2*ii,i-6),ii=1,nint)
           else
              write(unit = *, fmt = *) '[E2] Error in enkf_output_stddev'
              stop 1
           end if
        else
           if(i<7) then
              WRITE(18,'(a)') 'SCALARS '//variable_names_vtk(i)//' float 1'
              WRITE(18,'(a)') 'LOOKUP_TABLE default:mean'
              write(unit = 18, fmt = '(10e16.6)') (stdvar(4*ii-3,i),ii=1,nint)
           else if(i<13) then
              WRITE(18,'(a)') 'SCALARS '//variable_names_vtk(i)//' float 1'
              WRITE(18,'(a)') 'LOOKUP_TABLE default:mean'
              write(unit = 18, fmt = '(10e16.6)') (stdvar(4*ii-2,i-6),ii=1,nint)
           else
              write(unit = *, fmt = *) '[E3] Error in enkf_output_stddev'
              stop 1
           end if
        end if

     end do
     
     close(unit = 18)	      	      

  end if

end subroutine enkf_output_stddev
