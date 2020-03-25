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

!> @brief Residual output
!> @details
!> Residual output.
subroutine enkf_output_resid()

  use arrays, only:&
       project_sfx

  use mod_simul, only:&
       senkf_outdir

  use mod_genrl, only: &
       key_char,&
       tec_out

  use mod_enkf, only:&
       resvar,&
       obst,&
       nrobs_int,&
       lstate,&
       dual_enkf_switch,&
       nrens,&
       vtk_out_resid

  implicit none

  integer :: i, ii, ni
  character (len=12), dimension(26), parameter :: variable_names_vtk=&
               (/ "obsnr       ", "obstime     ", "res_head_bef", "rms_head_bef",&
               "res_temp_bef", "rms_temp_bef", "res_conc_bef", "rms_conc_bef",&
               "res_kz_bef  ", "rms_kz_bef  ", "res_lz_bef  ", "rms_lz_bef  ",&
               "res_por_bef ", "rms_por_bef ", "res_head_aft", "rms_head_aft",&
               "res_temp_aft", "rms_temp_aft", "res_conc_aft", "rms_conc_aft",&
               "res_kz_aft  ", "rms_kz_aft  ", "res_lz_aft  ", "rms_lz_aft  ",&
               "res_por_aft ", "rms_por_aft " /)
      

     !-------------------------------------------------------------------
     !-------------------------- TECPLOT --------------------------------
     !-------------------------------------------------------------------

  if(tec_out) then

     OPEN(unit=17,file=senkf_outdir//'residual'//trim(project_sfx(1))//'.plt')

     !----------------------------------------------------
     !-------------- WRITING THE OUTPUT ------------------
     !----------------------------------------------------

     WRITE(17,'(1A)') 'TITLE = "Residuals"'
     WRITE(17,'(7A)') &
          'VARIABLES = "obsnr", "obstime", &
          &"res_head_bef", ', '"rms_head_bef", "res_temp_bef", "rms_temp_bef", &
          &"res_conc_bef", "rms_conc_bef", ', '"res_kz_bef", "rms_kz_bef", "res_lz_bef", "rms_lz_bef", &
          &', '"res_por_bef", "rms_por_bef", &
          &"res_head_aft", ', '"rms_head_aft", "res_temp_aft", "rms_temp_aft", &
          &"res_conc_aft", "rms_conc_aft", ', '"res_kz_aft", "rms_kz_aft", "res_lz_aft", "rms_lz_aft", &
          &', '"res_por_aft", "rms_por_aft"'
112  FORMAT (2X,I5,6X,25(D11.5,4X))
     DO ii = 1, nrobs_int
        WRITE(17,112) ii, obst(ii), resvar(2*ii - 1,1), &
             dsqrt(resvar(2*ii - 1,1)/lstate), resvar(2*ii - 1,2), &
             dsqrt(resvar(2*ii - 1,2)/lstate), resvar(2*ii - 1,3), &
             dsqrt(resvar(2*ii - 1,3)/lstate), resvar(2*ii - 1,4), &
             dsqrt(resvar(2*ii - 1,4)/lstate), resvar(2*ii - 1,5), &
             dsqrt(resvar(2*ii - 1,5)/lstate), resvar(2*ii - 1,6), &
             dsqrt(resvar(2*ii - 1,6)/lstate), &
             resvar(2*ii,1), &
             dsqrt(resvar(2*ii,1)/lstate), resvar(2*ii,2), &
             dsqrt(resvar(2*ii,2)/lstate), resvar(2*ii,3), &
             dsqrt(resvar(2*ii,3)/lstate), resvar(2*ii,4), &
             dsqrt(resvar(2*ii,4)/lstate), resvar(2*ii,5), &
             dsqrt(resvar(2*ii,5)/lstate), resvar(2*ii,6), &
             dsqrt(resvar(2*ii,6)/lstate) 
     END DO
     CLOSE(17)

  end if

     !-------------------------------------------------------------------
     !-------------------------- PARAVIEW -------------------------------
     !-------------------------------------------------------------------

  if(vtk_out_resid) then

     !Variable names
     
     !Open the file
     open(unit = 17, file = senkf_outdir//'residual'//&
          trim(project_sfx(1))//'.vtk')

     !Number of outputs depending on dual enkf
     if(dual_enkf_switch) then
        ni = nrobs_int/2
     else
        ni = nrobs_int
     end if

     !----------------------------------------------------
     !-------------- WRITING THE OUTPUT ------------------
     !----------------------------------------------------
     
     !HEADER
     WRITE(unit = 17, fmt = '(a/a/a/a/a,3I7)') key_char//&
          ' vtk DataFile Version 2.0', &
          'residuals',&
          'ASCII',&
          'DATASET RECTILINEAR_GRID', &
          'DIMENSIONS  ', ni+1, 1, 1

     !X_COORDINATES
     WRITE(17,'(a,I7,a)') 'X_COORDINATES  ', ni+1, ' float'
     write(unit = 17, fmt = '(10e16.6)') (real(i,8),i=0,ni)

     !Y_COORDINATES
     WRITE(17,'(a,I7,a)') 'Y_COORDINATES  ', 1, ' float'
     write(unit = 17, fmt = '(10e16.6)') 1.0d0


     !Z_COORDINATES
     WRITE(17,'(a,I7,a)') 'Z_COORDINATES  ', 1, ' float'
     WRITE(17,'(10e16.6)') 1.0d0

     !VARIABLES
     WRITE(17,'(/a,i8)') 'CELL_DATA', ni
     do i=1,26

        ! When to write out residuals
        WRITE(17,'(a)') 'SCALARS '//variable_names_vtk(i)//' float 1'
        WRITE(17,'(a)') 'LOOKUP_TABLE default:mean'
        
        ! Input of variables 
        if(.not. dual_enkf_switch) then
           select case(i)
           case(1)
              write(unit = 17, fmt = '(10i5)') (ii,ii=1,ni)
           case(2)
              write(unit = 17, fmt = '(10e16.6)') (obst(ii),ii=1,ni)
           case(3)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii-1,1),ii=1,ni)
           case(4)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii-1,1)/lstate),ii=1,ni)
           case(5)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii-1,2),ii=1,ni)
           case(6)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii-1,2)/lstate),ii=1,ni)
           case(7)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii-1,3),ii=1,ni)
           case(8)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii-1,3)/lstate),ii=1,ni)
           case(9)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii-1,4),ii=1,ni)
           case(10)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii-1,4)/lstate),ii=1,ni)
           case(11)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii-1,5),ii=1,ni)
           case(12)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii-1,5)/lstate),ii=1,ni)
           case(13)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii-1,6),ii=1,ni)
           case(14)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii-1,6)/lstate),ii=1,ni)
           case(15)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii,1),ii=1,ni)
           case(16)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii,1)/lstate),ii=1,ni)
           case(17)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii,2),ii=1,ni)
           case(18)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii,2)/lstate),ii=1,ni)
           case(19)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii,3),ii=1,ni)
           case(20)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii,3)/lstate),ii=1,ni)
           case(21)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii,4),ii=1,ni)
           case(22)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii,4)/lstate),ii=1,ni)
           case(23)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii,5),ii=1,ni)
           case(24)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii,5)/lstate),ii=1,ni)
           case(25)
              write(unit = 17, fmt = '(10e16.6)') (resvar(2*ii,6),ii=1,ni)
           case(26)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(2*ii,6)/lstate),ii=1,ni)
           case default
              write(unit = *, fmt = *) '[E2] Error in enkf_output_resid'
           end select
        else
           select case(i)
           case(1)
              write(unit = 17, fmt = '(10i5)') (ii,ii=1,ni)
           case(2)
              write(unit = 17, fmt = '(10e16.6)') (obst(2*ii-1),ii=1,ni)
           case(3)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-3,1),ii=1,ni)
           case(4)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-3,1)/lstate),ii=1,ni)
           case(5)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-3,2),ii=1,ni)
           case(6)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-3,2)/lstate),ii=1,ni)
           case(7)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-3,3),ii=1,ni)
           case(8)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-3,3)/lstate),ii=1,ni)
           case(9)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-3,4),ii=1,ni)
           case(10)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-3,4)/lstate),ii=1,ni)
           case(11)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-3,5),ii=1,ni)
           case(12)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-3,5)/lstate),ii=1,ni)
           case(13)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-3,6),ii=1,ni)
           case(14)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-3,6)/lstate),ii=1,ni)
           case(15)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-2,1),ii=1,ni)
           case(16)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-2,1)/lstate),ii=1,ni)
           case(17)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-2,2),ii=1,ni)
           case(18)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-2,2)/lstate),ii=1,ni)
           case(19)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-2,3),ii=1,ni)
           case(20)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-2,3)/lstate),ii=1,ni)
           case(21)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-2,4),ii=1,ni)
           case(22)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-2,4)/lstate),ii=1,ni)
           case(23)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-2,5),ii=1,ni)
           case(24)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-2,5)/lstate),ii=1,ni)
           case(25)
              write(unit = 17, fmt = '(10e16.6)') (resvar(4*ii-2,6),ii=1,ni)
           case(26)
              write(unit = 17, fmt = '(10e16.6)') (dsqrt(resvar(4*ii-2,6)/lstate),ii=1,ni)
           case default
              write(unit = *, fmt = *) '[E2] Error in enkf_output_resid'
           end select
        end if
        
     end do
     
     close(unit = 17)	      	      

  end if

end subroutine enkf_output_resid
