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

!> @brief Write means and variances to Paraview-file
!> @param[in] afile switch before or after assimilation ('bef'/'aft')
!> @details
!> Cell centered works only for 2 dimensions.
!>
!> USE: Some dummy vectors, the active species of the state vector
!> and, most importantly, the ENSEMBLE averages and ENSEMBLE
!> variances for each state vector variable.
!>
!> SET: A Paraviewfile called
!> 'enkf_output/assim_variables_E1_100i.vtk' containing the
!> positions of the state vector variables in the grid and their
!> Ensemble average and Ensemble variances.
subroutine enkf_output_vtk(afile)

  use arrays, only: &
       project_sfx,&
       delx,&
       dely,&
       delz,&
       delxa,&
       delya,&
       delza

  use mod_genrl, only: &
       i0,&
       j0,&
       k0,&
       key_char

  use mod_simul, only: &
       senkf_outdir
       
  use mod_enkf, only: &
       irobs,&
       n_act_s,&
       lstate,&
       ave,&
       var,&
       act_s,&
       mem,&
       vtk_out_realisations,&
       num_out_realisations,&
       num_enkf_vars,&
       dual_enkf_switch,&
       enkf_variable_names,&
       cell_centered

  implicit none

  character (len=3), intent(in) :: afile

  double precision, dimension(max(i0,j0,k0)) :: val
  
  character (len=40), dimension(6) :: avt
  character (len=40), dimension(6) :: bvt
  character (len=40), dimension(6) :: dvt

  character (len=4) :: arobs

  integer :: lvar, i, j, k, ivar
  integer :: ifile

  integer :: irobs_vtk_out

  integer :: irealisation

  if(dual_enkf_switch) then
     irobs_vtk_out = (irobs + mod(irobs,2))/2
  else
     irobs_vtk_out = irobs
  end if

  if (afile == 'bef') then
     ifile = irobs_vtk_out
  else if (afile == 'aft') then
     ifile = irobs_vtk_out
  else
     write(unit = *, fmt = *) "[E1] Error in enkf_output_vtk"
     stop 1
  end if

  WRITE(arobs,'(i4.4)') ifile

  if (dual_enkf_switch) then
     if(mod(irobs,2) == 1) then
        OPEN(unit=42,file=senkf_outdir//'assim_variables'&
             //trim(project_sfx(1))//'_'//afile &
             //'_param_'//arobs//trim('.vtk'))
     else
        OPEN(unit=42,file=senkf_outdir//'assim_variables'&
             //trim(project_sfx(1))//'_'//afile &
             //'_'//arobs//trim('.vtk'))
     end if
  else
     OPEN(unit=42,file=senkf_outdir//'assim_variables'&
          //trim(project_sfx(1))//'_'//afile &
          //'_'//arobs//trim('.vtk'))
  end if

  !header
  if(cell_centered) then
     WRITE(42,'(a/a/a//a/a,3I7)') key_char//&
          ' vtk DataFile Version 2.0', &
          'assim_variables1000', 'ASCII',&
          'DATASET RECTILINEAR_GRID', &
          'DIMENSIONS  ', i0+1, j0+1, k0
  else
     WRITE(42,'(a/a/a//a/a,3I7)') key_char//&
          ' vtk DataFile Version 2.0', &
          'assim_variables1000', 'ASCII',&
          'DATASET RECTILINEAR_GRID', &
          'DIMENSIONS  ', i0, j0, k0
  end if

  ! grid structure
  if(cell_centered) then
     WRITE(42,'(a,I7,a)') 'X_COORDINATES  ', i0+1, ' float'
     WRITE(42,'(10e16.6)') (delxa(i)-0.5d0*delx(i),i=1,i0)
     write(unit = 42, fmt = '(e16.6)') delxa(i0)+0.5d0*delx(i0)

     WRITE(42,'(a,I7,a)') 'Y_COORDINATES  ', j0+1, ' float'
     WRITE(42,'(10e16.6)') (delya(i)-0.5d0*dely(i),i=1,j0)
     write(unit = 42, fmt = '(e16.6)') delya(j0)+0.5d0*dely(j0)

     WRITE(42,'(a,I7,a)') 'Z_COORDINATES  ', k0, ' float'
     WRITE(42,'(10e16.6)') (delza(i),i=1,k0)

     WRITE(42,'(/a,i8)') 'CELL_DATA', i0*j0*k0
  else
     WRITE(42,'(a,I7,a)') 'X_COORDINATES  ', i0, ' float'
     WRITE(42,'(10e16.6)') (delxa(i),i=1,i0)

     WRITE(42,'(a,I7,a)') 'Y_COORDINATES  ', j0, ' float'
     WRITE(42,'(10e16.6)') (delya(i),i=1,j0)

     WRITE(42,'(a,I7,a)') 'Z_COORDINATES  ', k0, ' float'
     WRITE(42,'(10e16.6)') (delza(i),i=1,k0)

     WRITE(42,'(/a,i8)') 'POINT_DATA', i0*j0*k0
  end if
  
  !state vector variables
  lvar = 1
    do ivar = 1, num_enkf_vars
     if (act_s(ivar)==1) then
        avt(lvar) = trim(enkf_variable_names(ivar))//'_mean'
        bvt(lvar) = trim(enkf_variable_names(ivar))//'_std'
        dvt(lvar) = trim(enkf_variable_names(ivar))//'_'
        lvar = lvar + 1
     end if
  end do
  IF ((lvar-1).ne.n_act_s) then
     PRINT *, '[E2] vtk mean output error'
     STOP 1
  END IF
  
  DO i=1, n_act_s
     k=(i-1)*lstate
     WRITE(42,'(a)') 'SCALARS '//avt(i)//' float 1'
     WRITE(42,'(a)') 'LOOKUP_TABLE default:mean'
     WRITE(42,'(100e16.6)') (ave(k+j),j=1,lstate)
     WRITE(42,'(a)') 'SCALARS '//bvt(i)//' float 1'
     WRITE(42,'(a)') 'LOOKUP_TABLE default'
     WRITE(42,'(100e16.6)') (dsqrt(var(k+j)+0.000000001D0), &
          j=1,lstate)                        
     if (vtk_out_realisations) then
        do irealisation = 1, num_out_realisations
           WRITE(42,'(a,I0.2,a)') 'SCALARS '//trim(dvt(i)),irealisation,' float 1'
           WRITE(42,'(a)') 'LOOKUP_TABLE default:mean'
           WRITE(42,'(100e16.6)') (mem(k+j,irealisation),j=1,lstate)
        end do
     end if
  END DO
  CLOSE(42)	      	      

end subroutine enkf_output_vtk
