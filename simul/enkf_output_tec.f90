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

!> @brief Write means and variances to Tecplot-file
!> @param[in] afile switch before or after assimilation ('bef'/'aft')
!> @details
!> USE: A lot of character variables for the Tecplotfile-Header, the
!> observation times, the active species of the state vector, a
!> zero and, most importantly, the ENSEMBLE averages and ENSEMBLE
!> variances for each state vector variable.
!>
!> SET: A Tecplotfile called 'enkf_output/assim_variables_E1_100i.plt'
!> containing the positions of the state vector variables in the
!> grid and their Ensemble average and Ensemble variances.
!>
!> Structure of the file:
!> 1. Preparations: Setting variables and strings for the header
!> and the node structure.
!> 2. Cell centered: Writing i, j, k, ave, var for cell centered
!> structure
!> 3. Node centered: Writing i, j, k, ave, var for node centered
!> structure
subroutine enkf_output_tec(afile)

  use arrays, only:&
       project_sfx

  use mod_genrl, only:&
       i0,&
       j0,&
       k0

  use mod_simul, only:&
       senkf_outdir

  use mod_enkf, only:&
       irobs,&
       obst,&
       nstate,&
       n_act_s,&
       lstate,&
       sysindx,&
       ave,&
       var,&
       act_s,&
       zero,&
       cell_centered,&
       dual_enkf_switch

  implicit none

  character (len=3), intent(in) :: afile

  integer :: strandid

  character (len=5) :: tecti
  character (len=15) :: tect
  character (len=160) :: tec
  character (len=5)  :: teci
  character (len=5)  :: tecj
  character (len=5)  :: teck
  character (len=45) :: tec1
  character (len=45) :: tec2
  character (len=22) :: tec2a
  character (len=22) :: tec2b
  character (len=22) :: tec2c
  character (len=22) :: tec2d
  character (len=22) :: tec2e
  character (len=22) :: tec2f
  character (len=22) :: tec3a
  character (len=22) :: tec3b
  character (len=45) :: tec3c
  character (len=55) :: tec4

  integer :: i, j, k, l, l0, l1, l2, l3, l4, l5
  integer :: ifile
  character (len=4) :: arobs

  integer :: line_break_cc

  ! Set ifile (part of filename) and strandid (tecplot time series-id)
  ! depending on before/after assimilation
  ! strandid Identification-integer for time-series output (before
  ! assimilation: 1, after assimilation:2)
  if(dual_enkf_switch) then
     if (afile == 'bef') then
        ifile = (irobs + mod(irobs,2))/2
        strandid = 1
     else if (afile == 'aft') then
        ifile = (irobs + mod(irobs,2))/2
        strandid = 2
     else
        write(unit = *, fmt = *) "[E1] Error in enkf_output_tec(dual)"
     end if
  else
     if (afile == 'bef') then
        ifile = irobs
        strandid = 1
     else if (afile == 'aft') then
        ifile = irobs
        strandid = 2
     else
        write(unit = *, fmt = *) "[E2] Error in enkf_output_tec"
     end if
  end if

  ! Number of Cells (cell- or node-centered)
  if(cell_centered) then
     WRITE(teci,'(i5)')i0+1
     WRITE(tecj,'(i5)')j0+1
     WRITE(teck,'(i5)')k0+1
  else
     WRITE(teci,'(i5)')i0
     WRITE(tecj,'(i5)')j0
     WRITE(teck,'(i5)')k0
  end if

  ! Strings for Tecplot file header
  tec1 = 'TITLE = ' // '"assim variables"'
  tec2 = 'VARIABLES = ' // '"i"' // ',' // '"j"' // ',' // '"k"'
  tec2a = ',' // '"ass_kz"' // ',' // '"std_kz"'
  tec2b = ',' // '"ass_temp"' // ',' // '"std_temp"'
  tec2c = ',' // '"ass_head"' // ',' // '"std_head"'
  tec2d = ',' // '"ass_conc"' // ',' // '"std_conc"'
  tec2e = ',' // '"ass_lz"' // ',' // '"std_lz"'
  tec2f = ',' // '"ass_por"' // ',' // '"std_por"'
  tec3a = 'ZONE  SolutionTime='
  tec3b = ', StrandID='

  if(cell_centered) then
     tec3c = ', I= '//teci//',J='//tecj//', K='//teck//', f=BLOCK'     
     if (n_act_s == 1) then
        tec4 = 'VARLOCATION=([1,2,3]=NODAL, [4,5]=CELLCENTERED)          '
     else if (n_act_s == 2) then
        tec4 = 'VARLOCATION=([1,2,3]=NODAL, [4,5,6,7]=CELLCENTERED)        '
     else if (n_act_s == 3) then
        tec4 = 'VARLOCATION=([1,2,3]=NODAL, [4,5,6,7,8,9]=CELLCENTERED)      '
     else if (n_act_s == 4) then
        tec4 = 'VARLOCATION=([1,2,3]=NODAL, [4,5,6,7,8,9,10,11]=CELLCENTERED)    '
     else if (n_act_s == 2) then
        tec4 = 'VARLOCATION=([1,2,3]=NODAL, [4,5,6,7,8,9,10,11,12,13]=CELLCENTERED)  '
     else if (n_act_s == 2) then
        tec4 = 'VARLOCATION=([1,2,3]=NODAL, [4,5,6,7,8,9,10,11,12,13,14,15]=CELLCENTERED)'
     end if
  else
     tec3c = ', I= '//teci//',J='//tecj//', K='//teck//', f=POINT'
  end if

  !--------------------------------------------------------------------------
  if(cell_centered) then

     WRITE(arobs,'(i4.4)') ifile

     if (dual_enkf_switch) then
        if(mod(irobs,2) == 1) then
           OPEN(unit=14,file=senkf_outdir//'assim_variables'//trim(project_sfx(1))//'_'//afile &
                //'_param_'//arobs//trim('.plt'))
        else
           OPEN(unit=14,file=senkf_outdir//'assim_variables'//trim(project_sfx(1))//'_'//afile &
             //'_state_'//arobs//trim('.plt'))
        end if
     else
        OPEN(unit=14,file=senkf_outdir//'assim_variables'//trim(project_sfx(1))//'_'//afile &
             //'_'//arobs//trim('.plt'))
     end if

     OPEN(unit=14,file=senkf_outdir//'assim_variables'//trim(project_sfx(1))//'_'//afile &
          //'_'//arobs//trim('.plt'))
     WRITE(tecti,'(i5)') strandid
     WRITE(tect,'(e15.8)') obst(irobs)
     WRITE(14,*) tec1
     tec = tec2

     IF (act_s(1)==1) tec = trim(tec) // tec2c 
     IF (act_s(2)==1) tec = trim(tec) // tec2b
     IF (act_s(3)==1) tec = trim(tec) // tec2d
     IF (act_s(4)==1) tec = trim(tec) // tec2a
     IF (act_s(5)==1) tec = trim(tec) // tec2e
     IF (act_s(6)==1) tec = trim(tec) // tec2f
     WRITE(14,'(a160)') tec
     write(unit = 14, fmt = '(a109)') tec3a//tect//tec3b//tecti//tec3c
     write(unit = 14, fmt = *) tec4

     line_break_cc = 0

     do k = 1, k0+1
        do j = 1, j0+1
           do i = 1, i0+1
              write(unit = 14, fmt = '(i6)', advance = 'no') i 
              ! Line Break after 100 entries
              line_break_cc = line_break_cc + 1
              if(line_break_cc >= 100) then
                 write(unit = 14, fmt = *)
                 line_break_cc = 0
              end if
           end do
        end do
     end do

     line_break_cc = 0
     write(unit = 14, fmt = *) 
     do k = 1, k0+1
        do j = 1, j0+1
           do i = 1, i0+1
              write(unit = 14, fmt = '(i6)', advance = 'no') j
              ! Line Break after 100 entries
              line_break_cc = line_break_cc + 1
              if(line_break_cc >= 100) then
                 write(unit = 14, fmt = *)
                 line_break_cc = 0
              end if
           end do
        end do
     end do

     line_break_cc = 0
     write(unit = 14, fmt = *) 
     do k = 1, k0+1
        do j = 1, j0+1
           do i = 1, i0+1
              write(unit = 14, fmt = '(i6)', advance = 'no') k 
              ! Line Break after 100 entries
              line_break_cc = line_break_cc + 1
              if(line_break_cc >= 100) then
                 write(unit = 14, fmt = *)
                 line_break_cc = 0
              end if
           end do
        end do
     end do

     line_break_cc = 0
     write(unit = 14, fmt = *) 
     write(unit = 14, fmt = *) 
     write(unit = 14, fmt = *)
     if( n_act_s >= 1) then 
        do k = 1, k0
           do j = 1, j0
              do i = 1, i0
                 l = (k-1)*i0*j0 + (j-1)*i0 + i
                 l0 = sysindx(l)
                 write(unit = 14, fmt = '(es16.8,1x)', advance = 'no') ave(l0) 
                 ! Line Break after 100 entries
                 line_break_cc = line_break_cc + 1
                 if(line_break_cc >= 100) then
                    write(unit = 14, fmt = *)
                    line_break_cc = 0
                 end if
              end do
           end do
        end do

        line_break_cc = 0
        write(unit = 14, fmt = *) 
        do k = 1, k0
           do j = 1, j0
              do i = 1, i0
                 l = (k-1)*i0*j0 + (j-1)*i0 + i
                 l0 = sysindx(l)
                 write(unit = 14, fmt = '(es16.8,1x)', advance = 'no') dsqrt(var(l0)+1.0d-6)
                 ! Line Break after 100 entries
                 line_break_cc = line_break_cc + 1
                 if(line_break_cc >= 100) then
                    write(unit = 14, fmt = *)
                    line_break_cc = 0
                 end if
              end do
           end do
        end do

        line_break_cc = 0
        write(unit = 14, fmt = *) 
        if( n_act_s >= 2) then
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    l = (k-1)*i0*j0 + (j-1)*i0 + i
                    l1 = lstate + sysindx(l)
                    write(unit = 14, fmt = '(es16.8,1x)', advance = 'no') ave(l1) 
                    ! Line Break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = 14, fmt = *)
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do

           line_break_cc = 0
           write(unit = 14, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    l = (k-1)*i0*j0 + (j-1)*i0 + i
                    l1 = lstate + sysindx(l)
                    write(unit = 14, fmt = '(es16.8,1x)', advance = 'no') dsqrt(var(l1)+1.0d-6)
                    ! Line Break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = 14, fmt = *)
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do

           line_break_cc = 0
           write(unit = 14, fmt = *)            
           if( n_act_s >= 3) then
              do k = 1, k0
                 do j = 1, j0
                    do i = 1, i0
                       l = (k-1)*i0*j0 + (j-1)*i0 + i
                       l2 = 2*lstate + sysindx(l)
                       write(unit = 14, fmt = '(es16.8,1x)', advance = 'no') ave(l2) 
                       ! Line Break after 100 entries
                       line_break_cc = line_break_cc + 1
                       if(line_break_cc >= 100) then
                          write(unit = 14, fmt = *)
                          line_break_cc = 0
                       end if
                    end do
                 end do
              end do

              line_break_cc = 0
              write(unit = 14, fmt = *) 
              do k = 1, k0
                 do j = 1, j0
                    do i = 1, i0
                       l = (k-1)*i0*j0 + (j-1)*i0 + i
                       l2 = 2*lstate + sysindx(l)
                       write(unit = 14, fmt = '(es16.8,1x)', advance = 'no') dsqrt(var(l2)+1.0d-6)
                       ! Line Break after 100 entries
                       line_break_cc = line_break_cc + 1
                       if(line_break_cc >= 100) then
                          write(unit = 14, fmt = *)
                          line_break_cc = 0
                       end if
                    end do
                 end do
              end do
              write(unit = 14, fmt = *) 

           end if
        end if
     end if

     !--------------------------------------------------------------------------
  else
     !--------------------------------------------------------------------------

     WRITE(arobs,'(i4.4)') ifile

     OPEN(unit=14,file=senkf_outdir//'assim_variables'//trim(project_sfx(1))//'_'//afile &
          //'_'//arobs//trim('.plt'))
     WRITE(tecti,'(i5)') strandid
     WRITE(tect,'(e9.2)') obst(irobs)
     WRITE(14,*) tec1
     tec = tec2

     IF (act_s(1)==1) tec = trim(tec) // tec2c 
     IF (act_s(2)==1) tec = trim(tec) // tec2b
     IF (act_s(3)==1) tec = trim(tec) // tec2d
     IF (act_s(4)==1) tec = trim(tec) // tec2a
     IF (act_s(5)==1) tec = trim(tec) // tec2e
     IF (act_s(6)==1) tec = trim(tec) // tec2f
     WRITE(14,'(a160)') tec
     WRITE(14,*) tec3a//tect//tec3b//tecti//tec3c

     DO k = 1, k0
        DO j = 1, j0
           DO i = 1, i0
              l = (k-1)*i0*j0 + (j-1)*i0 + i
              IF (sysindx(l) > 0) THEN
                 l0 = sysindx(l)
                 l1 = lstate + l0
                 l2 = 2*lstate + l0
                 l3 = 3*lstate + l0
                 l4 = 4*lstate + l0
                 l5 = 5*lstate + l0
                 IF (n_act_s==1) WRITE(14,*) i, j, k, ave(l0), &
                      dsqrt(var(l0)+0.000001D0)
                 IF (n_act_s==2) WRITE(14,*) i, j, k, ave(l0), &
                      dsqrt(var(l0)+0.000001D0), ave(l1), &
                      dsqrt(var(l1)+0.000001D0)
                 IF (n_act_s==3) WRITE(14,*) i, j, k, ave(l0), &
                      dsqrt(var(l0)+0.000001D0), ave(l1), &
                      dsqrt(var(l1)+0.000001D0), ave(l2), &
                      dsqrt(var(l2)+0.000001D0)
                 IF (n_act_s==4) WRITE(14,*) i, j, k, ave(l0), &
                      dsqrt(var(l0)+0.000001D0), ave(l1), &
                      dsqrt(var(l1)+0.000001D0), ave(l2), &
                      dsqrt(var(l2)+0.000001D0), ave(l3), &
                      dsqrt(var(l3)+0.000001D0)
                 IF (n_act_s==5) WRITE(14,*) i, j, k, ave(l0), &
                      dsqrt(var(l0)+0.000001D0), ave(l1), &
                      dsqrt(var(l1)+0.000001D0), ave(l2), &
                      dsqrt(var(l2)+0.000001D0), ave(l3), &
                      dsqrt(var(l3)+0.000001D0), ave(l4), &
                      dsqrt(var(l4)+0.000001D0)
                 IF (n_act_s==6) WRITE(14,*) i, j, k, ave(l0), &
                      dsqrt(var(l0)+0.000001D0), ave(l1), &
                      dsqrt(var(l1)+0.000001D0), ave(l2), &
                      dsqrt(var(l2)+0.000001D0), ave(l3), &
                      dsqrt(var(l3)+0.000001D0), ave(l4), &
                      dsqrt(var(l4)+0.000001D0), ave(l5), &
                      dsqrt(var(l5)+0.000001D0)
              ELSE
                 IF (n_act_s==1) WRITE(14,*) i, j, k, zero, zero
                 IF (n_act_s==2) WRITE(14,*) i, j, k, zero, zero, &
                      zero, ZERO 
                 IF (n_act_s==3) WRITE(14,*) i, j, k, zero, zero, &
                      zero, zero, zero, zero
                 IF (n_act_s==4) WRITE(14,*) i, j, k, zero, zero, &
                      zero, zero, zero, zero, zero, zero
                 IF (n_act_s==5) WRITE(14,*) i, j, k, zero, zero, &
                      zero, zero, zero, zero, zero, zero, zero, zero
                 IF (n_act_s==6) WRITE(14,*) i, j, k, zero, zero, &
                      zero, zero, zero, zero, zero, zero, zero, zero, &
                      zero, zero
              END IF
           END DO
        END DO
     END DO

  end if
  !--------------------------------------------------------------------------

  CLOSE(14)

end subroutine enkf_output_tec
