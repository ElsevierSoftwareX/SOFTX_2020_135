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

!> @brief Output correlation matrix
!> @param[in] a_before_after switch before or after assimilation ('bef'/'aft')
!> @details
!> This subroutine produces Tecplot or Paraview output for the
!> correlation matrix of all state vector variables with respect to
!> a certain state vector variable specified in the inputfile.
!>
!> __ATTENTION__: This can only be applied to trivial sysindx, i.e. no
!> excluded units.
subroutine enkf_output_covs(a_before_after)

  use arrays, only: &
       delxa,&
       delya,&
       delza,&
       delx,&
       dely,&
       delz

  use mod_simul, only: &
       senkf_outdir
  
  use mod_genrl, only: &
       i0,&
       j0,&
       k0,&
       tec_out,&
       key_char

  use mod_enkf, only:&
       corr_matrix,&
       irobs,&
       mat_cov_ref,&
       num_cov_ref,&
       nstate,&
       obst,&
       sysindx,&
       n_act_s,&
       act_s,&
       lstate,&
       num_enkf_vars,&
       vtk_out_covs
  
  implicit none

  character (len=3), intent (in) :: a_before_after

  integer :: ind_cov_ref
  character (len=4) :: arobs

  character (len=5) :: tecti
  character (len=15) :: tect
  character (len=5) :: teci
  character (len=5) :: tecj
  character (len=5) :: teck
  character (len=4) :: a_i_cov_ref, a_j_cov_ref, a_k_cov_ref, a_ivar_cov_ref
  
  integer :: i, j, k, l, imem
  integer :: strandid
  integer :: active_var
  integer :: ivar

  integer :: i_cov_ref
  integer :: j_cov_ref
  integer :: k_cov_ref
  integer :: ivar_cov_ref

  integer :: line_break_cc

  integer, external :: index_loc_to_mem
  
  ! Select different strandid for Tecplot output depending
  ! on the output time: before or after assimilation
  select case (a_before_after)
  case ('bef')
     strandid = 1
  case ('aft')
     strandid = 2
  case default
     write(unit = *, fmt = *) "[E1] Error in enkf_output_covs"
     stop
  end select

  do ind_cov_ref = 1, num_cov_ref
     
     i_cov_ref = mat_cov_ref(ind_cov_ref,1)
     j_cov_ref = mat_cov_ref(ind_cov_ref,2)
     k_cov_ref = mat_cov_ref(ind_cov_ref,3)
     ivar_cov_ref = mat_cov_ref(ind_cov_ref,4)


     !--------------------------------------------------------
     !------------------------ TEC ---------------------------
     !--------------------------------------------------------
     if(tec_out) then

        ! Characters used inside the file
        write(tecti,'(i5)') strandid
        write(tect, '(e15.8)') obst(irobs)
        !Cell_centered: plus one
        write(teci,'(i5)') i0+1
        write(tecj,'(i5)') j0+1
        write(teck,'(i5)') k0+1

        ! Characters used in filename
        write(arobs,'(i4.4)') irobs
        write(a_i_cov_ref,'(i4.4)') i_cov_ref
        write(a_j_cov_ref,'(i4.4)') j_cov_ref
        write(a_k_cov_ref,'(i4.4)') k_cov_ref
        write(a_ivar_cov_ref,'(i4.4)') ivar_cov_ref

        !Open the file
        open(unit = 32, file = senkf_outdir//'correlation'//&
             '_'//a_i_cov_ref//'_'//a_j_cov_ref//&
             '_'//a_k_cov_ref//'_'//a_ivar_cov_ref//&
             '_'//trim(a_before_after)//'_'//arobs//'.plt')

        ! Write output file
        write(unit = 32, fmt = *) 'TITLE = "correlation"'

        select case (n_act_s)
        case (1)
           write(unit = 32, fmt = *) 'VARIABLES = "i", "j", "k", "corr"'
        case (2)
           write(unit = 32, fmt = *) 'VARIABLES = "i", "j", "k", "corr", "corr2"'
        case (3)
           write(unit = 32, fmt = *) 'VARIABLES = "i", "j", "k", "corr", "corr2", "corr3"'
        case (4)
           write(unit = 32, fmt = *) 'VARIABLES = "i", "j", "k",',&
                '"corr", "corr2", "corr3", "corr4"'
        case (5)
           write(unit = 32, fmt = *) 'VARIABLES = "i", "j", "k",',&
                '"corr", "corr2", "corr3", "corr4", "corr5"'
        case (6)
           write(unit = 32, fmt = *) 'VARIABLES = "i", "j", "k",',&
                '"corr", "corr2", "corr3", "corr4", "corr5"',&
                '"corr6"'
        case default
           write(unit = *, fmt = *) '[E3] Error in enkf_output_covs'
           stop
        end select

        write(unit = 32, fmt = *) 'ZONE SolutionTime=', tect, ', StrandID=', tecti,&
             ', I=', teci, ', J=', tecj , ', K=', teck, ', f=BLOCK'

        select case (n_act_s)
        case(1)
           write(unit = 32, fmt = *) 'VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)'
        case(2)
           write(unit = 32, fmt = *) 'VARLOCATION=([1,2,3]=NODAL, [4,5]=CELLCENTERED)'
        case(3)
           write(unit = 32, fmt = *) 'VARLOCATION=([1,2,3]=NODAL, [4,5,6]=CELLCENTERED)'
        case(4)
           write(unit = 32, fmt = *) 'VARLOCATION=([1,2,3]=NODAL, [4,5,6,7]=CELLCENTERED)'
        case(5)
           write(unit = 32, fmt = *) 'VARLOCATION=([1,2,3]=NODAL, [4,5,6,7,8]=CELLCENTERED)'
        case(6)
           write(unit = 32, fmt = *) 'VARLOCATION=([1,2,3]=NODAL, [4,5,,6,7,8,9]=CELLCENTERED)'
        case default
           write(unit = *, fmt = *) '[E4] Error in enkf_output_covs'
           stop
        end select


        line_break_cc = 0

        do k = 1, k0+1
           do j = 1, j0+1
              do i = 1, i0+1
                 write(unit = 32, fmt = '(i6)', advance = 'no') i
                 ! Line Break after 100 entries
                 line_break_cc = line_break_cc + 1
                 if(line_break_cc >= 100) then
                    write(unit = 32, fmt = *)
                    line_break_cc = 0
                 end if
              end do
           end do
        end do

        write(unit = 32, fmt = *)
        line_break_cc = 0
        write(unit = 32, fmt = *)
        do k = 1, k0+1
           do j = 1, j0+1
              do i = 1, i0+1
                 write(unit = 32, fmt = '(i6)', advance = 'no') j
                 ! Line Break after 100 entries
                 line_break_cc = line_break_cc + 1
                 if(line_break_cc >= 100) then
                    write(unit = 32, fmt = *)
                    line_break_cc = 0
                 end if
              end do
           end do
        end do

        write(unit = 32, fmt = *)
        line_break_cc = 0
        write(unit = 32, fmt = *)
        do k = 1, k0+1
           do j = 1, j0+1
              do i = 1, i0+1
                 write(unit = 32, fmt = '(i6)', advance = 'no') k
                 ! Line Break after 100 entries
                 line_break_cc = line_break_cc + 1
                 if(line_break_cc >= 100) then
                    write(unit = 32, fmt = *)
                    line_break_cc = 0
                 end if
              end do
           end do
        end do

        write(unit = 32, fmt = *)
        line_break_cc = 0
        write(unit = 32, fmt = *)
        do ivar = 1, num_enkf_vars
           if(act_s(ivar) == 1) then

              do k = 1, k0
                 do j = 1, j0
                    do i = 1, i0
                       imem = index_loc_to_mem(i,j,k,ivar)
                       write(unit = 32, fmt = '(es16.8,1x)', advance = 'no') corr_matrix(imem, ind_cov_ref)

                       ! Line Break after 100 entries
                       line_break_cc = line_break_cc + 1
                       if(line_break_cc >= 100) then
                          write(unit = 32, fmt = *)
                          line_break_cc = 0
                       end if
                    end do
                 end do
              end do

           end if
        end do

        close(unit = 32)

     end if

     !------------------------------------------------------------------
     !----------------------------- VTK --------------------------------
     !------------------------------------------------------------------
     if(vtk_out_covs) then

        ! Characters used in filename
        write(arobs,'(i4.4)') irobs
        write(a_i_cov_ref,'(i4.4)') i_cov_ref
        write(a_j_cov_ref,'(i4.4)') j_cov_ref
        write(a_k_cov_ref,'(i4.4)') k_cov_ref
        write(a_ivar_cov_ref,'(i4.4)') ivar_cov_ref

        !Open the file
        open(unit = 32, file = senkf_outdir//'correlation'//&
             '_'//a_i_cov_ref//'_'//a_j_cov_ref//&
             '_'//a_k_cov_ref//'_'//a_ivar_cov_ref//&
             '_'//trim(a_before_after)//'_'//arobs//'.vtk')


        !HEADER
        WRITE(unit = 32, fmt = '(a/a/a/a/a,3I7)') key_char//&
             ' vtk DataFile Version 2.0', &
             'correlation'//&
             '_'//a_i_cov_ref//'_'//a_j_cov_ref//&
             '_'//a_k_cov_ref//'_'//a_ivar_cov_ref//&
             '_'//trim(a_before_after)//'_'//arobs,&
             'ASCII',&
             'DATASET RECTILINEAR_GRID', &
             'DIMENSIONS  ', i0+1, j0+1, k0

        !X_COORDINATES
        WRITE(32,'(a,I7,a)') 'X_COORDINATES  ', i0+1, ' float'
        !WRITE(32,'(10e16.6)') (delxa(i),i=1,i0)
        write(unit = 32, fmt = '(10e16.6)') (delxa(i)-0.5d0*delx(i),i=1,i0)
        write(unit = 32, fmt = '(e16.6)')  delxa(i0)+0.5d0*delx(i0)

        !Y_COORDINATES
        WRITE(32,'(a,I7,a)') 'Y_COORDINATES  ', j0+1, ' float'
        !WRITE(32,'(10e16.6)') (delya(i),i=1,j0)
        write(unit = 32, fmt = '(10e16.6)') (delya(i)-0.5d0*dely(i),i=1,j0)
        write(unit = 32, fmt = '(e16.6)')  delya(j0)+0.5d0*dely(j0)

        !Z_COORDINATES
        WRITE(32,'(a,I7,a)') 'Z_COORDINATES  ', k0, ' float'
        WRITE(32,'(10e16.6)') (delza(i),i=1,k0)
        !write(unit = 32, fmt = '(10e16.6)') (delza(i)-0.5d0*delz(i),i=1,k0)
        !write(unit = 32, fmt = '(e16.6)')  (delza(i)+0.5d0*delz(i),i=1,k0)

        !CORRELATIONS
        !WRITE(32,'(/a,i8)') 'POINT_DATA', i0*j0*k0
        WRITE(32,'(/a,i8)') 'CELL_DATA', i0*j0*k0

        do ivar = 1, num_enkf_vars
           if(act_s(ivar)==1) then

              write(unit = 32, fmt = *)
              write(unit = 32, fmt = '(a,i4.4,a)') 'SCALARS   correlations', ivar,'    float1'
              write(unit = 32, fmt = '(a)') 'LOOKUP_TABLE default'

              line_break_cc = 0
              do k = 1, k0
                 do j = 1, j0
                    do i = 1, i0
                       imem = index_loc_to_mem(i,j,k,ivar)
                       write(unit = 32, fmt = '(es16.8,1x)', advance = 'no') corr_matrix(imem, ind_cov_ref)
                       ! Line Break after 10 entries
                       line_break_cc = line_break_cc + 1
                       if(line_break_cc >= 10) then
                          write(unit = 32, fmt = *)
                          line_break_cc = 0
                       end if
                    end do
                 end do
              end do

           end if
        end do


        close(unit = 32)

     end if

  end do

end subroutine enkf_output_covs
