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

subroutine enkf_variable_dist()

  use arrays, only: &
       head,&
       temp,&
       conc,&
       propunit,&
       ibc_data,&
       dbc_data


  use mod_genrl, only: &
       i0,&
       j0,&
       k0

  use mod_enkf, only: &
       nrens,&
       h_dist_in,&
       h_dist_len,&
       num_ref_dist,&
       ref_dist_file_name,&
       ref_dist_var
  

  implicit none

  integer :: i_dist, j_dist
  double precision :: h_dist_dummy
  double precision, allocatable, dimension(:,:,:,:) :: ran_num_array
  double precision :: ran_num
  double precision :: ran_num_2

  integer :: i,j,k,irens, l, i_ref_dist
  integer :: ran_ind, ran_ind_same


  !h_dist_len = 1000

  
  ! Loop for input reference distributions
  do i_ref_dist = 1, num_ref_dist
     
     !Read the reference distribution
     open(unit = 23, file= ref_dist_file_name(i_ref_dist))

     read(unit = 23, fmt = *)
     read(unit = 23, fmt = *) h_dist_len
     read(unit = 23, fmt = *)



     allocate(h_dist_in(h_dist_len))
     allocate(ran_num_array(i0,j0,k0,nrens))


     do i_dist = 1, h_dist_len
        read(unit = 23, fmt = *) h_dist_in(i_dist)
     end do

     close(unit = 23)

     !Sort the reference distribution
     do i_dist = 1, h_dist_len-1

        do j_dist = i_dist+1, h_dist_len

           if(h_dist_in(j_dist) < h_dist_in(i_dist)) then
              h_dist_dummy = h_dist_in(i_dist)
              h_dist_in(i_dist) = h_dist_in(j_dist)
              h_dist_in(j_dist) = h_dist_dummy
           end if

        end do

     end do


     ! Array of random numbers between 0 and 1
     ! uniform distribution
     call random_number(ran_num_array)



     do irens = 1, nrens
        ! Same noise is added for all values of one ensemble
        ran_ind_same = int(ran_num_array(1,1,1,irens)*h_dist_len,4) + 1
        
        do k = 1, k0
           do j = 1, j0 
              do i = 1, i0

                 ! Different value for every grid point
                 ran_ind = int(ran_num_array(i,j,k,irens)*h_dist_len,4) + 1

                 if(ref_dist_var(i_ref_dist) > 3) then
                    l = (k-1)*i0*j0&
                         + (j-1)*i0&
                         + i&
                         + 1
                 else
                    l = 0
                 end if

                 select case (ref_dist_var(i_ref_dist))
                 case(1)
                    !ADDED NOISE TO INITIAL CONDITIONS
                    head(i,j,k,irens) = head(i,j,k,irens) + h_dist_in(ran_ind) 
                 case(2)
                    temp(i,j,k,irens) = h_dist_in(ran_ind) 
                 case(3)
                    ! ADDED NOISE TO INITIAL CONC
!!$                    conc(i,j,k,1,irens) = 10.0d0**(dlog10(conc(i,j,k,1,irens)) + h_dist_in(ran_ind))
                    conc(i,j,k,1,irens) = conc(i,j,k,1,irens) + h_dist_in(ran_ind)
                    if(conc(i,j,k,1,irens) < 0.0d0) then
                       write(unit = *, fmt = *) "[E1.1] negative concs in enkf_variable_dist()"
                       stop
                    end if
                 case(4)
                    !ADDED NOISE TO INITIAL PARAMTER kz
                    propunit(l,4,irens) = 10.0d0**(dlog10(propunit(l,4,irens)) + h_dist_in(ran_ind_same))
                 case(5)
                    propunit(l,8,irens) = h_dist_in(ran_ind) 
                 case(6)
                    propunit(l,1,irens) = h_dist_in(ran_ind) 
                 case default
                    write(unit = *, fmt = *) "[E2] Error in enkf_variable_dist()"
                    stop
                 end select
                 call random_number(ran_num)
                 call random_number(ran_num_2)
                 select case (ref_dist_var(i_ref_dist))
                 case(1)
                    head(i,j,k,irens) = head(i,j,k,irens)*(1 + 1.0d-6 * ( 1.0d0 - 2.0d0 * ran_num))
                 case(2)
                    temp(i,j,k,irens) = temp(i,j,k,irens)*(1 + 1.0d-6 * ( 1.0d0 - 2.0d0 * ran_num))
                 case(3)
                    conc(i,j,k,1,irens) = conc(i,j,k,1,irens)*(1 + 1.0d-6 * ( 1.0d0 - 2.0d0 * ran_num))
                 case(4)
                    propunit(l,4,irens) = propunit(l,4,irens)*(1 + 1.0d-6 * ( 1.0d0 - 2.0d0 * ran_num))
                    conc(i,j,k,1,irens) = conc(i,j,k,1,irens)*(1 + 1.0d-6 * ( 1.0d0 - 2.0d0 * ran_num))
                 case(5)
                    propunit(l,8,irens) = propunit(l,8,irens)*(1 + 1.0d-6 * ( 1.0d0 - 2.0d0 * ran_num))
                 case(6)
                    propunit(l,1,irens) = propunit(l,1,irens)*(1 + 1.0d-6 * ( 1.0d0 - 2.0d0 * ran_num))
                 case default
                    write(unit = *, fmt = *) "[E3] Error in enkf_variable_dist()"
                    stop
                 end select

              end do
           end do
        end do
     end do


     if(allocated(h_dist_in)) deallocate(h_dist_in)
     if(allocated(ran_num_array)) deallocate(ran_num_array)


  end do

end subroutine enkf_variable_dist
