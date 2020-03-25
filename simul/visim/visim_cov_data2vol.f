!>Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %\n
!>Junior University.  All rights reserved.                             %\n
!>
!>The programs in GSLIB are distributed in the hope that they will be  %\n
!>useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %\n
!>responsibility to anyone for the consequences of using them or for   %\n
!>whether they serve any particular purpose or work at all, unless he  %\n
!>says so in writing.  Everyone is granted permission to copy, modify  %\n
!>and redistribute the programs in GSLIB, but only under the condition %\n
!>that this notice and the above copyright notice remain intact.       %\n
 3    subroutine cov_data2vol(lindex,x1,y1,z1,ivol,vvcov)
c-----------------------------------------------------------------------
c     
c     Returns the covariance between to volume, ivol1 and ivol2
c     in case volume average data is present
c     *********************************************
c     
c     INPUT VARIABLES:
c  
c     lindex        index of data point considered
c     x1,y1,z1     location of data point
c     ivol         volume number
c     ivol2        number of volume 2 
c     
c     OUTPUT VARIABLES:
c     
c     vvcov       volume to volume covariance 
c     
c     ORIGINAL : Thomas Mejer nsen                       DATE: June 2004
c
c     TODO : Use either lookup table in RAM or a lookup table on disk 
c
c-----------------------------------------------------------------------
      use simul_arrays
        use mod_simul
      implicit none
      include 'visim.inc'
      double precision vvcov,ddcov, covsum
      integer :: i,j,k,ivol,lindex
      double precision x1,y1,z1
      double precision cov
      
      if (ivol.eq.0) then
         if (idbg.gt.0) then 
            write(*,*) 'Initializing data2vol covar lookup table' 
         endif
         k=0
         do i=1,(nx*ny*nz)
            do j=1,MAXVOLS
               cd2v(i,j)=UNEST
            enddo
         enddo
         vvcov=0
         return
      endif

c     UNCOMMENT NEXT LINE TO NOT USE LOOKUP TABLE
c      cd2v(lindex,ivol)=UNEST

      if (cd2v(lindex,ivol).eq.UNEST) then
c     CALCULATE THE VALUE

c
c HERE THE TESTING TAKES PLACE !!!!!!
c
         covsum=0
         do i=1,ndatainvol(ivol) 
            call cova3(volx(ivol,i),voly(ivol,i),volz(ivol,i),
     +           x1,y1,z1,1,nst,MAXNST,c0,it,
     +           cc,aa,1,MAXROT,rotmat,cmax,cov)
c
            covsum=voll(ivol,i)*dble(cov)+covsum
c
c            call cov_data2data(volx(ivol,i),voly(ivol,i),volz(ivol,i),
c     +           x1,y1,z1,ddcov)
c            covsum=voll(ivol,i)*dble(ddcov)+covsum
c            write(*,*) 'INIT ',cov,ddcov
         enddo
c
c
         vvcov=covsum
c     put the value inthe lookup table
c     comment this line out to disable the look table
c     in this case the cd2v variable should be removed from the visim.inc file.
        cd2v(lindex,ivol)=vvcov;
      else
         vvcov=cd2v(lindex,ivol);         
      end if
c
      return
      end
