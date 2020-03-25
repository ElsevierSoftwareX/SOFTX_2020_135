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
      subroutine rayrandpath(dummyorder)
!AW moved "order" to "dummyorder", because is it a global/common defined array
!-----------------------------------------------------------------------
!  Descriptions:
!     We can define both random and sequential path for 1 ray using this subroutine
!     The random path are store in the array called 'order'.
!     The random value and random path are writen to a temporary output file for check. 
!
!  parameters:
!
!      densitypr : Density priority : Give higher priority (sample early) 
!                  to data points sensitive to larger volumes 
!                  [2] : order by number of volumes at data
!                  [1] : order by sum of density at data
!                  [0] : Dont use density priority
!
!      shuffvol: [1] randomly shuffle volumes
!                [0] use volumes in the order they are read
!      
!      shuffinvol : [0] sort by distance from source (AS READ)
!                   [1] shuffle within volume
!                   [2] shuffle withine all volumes 
!                       This means a) random point i any volume
!                             then b) random point outside volume
!                       This option overrides 'shuffvol'
!
!
! ORIGINAL: Yongshe Liu Thomas Mejer Hansen   DATE: June, 2004
!-----------------------------------------------------------------------
      use simul_arrays
        use mod_simul
      implicit none
      include 'visim.inc'
      integer :: ind,ix,iy,iz,j,k,nvp, MXYZ,i,idum,idata,ivol
      double precision vvx(MAXGEOMDATA), dummyorder
!AW
      integer, allocatable, dimension (:) :: vorder
      integer, allocatable, dimension (:) :: ivoll
      double precision, allocatable :: tempsim(:),simrest(:)
      double precision, allocatable :: svoll(:),nvoll(:)
      double precision, allocatable :: svoll2(:),nvoll2(:)
!
      integer  varr(MAXVOLS)
!AW-end
      double precision tempvol(MAXVOLS)
      double precision p, c,d,e,f,g,h, acorni
      external acorni
!      integer  shuffvol,shuffinvol,densitypr
      character (len=80) :: tmpfl
      logical testind

!
! these next variables COULD be set in visim.par file      
! but, since there is a clear benefit setting shiffinvol=2,
! this is chosen as default.
! The defaults are chosen here : 

      shuffvol=1;
      shuffinvol=2;
!      densitypr=2;
      
      if (idbg.gt.0) then
         write(*,*) 'Random Path : densitypr=',densitypr,
     +        '  shuffvol=',shuffvol,
     +        '  shuffinvol=',shuffinvol
      endif
      
      nxy=nx*ny;

!     Classic independant path
      if (densitypr.eq.0) then         
         p=acorni(idum)
         call sortem(1,nxyz,sim,1,order,c,d,e,f,g,h)
         return
      endif
!
!AW
      MXYZ=MAXX*MAXY*MAXZ
      allocate(tempsim(MXYZ))
      allocate(simrest(MXYZ))
      allocate(svoll(MXYZ))
      allocate(nvoll(MXYZ))
      allocate(vorder(MXYZ))
      allocate(ivoll(MXYZ))
      allocate(svoll2(MXYZ))
      allocate(nvoll2(MXYZ))
!AW-end

!     SORT VOLUMES IF NEEDED
      if (shuffvol.eq.1) then
         do ivol=1,nvol
            tempvol(ivol) = acorni(idum)
            varr(ivol) = ivol

         enddo
         call isortem(1,nvol,tempvol,1,varr,c,d,e,f,g,h)
      else 
!     ELSE DONT SORT VOLUMES
         do ivol=1,nvol
            varr(ivol)=ivol
         enddo
      endif

      if (idbg.gt.3) then 
         do ivol=1,nvol
            write(*,*) 'varr(',ivol,')=',varr(ivol)
         enddo		  
      end if

!     INITIALIZE THE SORT OF ALL THE POINTS
!     ASSIGN A RANDOM VALUE BETWEEN 0 and 1 TO ALL DATA
      do i=1,nxyz
            tempsim(i)   = acorni(idum)
            order(i) = i
      enddo

!
!     the nvoll2 and svoll2 are only initialized since the sortem function
!     alters the values of nvoll and svoll when called !
!     

!     APPLY DENSITY PRIORITY IF NEEDED
      nvp=1
      if (densitypr.gt.1) then
!         write(*,*) 'Density Prioirity'
         do ind=1,nxyz
            svoll(nvp)=0;
            nvoll(nvp)=0;
            svoll2(nvp)=0;
            nvoll2(nvp)=0;
            do ivol=1,nvol
               do idata=1,ndatainvol(ivol)
                  if (voli(ivol,idata).eq.ind) then
                     svoll(nvp)=svoll(nvp)+voll(ivol,idata);
                     nvoll(nvp)=nvoll(nvp)+1
                     nvoll2(nvp)=nvoll(nvp)
                     svoll2(nvp)=svoll(nvp)
                     ivoll(nvp)=ind
                     vorder(nvp)=nvp;
                  end if
               end do
            end do
!     ONLY CONSIDER THIS DATA OF MORE THAN ONE VOLUME GOES THROUGH IT
            if (nvoll(nvp).gt.1) then 
               if (idbg.gt.-13) then 
!                  write(*,*) 'nv(',ind,')=',svoll(nvp),nvoll(nvp),nvp
               end if
               nvp=nvp+1
            end if
         end do         
         nvp=nvp-1
 
!         do i=1,20
!            write(*,*) 'ivp=',i,' vorder=',vorder(i),svoll(i),nvoll(i)
!         end do
 
!     NOW SORT THE nvp DATA USING EITHER OF TWO CRITERIA
         if (densitypr.eq.2) then 
!     SORT BY SUM OF DENSITY AT POINT
            if (idbg.gt.0) write(*,*) 'SORT BY DENSITY'
            call isortem(1,nvp,svoll,1,vorder,c,d,e,f,g,h)
         else
!     SORT BY NUMBER VOLUME DATA POINT
            if (idbg.gt.0) write(*,*) 'SORT BY MVOLS THROUGH POINT'
            call isortem(1,nvp,nvoll,1,vorder,c,d,e,f,g,h)
         end if
         

!         do i=1,20
!            j=vorder(i)
!            write(*,*) 'i=',i,' j=',j,' ',svoll(j),svoll2(j)
!            write(*,*) 'i=',i,' j=',j,' ',nvoll(j),nvoll2(j)
!         end do


         do i=1,nvp
            if (densitypr.eq.2) then 
               tempsim(ivoll(vorder(i))) = tempsim(ivoll(vorder(i))) - 
     +              (nvol + dble(i)/10000 + svoll2(vorder(i)) )
            else
               tempsim(ivoll(vorder(i))) = tempsim(ivoll(vorder(i))) - 
     +              (nvol + dble(i)/10000 + nvoll2(vorder(i)) )
            end if 
!            write(*,*) 'ivp=',i,' vorder=',vorder(i),svoll2(i),nvoll2(i),
         end do
         


      end if
   
   
!     GET INDEX OF DATA IN VOLUME
      i=0
      do ivol=1,nvol
         do idata=1,ndatainvol(varr(ivol))
            i=i+1;
            call getindx(nx,xmn,xsiz,volx(varr(ivol),idata),ix,testind)
            call getindx(ny,ymn,ysiz,voly(varr(ivol),idata),iy,testind)
            call getindx(nz,zmn,zsiz,volz(varr(ivol),idata),iz,testind)
            ind = ix + (iy-1)*nx + (iz-1)*nxy
!     ONLY CHANGE THE TEMPSIM FOR THE INDEX IF NOT PREVIOUSLY SAMPLED
            if (tempsim(ind).gt.0) then              
               if (shuffinvol.eq.0) then
                  tempsim(ind)   = dble(idata)/10000 - (nvol - ivol +1)
               elseif (shuffinvol.eq.1) then
                  tempsim(ind)   = tempsim(ind) - (nvol - ivol +1 )
               elseif (shuffinvol.eq.2) then
                  tempsim(ind)   = tempsim(ind) - 1
               end if 
!               write(*,*) 'ind,tempsim : ',ind,tempsim(ind),nvol,ivol
            else
!     DO NOTHING
            endif
         end do
      end do


! SORT THE DATA 
      call sortem(1,nxyz,tempsim,1,order,c,d,e,f,g,h)
!
      deallocate(tempsim)
      deallocate(simrest)
      deallocate(svoll)
      deallocate(nvoll)
      deallocate(vorder)
      deallocate(ivoll)
      deallocate(svoll2)
      deallocate(nvoll2)
!
      return
      end
