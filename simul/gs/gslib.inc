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
!
! GS-Library Common Blocks:
!
      integer        ccmin,ccmax,cint(4),cscl
      common /color/ ccmin,ccmax,cint   ,cscl
!
      integer, parameter ::  KORDEI=12, MAXOP1=KORDEI +1
      integer ixv(MAXOP1), ixv2(MAXOP1)
      common /iaco/  ixv
      common /iaco2/ ixv2
!
      integer lpsout
      double precision pscl,pxmin,pxmax,pymin,pymax
      double precision xmin,xmax,ymin,ymax
      common /psdata/ pscl,pxmin,pxmax,pymin,pymax,xmin
      common /psdata/ xmax,ymin,ymax,lpsout
!
      integer itab(55),n1,n2,gs_nseed
      common /unusual/ itab,n1,n2,gs_nseed
!
!     instead of 'save'
      double precision      p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
      common /gslib_gauinv/ p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
      double precision      lwold
      common /gslib_psline/ lwold
      integer ifold, izip
      double precision tsold
      character (len=32), dimension (10) :: fnnt
      common /gslib_pstext/ fnnt,tsold,ifold,izip
!AW
!     for OpenMP parallelisation, thread-save
!$OMP threadprivate(/color/,/psdata/,/unusual/)
!$OMP threadprivate(/gslib_psline/,/gslib_gauinv/,/gslib_pstext/)
!$OMP threadprivate(/iaco/,/iaco2/)
!AW-end
