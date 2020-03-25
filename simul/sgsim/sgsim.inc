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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
! The following Parameters control static dimensioning for SGSIM:
!
!   MAXX      maximum nodes in X
!   MAXY      maximum nodes in Y
!   MAXZ      maximum nodes in Z
!   MAXCTX    maximum X points in covariance table (ODD number)
!   MAXCTY    maximum Y points in covariance table (ODD number)
!   MAXCTZ    maximum Z points in covariance table (ODD number)
!   MAXSBX    maximum super block nodes in X direction
!   MAXSBY    maximum super block nodes in Y direction
!   MAXSBZ    maximum super block nodes in Z direction
!   MAXDAT    maximum number of data
!   MAXNOD    maximum number of previously simulated nodes to use
!   MAXSAM    maximum number of data for one simulation
!   MAXNST    maximum number of nested structures for variogram
!
!-----------------------------------------------------------------------
!
! User Adjustable Parameters:
!
      integer MAXCTX, MAXCTY, MAXCTZ
      integer MAXSBX, MAXSBY, MAXSBZ
      integer MAXDAT, MAXNOD, MAXSAM
      integer MAXNST
      parameter(MAXCTX =    71, MAXCTY =  71, MAXCTZ = 21)
      parameter(MAXSBX =    21, MAXSBY =  21, MAXSBZ = 11)
      parameter(MAXDAT = 10000, MAXNOD =  48, MAXSAM = 48)
      parameter(MAXNST =     4)
!
! Fixed Parameters:
!
      integer MAXKR1,MAXKR2
      integer MAXROT,MAXCXY
      integer MAXXYZ,MAXSB
      double precision UNEST,EPSLON,VERSION
      parameter(MAXKR1=MAXNOD+MAXSAM+1,MAXKR2=MAXKR1*MAXKR1)
      parameter(MAXROT=MAXNST+1,MAXCXY=MAXCTX*MAXCTY)
      parameter(MAXXYZ=MAXCTX*MAXCTY*MAXCTZ,MAXSB=MAXSBX*MAXSBY*MAXSBZ)
      parameter(UNEST=-99.0,EPSLON=1.0e-20)
      parameter(VERSION=2.000)
!
! ACORN parameters:
!
      integer KORDEI,MAXOP1,MAXINT
      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
!
! Variable Declaration:
!
      double precision x(MAXDAT),y(MAXDAT),z(MAXDAT),vr(MAXDAT)
      double precision vrtr(MAXDAT),vrgtr(MAXDAT),dclose(MAXDAT)
      double precision sec(MAXDAT),wt(MAXDAT)
      double precision c0(1),cc(MAXNST),aa(MAXNST),ang1(MAXNST)
      double precision ang3(MAXNST),anis1(MAXNST),anis2(MAXNST)
      double precision covtab(MAXCTX,MAXCTY,MAXCTZ),ang2(MAXNST)
      double precision cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      double precision cnodev(MAXNOD),vra(MAXKR1),vrea(MAXKR1),ltpar
      double precision rotmat(MAXROT,3,3),r(MAXKR1),rr(MAXKR1)
      double precision s(MAXKR1),a(MAXKR2)
      integer ltail,utail,sstrat,mults,nisb(MAXSB),nst(1)
      integer it(MAXNST),icnode(MAXNOD)
      integer ixnode(MAXXYZ),iynode(MAXXYZ),iznode(MAXXYZ)
      integer ixsbtosr(8*MAXSB),iysbtosr(8*MAXSB),izsbtosr(8*MAXSB)
!
! Common Blocks:
!
      common /iaco/   ixv(MAXOP1)

      common /generl/ xsiz,ysiz,zsiz,xmn,ymn,zmn,x,y,z,vr,wt,itrans
      common /generl/ ntr,vrtr,vrgtr,sec,idbg,lin,lout,ldbg,llvm
      common /generl/ icollvm,nvaril,ktype,colocorr,zmin,zmax,varred
      common /generl/ ltpar,utpar,nx,ny,nz,nxy,nxyz,nd,ltail,utail

      common /simula/ nsim

      common /search/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2
      common /search/ noct,nclose,dclose,ndmin,ndmax,sstrat,mults,nmult

      common /cova3d/ cmax,c0,cc,aa,ang1,ang2,ang3
      common /cova3d/ anis1,anis2,isrot,nst,it

      common /clooku/ covtab,nlooku,ncnode,icnode
      common /clooku/ cnodex,cnodey,cnodez,cnodev,nodmax,ixnode,iynode
      common /clooku/ iznode,nctx,ncty,nctz

      common /krigev/ rotmat,r,rr,s,a,vra,cbb

!aw
!$OMP threadprivate(/iaco/,/generl/,/simula/,/search/,/cova3d/)
!$OMP threadprivate(/clooku/,/krigev/)
!aw
