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
!   MAXVOLS   maximum number of volumes type data to include 
!   MAXDINVOL MAXimum number of Data IN a VOLumes 
!
!-----------------------------------------------------------------------
!
! User Adjustable Parameters:
!
      integer MAXCTX, MAXCTY, MAXCTZ
      integer MAXSBX, MAXSBY, MAXSBZ, MAXDAT, MAXNOD, MAXSAM
      integer MAXNST, MAXREF, MAXCAT, MV, MAXDINVOL

      parameter(MAXCTX =    71, MAXCTY =  71, MAXCTZ = 21)
      parameter(MAXSBX =    21, MAXSBY =  21, MAXSBZ = 11)
      parameter(MAXDAT = 50000, MAXNOD =  1448, MAXSAM = 148)
      parameter(MAXNST =     4, MAXREF=10000,MAXCAT=24)
      parameter(MV=20         , MAXDINVOL=199)
!
! Fixed Parameters:
!
      integer MAXKR1, MAXKR2, MAXROT, MAXCXY, MAXXYZ, MAXSB
      integer MAXGEOMDATA
      double precision VERSION, UNEST, EPSLON
      parameter(MAXKR1=MAXNOD+MAXSAM+1,MAXKR2=MAXKR1*MAXKR1)
      parameter(MAXROT=MAXNST+1,MAXCXY=MAXCTX*MAXCTY)
      parameter(MAXXYZ=MAXCTX*MAXCTY*MAXCTZ,MAXSB=MAXSBX*MAXSBY*MAXSBZ)
      parameter(UNEST=-99999999.0d0,EPSLON=1.0d-20)
      parameter(MAXGEOMDATA=MAXVOLS*MAXDINVOL)
      parameter(VERSION=1.2008d0)
!
! ACORN parameters:
!
      integer MAXINT
      parameter(MAXINT=2**30)
!
! Variable Declaration:
!
      double precision dclose(MAXDAT)
      double precision x(MAXDAT),y(MAXDAT),z(MAXDAT),vr(MAXDAT)
      double precision vrtr(MAXDAT),vrgtr(MAXDAT),sec(MAXDAT)
      double precision c0(1),cc(MAXNST),aa(MAXNST)
      double precision ang1(MAXNST),ang2(MAXNST)
      double precision ang3(MAXNST),anis1(MAXNST),anis2(MAXNST)
      double precision wt(MAXDAT)
      double precision covtab(MAXCTX,MAXCTY,MAXCTZ)
      double precision datacov(MAXVOLS,MAXVOLS)
      double precision cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD)
      double precision cnodev(MAXNOD),vra(MAXKR1),vrea(MAXKR1),ltpar
      double precision utpar, pkr
      integer   cnodeindex(MAXNOD)
      double precision rcdf(MAXREF),rvr(MAXREF),dcdf(MAXDAT)
      double precision indx(MAXDAT),fuzzcat(MAXDAT),catcdf(MAXCAT)

      double precision rotmat(MAXROT,3,3),r(MAXKR1),rr(MAXKR1)
      double precision s(MAXKR1),a(MAXKR2)
      integer   ltail,utail,sstrat,mults,nisb(MAXSB),nst(1)
      integer   it(MAXNST),icnode(MAXNOD)
      integer   category(MAXDAT), ncut, inovar
      integer   ixnode(MAXXYZ),iynode(MAXXYZ),iznode(MAXXYZ)
      integer   ixsbtosr(8*MAXSB),iysbtosr(8*MAXSB),izsbtosr(8*MAXSB)
      double precision      wtt, vrt, wtd, vrd
      double precision	ckurt, colocorr, zmaxvar
      integer   wx, wy, wz,  isim, idrawopt   
      character (len=40) :: distin
      character (len=40) :: localfl
      character (len=40) :: outfl
      character (len=40) :: transoutfl      
      character (len=40) :: lambdafl

      integer   ibt,ibtw,nbt
      double precision	bootvar(MAXDAT),bootwt(MAXDAT),btmean,btvar
      double precision	bootcdf(MAXDAT)	
! FOR VOLUMES
      integer   nvol,ndatainvol(MAXVOLS)
      double precision	volobs(MAXVOLS)
      double precision	volvar(MAXVOLS)
      double precision	volobs_ref(MAXVOLS)  
      double precision	volx(MAXVOLS,MAXDINVOL)	
      double precision	voly(MAXVOLS,MAXDINVOL)	
      double precision	volz(MAXVOLS,MAXDINVOL)	
      double precision	voll(MAXVOLS,MAXDINVOL)	
      integer	voli(MAXVOLS,MAXDINVOL)	
      integer   nusev,usev(MAXVOLS)	
      integer   musevols,nusevols
      double precision 	accept_fract
      integer   shuffvol,shuffinvol,densitypr
      integer   lout_mean,lout_std,lout_krig
	
! FOR Histrogram Reproduction ..
      double precision 	condlookup_mean(400,400)
      double precision 	condlookup_var(400,400)
      double precision 	condlookup_cpdf(400,400,300)
      integer n_Gmean,n_Gvar,n_q,n_monte
      double precision min_Gmean,max_Gmean
      double precision min_Gvar,max_Gvar

      double precision x_quan(500)

! 	estimation/simulation
      integer doestimation

!
! The matrices below are covariance lookup tables. 
! they can optionally be excluded
!
      double precision  cv2v(MAXVOLS,MAXVOLS)
      integer read_covtable, read_lambda
      integer n_mean, n_var
      integer lin, ldbg, llvm, lkv, lout, nctx,ncty,nctz, nxyz, nd
      integer nsim, ntr, nvaril, nx,ny,nz, nxy
      double precision radius, sang1,sang2,sang3,sanis1,sanis2,rn,radsqd
      double precision tmin, tmax, skgmean, xmn,ymn,zmn,varred,zmin,zmax
      double precision wtfac, xsiz,ysiz,zsiz
      double precision gvar, cbb, cmax
      integer noct, icoll, icollvm, idbg, icond, isrot
      integer ivar, itr, itrans, iwtr, ivrd, ivrr, iwtd, ktype
      integer nclose, ndmin, ndmax, nmult, nlooku, ncnode, nodmax

!
! Common Blocks:
!
      common /kurto/  ckurt,pkr,idrawopt

      common /generl/ xsiz,ysiz,zsiz,xmn,ymn,zmn,nxy,nxyz, &
      	              x,y,z,vr,wt,vrtr,vrgtr,sec,idbg,ivar, &
	              lin,lout,ldbg,llvm,lkv,icollvm,nvaril,ktype, &
                      ltail,utail,lout_mean,lout_std,lout_krig, &
                      musevols,nusevols,nusev,usev,icond,ltpar,gvar, &
                      shuffvol,shuffinvol,zmin,zmax, &
                      accept_fract,varred,skgmean,utpar,colocorr, &
                      nx,ny,nz,nd,ntr,densitypr

      common /volume/ volx,voly,volz,volobs, &
      	volvar,volobs_ref,nvol,ndatainvol,voll,voli

      common /simula/ doestimation,isim, nsim

      common /search/ noct,nclose,ndmin,ndmax,dclose,sstrat,mults, &
                      radius,radsqd,sang1,sang2,sang3,sanis1,sanis2,nmult

      common /cova3d/ cmax,c0,cc,aa,ang1,ang2,ang3, &
                      anis1,anis2,cv2v,read_covtable, &
     		      read_lambda,nst,it,isrot

      common /cd/     datacov


      common /clooku/ covtab,cnodex,nodmax,ixnode,iynode,iznode, &
                      cnodey,cnodez,cnodev,cnodeindex, &
                      nctx,ncty,nctz,nlooku,ncnode,icnode

      common /krigev/ rotmat,r,rr,s,a,vra,cbb, &
                      zmaxvar,lambdafl,inovar
      common /transcon/  ivrr, iwtr, wx,wy, wz, &
     		      icoll,wtfac,rn,distin, &
     		      localfl,outfl,tmin,tmax,transoutfl, &
                      rcdf,rvr,dcdf,indx,fuzzcat, &
                      catcdf,category,wtt,vrt,wtd,vrd,ncut, &
		      ivrd,iwtd,itrans

      common /bt/   bootvar,bootwt,btmean,btvar,bootcdf,ibt,ibtw,nbt
	
      common /hrdssim/   condlookup_mean,condlookup_var, &
                    condlookup_cpdf,n_mean,n_var,n_q,n_monte, &
      	            x_quan, &
                    min_Gmean,max_Gmean, &
                    min_Gvar,max_Gvar,n_Gvar,n_Gmean

!AW
!     for OpenMP parallelisation, thread-save
!$OMP threadprivate(/kurto/,/generl/,/volume/,/simula/,/search/)
!$OMP threadprivate(/cd/,/bt/,/clooku/,/krigev/)
!$OMP threadprivate(/transcon/,/hrdssim/,/cova3d/)
!AW-end
