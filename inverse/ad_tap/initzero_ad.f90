!>    @brief parallelisation wrapper for "omp_initzero_ad"
!>    @param[in] ismpl local sample index
      SUBROUTINE initzero_ad(ismpl)
        USE mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_initzero_ad(ismpl)
#ifdef fOMP
!$OMP end parallel
#endif

        RETURN
      END

!>    @brief to initialize all derived-variables with zero
!>    @param[in] ismpl local sample index
!>    @details
!>    Initialize all derived-variables with zero needed before seeding\n
      SUBROUTINE omp_initzero_ad(ismpl)
        use arrays
        USE arrays_ad
        USE mod_genrl
        USE mod_conc
        USE mod_time
        USE mod_data
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l, m
!     OpenMP stuff
        INTEGER tpos, tanz

!     OpenMP stuff
        CALL omp_part(i0*j0*k0,tpos,tanz)
        CALL ijk_m(tpos,i,j,k)
!
!       Set all Head, Temp, Pres, Conc, Epot to 0
        CALL set_dval(tanz,0.D0,head_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,temp_ad(i,j,k,ismpl))
        DO l = 1, ntrans
          CALL set_dval(tanz,0.D0,conc_ad(i,j,k,l,ismpl))
        END DO
        CALL set_dval(tanz,0.D0,pres_ad(i,j,k,ismpl))
!       Set all matrix diagonal elements to 0
        CALL set_dval(tanz,0.D0,a_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,b_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,c_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,d_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,e_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,f_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_ad(i,j,k,ismpl))
!       Set all left/right hand side to 0
        CALL set_dval(tanz,0.D0,w_ad(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,x_ad(i,j,k,ismpl))
!
!       Set all *old values to 0
        DO m = 1, ncgen
          CALL set_dval(tanz,0.D0,headold_ad(tpos,m,ismpl))
          CALL set_dval(tanz,0.D0,tempold_ad(tpos,m,ismpl))
          DO l = 1, ntrans
            CALL set_dval(tanz,0.D0,concold_ad(tpos,l,m,ismpl))
          END DO

          CALL set_dval(tanz,0.D0,presold_ad(tpos,m,ismpl))
        END DO
!
!$OMP master
        CALL set_dval(nunits*nprop,0.D0,propunit_ad(1,1,ismpl))
        CALL set_dval(ngsmax*3*nbctp,0.D0,bcperiod_ad(1,1,1,ismpl))
!        simtime_ad=0.0d0
        CALL set_dval(nsmpl,0.D0,simtime_ad(1))
        IF (ndata>0) THEN
          CALL set_dval(ndata,0.D0,sdata(1,ismpl))
          CALL set_dval(ndata,0.D0,sdata_ad(1,ismpl))
        END IF
!$OMP end master
!
        RETURN
      END SUBROUTINE
