!>    @brief free memory of derivative objects
!>    @param[in] ismpl local sample index
      SUBROUTINE g_dealloc_arrays(ismpl)
        use arrays
        use g_arrays
        use mod_genrl
        use mod_genrlc
        use mod_inverse
        use mod_data
        use mod_conc
        use mod_time
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INTRINSIC max

        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... g_dealloc_arrays'
! deallocating g_prad
        DEALLOCATE(g_prad)
        memory = memory-i0*j0*k0*nsmpl
! deallocating g_conc_conv - Temporary conversions list
        DEALLOCATE(g_conc_conv)
        memory = memory-ntrac*nsmpl
! deallocating g_omp_dglobal
        DEALLOCATE(g_omp_dglobal)
        memory = memory-tlevel_1*9*nsmpl

! hydraulic potential,pressure
        DEALLOCATE(g_head)
        memory = memory - i0*j0*k0*nsmpl
! temperature
        DEALLOCATE(g_temp)
        memory = memory - i0*j0*k0*nsmpl

        DEALLOCATE(g_pres)
        memory = memory - i0*j0*k0*nsmpl

        DEALLOCATE(g_conc)
        memory = memory - i0*j0*k0*max(ntrans,1)*nsmpl

        DEALLOCATE(g_tsal)
        memory = memory - i0*j0*k0*nsmpl

!     inversion
        DEALLOCATE(g_propunit)
        memory = memory - nunits*nprop*nsmpl

        DEALLOCATE(g_a)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_b)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_c)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_d)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_e)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_f)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_g)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_w)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_x)
        memory = memory - i0*j0*k0*nsmpl

        DEALLOCATE(g_headold)
        DEALLOCATE(g_tempold)
        DEALLOCATE(g_presold)
        DEALLOCATE(g_concold)
        memory = memory - i0*j0*k0*max(ntrans,1)*4*nsmpl

!     boundary conditions
!      value
        DEALLOCATE(g_dbc_data)
        memory = memory - nbc_data*ndbc*nsmpl
        DEALLOCATE(g_sdata)
        memory = memory - max(ndata,1)*nsmpl

        DEALLOCATE(g_diff_c)
        memory = memory - max(ntrans,1)
        DEALLOCATE(g_mmas_c)
        memory = memory - max(ntrans,1)

        DEALLOCATE(g_bcperiod)
        memory = memory - ngsmax*3*max(nbctp,1)*nsmpl
        DEALLOCATE(g_delta_time)
        memory = memory - ntimestep
        DEALLOCATE(g_simtime)
        memory = memory - nsmpl

        write(*,*) "Cleared allocated memory for AD memory usage is now ",memory*8, " Byte"
        RETURN
      END
