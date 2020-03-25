!>    @brief free memory of derivative objects (used for AD)
!>    @param[in] ismpl local sample index
      SUBROUTINE dealloc_arrays_ad(ismpl)
        use arrays
        USE arrays_ad
        USE mod_genrl
        USE mod_data
        USE mod_time
        USE mod_conc
        IMPLICIT NONE
        integer :: ismpl
        INTRINSIC max

        DEALLOCATE(a_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(b_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(c_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(d_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(e_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(f_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(g_ad)
        memory = memory - i0*j0*k0*nsmpl
!
        DEALLOCATE(w_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(x_ad)
        memory = memory - i0*j0*k0*nsmpl
!
        DEALLOCATE(head_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(temp_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(conc_ad)
        memory = memory - i0*j0*k0*max(ntrans,1)*nsmpl
        DEALLOCATE(pres_ad)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(tsal_ad)
        memory = memory - i0*j0*k0*nsmpl

!
        DEALLOCATE(headold_ad)
        memory = memory - i0*j0*k0*ncgen*nsmpl
        DEALLOCATE(tempold_ad)
        memory = memory - i0*j0*k0*ncgen*nsmpl
        DEALLOCATE(concold_ad)
        memory = memory - i0*j0*k0*max(ntrans,1)*ncgen*nsmpl
        DEALLOCATE(presold_ad)
        memory = memory - i0*j0*k0*ncgen*nsmpl
!
        DEALLOCATE(propunit_ad)
        memory = memory - nunits*nprop*nsmpl
        DEALLOCATE(bcperiod_ad)
        memory = memory - ngsmax*3*max(nbctp,1)*nsmpl
        DEALLOCATE(sdata_ad)
        memory = memory - max(ndata,1)*nsmpl
        DEALLOCATE(dbc_data_ad)
        memory = memory - nbc_data*ndbc*nsmpl
!
        deallocate(adm_ad)
        deallocate(dms_ad)
        deallocate(ctgt_ad)
        deallocate(gravm_ad)
        deallocate(gram_ad)     

        RETURN
      END SUBROUTINE dealloc_arrays_ad
