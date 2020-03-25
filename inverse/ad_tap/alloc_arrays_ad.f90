!>    @brief allocate derivative objects (used for AD)
!>    @param[in] ismpl local sample index
      SUBROUTINE alloc_arrays_ad(ismpl)
        use arrays
        USE arrays_ad
        USE mod_genrl
        USE mod_data
        USE mod_time
        USE mod_conc
        IMPLICIT NONE
        integer :: ismpl
        INTRINSIC max


        write(*,*) "Initzero_ad called"

        ALLOCATE(a_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(b_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(c_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(d_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(e_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(f_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
!
        ALLOCATE(w_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(x_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
!
        ALLOCATE(head_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(temp_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(conc_ad(I0,J0,K0,max(ntrans,1),nsmpl))
        memory = memory + i0*j0*k0*max(ntrans,1)*nsmpl
        ALLOCATE(pres_ad(I0,J0,K0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(tsal_ad(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl


!
        ALLOCATE(headold_ad(I0*J0*K0,ncgen,nsmpl))
        memory = memory + i0*j0*k0*ncgen*nsmpl
        ALLOCATE(tempold_ad(I0*J0*K0,ncgen,nsmpl))
        memory = memory + i0*j0*k0*ncgen*nsmpl
        ALLOCATE(concold_ad(I0*J0*K0,max(ntrans,1),ncgen,nsmpl))
        memory = memory + i0*j0*k0*max(ntrans,1)*ncgen*nsmpl
        ALLOCATE(presold_ad(I0*J0*K0,ncgen,nsmpl))
        memory = memory + i0*j0*k0*ncgen*nsmpl
!
        ALLOCATE(propunit_ad(nunits,nprop,nsmpl))
        memory = memory + nunits*nprop*nsmpl
        ALLOCATE(bcperiod_ad(ngsmax,3,max(nbctp,1),nsmpl))
        memory = memory + ngsmax*3*max(nbctp,1)*nsmpl
        ALLOCATE(sdata_ad(max(ndata,1),nsmpl))
        memory = memory + max(ndata,1)*nsmpl
        ALLOCATE(dbc_data_ad(nbc_data,ndbc,nsmpl))
        memory = memory + nbc_data*ndbc*nsmpl
!
        allocate(simtime_ad(nsmpl))
        memory = memory +nsmpl
        allocate(adm_ad(i0,j0, k0, nsmpl))
        allocate(dms_ad(i0,j0, k0, nsmpl))
        allocate(ctgt_ad(i0,j0, k0, nsmpl))
        allocate(gravm_ad(i0,j0, k0, nsmpl))
        allocate(gram_ad(i0,j0, k0, nsmpl))


        RETURN
      END SUBROUTINE alloc_arrays_ad
