!>    @brief check changes between vectors [g_new] and [g_old]
!>    @param[in] mode switch absolute/relative
!>    @param[out] rms return value
!>    @param[out] difmax maximal difference
!>    @param[out] g_difmax (ignored) maximal difference
!>    @param[in] ni I-dimension for vectors [g_new], [g_old]
!>    @param[in] nj J-dimension for vectors [g_new], [g_old]
!>    @param[in] nk K-dimension for vectors [g_new], [g_old]
!>    @param[in] new (ignored) vector with new values
!>    @param[in] g_new vector with new derivative values
!>    @param[in] old (ignored) vector with old values
!>    @param[in] g_old vector with old derivative values
!>    @param[in] pv_idx index number (physical value), only needed for AD code generation/modification
!>    @param[out] loc_nltol tolerance criteria, only needed for AD code generation/modification
!>    @param[in] ismpl local sample index
!>    @details
!>      see subroutine "check_change"\n
      SUBROUTINE g_check_change( mode, pv_idx, loc_nltol, rms, difmax, g_difmax, ni, nj, nk, new, g_new, old, g_old, ismpl )
        USE arrays
        IMPLICIT NONE
        DOUBLE PRECISION difmax
        DOUBLE PRECISION g_difmax
        INTEGER ni, nj, nk, mode, pv_idx, ismpl
        DOUBLE PRECISION g_new(ni,nj,nk)
        DOUBLE PRECISION g_old(ni,nj,nk)
        DOUBLE PRECISION new(ni,nj,nk)
        DOUBLE PRECISION old(ni,nj,nk)
        DOUBLE PRECISION rms, loc_nltol

!
!       exchange original with derivative values
        CALL check_change( mode, pv_idx, loc_nltol, rms, difmax, ni, nj, nk, g_new, g_old, ismpl )
!       modify relative tolerance criteria
        CALL s_damax( ni*nj*nk, g_new, loc_nltol )
        loc_nltol = max( loc_nltol*nltol_g(pv_idx), const_dble(2) )
!
        RETURN
      END SUBROUTINE g_check_change
