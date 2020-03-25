!>    @brief declaration for all AD gradient objects
      MODULE g_arrays
        IMPLICIT NONE

!       hydraulic potential
        DOUBLE PRECISION, ALLOCATABLE :: g_head(:,:,:,:)
!       hydraulic pressure
        DOUBLE PRECISION, ALLOCATABLE :: g_pres(:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE :: g_conc(:,:,:,:,:)

!       temperature
        DOUBLE PRECISION, ALLOCATABLE :: g_temp(:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE :: g_tsal(:,:,:,:)

!       inverse
        DOUBLE PRECISION, ALLOCATABLE :: g_propunit(:,:,:)

!       coefficients for linear system solver
        DOUBLE PRECISION, ALLOCATABLE :: g_a(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_b(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_c(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_d(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_e(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_f(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_g(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_w(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_x(:,:,:,:)

!       storing old temp,head,conc for iteration
        DOUBLE PRECISION, ALLOCATABLE :: g_headold(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_tempold(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_presold(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_concold(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_prad(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_conc_conv(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_omp_dglobal(:,:,:)

!       boundary conditions
        DOUBLE PRECISION, ALLOCATABLE :: g_dbc_data(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_sdata(:,:)

        DOUBLE PRECISION, ALLOCATABLE :: g_diff_c(:)
        DOUBLE PRECISION, ALLOCATABLE :: g_mmas_c(:)

        DOUBLE PRECISION, ALLOCATABLE :: g_bcperiod(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: g_delta_time(:)
        DOUBLE PRECISION, ALLOCATABLE :: g_simtime(:)

        CONTAINS

!>    @brief dummy routine for TAF (content moved to "g_alloc_arrays")
        SUBROUTINE g_arrays_constructor
        RETURN
        END SUBROUTINE g_arrays_constructor

!>    @brief dummy routine for TAF (content moved to "g_dealloc_arrays")
        SUBROUTINE g_arrays_destructor
        RETURN
        END SUBROUTINE g_arrays_destructor

      END MODULE g_arrays
