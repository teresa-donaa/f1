MODULE observations
!
! Collects the data
!
USE constants
!
IMPLICIT NONE
!
SAVE
!
! Data
!
INTEGER :: i_trial                                  ! Loop indexes
INTEGER :: seedrn(4)								! Starting seed for the random number generator
INTEGER :: seedini(4)								! Starting seed for the starting values
CHARACTER(len=3) :: numerical_error_flag			! Numerical error flag 
CHARACTER(len=1) :: qqq                             ! Reading string
!
! Data
!
INTEGER :: num_C(num_T)                             ! Number of observed options per day
INTEGER :: num_tau(num_T)                           ! Number of observed maturities per day
INTEGER :: date(num_T)                              ! Observation date
REAL(8) :: S(num_T)                                 ! Observed prices
REAL(8) :: RV(num_T)                                ! Observed realized volatility
REAL(8) :: r                                        ! Observed average risk free rate 
REAL(8) :: q                                        ! Observed average dividend flow rate 
REAL(8) :: rmq                                      ! Observed r-q
REAL(8) :: y(num_T)                                 ! Observed return on the underlying asset
                                                    ! NB: Available only in t >= 2
!
REAL(8) :: K(num_max_C,num_T)                       ! Observed strike prices with ME
INTEGER :: tau(num_max_C,num_T)                     ! Observed times to maturity with ME
INTEGER :: num_K(num_max_C,num_T)                   ! Observed options with the same tau
REAL(8) :: C(num_max_C,num_T)                       ! Observed option midprice with ME
REAL(8) :: X(num_max_C,num_T)					    ! Observed discounted moneyness with ME
REAL(8) :: H(num_max_C,num_T)					    ! Observed normalized option prices
!
! Loglikelihood evaluation
!
REAL(8) :: K_pdf_gamma(0:pdf_gamma_trunc)
!
! Option pricing
!
! ffo
!
COMPLEX(8) :: z_ffo(num_N)                  
REAL(8) :: vk_ffo(num_N,num_max_C)
COMPLEX(8) :: zmat_ffo(num_N,num_max_C)   
!
! bsp / lin
!
REAL(8), DIMENSION(num_grid_x) :: Xgrid
INTEGER, DIMENSION(num_grid_tau) :: taugrid
REAL(8), DIMENSION(num_grid_f) :: fgrid
REAL(8) :: xstep, fstep
INTEGER :: taustep
REAL(8), DIMENSION(num_N,num_grid_X) :: cosdx, sindx
!
! lin
!
INTEGER, DIMENSION(num_max_C,num_T) :: jx_mat, jtau_mat 
REAL(8), DIMENSION(num_max_C,num_T) :: dx_mat, dtau_mat
!
! Identity matrices
!
REAL(8), DIMENSION(num_theta,num_theta) :: eye_theta
REAL(8), DIMENSION(num_M,num_M) :: eye_M
!
! Ending module
!
END MODULE observations