MODULE constants
!
! Declare constants used throughout the project
!
IMPLICIT NONE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Choosing the task
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER, PARAMETER :: switch_timing = 1         ! Timing the code
INTEGER, PARAMETER :: num_timing = 10           ! Number of likelihood evaluations in timing the code
INTEGER, PARAMETER :: switch_stima1 = 1         ! Compute first round estimates
INTEGER, PARAMETER :: switch_stima2 = 0         ! Compute second round estimates 
INTEGER, PARAMETER :: switch_asvar = 0          ! Compute asymptotic variance matrix
!
! Number of threads to be used in openmp
!
INTEGER, PARAMETER :: num_threads = 20
!
! Switches to select "special" starting values
!
INTEGER, PARAMETER :: switch_readtheta = 1      ! Start from theta read from file
!
! Switches to show something on screen
!
INTEGER, PARAMETER :: switch_show_tt = 0        ! = 1: Show date index on screen
INTEGER, PARAMETER :: switch_show_tt_ld = 0     ! = 1: Show date index in loading data
!
! Switches to print something to file
!
INTEGER, PARAMETER :: switch_print_filter_control = 0
                                                ! = 1: Print output filter control file
INTEGER, PARAMETER :: switch_print_theta = 0    ! = 1: Print theta to file
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Model parameters
!
INTEGER, PARAMETER :: num_psi = 4               ! (lambda, nu, beta, mu)
INTEGER, PARAMETER :: num_phi = 1               ! phi = 1/(1 - ystar*mu)
INTEGER, PARAMETER :: num_eta = 1               ! Total number of parameters in expected RV 
INTEGER, PARAMETER :: num_chi = 1               ! Total number of parameters in sigma_RV
INTEGER, PARAMETER :: num_xi = 6                ! Total number of parameters in sigma_H
INTEGER, PARAMETER :: num_theta = num_psi+num_phi+num_eta+num_chi+num_xi
                                                ! Total number of parameters
!
! True parameters 
!
INTEGER, PARAMETER :: num_param = &
    5+ &                                        ! (lambdaP, nuP, betaP, muP, rhoP)
    5+ &                                        ! (lambdaQ, nuQ, betaQ, muQ, rhoQ)
    4+ &                                        ! phi, ystar, delta1, delta2
    2+ &                                        ! Parameters in expected RV 
    2+ &                                        ! Parameters in sigma_RV
    6                                           ! Parameters in sigma_H
!
! Bounds on initial parameters generation
!
REAL(8), PARAMETER :: buini = 1.d0
!
! Number of estimation trials
!
INTEGER, PARAMETER :: num_trial = 100*(switch_stima1)+1*switch_stima2
                                                ! Total number of estimation trials
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Parameters about the observed sample
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Parameters about the total sample size
!
CHARACTER(len=30), PARAMETER :: file_data = 'optsp500_crit1.txt'    
                                                ! Name of the data file containing sample observations
INTEGER, PARAMETER :: num_obs_opt = 80979       ! Total number of options used for estimation
INTEGER, PARAMETER :: num_T = 2517              ! Total number of daily observationsed dates *** MINUS 1 ***
!
! Switches to select the options sample
!
INTEGER, PARAMETER :: min_volume = 5            ! Minimum admissible transaction volume
INTEGER, PARAMETER :: min_tau = 15              ! Minimum admissible time to maturity
INTEGER, PARAMETER :: max_tau = 365             ! Maximum admissible time to maturity
REAL(8), PARAMETER :: min_x = -1.d-1            ! Minimum admissible moneyness (in terms of discounted moneyness)
REAL(8), PARAMETER :: max_x =  1.d-1            ! Maximum admissible moneyness (in terms of discounted moneyness)
!
! General sample features
!        
INTEGER, PARAMETER :: num_max_C = 75            ! Maximum number of options per day used in estimation
INTEGER, PARAMETER :: num_max_C_ds = 150        ! Maximum number of options per day in dataset
INTEGER, PARAMETER :: num_max_C_Ktau = 35       ! Maximum number of strike prices per day/maturity
INTEGER, PARAMETER :: num_open_days = 365       ! Agerage number of open market days
REAL(8), PARAMETER :: P_norm = 1.d2             ! Log(S) normalizing constant
REAL(8), PARAMETER :: V_norm = P_norm**2        ! Volatility normalizing constant
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Parameters about loglikelihood evaluation
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER, PARAMETER :: pdf_gamma_trunc = 100     ! Truncation order for the noncentral Gamma pdf
REAL(8), PARAMETER :: bulletproof = EXP(-5.d2)  ! see Tauchen
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Parameters about discretization
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER, PARAMETER :: num_M = 8                 ! Size of discretization grid
INTEGER, PARAMETER :: max_num_L = 2             ! Maximum number of conditional moments to match
INTEGER, PARAMETER :: num_L = 2                 ! Number of conditional moments to match
INTEGER, PARAMETER :: X_method = 1              ! = 1: quantile based on the Gamma marginal pdf
INTEGER, PARAMETER :: Q_method = 4              ! = 1: initial Q probabilities set using the modified gamma pdf
                                                ! = 2: initial Q probabilities set using the true gamma pdf
                                                ! = 3: initial Q probabilities set to 1/M
                                                ! = 4: initial Q probabilities set using the true gamma cdf
!
INTEGER, PARAMETER :: lambda_method = 1         ! = 1: L-BFGS; = 2: Newton-Raphson
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Parameters used the numerical option pricing method
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Parameters used in Fang and Oosterlee's COS method
!
INTEGER, PARAMETER :: num_N = 2000              ! Number of Fourier-cosine series approximation points
REAL(8), PARAMETER :: a_fs = -4.d0              ! Lower bound of Fourier cosine series integration
REAL(8), PARAMETER :: b_fs = 4.d0               ! Upper bound of Fourier cosine series integration
!
! Parameters used in computing option prices on the grid (splines or linear interpolation)
!
INTEGER, PARAMETER :: num_grid_X = 32
INTEGER, PARAMETER :: num_grid_tau = 96
!INTEGER, PARAMETER :: num_grid_tau = 64*switch_stima1 &
!                            + 96*(switch_stima2+switch_asvar) 
INTEGER, PARAMETER :: num_grid_f = 48
!INTEGER, PARAMETER :: num_grid_f = 16*switch_stima1 &
!                            + 40*(switch_stima2+switch_asvar)
!
INTEGER, PARAMETER :: switch_log_fgrid = 1      ! = 1: fgrid in logs; 0 = fgrid in levels
INTEGER, PARAMETER :: switch_min_f_level = 0    ! = 1: min(fgrid) fixed by min_f; = 0: fixed in percent
REAL(8), PARAMETER :: min_f = 0.d0              ! Smallest value in fgrid
REAL(8), PARAMETER :: min_f_frac = 1.d-1        ! Coefficient to multiply MINVAL(LOG(RV)) to define first point in fgrid
REAL(8), PARAMETER :: max_f = 2.d-3             ! Largest value in fgrid
REAL(8), PARAMETER :: min_l = 0.d0              ! Smallest value in lgrid
REAL(8), PARAMETER :: max_l = 1.d1              ! Largest value in lgrid (l approx chi2(1))
!
! Parameters used in linear interpolation
!
INTEGER, PARAMETER :: Lmat(8,3) = RESHAPE( (/ & ! NOTE: LMAT is written transposed!
    0, 0, 0, 0, 1, 1, 1, 1, &
    0, 0, 1, 1, 0, 0, 1, 1, &
    0, 1, 0, 1, 0, 1, 0, 1  &
    /), (/ 8, 3 /) )
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters for POLITOPE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
INTEGER, PARAMETER :: max_iters_politope = (25*num_theta)*switch_stima1+(100*num_theta)*switch_stima2
                                                ! Maximum number of iterations in politope
REAL(8), PARAMETER :: pert_theta(num_theta) = 0.5d0     ! Relative perturbation of the parameters
INTEGER, PARAMETER :: max_repeat_politope = num_theta*switch_stima1+8*num_theta*switch_stima2
                                                ! Maximum number of politope restarts in second round optimization
REAL(8), PARAMETER :: tol = 1.d-7*switch_stima1+1.d-9*switch_stima2
                                                ! Tolerance in second round optimization
REAL(8), PARAMETER :: tol_conv = 1.d-6*switch_stima1+1.d-8*switch_stima2
                                                ! Convergence tolerance in second round optimization
REAL(8), PARAMETER :: tol_politope_p = tol
REAL(8), PARAMETER :: crit_politope_conv_p = tol_conv
REAL(8), PARAMETER :: tol_politope_y = tol
REAL(8), PARAMETER :: crit_politope_conv_y = tol_conv
INTEGER, PARAMETER :: rtol_formula = 2          ! Chooses the formula used to compute rtol
                                                ! See lines 150-190 in simplex_M.f90
INTEGER, PARAMETER :: crit_conv_formula = 1     ! Politope minimizations are restarted looking 
                                                ! at improvements in: 
                                                ! 1 = y
                                                ! 2 = p
REAL(8), PARAMETER :: loglik_penalty = 1.d25                                                
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters about the output files
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
CHARACTER(len=30), PARAMETER :: file_res_stima = 'res_stima.txt'
                                                ! Name of the stima results file
CHARACTER(len=30), PARAMETER :: file_res_stima1 = 'res_stima1.txt'
CHARACTER(len=30), PARAMETER :: file_res_stima2 = 'res_stima2.txt'
CHARACTER(len=30), PARAMETER :: file_ll_params = 'll_params.txt'
                                                ! Name of the ll_params results file
CHARACTER(len=30), PARAMETER :: file_politope = 'politope.txt'
                                                ! Name of the politope output file
CHARACTER(len=30), PARAMETER :: file_Imat = 'Imat.txt'
                                                ! Name of the Imat output file
CHARACTER(len=30), PARAMETER :: file_Jmat = 'Jmat.txt'
                                                ! Name of the Jmat output file
CHARACTER(len=30), PARAMETER :: file_invJmat = 'invJmat.txt'
                                                ! Name of the invJmat output file
CHARACTER(len=30), PARAMETER :: file_dl = 'dl.txt'
                                                ! Name of the dl output file
CHARACTER(len=30), PARAMETER :: file_dl2 = 'dl2.txt'
                                                ! Name of the dl2 output file
CHARACTER(len=30), PARAMETER :: file_finres = 'finres.txt'
                                                ! Name of the finres output file
CHARACTER(len=30), PARAMETER :: file_filter_control = 'filter_control.txt'
                                                ! Name of the filter control file
CHARACTER(len=30), PARAMETER :: file_discretization_control = 'filter_discretization.txt'
                                                ! Name of the filter control file
CHARACTER(len=30), PARAMETER :: file_theta = 'theta_M16.txt'
                                                ! Name of the theta file
CHARACTER(len=30), PARAMETER :: file_dparam_dtheta = 'dparam_dtheta.txt'
                                                ! Name of the dparam_dtheta file
CHARACTER(len=30), PARAMETER :: file_param = 'param.txt'
                                                ! Name of the param file
CHARACTER(len=30), PARAMETER :: file_var_param = 'var_param.txt'
                                                ! Name of the var_param file
!
INTEGER, PARAMETER :: unit_data = 13            ! Number of unit of the data file containing sample observations
INTEGER, PARAMETER :: unit_data_sim = 131       ! Number of unit of the simulated data file 
INTEGER, PARAMETER :: unit_res_stima = 22       ! Number of the stima results file unit
INTEGER, PARAMETER :: unit_res_stima1 = 221     ! Number of the stima1 results file unit
INTEGER, PARAMETER :: unit_res_stima2 = 222     ! Number of the stima2 results file unit
INTEGER, PARAMETER :: unit_ll_params = 32       ! Number of the ll_params results file unit
INTEGER, PARAMETER :: unit_politope = 5         ! Number of politope  output file unit
INTEGER, PARAMETER :: unit_Imat = 53            ! Number of Imat output file unit
INTEGER, PARAMETER :: unit_Jmat = 531           ! Number of Jmat output file unit
INTEGER, PARAMETER :: unit_invJmat = 54         ! Number of invJmat output file unit
INTEGER, PARAMETER :: unit_dl = 55              ! Number of dl output file unit
INTEGER, PARAMETER :: unit_dl2 = 56             ! Number of dl2 output file unit
INTEGER, PARAMETER :: unit_finres = 61          ! Number of finres output file unit
INTEGER, PARAMETER :: unit_filter_control = 94  ! Number of filter control file unit
INTEGER, PARAMETER :: unit_discretization_control = 95  
                                                ! Number of filter control file unit
INTEGER, PARAMETER :: unit_theta = 98           ! Number of theta file unit
INTEGER, PARAMETER :: unit_dparam_dtheta = 71   ! Number of dparam_dtheta file unit
INTEGER, PARAMETER :: unit_param = 72           ! Number of param file unit
INTEGER, PARAMETER :: unit_var_param = 73       ! Number of var_param file unit
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring constants
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
REAL(8), PARAMETER :: greek_pi = 3.14159265358979323846264338328d0          ! Greek pi
REAL(8), PARAMETER :: invsqrt2pi = 0.398942280401432677939946059934d0       ! 1/SQRT(2*pi)
REAL(8), PARAMETER :: dueinvpi = 0.636619772367581343075535053490d0         ! 2/Pi
REAL(8), PARAMETER :: log1oversqrt2pi = -0.918938533204672741780329736406d0 ! log(1/sqrt(2*pi))
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE constants
