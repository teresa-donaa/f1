PROGRAM main
!
USE constants
USE observations
USE load_data
USE io_stuff
USE option_prices
USE starting_points
USE estimates
!USE asymptotic_variance
!
IMPLICIT NONE
!
! Declaring local variables
!
REAL(8) :: theta(num_theta), llcheck, time0, time1, objf, ll, chiM(num_M), PthetaM(num_M,num_M)
CHARACTER(len=4) :: ichar
INTEGER :: seed(2), i
REAL(8), DIMENSION(num_M,0:num_T) :: zetap, zetau
INTEGER(8) :: tclock1, tclock2, clock_rate
REAL(8) :: elapsed_time
!
! Beginning execution
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Setup
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Reading data
!
CALL load_obs(unit_data,file_data)
!
! Picking the setup routine to compute normalized option prices
!
CALL option_prices_bsp_setup()
CALL option_prices_lin_setup()
CALL option_prices_ffo_setup(num_N,a_fs,b_fs)
!
! Setting the number of threads
!
!$ CALL OMP_SET_NUM_THREADS(num_threads)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Estimation: timing, stima1 or stima2
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
IF ((switch_timing .EQ. 1) .OR. (switch_stima1 .EQ. 1) .OR. (switch_stima2 .EQ. 1)) THEN
    !
    ! Initialising random number generator
    !
    IF ((switch_stima1 .EQ. 1) .OR. (switch_timing .EQ. 1)) THEN
        !
        CALL GETARG(1,ichar)                ! seed for initial values generation
        seed = INUM(ichar) 
        IF (ALL(seed .LE. 0)) seed = num_threads
        !
    END IF
    !
    ! Opening control and files, if required
    !
    IF (switch_print_filter_control .EQ. 1)  THEN
        !
        CALL open_filter_control()
        CALL open_discretization_control()
        !
    END IF
    IF (switch_print_theta .EQ. 1) CALL open_write_file(unit_theta,file_theta)
    !
    ! Opening input and output files
    !
    CALL open_write_file(unit_res_stima,file_res_stima)
    CALL open_write_file(unit_ll_params,file_ll_params)
    CALL head_ll_params()
    !
    ! Beginning loop over number of trials
    !
    i_trial = 0
    DO 
        !
        i_trial = i_trial+1
        IF (i_trial .GT. num_trial) EXIT
        WRITE(*,1976) i_trial, num_trial
        1976 FORMAT('trial: ', I5, '/', I5)
        !
        ! Initialising starting points
        !
        CALL admissible_starting_point(seed,theta,llcheck)
        !
        ! Timing segment
        !
        IF (switch_timing .EQ. 1) THEN
            !
            CALL SYSTEM_CLOCK(tclock1)
            PRINT*, 'time0 = ', time0
            DO i = 1, num_timing
                WRITE(*, 3) i
3               FORMAT('i    = ', I3)
                objf = loglik_politope(theta)
                WRITE(*, 4) objf
4               FORMAT('objf = ', ES16.9)
                CALL SYSTEM_CLOCK(tclock2,clock_rate)
                elapsed_time = FLOAT(tclock2-tclock1)/FLOAT(clock_rate)
                WRITE(*,5) elapsed_time/i
            END DO
            CALL SYSTEM_CLOCK(tclock2,clock_rate)
            elapsed_time = FLOAT(tclock2-tclock1)/FLOAT(clock_rate)
            WRITE(*,5) elapsed_time/num_timing
5           FORMAT('ttime = ', F16.4)            
            STOP
            !
        END IF
        !
        ! SML estimation
        ! 
        CALL compute_estimate(theta,objf)
        !
        ! Retrieving eiscount statistics
        !
        CALL loglik(theta,ll,chiM,PthetaM,zetap,zetau)
     	!
        ! Printing intermediate trace output 
        !
        CALL print_res(i_trial,-ll,theta)
        !
    END DO
    !
    ! Closing output files
    !
    CLOSE(UNIT=unit_res_stima)
    CLOSE(UNIT=unit_ll_params)
    !
    ! Closing control files, if required
    !
    IF (switch_print_filter_control .EQ. 1) THEN
        !
        CALL close_filter_control()
        CALL close_discretization_control()
        !
    END IF
    IF (switch_print_theta .EQ. 1) CLOSE(unit_theta)
    !
END IF    
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Compute asymptotic variance matrix
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!IF (switch_asvar .EQ. 1) THEN
!    !
!    ! Reading parameters estimates 
!    !
!    CALL admissible_starting_point(seed,theta,llcheck)
!    !
!    ! Random number generation
!    !
!    seed = 1234567890                   ! Seed for initial r.n. generation
!    CALL generate_rn(seed)              ! r.n. generation
!    !
!    ! Computing and writing asymptotic standard errors
!    !
!    CALL compute_asymptotic_variance(theta,llcheck,objf,grad,param,stderr)
!    CALL print_final_results(param,stderr,objf,grad)
!    !
!END IF
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! End of computation of the asymptotic variance matrix 
! Halting execution
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END PROGRAM main
