MODULE load_data
!
USE constants
USE observations
USE io_stuff
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE load_obs ( unit_number, file_name )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit_number
    CHARACTER(len=30), INTENT(IN) :: file_name
    !
    ! Declaring local variables
    !
    INTEGER :: tot_num_C, tt, cc, flag_eof, im
    INTEGER :: date_temp, tau_temp, NC_temp, NT_temp, NK_temp
    REAL(8) :: S_temp, C_temp, K_temp, X_temp, r_temp, q_temp, RV_temp
    !
    ! Beginning execution
    !
    CALL open_read_file(unit_number,file_name)   
    !
    ! Skip the first row of the data files (containing names)
    !
	READ(unit_number,101)
101 FORMAT()
    !
    ! Beginning loop over observation dates
    !
    tot_num_C = 0
    DO tt = 1, num_T
        !
        ! Beginning loop over observations records at date tt
        !
        DO cc = 1, num_max_C_ds
            !
            ! Reading new record
            !
            READ(unit_number,*,IOSTAT=flag_eof) date_temp, S_temp, C_temp, &
                K_temp, X_temp, tau_temp, NC_temp, NT_temp, NK_temp, &
                r_temp, q_temp, RV_temp
            !
            ! Check if end of file is attained
            !
            IF (flag_eof < 0) THEN
                EXIT
            END IF
            !
            ! Reading variables constant over dates and options
            !
            IF ((tt .EQ. 1) .AND. (cc == 1)) THEN
                !
                r = r_temp
                q = q_temp
                rmq = r-q
                !
            END IF
            !
            ! Reading variables constant across options within the same day
            !
            IF (cc == 1) THEN
                !
                date(tt)    = date_temp
                S(tt)       = S_temp
                num_C(tt)   = NC_temp
                num_tau(tt) = NT_temp
                RV(tt)      = RV_temp
                !
                IF (tt .EQ. 1) y(tt) = 0.d0
                IF (tt .GT. 1) y(tt) = LOG(S(tt)/S(tt-1))
                !
            END IF
            !
            ! If date has changed, backspace one record and exit inner cycle
            !
            IF (date_temp /= date(tt)) THEN
                BACKSPACE(unit_number)
                EXIT
            END IF
            !
            ! Reading variables different across options within the same day
            !
            C(cc,tt)     = C_temp
            K(cc,tt)     = K_temp
            X(cc,tt)     = X_temp
            tau(cc,tt)   = tau_temp
            num_K(cc,tt) = NK_temp
            H(cc,tt)     = C_temp/(K_temp*EXP(-r_temp*tau_temp))
            !
        END DO
        !
        ! Screen output
        !
        tot_num_C = tot_num_C+num_C(tt)
        IF (switch_show_tt_ld .EQ. 1) THEN
            !
            WRITE(*,3) date(tt), tt, num_C(tt), tot_num_C
3           FORMAT('date: ', I8, ' / i_T : ', I5, ' / num_C : ', I3, ' / tot_num_C : ', I7)
            !
        END IF
        !
    END DO
    !
    ! Closing data files
    !
    CLOSE(UNIT=unit_number)
    !
    ! Constructing identity matrices of orders num_theta and num_M
    !
    eye_theta = 0.d0
    DO im = 1, num_theta
        eye_theta(im,im) = 1.d0
    END DO
    eye_M = 0.d0
    DO im = 1, num_M
        eye_M(im,im) = 1.d0
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE load_obs 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE load_data
