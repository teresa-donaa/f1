MODULE starting_points
!
USE constants
USE io_stuff
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	SUBROUTINE admissible_starting_point ( seed, theta, llcheck )
	!
	! Computes an admissible starting parameters vector for the PPR iterations
	!
	IMPLICIT NONE
	!
	! Declaring dummy variables
	!
	INTEGER, INTENT(INOUT) :: seed(2)				! Seed for r.n. generation
	REAL(8), INTENT(OUT) :: theta(num_theta)		! Admissible value of the parameters vector
	REAL(8), INTENT(OUT) :: llcheck
    !
    ! Declaring external routine
    !
    REAL(8), EXTERNAL :: r8_uniform_01
    !
	! Declaring local variables
	!
	INTEGER :: i									! Index
    CHARACTER(len=2) :: ichar
	!
	! Beginning execution
    !
    ! Start from theta read from file (only first trial)
    !
    IF ((switch_readtheta .EQ. 1) .AND. (i_trial .EQ. 1)) THEN
        !
        CALL open_read_file(unit_theta,file_theta)
        READ(unit_theta,*) theta
        CLOSE(UNIT = unit_theta)
        RETURN
        !
    END IF
	!
	! Stima1: Starting points chosen randomly
	!
    IF (switch_stima1 .EQ. 1) THEN
		!
        CALL initialize()
        CALL set_initial_seed(seed(1),seed(2))
		DO i = 1, num_theta
            theta(i) = -buini+2*buini*r8_uniform_01()
		END DO 
        RETURN
		!
    END IF
	!
	! Stima2: Starting point read from the res_stima1.txt file
	!
	IF (switch_stima2 .EQ. 1) THEN
        !
        CALL GETARG(1,ichar)
        i = INUM(ichar)
        IF (i .EQ. 0) i = 1
        !
        CALL open_read_file(unit_res_stima1,file_res_stima1)
        IF (i .GT. 1) THEN
            READ(unit_res_stima1, 2534) 
2534        FORMAT(<i-1>(/))
            BACKSPACE(UNIT=unit_res_stima1)
        END IF
        READ(unit_res_stima1,*) llcheck, theta
        CLOSE(UNIT=unit_res_stima1)
        RETURN
        !
	END IF 
	!
	! Asvar: Starting point read from the res_stima2.txt file
	!
	IF (switch_asvar .EQ. 1) THEN
        !
        CALL GETARG(1,ichar)
        i = INUM(ichar)
        IF (i .EQ. 0) i = 1
        !
        CALL open_read_file(unit_res_stima2,file_res_stima2)
        IF (i .GT. 1) THEN
            READ(unit_res_stima2, 2535) 
2535        FORMAT(<i-1>(/))
            BACKSPACE(UNIT=unit_res_stima2)
        END IF
        READ(unit_res_stima2,11) llcheck, theta
11      FORMAT(34X, ES25.18, 2X, <num_theta>(1X, ES25.18))
        CLOSE(UNIT=unit_res_stima2)
        RETURN
        !
    END IF 
    !
END SUBROUTINE admissible_starting_point
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE starting_points