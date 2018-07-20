MODULE estimates
!
USE constants
USE loglikelihood
USE observations
USE simplex
USE io_stuff
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE compute_estimate ( theta, objf )
    !
    ! Computes the sml estimate
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
	REAL(8), INTENT(INOUT) :: theta(num_theta)          ! On input, the starting values of the parameters;
                                                        ! On output, the estimated values of the parameters
    REAL(8), INTENT(OUT) :: objf                        ! Loglikelihood function at the optimum
    !
    ! Declaring local variables
    !
    INTEGER :: iter_opt									! Number of iterations performed in the optimization
    REAL(8) :: theta_ini(num_theta)						! Starting values of the parameters
	!
    ! Beginning execution
    !
    CALL open_write_file(unit_politope,file_politope)
	!
    ! Call to politope
    !
    CALL politope(i_trial,loglik_politope,theta,pert_theta,1,10,objf,iter_opt)
	!
    CLOSE(UNIT=unit_politope)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE compute_estimate
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION loglik_politope ( theta ) 
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)   
    !
    ! Declaring function's type
    !
    REAL(8) :: loglik_politope
    !
    ! Declaring local variables
    !
    REAL(8) :: ll
    REAL(8) :: chiM(num_M)
    REAL(8) :: PthetaM(num_M,num_M)
    REAL(8), DIMENSION(num_M,0:num_T) :: zetap, zetau
    !
    ! Beginning execution
    !
    ! Call to loglikelihood subroutine
    !
    CALL loglik(theta,ll,chiM,PthetaM,zetap,zetau)
    !
    IF (numerical_error_flag .NE. '   ') THEN
        !
        loglik_politope = loglik_penalty
        PRINT*, 'numerical_error_flag = ', numerical_error_flag
!        PAUSE
        !
    ELSE
        !
        loglik_politope = -ll
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END FUNCTION loglik_politope
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE estimates