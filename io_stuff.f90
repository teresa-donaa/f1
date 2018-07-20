MODULE io_stuff
!
USE constants
USE observations
USE loglikelihood
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_write_file ( unit_number, file_name )
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
    INTEGER :: open_err,i                   ! Open file error code
    CHARACTER(len=1) :: q
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='replace', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        READ*, q
        STOP
    END IF
    1234 FORMAT ('Unable to open output file ', A30, /, &
                 'Press any key to stop')
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_write_file 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_read_file ( unit_number, file_name )
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
    INTEGER :: open_err,i                   ! Open file error code
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='old', ACTION='read', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        STOP
    END IF
    1234 FORMAT ('Unable to open output file ', A80)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_read_file 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_screen ( theta )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)
    !
    ! Declaring local variables
    !
    REAL(8) :: lambda, nu, mu, beta
    REAL(8) :: lambdaQ, nuQ, muQ, betaQ
    REAL(8) :: eta(num_eta), chi(num_chi), xi(num_xi)
    !
    ! Beginning execution
    !
    CALL psi_params(theta,'P',lambda,nu,mu,beta)
    CALL psi_params(theta,'Q',lambdaQ,nuQ,muQ,betaQ)
    CALL sigma_params(theta,eta,chi,xi)
    !
    WRITE(*,1234) lambda, nu, mu, beta, mu*beta, &
        lambdaQ, nuQ, muQ, betaQ, muQ*betaQ, &
        eta, chi, xi
1234    FORMAT('Parameters               : ', /, &
        <5>(ES9.2,1X), /, &
        <5>(ES9.2,1X), /, &
        <num_eta>(ES9.2,1X), /, &
        <num_chi>(ES9.2,1X), /, &
        <num_xi>(ES9.2,1X))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_screen 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_res ( i_stime, objf, theta )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime      ! Estimation trial number
    REAL(8), INTENT(IN) :: objf         ! Latent criterion function at the optimum
    REAL(8), INTENT(IN) :: theta(num_theta)
    !
    ! Declaring local variables
    !
    INTEGER :: i
    REAL(8) :: lambda, nu, mu, beta
    REAL(8) :: lambdaQ, nuQ, muQ, betaQ
    REAL(8) :: eta(num_eta), chi(num_chi), xi(num_xi)    
    !
    ! Beginning execution
    !
    ! Print stima_res file
    !
    WRITE (unit_res_stima,35) i_stime, objf, theta
35  FORMAT ( I3, " # ", ES25.18, " #", <num_theta>(1X, ES25.18))
    !
    ! Print ll_params file
    !
    CALL psi_params(theta,'P',lambda,nu,mu,beta)
    CALL psi_params(theta,'Q',lambdaQ,nuQ,muQ,betaQ)
    CALL sigma_params(theta,eta,chi,xi)
    WRITE (unit_ll_params,36) i_stime, objf, &
        lambda, nu, mu, beta, mu*beta, &
        lambdaQ, nuQ, muQ, betaQ, muQ*betaQ, &
        eta, chi, xi
36  FORMAT ( I3, " # ", ES25.18, &
        " #", <5>(1X, ES25.18), " #", <5>(1X, ES25.18), " #", <num_eta>(1X, ES25.18), " #", &
        <num_chi>(1X, ES25.18), " #", <num_xi>(1X, ES25.18))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_res
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE head_ll_params ( )
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: i
    !
    ! Beginning execution
    !
    WRITE(unit_ll_params,11) (i, i = 1, num_eta), (i, i = 1, num_chi), (i, i = 1, num_xi)
11  FORMAT("  i # ", &
        "                   '-ll/T # ", &
        "                   lambda ", &
        "                       nu ", &
        "                       mu ", &
        "                     beta ", &
        "                      rho # ", &
        "                  lambdaQ ", &
        "                      nuQ ", &
        "                      muQ ", &
        "                    betaQ ", &
        "                     rhoQ # ", &
        <num_eta>("                   eta(", I1, ") # "), &
        <num_chi>("                   chi(", I1, ") # "), &
        <num_xi>("                    xi(", I1, ") "))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE head_ll_params
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_filter_control ( )
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: i
    !
    ! Beginning execution
    !
    CALL open_write_file(unit_filter_control,file_filter_control)
    !
    WRITE(unit_filter_control, 1) (i, i = 1, num_M), (i, i = 1, num_M), (i, i = 1, num_M), &
        (i, i = 1, num_M), (i, i = 1, num_M), (i, i = 1, num_M)
1   FORMAT('  tt #               lik_tt #            llvec(tt) #', &
        <num_M>('            zetap(', I2, ')'), ' #', &
        <num_M>('            zetau(', I2, ')'), ' #', &
        <num_M>('             g_tt(', I2, ')'), ' #', &
        <num_M>('        logpdf_RV(', I2, ')'), ' #', &
        <num_M>('         logpdf_H(', I2, ')'), ' #', &
        <num_M>('         logpdf_y(', I2, ')'))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_filter_control
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_discretization_control ( )
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: i
    !
    ! Beginning execution
    !
    CALL open_write_file(unit_discretization_control,file_discretization_control)
    !
    WRITE(unit_discretization_control, 1) (i, i = 1, num_M), &
        (i, i = 1, num_L), (i, i = 1, num_L), (i, i = 1, num_M)
1   FORMAT('  mm #             chiM(mm) # ', <num_M>('      PthetaM(mm,', I2, ') '), '#         pithetaM(mm)', &
        ' #        num_L #  match_L(mm) # ', &
        <num_L>('        check(mm,', I2, ') '), '# ', <num_L>('  lambda_star(mm,', I2, ') '), '# ', <num_M>('      QthetaM(mm,', I2, ') '))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_discretization_control
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE close_filter_control ( )
    !
    IMPLICIT NONE
    !
    ! Beginning execution
    !
    CLOSE(UNIT=unit_filter_control)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE close_filter_control
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE close_discretization_control ( )
    !
    IMPLICIT NONE
    !
    ! Beginning execution
    !
    CLOSE(UNIT=unit_discretization_control)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE close_discretization_control
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    !SUBROUTINE print_final_results ( param, stderr, objf, grad )
    !!
    !! Computes final results and writes them on file
    !!
    !IMPLICIT NONE
    !!
    !! Declaring dummy variables
    !!
    !REAL(8), INTENT(IN) :: param(num_param)             ! Estimates of parameters
    !REAL(8), INTENT(IN) :: stderr(num_param)            ! Asymptotic standard errors
    !REAL(8), INTENT(IN) :: objf                         ! Optimized latent criterion 
    !REAL(8), INTENT(IN) :: grad(num_param)              ! Gradient at theta_star
    !!
    !! Declaring local variables
    !!
    !INTEGER :: i
    !CHARACTER(len=12), DIMENSION(num_param) :: names
    !!
    !! Beginning execution
    !!
    !names = (/ 'lambda', 'nu', 'mu', 'beta', 'rho', &
    !    'lambdaQ', 'nuQ', 'muQ', 'betaQ', 'rhoQ', &
    !    'phi', 'ystar', 'delta1', 'delta2', &
    !    'eta', 'etastar', 'chi', 'sigma', 'xi0', 'xiX', 'xiX2', 'xiTY', 'xiTY2', 'xiXTY' /)
    !CALL open_write_file(unit_finres,file_finres)
    !WRITE (unit_finres,3525) 
    !3525 FORMAT ( /, &
    !    '----------------------------------------------------', /, &
    !    '    Parameter     Estimate       AsySEI     Gradient', /, &
    !    '----------------------------------------------------' )
    !!
    !DO i = 1, num_param
    !    WRITE (unit_finres, 9999) names(i), param(i), stderr(i), grad(i)
    !    9999 FORMAT ( 1X, A12, 3(1X, ES12.5) )
    !END DO
    !!
    !WRITE (unit_finres,3524) 
    !!
    !WRITE (unit_finres, 9997) -objf*num_T
    !9997 FORMAT ( 1X, '  Total objf', 1X, ES12.5 )
    !!
    !WRITE (unit_finres, 9992) -objf
    !9992 FORMAT ( 1X, '   Avg. objf', 1X, ES12.5 )
    !!
    !WRITE (unit_finres,3524) 
    !!        
    !3524 FORMAT ( &
    !     '---------------------------------------------------' )
    !!         
    !CLOSE(UNIT=unit_finres)
    !!
    !! Ending execution and returning control
    !!
    !END SUBROUTINE print_final_results
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE io_stuff
