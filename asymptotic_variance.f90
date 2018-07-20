MODULE asymptotic_variance
!
USE constants
USE observations
USE loglikelihood
USE estimates
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE compute_asymptotic_variance ( theta, llcheck, objf, grad, param, stderr )
    !
    ! Computes the asymptotic variance matrix of the maximum likelihood estimator 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)                     ! Estimates
    REAL(8), INTENT(IN) :: llcheck  
    REAL(8), INTENT(OUT) :: objf                                ! Optimized loglikelihood 
    REAL(8), INTENT(OUT) :: param(num_param)                    ! Estimates of true parameters
    REAL(8), INTENT(OUT) :: grad(num_param)                     ! Gradient at psi
    REAL(8), INTENT(OUT) :: stderr(num_param)                   ! Estimated asymptotic standard errors
    !
    ! Declaring local variables
    !
    REAL(8), DIMENSION(num_S,num_T) :: f_up, w_up
    REAL(8) :: eiscount_mean
    INTEGER :: tt, i, eiscount_min, eiscount_max, eiscount_100
    REAL(8) :: llvec(num_T), llvecd(num_theta,num_T), dl(num_T,num_theta), dobjf(num_theta)
    REAL(8) :: llvecdd(num_theta,num_theta,num_T)
    REAL(8), DIMENSION(num_theta,num_theta) :: Imat, Jmat, invJmat
    INTEGER :: ipiv(num_theta), info                            ! Used in DGESV
    REAL(8) :: dparam_dtheta(num_param,num_theta), dparam(num_theta,num_param)
    REAL(8) :: var_theta(num_theta,num_theta), var_param(num_param,num_param)
    !
    ! Beginning execution
    !
    WRITE (*,11)
11  FORMAT (//, 'Computing asymptotic variance matrix', //)
    !
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! Reading parameters
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !
    objf = loglik_politope(theta)
    WRITE (*,12) theta, objf, llcheck
12  FORMAT('Parameters on input:', /, <num_theta>(ES9.2,1X), /, &
    'With these parameters, the average loglikelihood equals ', ES20.13, /, &
    'Should have been                                        ', ES20.13)
    IF (ABS(objf-llcheck) .GT. 1.d-10) THEN 
        PRINT*, 'Problem:'
        PRINT*, 'The loglikelihood value of the starting parameters'
        PRINT*, 'is not equal to the previously computed value'
        PRINT*, 'Type any key to stop'
    END IF
    objf = -objf
!    !
!    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    ! Computing the I matrix (outer product of gradients) and 
!    ! Computing the J matrix (-Hessian/N) and its inverse
!    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    !
!    WRITE (*,*) 'Computing gradients and hessian ...'
!    !
!    CALL LOGLIK_VEC_DV_DV(theta,eye_theta,eye_theta,llvec,llvecd, &
!        llvecdd,f_up,w_up,eiscount_mean,eiscount_min,eiscount_max, &
!        eiscount_100,num_theta,num_theta)
!    !
!    ! Gradients
!    !
!    WRITE (*,*) 'Computing Imat ...'
!    !
!    dl = TRANSPOSE(llvecd)
!    CALL open_write_file(unit_dl,file_dl)
!    WRITE(unit_dl,15) (llvec(tt), dl(tt,:), tt = 1, num_T)
!15  FORMAT (<num_T>(ES20.13, ' #', <num_theta>(1X, ES20.13), / ))
!    CLOSE(UNIT=unit_dl)
!    dobjf = -SUM(llvecd,DIM=2)/num_T
!    !
!    dl = dl-SPREAD(SUM(dl,DIM=1)/num_T,DIM=1,NCOPIES=num_T)
!    Imat = MATMUL(TRANSPOSE(dl),dl)/num_T
!    CALL open_write_file(unit_Imat,file_Imat)
!    WRITE(unit_Imat,14) (Imat(i,:), i = 1, num_theta)
!14  FORMAT(<num_theta>(<num_theta>(ES20.13,1X), /))
!    CLOSE(UNIT=unit_Imat)
!    !
!    WRITE (*,*) 'Computing Jmat ...'
!    !
!    CALL open_write_file(unit_dl2,file_dl2)
!    WRITE(unit_dl2,18) (llvecdd(:,:,tt), tt = 1, num_T)
!18  FORMAT (<num_T>(<num_theta**2>(1X, ES20.13), / ))
!    CLOSE(UNIT=unit_dl2)
!    !
!    Jmat = -SUM(llvecdd,DIM=3)/num_T
!    CALL open_write_file(unit_Jmat,file_Jmat)
!    WRITE(unit_Jmat,14) (Jmat(i,:), i = 1, num_theta)
!    CLOSE(UNIT=unit_Jmat)
!    !
!    ! Inverting the J matrix
!    !
!    invJmat = eye_theta
!    CALL DGESV(num_theta,num_theta,Jmat,num_theta,ipiv,invJmat,num_theta,info)
!    IF (info .NE. 0) THEN
!        PRINT*, "dgesv can't be applied to Jmat"
!        STOP
!    END IF
!    invJmat = (invJmat+TRANSPOSE(invJmat))/2.d0      ! Symmetrification
!    !
!    CALL open_write_file(unit_invJmat,file_invJmat)
!    WRITE(unit_invjmat,14) (invJmat(i,:), i = 1, num_theta)
!    CLOSE(UNIT=unit_invJmat)
!    !
!    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    ! Computing the param array and the dparam_dtheta matrix
!    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    !
!    WRITE (*,*) 'Computing param and dparam_dtheta ...'
!    !
!    CALL VEC_THETA_TO_PARAM_DV(theta,eye_theta,param,dparam,num_theta)
!    dparam_dtheta = TRANSPOSE(dparam)
!    grad = MATMUL(dparam_dtheta,dobjf)
!    CALL open_write_file(unit_dparam_dtheta,file_dparam_dtheta)
!    WRITE(unit_dparam_dtheta,14) (dparam_dtheta(i,:), i = 1, num_param)
!    CLOSE(UNIT=unit_dparam_dtheta)
!    CALL open_write_file(unit_param,file_param)
!    WRITE(unit_param,17) (param(i), i = 1, num_param)
!17  FORMAT(<num_param>(ES20.13, 1X))
!    CLOSE(UNIT=unit_param)
!    !
!    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    ! Computing the variance matrix V and the associated standard errors
!    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    !
!    var_theta = MATMUL(MATMUL(invJmat,Imat),invJmat)                                ! Var(theta)
!    var_param = MATMUL(MATMUL(dparam_dtheta,var_theta),TRANSPOSE(dparam_dtheta))    ! Var(param)
!    !
!    CALL open_write_file(unit_var_param,file_var_param)
!    WRITE(unit_var_param,16) (var_param(i,:), i = 1, num_param)
!16  FORMAT(<num_param>(<num_param>(ES20.13, 1X), /))
!    CLOSE(UNIT=unit_var_param)
!    !
!    DO i = 1, num_param
!        IF (var_param(i,i) .GE. 0.d0) stderr(i) = SQRT(var_param(i,i)/REAL(num_T))
!    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE compute_asymptotic_variance
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE asymptotic_variance
