MODULE minEntropy
!
USE constants
USE observations
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE estimate_lambda ( L, M, lambda, Tdiff, q, obj, grad, task )
    !
    IMPLICIT NONE
    ! 
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: L, M
    REAL(8), DIMENSION(M,L), INTENT(IN) :: Tdiff
    REAL(8), DIMENSION(M), INTENT(IN) :: q
    REAL(8), INTENT(INOUT) :: lambda(L)
    REAL(8), INTENT(OUT) :: obj
    REAL(8), INTENT(OUT) :: grad(L)
    CHARACTER(len=60), INTENT(OUT) :: task
    !
    ! Beginning execution
    !
    IF (lambda_method .EQ. 1) THEN
        !
        CALL estimate_lambda_lbfgs(L,M,lambda,Tdiff,q,obj,grad,task)
        !
    END IF
    IF (lambda_method .EQ. 2) THEN
        !
        CALL estimate_lambda_NR(L,M,lambda,Tdiff,q,obj,grad)
        task = '                                                            '
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE estimate_lambda
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE estimate_lambda_lbfgs ( L, M, lambda, Tdiff, q, obj, grad, task )
    !
    IMPLICIT NONE
    ! 
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: L, M
    REAL(8), DIMENSION(M,L), INTENT(IN) :: Tdiff
    REAL(8), DIMENSION(M), INTENT(IN) :: q
    REAL(8), INTENT(INOUT) :: lambda(L)
    REAL(8), INTENT(OUT) :: obj
    REAL(8), INTENT(OUT) :: grad(L)
    CHARACTER(len=60), INTENT(OUT) :: task
    !
    ! Declaring local variables
    !
    INTEGER :: isave(44), iter, i
    INTEGER, ALLOCATABLE :: nbd(:), iwa(:)
    REAL(8) :: dsave(29)
    REAL(8), ALLOCATABLE :: wa(:), low(:), u(:)
    CHARACTER(len=60) :: csave
    LOGICAL :: lsave(4)
    REAL(8) :: factr = 1.d+1
    REAL(8) :: pgtol = 1.d-12
    !
    ! Beginning execution
    !
    ALLOCATE ( nbd(L), iwa(3*L) )
    ALLOCATE ( wa(2*20*L + 5*L + 11*20*20 + 8*20), low(L), u(L) )
    nbd = 0
    task = 'START'
    iter = 0
    DO WHILE ((task(1:2) .EQ. 'FG') .OR. (task .EQ. 'NEW_X') .OR. (task .EQ. 'START')) 
        !
        CALL setulb ( L, 20, lambda, low, u, nbd, obj, grad, factr, pgtol, &
            wa, iwa, task, -1, csave, lsave, isave, dsave )
        ! 
        IF (task(1:2) .eq. 'FG') CALL objectiveFunction_lbfgs(L,M,lambda,Tdiff,q,obj,grad) 
        !
    END DO
    DEALLOCATE ( nbd, iwa )
    DEALLOCATE ( wa, low, u )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE estimate_lambda_lbfgs
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE objectiveFunction_lbfgs ( L, M, lambda, Tdiff, q, obj, grad ) 
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: L, M
    REAL(8), DIMENSION(L), INTENT(IN) :: lambda
    REAL(8), DIMENSION(M,L), INTENT(IN) :: Tdiff
    REAL(8), DIMENSION(M), INTENT(IN) :: q
    REAL(8), INTENT(OUT) :: obj
    REAL(8), DIMENSION(L), INTENT(OUT) :: grad
    !
    ! Declaring local variables
    !
    REAL(8), DIMENSION(M) :: temp
    REAL(8), DIMENSION(M,L) :: temp2
    !
    ! Beginning execution
    !
    ! Compute objective function
    !
    temp = q*EXP(MATMUL(Tdiff,lambda))
    obj = SUM(temp)
    !
    ! Compute gradient of objective function
    !
    temp2 = SPREAD(temp,DIM = 2,NCOPIES = L)*Tdiff
    grad = SUM(temp2,DIM = 1)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE objectiveFunction_lbfgs
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE estimate_lambda_NR ( L, M, lambda, Tdiff, q, obj, grad )
    !
    IMPLICIT NONE
    ! 
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: L, M
    REAL(8), DIMENSION(M,L), INTENT(IN) :: Tdiff
    REAL(8), DIMENSION(M), INTENT(IN) :: q
    REAL(8), INTENT(INOUT) :: lambda(L)
    REAL(8), INTENT(OUT) :: obj
    REAL(8), INTENT(OUT) :: grad(L)
    !
    ! Declaring local variables
    !
    REAL(8), DIMENSION(L) :: newlambda, oldlambda
    REAL(8) :: hess(L,L), slambda(L,1)
    INTEGER :: ipiv(L), info
    !
    ! Beginning execution
    !
    newlambda = lambda
    DO
        !
        oldlambda = newlambda
        CALL objectiveFunction_NR(L,M,oldlambda,Tdiff,q,obj,grad,hess) 
        slambda(:,1) = grad
        CALL DGESV(L,1,hess,L,ipiv,slambda,L,info)
        newlambda = oldlambda-slambda(:,1)
        IF (MAXVAL(ABS(slambda)) .LT. 1.D-11) EXIT
        !
    END DO
    lambda = newlambda
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE estimate_lambda_NR
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE objectiveFunction_NR ( L, M, lambda, Tdiff, q, obj, grad, hess ) 
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: L, M
    REAL(8), DIMENSION(L), INTENT(IN) :: lambda
    REAL(8), DIMENSION(M,L), INTENT(IN) :: Tdiff
    REAL(8), DIMENSION(M), INTENT(IN) :: q
    REAL(8), INTENT(OUT) :: obj
    REAL(8), DIMENSION(L), INTENT(OUT) :: grad
    REAL(8), DIMENSION(L,L), INTENT(OUT) :: hess
    !
    ! Declaring local variables
    !
    REAL(8), DIMENSION(M) :: temp
    REAL(8), DIMENSION(M,L) :: temp2
    !
    ! Beginning execution
    !
    ! Compute objective function
    !
    temp = q*EXP(MATMUL(Tdiff,lambda))
    obj = SUM(temp)
    !
    ! Compute gradient of objective function
    !
    temp2 = SPREAD(temp,DIM = 2,NCOPIES = L)*Tdiff
    grad = SUM(temp2,DIM = 1)
    !
    ! Compute hessian of objective function
    !
    hess = MATMUL(TRANSPOSE(temp2),Tdiff)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE objectiveFunction_NR
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE minEntropy
