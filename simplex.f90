MODULE simplex
!
USE constants
USE observations
USE io_stuff
!
IMPLICIT NONE
!
! Minimization of a function (passed through module latent_criterion) in N dimensions 
! by the downhill simplex method of Nelder and Mead. 
! Inputs are:
! - theta, the (num_theta x 1) vector of starting values of the variables
! - perttheta, the (num_theta x 1) vector of perturbations of theta. 
!   This is used to compute the vertices of the starting simplex
! - view, the control variable to switch between types of output
! - iuser, vector of integers passed to the objective functions
! - user, vector of reals passed to the objective functions
! On output:
! - theta holds the parameters' values associated to the minimum of the objective function
! - fmin holds the minimum value of the objective function
! - iter is the number of iterations made
!
! Version date: 2002/03/07
!
CONTAINS
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE politope ( i_stime, loglik_politope, theta, perttheta, &
        view, iter_print, fmin, iter )
    !
    REAL(8), EXTERNAL :: loglik_politope
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime              ! Index of estimation trial
    REAL(8), INTENT(INOUT) :: theta(num_theta)
                                                ! Vector of starting values of the variables
    REAL(8), INTENT(IN) :: perttheta(num_theta)
                                                ! Vector of perturbations of theta
    INTEGER, INTENT(IN) :: view                 ! Controls printing of intermediate results:
                                                ! if == 1, prints; if /=1, does not print
    INTEGER, INTENT(IN) :: iter_print           ! Controls number of iterations between intermediate output
    REAL(8), INTENT(OUT) :: fmin                ! Minimum value of the objective function
    INTEGER, INTENT(OUT) :: iter                ! Number of iterations
    !
    ! Declaring local variables
    !
    INTEGER :: i_repeat_politope                ! Loop control variable
    REAL(8) :: theta_old(num_theta)             ! Temporary parameters values
    REAL(8) :: p(num_theta+1,num_theta)         ! Coordinates of the vertices of the simplex 
    REAL(8) :: y(num_theta+1)                   ! Function values at the vertices of the simplex
    INTEGER :: i                                ! Loop control variable
    REAL(8) :: fmin_old,crit,crit_politope_conv
    !
    ! Beginning execution
    !
    ! Beginning loop of politope minimization trials
    !
    fmin = loglik_politope(theta)
    politope_loop: DO i_repeat_politope = 1, max_repeat_politope
        !
        WRITE(unit_politope, 9873) 
        9873 FORMAT (/, '(Re)Entering politope ...')
        !
        ! Computing the vertices of the starting simplex and the associated function values
        !
        IF (iter_print .GE. 1) THEN
            WRITE (unit_politope,7588) i_stime, i_repeat_politope
            7588 FORMAT ( 1X, /, &
                          ' Trial       : ', I5, /, &
                          ' Restart     : ', I5, /, &
                          ' simpl               loglik         rtol                theta ')
        END IF
        theta_old = theta
        fmin_old = fmin
        p(1,:) = theta_old
        y(1) = fmin_old
        WRITE (unit_politope,7586) y(1), p(1,:)
        7586 FORMAT ( 5X, '1', 1X, ES20.13, 13X, <num_theta>(1X, E20.13) )        
        rows_p_loop: DO i=2,num_theta+1
            p(i,:) = theta_old
            IF (ABS(p(i,i-1)) .NE. 0.d0) THEN
                p(i,i-1) = p(i,i-1)*(1.d0+perttheta(i-1))
            ELSE
               p(i,i-1) = 1.d-1
            END IF
            y(i) = loglik_politope(p(i,:))
            WRITE (unit_politope,7587) i, y(i), p(i,:)
        7587 FORMAT ( 1X, I5, 1X, ES20.13, 13X, <num_theta>(1X, E20.13) )        
        END DO rows_p_loop
        !
        ! Call to the minimization subroutine
        !
        IF (iter_print .GE. 1) THEN
            WRITE (unit_politope,7589) 
            7589 FORMAT ( 1X, /, &
                          '  iter               loglik         rtol                theta ')
        END IF
        CALL amoeba_private(i_stime,i_repeat_politope,loglik_politope,view,iter_print,p,y,perttheta,iter)
        !
        ! Recovering the minimum of func and its coordinates
        !
        theta = p(SUM(MINLOC(y)),:)
        fmin = MINVAL(y)
        !
        ! Writing message about the results of the current trial
        !
        IF (crit_conv_formula .EQ. 1) THEN
            crit = ABS(fmin-fmin_old)
            crit_politope_conv = crit_politope_conv_y

        ELSE IF (crit_conv_formula .EQ. 2) THEN
            crit = MAXVAL(ABS(theta-theta_old))
            crit_politope_conv = crit_politope_conv_p
        END IF
        !
        !IF (view .EQ. 1) THEN
        !   WRITE (*,7528) i_repeat_politope, iter, theta, fmin, &
        !       crit, crit_politope_conv
        !END IF
        WRITE (unit_politope,7528) i_repeat_politope, iter, theta, fmin, &
            crit, crit_politope_conv
        7528 FORMAT (' Restart:     ', I5, /, &
                     ' Iterations:  ', I5, /, &
                     ' Parameters:  ', <num_theta>(ES12.5, 1X), /, &
                     ' Objective:   ', ES30.23, /, &
                     ' Improvement: ', ES12.5, /, &
                     ' Crit. Thres.:', ES12.5)
        !
        ! Check if this politope trial has significantly improved the solution
        !
        IF (crit .LE. crit_politope_conv)  EXIT
        !
    END DO politope_loop
    !
    END SUBROUTINE politope
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE amoeba_private ( i_stime, i_repeat_politope, loglik_politope, view, iter_print, p, y, perttheta, iter )
    !
    EXTERNAL loglik_politope
    REAL(8) :: loglik_politope
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime              ! Index of estimation trial
    INTEGER, INTENT(IN) :: i_repeat_politope    ! Loop control variable
    INTEGER, INTENT(IN) :: view
    INTEGER, INTENT(IN) :: iter_print
    REAL(8), INTENT(INOUT) :: p(num_theta+1,num_theta) 
                                                ! Coordinates of the vertices of the simplex 
    REAL(8), INTENT(INOUT) :: y(num_theta+1)    ! Function values at the vertices of the simplex
    REAL(8), INTENT(IN) :: perttheta(num_theta) ! Vector of perturbations of theta
    INTEGER, INTENT(OUT) :: iter                ! Number of iterations
    !
    ! Declaring local variables 
    !
    INTEGER :: iter_old
    INTEGER :: i,ilo,inhi
    REAL(8) :: rtol,ysave,ytry,ytmp
    REAL(8) :: ptry(num_theta)
    CHARACTER(len=3) :: fac
    REAL(8) :: dum
    REAL(8) :: dumv(num_theta)
    INTEGER :: ihi                              ! Index of best point in the simplex
    REAL(8) :: psum(num_theta)                  ! Sum of elements of p along the rows
    REAL(8) :: tol_politope
    !
    ! Beginning execution
    !
    iter=0
    iter_old=0
    psum(:)=SUM(p(:,:),dim=1)
    DO                                          ! Iteration loop
    !
        ilo=SUM(MINLOC(y))                      ! Determines which point is the highest (worst),
        ihi=SUM(MAXLOC(y))                      ! next-highest, and lowest (best)
        !
        ytmp=y(ihi)
        y(ihi)=y(ilo)
        inhi=SUM(MAXLOC(y))
        y(ihi)=ytmp
        !
        ! Computing contraction measure. Several alternatives are available:
        !
        rtol_formula_choice: SELECT CASE (rtol_formula)
            !
            ! 1) The original contraction measure in Press et alii:
            !
            CASE (1)
                rtol=2.0d0*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo))+1.0d-10)
                tol_politope = tol_politope_y
            !
            ! 2) The contraction measure used in NAG subroutine e04ccf, which
            ! implements the Nelder and Mead algorithm:
            !
            CASE (2)
                rtol = SQRT(SUM((y-SUM(y)/REAL(num_theta+1))**2)/REAL(num_theta+1))
                tol_politope = tol_politope_y
            !
            ! 3) The contraction measure adopted in a GAUSS subroutine by J. H. McCulloch.
            ! (Notice that, contrary to the two measures above, it checks contraction in p,
            ! rather than y.)
            !
            CASE (3)
                rtol = MAXVAL((MAXVAL(p, DIM=1)-MINVAL(p, DIM=1))/perttheta)
                tol_politope = tol_politope_p
            !
            ! 4) The contraction measure adopted by R. Schoenberg in its implementation 
            ! of the Nelder and Mead algorithm. (Again, it checks contraction in p. I don't
            ! like it because it seems meaningless -- it doesn't work if y(ilo) is small or
            ! all the parameters are negative.)
            !
            CASE (4)
                rtol = MAXVAL(p(:,ilo)-p(:,ihi))/y(ilo)
                tol_politope = tol_politope_p
            !
            ! 5) The contraction measure equivalent to (2) above, but applied to the columns
            ! of p instead of y. (This one seems interesting and less restrictive than (3).)
            !
            CASE (5)
                rtol = MAXVAL( SQRT( SUM( ( &
                    p-SPREAD(SUM(p, DIM=1)/REAL(num_theta+1), DIM=1, NCOPIES=num_theta+1) &
                    )**2, DIM=1 )/REAL(num_theta+1) ) )
                tol_politope = tol_politope_p
            !
            ! Notice that the relevant critical value to gauge convergence depends on
            ! which alternative is chosen to compute rtol. Generally speaking, if the measure
            ! is based on y, then tol_politope_y is to be used; otherwise, tol_politope_p is
            ! to be preferred.
        END SELECT rtol_formula_choice
        !
        ! Print intermediate results
        !
        IF (((iter-iter_old) .GT. iter_print) .OR. (iter .EQ. 0)) THEN
            IF (iter .NE. 0) iter_old = iter_old+iter_print
            IF (view .EQ. 1) THEN
                !
                WRITE (*,7588) i_stime, i_repeat_politope, max_repeat_politope, iter, max_iters_politope, &
                    y(ilo), rtol
                7588 FORMAT ( ///, &
                              'Estimation trial    : ', I5, /, & 
                              'Restart number      : ', I5, ' / ', I5, /, & 
                              'At iteration number : ', I5, ' / ', I5, /, & 
                              'Best function value : ', ES30.23, /, &
                              'Criterion           : ', ES30.23)
                CALL print_screen(p(ilo,:))
                !
            END IF
            WRITE (unit_politope,7589) iter, y(ilo), rtol, p(ilo,:)
            7589 FORMAT ( 1X, I5, 1X, ES20.13, 1X, ES12.5, <num_theta>(1X, E20.13) )        
        END IF
        !
        ! Compute the fractional range from highest to lowest and return if satisfactory
        !   
        IF ((rtol .LT. tol_politope) .OR. (iter .GE. max_iters_politope)) THEN          
                                        ! If returning, put best point and value in slot 1
            dum = y(1)                  ! 
            y(1) = y(ilo)               ! Swaps y(1) with y(ilo)    
            y(ilo) = dum                !
            dumv = p(1,:)               !
            p(1,:) = p(ilo,:)           ! Swaps p(1,:) with p(ilo,:)
            p(ilo,:) = dumv             !
            RETURN
        END IF
        !
        ! Begin a new iteration. 
        ! First, extrapolate by a factor -1 through the face of the simplex
        ! across from the high point, i.e., reflect the simplex from the high point
        !
        fac = 'ref'
        CALL amotry ( loglik_politope, -1.0d0, ytry, ptry, ihi, psum, p, y )
        iter=iter+1
        IF (ytry <= y(ilo)) THEN        ! Gives a result better than the best point, so
            fac = 'exp'                 ! try an additional extrapolation by a factor of 2
            CALL amotry ( loglik_politope, 2.0d0, ytry, ptry, ihi, psum, p, y )
            iter=iter+1
        ELSE IF (ytry >= y(inhi)) THEN  ! The reflected point is worse than the second
            ysave=y(ihi)                ! highest, so look for an intermediate lower point,
            fac = 'hal'                 ! i.e., do a one-dimensional contraction
            CALL amotry ( loglik_politope, 0.5d0, ytry, ptry, ihi, psum, p, y )
            iter=iter+1
            IF (ytry >= ysave) THEN
            !
            ! Can't seem to get rid of that high point. 
            ! Better contract around the lowest (best) point
            !
                p(:,:)=0.5d0*(p(:,:)+SPREAD(p(ilo,:),1,SIZE(p,1)))
                DO i=1,num_theta+1
                    IF (i .NE. ilo) THEN
                        y(i)=loglik_politope(p(i,:))
                    END IF
                END DO
                iter=iter+num_theta
                psum(:)=sum(p(:,:),dim=1)
                fac = 'con'
            END IF
        END IF
        !
        ! Prints the current state of the minimization procedure
        !
        IF (view .EQ. 2) THEN
            WRITE (*,1211) fac, iter, ptry, ytry, y(ilo), rtol
            WRITE (unit_politope,1211) fac, iter, ptry, ytry, y(ilo), rtol
            1211 FORMAT (A3, 1X, I5, 1X, <num_theta>(F7.4, 1X), 2(F20.10, 1X), F7.4)
        END IF  
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE amoeba_private
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    SUBROUTINE amotry( loglik_politope, fac, ytry, ptry, ihi, psum, p, y )
    !
    ! Extrapolates by a factor FAC through the face of the simplex across from the high point,
    ! tries it, and replaces the high point if the new point is better
    !
    EXTERNAL loglik_politope
    REAL(8) :: loglik_politope
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: fac
    REAL(8), INTENT(OUT) :: ytry
    REAL(8), INTENT(OUT) :: ptry(num_theta)
    INTEGER, INTENT(IN) :: ihi                      ! Index of best point in the simplex
    REAL(8), INTENT(INOUT) :: psum(num_theta)       ! Sum of elements of p along the rows
    REAL(8), INTENT(INOUT) :: p(num_theta+1,num_theta) 
                                                    ! Coordinates of the vertices of the simplex 
    REAL(8), INTENT(INOUT) :: y(num_theta+1)        ! Function values at the vertices of the simplex
    !
    ! Declaring local variables 
    !
    REAL(8) :: fac1,fac2
    !
    ! Beginning execution
    !
    fac1=(1.0d0-fac)/num_theta
    fac2=fac1-fac
    ptry(:)=psum(:)*fac1-p(ihi,:)*fac2          
    ytry=loglik_politope(ptry)                      ! Evaluate the function at the trial point.
    IF (ytry < y(ihi)) THEN                         ! If it is better than the highest, then replace
        y(ihi)=ytry                                 ! the highest
        psum(:)=psum(:)-p(ihi,:)+ptry(:)
        p(ihi,:)=ptry(:)
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE amotry
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
END MODULE simplex
