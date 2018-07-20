MODULE loglikelihood
!
USE constants
USE observations
USE option_prices
USE cdflib
USE minEntropy
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE loglik ( theta, ll, chiM, PthetaM, zetap, zetau )
    !
    ! Modified Kalman filter - loglikelihood contributions
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN)  :: theta(num_theta)
    REAL(8), INTENT(OUT) :: ll
    REAL(8), INTENT(OUT) :: chiM(num_M)
    REAL(8), INTENT(OUT) :: PthetaM(num_M,num_M)
    REAL(8), INTENT(OUT), DIMENSION(num_M,0:num_T) :: zetap, zetau
    !
    ! Declaring local variables
    !
    REAL(8) :: llvec(num_T)
    !
    ! Beginning execution
    !
    ! Computing loglikelihood
    !
    CALL loglik_vec(theta,llvec,chiM,PthetaM,zetap,zetau)
    IF (numerical_error_flag .NE. '   ') RETURN
    !
    ! Average loglikelihood
    !
    ll = SUM(llvec)/num_T
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE loglik
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE loglik_vec ( theta, llvec, chiM, PthetaM, zetap, zetau )
    !
    IMPLICIT NONE
    !
    ! Dummy variables
    !
    REAL(8), INTENT(IN)  :: theta(num_theta)			    ! Total vector of parameters
    REAL(8), INTENT(OUT) :: llvec(num_T)
    REAL(8), INTENT(OUT) :: chiM(num_M)
    REAL(8), INTENT(OUT) :: PthetaM(num_M,num_M)
    REAL(8), INTENT(OUT), DIMENSION(num_M,0:num_T) :: zetap, zetau
    !
    ! Local variables 
    !
    REAL(8) :: lambda, nu, mu, beta, lambdaQ, muQ, betaQ
    REAL(8) :: eta(num_eta), chi(num_chi), xi(num_xi)
    REAL(8) :: shape, scale, pgam, qgam, xgam, bound
    INTEGER :: i, mm, status, kk, match_L(num_M), tt, NC_tt, info
    REAL(8) :: QthetaM(num_M,num_M), Qvec(num_M), obj
    REAL(8) :: lambda_star(num_M,num_L), lambda0(num_L), Tdiff(num_M,num_L), grad(num_L), check(num_M,num_L)
    REAL(8) :: lambda_star1(num_M,1), lambda01(1), Tdiff1(num_M,1), grad1(1)
    REAL(8) :: lambda_star2(num_M,2), lambda02(2), Tdiff2(num_M,2), grad2(2), check2(2)
    CHARACTER(len=60) :: task
    REAL(8), DIMENSION(num_M,num_M) :: TPthetaM, VL, VR, B_rhs(num_M,1)
    INTEGER :: ipiv(num_M)
    REAL(8), DIMENSION(num_M) :: lchiM, wr, wi, pithetaM, g_tt, logpdf_RV, logpdf_H, logpdf_y
    REAL(8) :: work(10*num_M)
    REAL(8) :: lRV_tt, y_tt, lik_tt
    REAL(8), DIMENSION(num_grid_x,num_grid_tau,num_grid_f) :: Hgrid 
    COMPLEX(8), DIMENSION(num_N,min_tau:max_tau) :: Acf_ffo, Bcf_ffo
    !
    ! Beginning execution
    !
    numerical_error_flag = '   '
    IF (switch_print_theta .EQ. 1) WRITE(unit_theta,11) theta
    !
    ! P parameters
    !
    CALL psi_params(theta,'P',lambda,nu,mu,beta)
    !
    ! Q parameters
    !
    CALL psi_params(theta,'Q',lambdaQ,nu,muQ,betaQ)
    !
    ! Measurement errors variance parameters
    !
    CALL sigma_params(theta,eta,chi,xi)
    !
    ! Picking the loglikelihood routine to compute normalized option prices 
    ! and checking for NaNs
    !
    CALL option_prices_ffo_loglik(num_N,lambdaQ,nu,muQ,betaQ,Acf_ffo,Bcf_ffo)
    CALL option_prices_lin_loglik(lambdaQ,nu,muQ,betaQ,Hgrid,Acf_ffo,Bcf_ffo)
    IF (numerical_error_flag .NE. '   ') RETURN
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Computing discretization grid 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Setting the state discrete grid
    !
    IF (X_method .EQ. 1) THEN
        !
        ! Quantile, based on the Gamma marginal pdf
        !
        shape = nu
        scale = mu/(1.d0-beta*mu)
        DO mm = 1, num_M
            !
            pgam = (DBLE(mm)-0.5d0)/DBLE(num_M)
            qgam = 1.d0-pgam
            CALL cdfgam(2,pgam,qgam,xgam,shape,1.d0/scale,status,bound)
            IF (status .NE. 0) THEN
                !
                numerical_error_flag = 'gam' 
                RETURN
                !
            END IF
            chiM(mm) = xgam
            !
        END DO
        lchiM = LOG(chiM)
        !
    END IF
    !
    ! Computing Farmer and Toda's transition probabilities
    !
    lambda_star = 0.d0
    DO mm = 1, num_M
        !
        ! Setting the initial probability distribution
        !
        IF (Q_method .EQ. 1) THEN
            !
            ! Modified Gamma pdf
            !
            CALL pdf_gamma(1.d0,mu/(1.d0-beta*mu),num_M,chiM,'exp',Qvec)
            !
        END IF
        IF (Q_method .EQ. 2) THEN
            !
            ! Modified Gamma pdf
            !
            CALL pdf_gamma(mu,mu/(1.d0-beta*mu),num_M,chiM,'exp',Qvec)
            !
        END IF
        IF (Q_method .EQ. 3) THEN
            !
            ! Equal
            !
            Qvec = 1.d0
            !
        END IF
        IF (Q_method .EQ. 4) THEN
            !
            ! True Gamma cdf
            !
            Qvec = 1.d0
            CALL cdf_gamma(nu,mu/(1.d0-beta*mu),num_M-1,0.5d0*(chiM(2:)+chiM(:num_M-1)),'cdf',Qvec(:num_M-1))
            Qvec(2:) = Qvec(2:)-Qvec(:num_M-1)
            !
        END IF
        Qvec = Qvec/SUM(Qvec)
        QthetaM(mm,:) = Qvec
        !
        ! Solving the dual problem, as in Farmer and Toda code 
        !
        IF (num_L .EQ. 1) THEN      ! Match first moments only
            !
            Tdiff1 = T(1,chiM,chiM(mm),nu,beta,mu)-Tbar(1,chiM,chiM(mm),nu,beta,mu)
            lambda01 = 0.d0
            CALL estimate_lambda(1,num_M,lambda01,Tdiff1,Qvec,obj,grad1,task)
            lambda_star(mm,1) = lambda01(1)
            PthetaM(mm,:) = Qvec*EXP(MATMUL(Tdiff1,lambda01))
            PthetaM(mm,:) = PthetaM(mm,:)/SUM(PthetaM(mm,:))
            check(mm,1) = grad1(1)/obj
            match_L(mm) = 1
            !
        ELSE        ! Match first two moments 
            !
            Tdiff2 = T(2,chiM,chiM(mm),nu,beta,mu)-Tbar(2,chiM,chiM(mm),nu,beta,mu)
            IF (mm .EQ. 1) lambda02 = lambda_star(mm,:2)
            IF (mm .GT. 1) lambda02 = lambda02
            CALL estimate_lambda(2,num_M,lambda02,Tdiff2,Qvec,obj,grad2,task)
            check2 = grad2/obj
            !
            IF (SQRT(SUM(check2**2)) .GT. 1.d-5) THEN       ! If first 2 moments fail, then match first moment only
                !
                Tdiff1 = T(1,chiM,chiM(mm),nu,beta,mu)-Tbar(1,chiM,chiM(mm),nu,beta,mu)
                lambda01 = 0.d0
                CALL estimate_lambda(1,num_M,lambda01,Tdiff1,Qvec,obj,grad1,task)
                lambda_star(mm,1) = lambda01(1)
                lambda02 = 0.d0
                PthetaM(mm,:) = Qvec*EXP(MATMUL(Tdiff1,lambda01))
                PthetaM(mm,:) = PthetaM(mm,:)/SUM(PthetaM(mm,:))
                check(mm,1) = grad1(1)/obj
                match_L(mm) = 1
                !
            ELSE 
                !
                lambda_star(mm,:2) = lambda02
                PthetaM(mm,:) = Qvec*EXP(MATMUL(Tdiff2,lambda02))
                PthetaM(mm,:) = PthetaM(mm,:)/SUM(PthetaM(mm,:))
                check(mm,:2) = check2
                match_L(mm) = 2
                !
            END IF
            !
        END IF
        !
    END DO
    !
    ! Compute stationary distribution 
    ! (Normalized) Eigenvector associated to the unit eigenvalue of PthetaM'
    !
    TPthetaM = TRANSPOSE(PthetaM)
    CALL DGEEV('N','V',num_M,TPthetaM,num_M,wr,wi,VL,num_M,VR,num_M,work,10*num_M,info)
    i = SUM(MAXLOC(wr))
    pithetaM = VR(:,i)/SUM(VR(:,i))
	!
    ! Write control files if required
    !
    IF (switch_print_filter_control .EQ. 1) THEN
        !
        WRITE(unit_discretization_control, 2) &
            (mm, chiM(mm), PthetaM(mm,:), pithetaM(mm), num_L, match_L(mm), check(mm,:), &
            lambda_star(mm,:), QthetaM(mm,:), mm = 1, num_M)
2       FORMAT(<num_M>(I4, ' # ', ES20.8, ' # ', <num_M>(ES20.8, 1X), '# ', ES20.8, ' # ', I12, ' # ', I12, ' # ', <num_L>(ES20.8, 1X), '# ', & 
            <num_L>(ES20.8, 1X), '# ', <num_M>(ES20.8, 1X), /))
        !
    ENDIF
!@SP
!CALL CPU_TIME(time1)
!PRINT*, 'after discretization: ', time1-time0
!time0 = time1
!@SP
    !
    ! Beginning cycle over observations
    !
	llvec = 0.d0
    zetap = 0.d0
    zetau = 0.d0
    zetau(:,0) = pithetaM
	DO tt = 1, num_T
        !
        IF (switch_show_tt .EQ. 1) WRITE(*,10) tt
		!
!@SP
!IF (tt .EQ. 813) THEN
!    PRINT*, 'Date: ', date(tt)
!    PRINT*, 'Halt here!'
!ENDIF
!@SP
        !
        lRV_tt = LOG(RV(tt))
        y_tt = y(tt)
        NC_tt = num_C(tt)
        !
        ! Compute predicted probabilities: zetap(:,tt)
        !
        zetap(:,tt) = MATMUL(TRANSPOSE(PthetaM),zetau(:,tt-1))
!@SP
!PRINT*, 'tt = ', tt
!PRINT*, 'zetau(:,tt-1) = ', zetau(:,tt-1)
!PRINT*, 'pithetaM = ', pithetaM
!PAUSE
!@SP
        !
        ! Compute observables' contributions: g(tt)
        !
        CALL eval_g_tt(tt,lambda,nu,mu,beta,eta,chi,xi,lRV_tt,y_tt,Hgrid,NC_tt,lchiM,chiM, &
            logpdf_RV,logpdf_H,logpdf_y) 
        g_tt = EXP(logpdf_H+logpdf_y+logpdf_RV)
        !
        ! Compute the loglikelihood contribution
        !
	    lik_tt = SUM(g_tt*zetap(:,tt))+bulletproof
	    llvec(tt) = LOG(lik_tt)
        IF (llvec(tt) .NE. llvec(tt)) THEN
            !
            numerical_error_flag = 'llt'
!@SP
!PRINT*, 'tt = ', tt
!PRINT*, 'g_tt = ', g_tt
!PRINT*, 'zetap(:,tt) = ', zetap(:,tt)
!PAUSE
!@SP
            RETURN
            !
        END IF
        !
        ! Compute updated probabilities: zetau(:,tt)
        !
        zetau(:,tt) = (g_tt*zetap(:,tt))/lik_tt
        !
        ! Write control files if required
        !
        IF (switch_print_filter_control .EQ. 1) THEN
            !
            WRITE(*, 1) tt
1           FORMAT('tt = ', I5)  
            !
            WRITE(unit_filter_control, 3) tt, lik_tt, llvec(tt), zetap(:,tt), zetau(:,tt), &
                g_tt, logpdf_RV, logpdf_H, logpdf_y
3           FORMAT(I4, ' # ', ES20.8E3, ' # ', ES20.8, ' # ', <num_M>(ES20.8E3, 1X), '# ', & 
                <num_M>(ES20.8E3, 1X), '# ', <num_M>(ES20.8E3, 1X), '# ', <num_M>(ES20.8E3, 1X), '# ', &
                <num_M>(ES20.8E3, 1X), '# ', <num_M>(ES20.8E3, 1X))
            !
        ENDIF
		!
	END DO
    !
10  FORMAT('tt = ', I5)        
11  FORMAT(<num_theta>(1X,ES30.23))
!@SP
!CALL CPU_TIME(time1)
!PRINT*, 'after loop over dates: ', time1-time0
!@SP
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE loglik_vec
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE eval_g_tt ( tt, lambda, nu, mu, beta, eta, chi, xi, &
        lRV_tt, y_tt, Hgrid, NC_tt, lf_tt, f_tt, logpdf_RV, logpdf_H, logpdf_y ) 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: tt
    REAL(8), INTENT(IN) :: lambda, nu, mu, beta
    REAL(8), INTENT(IN) :: eta(num_eta), chi(num_chi), xi(num_xi)
    REAL(8), INTENT(IN) :: lRV_tt, y_tt
    REAL(8), INTENT(IN), DIMENSION(num_grid_x,num_grid_tau,num_grid_f) :: Hgrid 
    INTEGER, INTENT(IN) :: NC_tt
    REAL(8), INTENT(IN), DIMENSION(num_M) :: lf_tt, f_tt
    REAL(8), DIMENSION(num_M), INTENT(OUT) :: logpdf_RV, logpdf_H, logpdf_y
    !
    ! Declaring local variables
    !
    REAL(8), DIMENSION(num_M) :: E_lRV, sd_lRV, zz_lRV, E_y, sd_y, zz_y
    INTEGER :: ss, j, tau_tt(NC_tt)
    REAL(8), DIMENSION(NC_tt) :: E_lH, sd_lH, zz_C, lH, X_tt, lpdf
    REAL(8), DIMENSION(NC_tt,num_M) :: H_hat
    !
    ! Beginning execution
    !
    ! Evaluate RV logdensity
    ! ===============================================================================
    !
    ! Computing moments
    E_lRV  = eta(1)+lf_tt
    sd_lRV = EXP(chi(1))
    !
    ! Computing log density
    zz_lRV = (lRV_tt-E_lRV)/sd_lRV
    CALL eval_logpdf_N01(num_M,zz_lRV,logpdf_RV)
    logpdf_RV = logpdf_RV-LOG(sd_lRV)
    !
    ! Check for NaNs
    IF (ANY(logpdf_RV .NE. logpdf_RV)) THEN
        !
        numerical_error_flag = 'pRV'
        RETURN
        !
    END IF
    !
    ! Evaluate H logdensity
    ! ===============================================================================
    !
    ! Computing standard errors
    X_tt = X(:NC_tt,tt)
    tau_tt = tau(:NC_tt,tt)
    sd_lH = xi(1)+ &
        xi(2)*X_tt+xi(3)*X_tt**2+ &
        xi(4)*tau_tt/365+xi(5)*(tau_tt/365)**2+ &
        xi(6)*X_tt*(tau_tt/365)
    sd_lH = EXP(sd_lH)
    !
    IF (ANY(sd_lH .EQ. 0.d0)) THEN
        !
        numerical_error_flag = 'vme'
        RETURN
        !
    END IF
    !
    ! Evaluate fitted normalized option prices
    H_hat = 0.d0
    CALL option_prices_lin_date(tt,NC_tt,f_tt,Hgrid,H_hat)
    WHERE (H_hat .LE. 0.d0)
        H_hat = bulletproof
    END WHERE
    H_hat = LOG(H_hat)
	!
    ! Evaluate log density
    logpdf_H = 0.d0
    lH = LOG(H(:NC_tt,tt))
	DO ss = 1, num_M
        !
        E_lH         = -0.5d0*sd_lH**2+H_hat(:,ss)
        zz_C         = (lH-E_lH)/sd_lH
		CALL eval_logpdf_N01(NC_tt,zz_C,lpdf)
		lpdf = lpdf-LOG(sd_lH)
		logpdf_H(ss) = SUM(lpdf)		
        !
	END DO	
    !
    ! Check for NaNs
    IF (ANY(logpdf_H .NE. logpdf_H)) THEN
        !
        numerical_error_flag = 'pdH'
        RETURN
        !
    END IF
    !
    ! Evaluate y logdensity
    ! ===============================================================================
    !
    ! Computing moments
    !
    E_y  = rmq+lambda*f_tt
    sd_y = SQRT(f_tt)
    !
    ! Computing log density
    !
    zz_y = (y_tt-E_y)/sd_y
    CALL eval_logpdf_N01(num_M,zz_y,logpdf_y)
    logpdf_y = logpdf_y-LOG(sd_y)
    !
    ! Check for NaNs
    IF (ANY(logpdf_y .NE. logpdf_y)) THEN
        !
        numerical_error_flag = 'pdy'
        RETURN
        !
    END IF
	!
    ! Ending execution and returning control
    !
    END SUBROUTINE eval_g_tt
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE psi_params ( theta, distrib, lambda, nu, mu, beta )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)
    CHARACTER(len=1), INTENT(IN) :: distrib
    REAL(8), INTENT(OUT) :: lambda
    REAL(8), INTENT(OUT) :: nu
    REAL(8), INTENT(OUT) :: mu
    REAL(8), INTENT(OUT) :: beta
    !
    ! Declaring local variables
    !
    REAL(8) :: rho, rhoQ, phi
    !
    ! Beginning execution
    !
    lambda = theta(1)
    nu     = EXP(theta(2))
    mu     = EXP(theta(3))/V_norm
    rho    = 0.5d0*(dueinvpi*ATAN(theta(4))+1.d0)
    beta   = rho/mu
    !
    IF (distrib .EQ. 'Q') THEN
        !
        lambda  = -0.5d0
        rhoQ    = 0.5d0*(dueinvpi*ATAN(theta(5))+1.d0)
        phi     = SQRT(rhoQ/rho)
        mu      = mu*phi
        beta    = beta*phi
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE psi_params
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE sigma_params ( theta, eta, chi, xi )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)
    REAL(8), INTENT(OUT) :: eta(num_eta)
    REAL(8), INTENT(OUT) :: chi(num_chi)
    REAL(8), INTENT(OUT) :: xi(num_xi)
    !
    ! Beginning execution
    !
    eta = theta(num_psi+num_phi+1:num_psi+num_phi+num_eta)
    chi = theta(num_psi+num_phi+num_eta+1:num_psi+num_phi+num_eta+num_chi)
    xi = theta(num_psi+num_phi+num_eta+num_chi+1:)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE sigma_params
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION Tbar ( L, chiM, x0, nu, beta, mu )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: L
    REAL(8), INTENT(IN) :: chiM(num_M)
    REAL(8), INTENT(IN) :: x0
    REAL(8), INTENT(IN) :: nu, beta, mu
    !
    ! Declaring function's type
    !
    REAL(8) :: Tbar(num_M,L)
    !
    ! Declaring local variables
    !
    REAL(8) :: betax0
    REAL(8) :: Tbar_tmp(num_M,max_num_L)
    !
    ! Beginning execution
    !
    betax0 = beta*x0
    IF (L .GE. 1) Tbar_tmp(:,1) = mu*(nu+betax0)/chiM(num_M)
    IF (L .GE. 2) Tbar_tmp(:,2) = mu**2*(2.d0*betax0+nu)/chiM(num_M)**2
    Tbar = Tbar_tmp(:,:L)
    !
    ! Ending execution and returning control
    !
    END FUNCTION Tbar
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION T ( L, chiM, x0, nu, beta, mu )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: L
    REAL(8), INTENT(IN) :: chiM(num_M)
    REAL(8), INTENT(IN) :: x0
    REAL(8), INTENT(IN) :: nu, beta, mu
    !
    ! Declaring function's type
    !
    REAL(8) :: T(num_M,L)
    !
    ! Declaring local variables
    !
    REAL(8) :: betax0
    REAL(8) :: T_tmp(num_M,max_num_L)
    !
    ! Beginning execution
    !
    betax0 = beta*x0
    IF (L .GE. 1) T_tmp(:,1) = chiM/chiM(num_M)
    IF (L .GE. 2) T_tmp(:,2) = (chiM-mu*(nu+betax0))**2/chiM(num_M)**2
    T = T_tmp(:,:L)
    !
    ! Ending execution and returning control
    !
    END FUNCTION T
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE eval_logpdf_N01 ( n, x, logpdf )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x(n)
    REAL(8), INTENT(OUT) :: logpdf(n)
    !
    ! Beginning execution
    !
    logpdf = log1oversqrt2pi-x**2/2.d0
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE eval_logpdf_N01
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE pdf_gamma ( shape, scale, n, x, dolog, output ) 
    !
	! Evaluate the density of a Gamma variable with parameters shape and scale at x
	! gamma(x; shape, scale) \propto x^{shape-1}*\exp(-x/scale)
	!
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: shape
    REAL(8), INTENT(IN) :: scale
	INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x(n)
    CHARACTER(len = 3), INTENT(IN) :: dolog
    REAL(8), INTENT(OUT) :: output(n)
    !
    ! Declaring local variables
    !
    INTEGER :: ifault
    REAL(8) :: lgamma
    !
    ! Declaring local variables
    !
    CALL alogam(shape,ifault,lgamma)
    output = -x/scale+(shape-1.d0)*LOG(x)-shape*LOG(scale)-lgamma
    IF (dolog .NE. 'log') output = EXP(output)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE pdf_gamma
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE cdf_gamma ( shape, scale, n, x, dolog, output ) 
    !
	! Evaluate the cdf of a Gamma variable with parameters shape and scale at x
	! gamma(x; shape, scale) \propto x^{shape-1}*\exp(-x/scale)
	!
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: shape
    REAL(8), INTENT(IN) :: scale
	INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x(n)
    CHARACTER(len = 3), INTENT(IN) :: dolog
    REAL(8), INTENT(OUT) :: output(n)
    !
    ! Declaring local variables
    !
    INTEGER :: i, status
    REAL(8) :: q, bound
    !
    ! Declaring local variables
    !
    DO i = 1, n
        !
        CALL cdfgam(1,output(i),q,x(i),shape,1.d0/scale,status,bound)
        !
    END DO
    IF (dolog .EQ. 'log') output = LOG(output)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE cdf_gamma
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE pdf_ncgamma ( nu, mu, beta, f_t, f_tm1, num_sim, pdf ) 
!    !
!    ! p(f_t|Y_{t-1})
!	! Evaluate the density of a noncentral Gamma variable with parameters mu and nu at ft
!	! p(f_t|f_{t-1}; \mu, \nu, \beta) where 
!	! beta*f_{t-1} is the Poisson intensity
!	! mu is the shape parameter
!	! nu are the degrees of freedom
!	!
!	! See Gourieroux and Jasiak J. of Forecasting p. 131 for details
!	!
!    IMPLICIT NONE
!    !
!    ! Declaring dummy variables
!    !
!    REAL(8), INTENT(IN) :: nu, mu, beta
!    REAL(8), INTENT(IN) :: f_t(1)
!    INTEGER, INTENT(IN) :: num_sim
!    REAL(8), INTENT(IN), DIMENSION(num_sim) :: f_tm1
!    REAL(8), INTENT(OUT) :: pdf(num_sim)
!    !
!    ! Declaring local variables
!    !
!    REAL(8), DIMENSION(num_sim) :: z0
!    REAL(8), DIMENSION(0:pdf_gamma_trunc,num_sim) :: mat_pdf_poisson
!    REAL(8), DIMENSION(0:pdf_gamma_trunc) :: vec_pdf_gamma
!    INTEGER :: kk 	
!    REAL(8) :: logpdf(1)
!    !
!    ! Beginning execution
!    !
!    z0 = f_tm1*beta
!	!
!    vec_pdf_gamma = 0.d0
!    CALL pdf_gamma(1.d0/mu,nu,1,f_t,'log',logpdf)
!    vec_pdf_gamma(0) = EXP(SUM(logpdf))
!    mat_pdf_poisson = 0.d0
!    mat_pdf_poisson(0,:) = EXP(-z0)
!    pdf = vec_pdf_gamma(0)*mat_pdf_poisson(0,:)
!    !
!	DO kk = 1, pdf_gamma_trunc
!        !
!        CALL pdf_gamma(1.d0/mu,nu+kk,1,f_t,'log',logpdf)
!        vec_pdf_gamma(kk) = EXP(SUM(logpdf))
!        mat_pdf_poisson(kk,:) = mat_pdf_poisson(kk-1,:)*z0/kk
!        pdf = pdf+vec_pdf_gamma(kk)*mat_pdf_poisson(kk,:)
!        !
!	END DO
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE pdf_ncgamma
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE alogam ( x, ifault, lgamma )
    !
    !*****************************************************************************80
    !
    !! ALOGAM computes the logarithm of the Gamma function.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 March 1999
    !
    !  Author:
    !
    !    Malcolm Pike, David Hill
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    Malcolm Pike, David Hill,
    !    Algorithm 291: 
    !    Logarithm of Gamma Function,
    !    Communications of the ACM,
    !    Volume 9, Number 9, September 1966, page 684.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the Gamma function.
    !    X should be greater than 0.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error flag.
    !    0, no error.
    !    1, X <= 0.
    !
    !    Output, real ( kind = 8 ) ALOGAM, the logarithm of the Gamma 
    !    function of X.
    !
    IMPLICIT NONE
    !
    REAL(8), INTENT(IN) :: x
    INTEGER, INTENT(OUT) :: ifault
    REAL(8), INTENT(OUT) :: lgamma
    !
    REAL(8) :: f
    REAL(8) :: y
    REAL(8) :: z
    !
    IF (x .LE. 0.0d0) THEN
        !
        ifault = 1
        lgamma = 0.0d0
        RETURN
        !
    END IF
    !
    ifault = 0
    y = x
    !
    IF (x .LT. 7.0d0) THEN
        !
        f = 1.0d0
        z = y
        !
        DO WHILE ( z < 7.0d0 )
            !
            f = f * z
            z = z + 1.0d0
            !
        END DO
        !
        y = z
        f = -LOG(f)
        !
    ELSE
        !
        f = 0.0d0
        !
    END IF
    !
    z = 1.0d0 / y / y
    !
    lgamma = f + ( y - 0.5d0 ) * LOG ( y ) - y &
    + 0.918938533204673d0 + &
    ((( &
    - 0.000595238095238d0   * z &
    + 0.000793650793651d0 ) * z &
    - 0.002777777777778d0 ) * z &
    + 0.083333333333333d0 ) / y
    !
    RETURN
    !
    END SUBROUTINE alogam
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE loglikelihood

