MODULE option_prices
!
USE constants
USE observations
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE option_prices_bsp_setup ( )
    !
    ! Creates grid coordinates
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: i
    REAL(8) :: f(num_T), meanf, sdf
    !
    ! Beginning execution
    !
    ! 1) Moneyness
    ! Equally spaced grid
    !
    Xgrid(1) = min_x
    Xgrid(num_grid_x) = max_x
    xstep = (Xgrid(num_grid_x)-Xgrid(1))/(num_grid_x-1)
    DO i = 2, num_grid_x-1
        Xgrid(i) = Xgrid(i-1)+Xstep
    END DO
    !
    ! 2) Time to maturity
    ! Equally spaced grid
    !
    taugrid(1) = min_tau
    taugrid(num_grid_tau) = max_tau
    taustep = (taugrid(num_grid_tau)-taugrid(1))/(num_grid_tau-1)
    DO i = 2, num_grid_tau-1
        taugrid(i) = taugrid(i-1)+taustep
    END DO
    !
    ! 3) Volatility
    ! Grid bounds: min = min_f_frac*MINVAL(RV), max = max_f
    !
    IF (switch_min_f_level .EQ. 1) fgrid(1) = min_f
    IF (switch_min_f_level .EQ. 0) fgrid(1) = min_f_frac*MINVAL(RV)
    fgrid(num_grid_f) = max_f
    !
    IF (switch_log_fgrid .EQ. 0) THEN
        !
        fstep = (fgrid(num_grid_f)-fgrid(1))/(num_grid_f-1)
        DO i = 2, num_grid_f-1
            fgrid(i) = fgrid(i-1)+fstep
        END DO
        !
    END IF
    IF (switch_log_fgrid .EQ. 1) THEN
        !
        fstep = (LOG(fgrid(num_grid_f))-LOG(fgrid(1)))/(num_grid_f-1)
        DO i = 2, num_grid_f-1
            fgrid(i) = EXP(LOG(fgrid(i-1))+fstep)
        END DO
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE option_prices_bsp_setup 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE option_prices_ffo_setup ( NN, aint, bint )
    !
    ! Setup routine for fast Fang and Oosterlee COS method
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: NN           ! Order of approximation
    REAL(8), INTENT(IN) :: aint         ! Lower integration bound
    REAL(8), INTENT(IN) :: bint         ! Upper integration bound
    !
    ! Declaring local variables
    !
    REAL(8) :: cint, dint, tmp
    REAL(8), DIMENSION(NN) :: u, tmp1, tmp2, CHI, tmp3, PSI, CHIPSIdiff
    INTEGER :: kvec(NN), i, j
    !
    ! Beginning execution
    !
    cint = 0.d0
    dint = bint
    !
    ! Approximation indexes
    !
    kvec = (/ (i, i = 0, NN-1) /)
    u = (kvec*greek_pi)/(bint-aint)
    z_ffo = CMPLX(0.d0,u,8)
    !
    ! Chi_k
    !
    tmp1 = COS((dint-aint)*u)*EXP(dint)-COS((cint-aint)*u)*EXP(cint)
    tmp2 = u*SIN((dint-aint)*u)*EXP(dint)-u*SIN((cint-aint)*u)*EXP(cint)
    CHI = (1.d0/(1.d0+u**2))*(tmp1+tmp2)
    !
    ! Psi_k
    !
    tmp3(1) = dint-cint
    tmp3(2:) = (SIN((dint-aint)*u(2:))-SIN((cint-aint)*u(2:)))*(1.d0/u(2:))
    PSI = tmp3
    !
    ! ChiPsidiff and vk [a (NNo x NKo) matrix]
    !
    CHIPSIdiff = CHI-PSI
!    vk_ffo = SPREAD(CHIPSIdiff*(2.d0/(bint-aint)),DIM=2,NCOPIES=num_max_C)
    DO i = 1, num_max_C
        !
        vk_ffo(:,i) = CHIPSIdiff*(2.d0/(bint-aint))
        !
    END DO
    !
    ! Computing component of option prices that does not depend on theta
    !
    DO i = 1, NN
        !
        DO j = 1, num_grid_X
            !
            tmp = u(i)*(Xgrid(j)-a_fs)
            cosdx(i,j) = COS(tmp)
            sindx(i,j) = SIN(tmp)
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE option_prices_ffo_setup 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE option_prices_ffo_loglik ( NN, lambdaQ, nuQ, muQ, betaQ, &
        Acf_ffo, Bcf_ffo )
    !
    ! Computes characteristic function coefficients for the fast COS method
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: NN
    REAL(8), INTENT(IN) :: lambdaQ, nuQ, muQ, betaQ
    COMPLEX(8), INTENT(OUT), DIMENSION(num_N,min_tau:max_tau) :: Acf_ffo, Bcf_ffo
    !
    ! Declaring local variables
    !
    INTEGER :: i, j, iflag
    COMPLEX(8), DIMENSION(NN) :: xx, Acf, Bcf
    !
    ! Beginning execution
    !
    ! Loop over the approximation order
    !
    Acf = (0.d0, 0.d0)
    Bcf = (0.d0, 0.d0)
    !
    DO i = 1, max_tau
        !
        xx = Bcf+z_ffo*lambdaQ+z_ffo**2/2.d0
        Acf = Acf+z_ffo*rmq-nuQ*LOG(1.d0-xx*muQ)
        Bcf = (muQ*xx)/(1.d0-muQ*xx)*betaQ
        !
        IF (i .GE. min_tau) THEN
            !
            Acf_ffo(:,i) = Acf
            Bcf_ffo(:,i) = Bcf
            !
        END IF
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE option_prices_ffo_loglik 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE gen_option_price_ffo ( NK, X, lambdaQ, nuQ, muQ, betaQ, &
        f, tau, NN, aint, bint, Acf_ffo, Bcf_ffo, opt_price )
    !
    ! Computes normalized option prices date by date using the fast COS method
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: NK           ! Number of strike prices
    REAL(8), INTENT(IN) :: X(NK)
    REAL(8), INTENT(IN) :: lambdaQ, nuQ, muQ, betaQ
    REAL(8), INTENT(IN) :: f            ! Latent volatility
    INTEGER, INTENT(IN) :: tau          ! Maturity
    INTEGER, INTENT(IN) :: NN           ! Order of approximation
    REAL(8), INTENT(IN) :: aint         ! Lower integration bound
    REAL(8), INTENT(IN) :: bint         ! Upper integration bound
    COMPLEX(8), INTENT(IN), DIMENSION(num_N,min_tau:max_tau) :: Acf_ffo, Bcf_ffo
    REAL(8), INTENT(OUT) :: opt_price(NK)
    !
    ! Declaring local variables
    !
    INTEGER :: i, j
    REAL(8) :: lbH, ubH
    REAL(8), DIMENSION(NN,NK) :: fk
    COMPLEX(8), DIMENSION(NN) :: cf
    REAL(8), DIMENSION(NN) :: cfr, cfi
    !
    ! Beginning execution
    !
    ! Characteristic function of log(S_T)
    !
    cf = EXP(Acf_ffo(:,tau)+Bcf_ffo(:,tau)*f)
    cfr = DBLE(cf)
    cfi = IMAG(cf)
    !
    ! Option price
    !
    DO i = 1, NN
        !
        DO j = 1, NK
            !
            fk(i,j) = cfr(i)*cosdx(i,j)-cfi(i)*sindx(i,j)
            !
        END DO
        !
    END DO
    fk(1,:) = 0.5d0*fk(1,:)
    !
    opt_price = 0.d0
    DO i = 1, NK
        !
        DO j = 1, NN
            !
            opt_price(i) = opt_price(i)+fk(j,i)*vk_ffo(j,i)
            !
        END DO
        !
    END DO
    !
    ! Imposing no arbitrage bounds
    !
    DO i = 1, NK
        !
        ubH = EXP(X(i)+rmq*tau)
        lbH = MAX(0.d0,EXP(X(i)+rmq*tau)-1.d0)
        IF (opt_price(i) .GT. ubH) opt_price(i) = ubH
        IF (opt_price(i) .LT. lbH) opt_price(i) = lbH
        !
    END DO
    !
    ! End execution and returning control
    !
    END SUBROUTINE gen_option_price_ffo
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE option_prices_lin_setup ( )
    !
    ! Creates grid coordinates
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: tt, firstX, lastX, NT, NK, taucc, ktau, kx, htau, hx
    REAL(8) :: xcc
    !
    ! Beginning execution
    !
    ! Beginning loop over dates
    !
    jx_mat = 0
    jtau_mat = 0
    DO tt = 1, num_T
        !
        lastX = 0
        NT = num_tau(tt)
        !
        ! Beginning loop over maturities
        !
        DO ktau = 1, NT
            !
            ! Extracting indexes of first and last moneynesses at maturity taucc
            !
            firstX = lastX+1
            NK = num_K(firstX,tt)
            lastX = lastX+NK
            !
            ! Extracting time to maturity
            !
            taucc = tau(firstX,tt)
            !
            ! Compute jtau and dtau 
            !
            htau = (taucc-min_tau)/taustep+1           ! Integer division
            IF (htau .LT. (num_grid_tau-1)) THEN
                !
                jtau_mat(firstX:lastX,tt) = htau
                dtau_mat(firstX:lastX,tt) = DBLE(MOD(taucc-min_tau,taustep))/taustep
                !
            ELSE IF (htau .GE. (num_grid_tau-1)) THEN
                !
                jtau_mat(firstX:lastX,tt) = num_grid_tau-1
                dtau_mat(firstX:lastX,tt) = &
                    DBLE((taucc-taugrid(num_grid_tau-1)))/(taugrid(num_grid_tau)-taugrid(num_grid_tau-1))
                !
            END IF
            !
            ! Beginning loop over moneynesses
            !
            DO kx = firstX, lastX
                !
                ! Extracting strike price xcc
                !
                xcc = X(kx,tt)
                !
                ! Compute jx and dx
                !
                hx = INT((xcc-min_x)/xstep+1)           ! Integer division
                IF (hx .LT. num_grid_x) THEN
                    !
                    jx_mat(kx,tt) = hx
                    dx_mat(kx,tt) = MOD(xcc-min_x,xstep)/xstep
                    !
                ELSE IF (hx .EQ. num_grid_x) THEN
                    !
                    jx_mat(kx,tt) = hx-1
                    dx_mat(kx,tt) = 1.d0
                    !
                END IF
                !
            END DO
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE option_prices_lin_setup 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE option_prices_lin_loglik ( lambdaQ, nuQ, muQ, betaQ, &
        Hgrid, Acf_ffo, Bcf_ffo )
    !
    ! Computes normalized option prices on the grid
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: lambdaQ, nuQ, muQ, betaQ
    REAL(8), INTENT(OUT), DIMENSION(num_grid_x,num_grid_tau,num_grid_f) :: Hgrid 
    COMPLEX(8), INTENT(IN), DIMENSION(num_N,min_tau:max_tau) :: Acf_ffo, Bcf_ffo
    !
    ! Declaring local variables
    !
    INTEGER :: i, j, k
    !
    ! Beginning execution
    !
    ! Computing normalized option prices on the grid
    !
    !$omp parallel do private(j)
    DO i = 1, num_grid_tau
        !
        DO j = 1, num_grid_f
            !
            CALL gen_option_price_ffo(num_grid_X,Xgrid,lambdaQ,nuQ,muQ,betaQ, &
                fgrid(j),taugrid(i),num_N,a_fs,b_fs,Acf_ffo,Bcf_ffo, &
                Hgrid(:,i,j))
            !
        END DO
        !
    END DO
    !$omp end parallel do
    !
    ! Checking for NaNs
    !
    IF (ANY(Hgrid .NE. Hgrid)) THEN
        !
        numerical_error_flag = 'Hgr'
        RETURN
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE option_prices_lin_loglik 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE option_prices_lin_date ( tt, NC, f1, Hgrid, H_hat )
    !
    ! Computes normalized option prices date by date
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: tt, NC
    REAL(8), INTENT(IN) :: f1(num_M)
    REAL(8), INTENT(IN), DIMENSION(num_grid_x,num_grid_tau,num_grid_f) :: Hgrid 
    REAL(8), INTENT(OUT) :: H_hat(NC,num_M)
    !
    ! Declaring local variables
    !
    INTEGER :: kf, hf, jf(num_M), kC, jxC, jtauC, jfC, iw, iw1, iw2
    REAL(8) :: fcc, df(num_M), dxC, dtauC, dfC, Wmat(8,4), prodWmat(8)
    REAL(8) :: lbH, ubH, taucc, Xcc
    !
    ! Beginning execution
    !
    ! Beginning loop over volatilities to compute f indexes and weights
    !
    jf = 0
    DO kf = 1, num_M
        !
        ! Extracting volatility fC
        !
        fcc = f1(kf)
        !
        ! Compute jf and df 
        !
        IF (switch_log_fgrid .EQ. 0) THEN
            !
            hf = INT((fcc-fgrid(1))/fstep)+1                    ! Integer division
            !
            IF (hf .LE. 0) THEN                                 ! fcc smaller than the first point in the grid
                !
                jf(kf) = 1
                df(kf) = (fcc-fgrid(1))/fstep
                !
            END IF
            IF ((hf .GT. 0) .AND. (hf .LT. (num_grid_f-1))) THEN
                !
                jf(kf) = hf
                df(kf) = MOD(fcc-fgrid(1),fstep)/fstep
                !
            END IF
            IF (fcc .GE. fgrid(num_grid_f-1)) THEN              ! Extrapolation
                !
                jf(kf) = num_grid_f-1
                df(kf) = (fcc-fgrid(num_grid_f-1))/(fgrid(num_grid_f)-fgrid(num_grid_f-1))
                !
            END IF
            !
        END IF
        IF (switch_log_fgrid .EQ. 1) THEN
            !
            hf = INT((LOG(fcc)-LOG(fgrid(1)))/fstep)+1                    ! Integer division
            !
            IF (hf .LE. 0) THEN                                 ! fcc smaller than the first point in the grid
                !
                jf(kf) = 1
                df(kf) = (fcc-fgrid(1))/(fgrid(2)-fgrid(1))
                !
            END IF
            IF ((hf .GT. 0) .AND. (hf .LT. (num_grid_f-1))) THEN
                !
                jf(kf) = hf
                df(kf) = (fcc-fgrid(hf))/(fgrid(hf+1)-fgrid(hf))
                !
            END IF
            IF (fcc .GE. fgrid(num_grid_f-1)) THEN              ! Extrapolation
                !
                jf(kf) = num_grid_f-1
                df(kf) = (fcc-fgrid(num_grid_f-1))/(fgrid(num_grid_f)-fgrid(num_grid_f-1))
                !
            END IF
            !
        END IF
        !
    END DO
    !
    ! Beginning loop over options of date tt
    !
    H_hat = 0.d0
    DO kC = 1, NC
        !
        ! Retrieving x and tau indexes and weights
        !
        jxC = jx_mat(kC,tt)
        jtauC = jtau_mat(kC,tt)
        dxC = dx_mat(kC,tt)
        dtauC = dtau_mat(kC,tt)
        !
        ! Computing first two columns of the Wmat matrix
        !
        Wmat(:,1) = (1.d0-dxC)*(1.d0-Lmat(:,1))+dxC*Lmat(:,1)
        Wmat(:,2) = (1.d0-dtauC)*(1.d0-Lmat(:,2))+dtauC*Lmat(:,2)
        !
        ! Computing no arbitrage bounds
        !
        taucc = DBLE(tau(kC,tt))
        Xcc = X(kC,tt)
        ubH = EXP(Xcc+rmq*taucc)
        lbH = MAX(0.d0,EXP(Xcc+rmq*taucc)-1.d0)
        !
        ! Beginning loop over volatilities
        !
        DO kf = 1, num_M
            !
            ! Retrieving f and l indexes and weights
            !
            jfC = jf(kf)
            dfC = df(kf)
            !
            ! Computing third and fourth columns of the Wmat matrix
            !
            Wmat(:,3) = (1.d0-dfC)*(1.d0-Lmat(:,3))+dfC*Lmat(:,3)
            !
            ! Computing fourth line of weight matrix
            !
            DO iw = 1, 8
                !
                Wmat(iw,4) = Hgrid(jxC+Lmat(iw,1),jtauC+Lmat(iw,2),jfC+Lmat(iw,3))
                !
            END DO
            !
            ! Computing predicted price
            ! 
            H_hat(kC,kf) = 0.d0
            DO iw1 = 1, 8
                ! 
                prodWmat(iw1) = 1.d0
                DO iw2 = 1, 4
                    !
                    prodWmat(iw1) = prodWmat(iw1)*Wmat(iw1,iw2)
                    !
                END DO
                H_hat(kC,kf) = H_hat(kC,kf)+prodWmat(iw1)
                !
            END DO
            !
            ! Imposing no arbitrage bounds
            !
            IF (H_hat(kC,kf) .GT. ubH) H_hat(kC,kf) = ubH
            IF (H_hat(kC,kf) .LT. lbH) H_hat(kC,kf) = lbH
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE option_prices_lin_date
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE option_prices
