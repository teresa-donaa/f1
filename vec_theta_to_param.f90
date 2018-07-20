SUBROUTINE vec_theta_to_param ( theta, param )
!
USE constants
USE observations
USE loglikelihood
!
! Transforms theta to real parameters 
!
IMPLICIT NONE
!
! Declaring dummy variables
!
REAL(8), INTENT(IN) :: theta(num_theta)
REAL(8), INTENT(OUT) :: param(num_param)
!
! Declaring local variables
!
REAL(8) :: lambda, nu, mu, beta, rho
REAL(8) :: lambdaQ, nuQ, muQ, betaQ, rhoQ
REAL(8) :: phi, ystar, delta1, delta2
REAL(8) :: eta(num_eta), chi(num_chi), xi(num_xi), sigma, etastar
INTEGER :: il, iu, ir
! 
! Beginning execution
!
! P parameters
!
CALL PSI_PARAMS(theta,'P',lambda,nu,mu,beta)
rho = mu*beta
!
! Q parameters
!
CALL PSI_PARAMS(theta,'Q',lambdaQ,nuQ,muQ,betaQ)
rhoQ   = muQ*betaQ
!
! Risk premia parameters
!
phi    = SQRT(rhoQ/rho)
ystar  = (1.d0-1.d0/phi)/mu
delta2 = lambda+0.5d0
delta1 = -ystar-delta2*lambda+0.5d0*delta2**2
!
! Measurement errors variance parameters
!
CALL SIGMA_PARAMS(theta,eta,chi,xi)
sigma = EXP(chi(1))
etastar = EXP(eta(1)+sigma**2/2.d0)
!
! Defining PARAM
!
! P parameters
!
param(1)  = lambda
param(2)  = nu
param(3)  = mu
param(4)  = beta
param(5)  = rho
!
! Q parameters
!
param(6)  = lambdaQ
param(7)  = nuQ
param(8) = muQ
param(9) = betaQ
param(10) = rhoQ
!
! Risk premia parameters
!
param(11) = phi
param(12) = ystar
param(13) = delta1
param(14) = delta2
!
! Other parameters
!
param(15) = eta(1)
param(16) = etastar
param(17) = chi(1)
param(18) = sigma
!
! xi
!
iu = 18
il = iu
DO ir = 1, num_xi
    param(il+ir) = xi(ir)
END DO    
iu = iu+num_xi
!
! End execution and returning control
!
END SUBROUTINE vec_theta_to_param
