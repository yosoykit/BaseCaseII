######################################################################################################
##  EXAMPLE PARAMETERS PROVIDED. Utilized for Screening paper 2014, fit to SEER EAC incidence data ###
##                			  all MALE PARAMETERS: AGE-COHORT MODEL									##
##  			birth cohort variables for sigmoidal g_P and g_M: See Supplementary 				##
######################################################################################################

param_list <- function(X, nu, alphaP, alphaM, rho, gP,gM,mu0,mu1,mu2){
	betaM <- alphaM-gM-rho
	mu2eff <- mu2*(1-betaM/alphaM)				## effective rate to pre-clinical cancer 
	betaeff<- alphaP-gP-mu2eff
	betaP<- alphaP-gP-mu2
	peff <- (1/2)*(-alphaP+betaeff + mu2eff - sqrt((alphaP+betaeff+mu2eff)^2-4*alphaP*betaeff))
	qeff <- (1/2)*(-alphaP+betaeff + mu2eff + sqrt((alphaP+betaeff+mu2eff)^2-4*alphaP*betaeff))
	pP <- (1/2)*(-alphaP+betaP + mu2 - sqrt((alphaP+betaP+mu2)^2-4*alphaP*betaP))
	qP <- (1/2)*(-alphaP+betaP + mu2 + sqrt((alphaP+betaP+mu2)^2-4*alphaP*betaP))
	pM <- (1/2)*(-alphaM+betaM + rho - sqrt((alphaM+betaM+rho)^2-4*alphaM*betaM))
	qM <- (1/2)*(-alphaM+betaM + rho + sqrt((alphaM+betaM+rho)^2-4*alphaM*betaM))
	zeta=-qM/pM 
	tlag <- log((zeta/(1+zeta)))/pM     				## exact T_2
	return(list(X=X,nu=nu,mu0=mu0,mu1=mu1,mu2=mu2, rho=rho,alphaP=alphaP,alphaM=alphaM, betaP=betaP,betaM=betaM, gP=gP,gM=gM, tlag=tlag, mu2eff=mu2eff,betaeff=betaeff,pP=pP, qP=qP,pM=pM,qM=qM,peff=peff,qeff=qeff))
}

am_g0<- 0.99061E-01        
am_g1<- 0.50880E+00       
am_g2<- 0.53768E-01    
am_trefc<- 0.19125E+04 
am_gC0<- 0.75000E+00     

gerd1<- 0.00061422    
gerd2<- 0.0070447     
gerd3<- 26.002   
nu0 <-amnu0<- 0.36494E-03 

ammu0<- ammu1 <-  0.79942E-03
ammu2<-    0.45439E-04 

kstem<- 4    				## stem cells/crypt
X<- 250000*kstem 			## total stem cells in a 5 cm BE segment assuming 250,000 crypts/5cm

amrho<- 1.0000E-09       	## detection rate of clinical cancers
amalphaP0<- 10				## premalignant cell division rate
amalphaM0 <- 150			## pre-clinical malignant cell division rate
birth_cohort<- 1950

am_gP<- am_g0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
am_gM<- am_gC0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
#am_gM=am_gM*1.25
amalphaP<- amalphaP0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
amalphaM<- amalphaM0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
allmales<- param_list(X, amnu0, amalphaP, amalphaM, amrho, am_gP,am_gM, ammu0, ammu1, ammu2)
amparams<- c(allmales$mu0,allmales$mu1,allmales$betaP,allmales$mu2,allmales$alphaP,allmales$betaM,allmales$rho,allmales$alphaM,allmales$pM,allmales$qM,allmales$tlag)

