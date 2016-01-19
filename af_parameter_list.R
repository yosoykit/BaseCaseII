## BASE CASE 2 PARAMETERS:

## ALL FEMALE PARAMETERS BIRTH COHORT fixed at 1950- AGE-COHORT MODEL
## birth cohort varaibles for sigmoidal g_P and g_M
param_list <- function(X, nu, alphaP, alphaM, rho, gP,gM,mu0,mu1,mu2){
	betaM <- alphaM-gM-rho
	mu2eff <- mu2*(1-betaM/alphaM)						## effective rate to pre-clinical cancer 
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

af_g0<- 0.12323E+00   
af_g1<- 0.63999E+00   
af_g2<- 0.29781E-01    
af_trefc<- 0.19453E+04   
af_gC0<- 0.75000E+00     

gerd1<- 0.18012E-03     
gerd2<- 0.71326E-02     
gerd3<- 0.25808E+02  
nu0<-afnu0<- 0.74828E-04 

## ADJUSTED VALUES FOR BASE CASE 2 need a factor change for BE segment length distribution, mean length of beta(.5,4)
crypts_per_avBE<- (1+15*(.5/4.5))*10*75/15*1000
afmu0<- 0.70504E-03*250000/crypts_per_avBE
afmu1 <-  0.70504E-03
afmu2<-   0.68903E-04 

kstem<- 4    				## stem cells/crypt
#X<- 250000*kstem 			## total stem cells in a 5 cm BE segment assuming 250,000 crypts/5cm
X<- crypts_per_avBE*kstem 			

afrho<- 1.0000E-09       ## detection rate of clinical cancers
afalphaP0<- 10				
afalphaM0 <- 150			## pre-clinical malignant cell division rate

birth_cohort<- 1950
#birth_cohort = 1925  

af_gP<- af_g0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
af_gM<- af_gC0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
afalphaP<- afalphaP0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
afalphaM<- afalphaM0*(af_g1+2/(1+exp(-af_g2*(birth_cohort-af_trefc))))
allfemales<- param_list(X, afnu0, afalphaP, afalphaM, afrho, af_gP,af_gM, afmu0, afmu1, afmu2)
afparams<- c(allfemales$mu0,allfemales$mu1,allfemales$betaP,allfemales$mu2,allfemales$alphaP,allfemales$betaM,allfemales$rho,allfemales$alphaM,allfemales$pM,allfemales$qM,allfemales$tlag)

