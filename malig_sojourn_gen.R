 #amparams<- c(allmales$mu0,allmales$mu1,allmales$betaP,allmales$mu2,allmales$alphaP,allmales$betaM,allmales$rho,allmales$alphaM,allmales$pM,allmales$qM,allmales$tlag)


## Samples sojourn times of EAC clinical detection beginning with 1 malignant cell at time tau, conditional on non-extinction
 malig_sojourn <- function(params, tau=0){
 	alphaM = params[8]
 	betaM = params[6]
 	pM = params[9]
 	qM = params[10]
 	# conditional on non-extinction
 	#r<-runif(1,0,(1-betaM/alphaM))
 	r <- runif(1,0,1)
 	if ( r >= (1-betaM/alphaM)){
 		u =Inf
 		#print(r)
 	}
 	else{
 	u <- -(1/(pM-qM))*log((pM*(alphaM*r+qM))/(qM*(alphaM*r+pM)))
 	}
 	return(list(sojourn = u, age_EAC = u+tau))
 }


  malig_sojourn_m <- function(params, tau=0,mtotal){
 	alphaM = params[8]
 	betaM = params[6]
 	pM = params[9]
 	qM = params[10]
 	# conditional on non-extinction
 	#r<-runif(1,0,(1-betaM/alphaM))
 	print((1-(betaM/alphaM)^mtotal))
 	r <- runif(1,0,1)
 	if( r >= (1-(betaM/alphaM)^mtotal)){
 		u =Inf
 		#print(r)
 	}
 	else{
 	u <- -(1/(pM-qM))*log((pM*(alphaM*(-r+1)^(1/mtotal)-alphaM-qM))/(qM*(alphaM*(-r+1)^(1/mtotal)-alphaM-pM)))
 	}
 	return(list(sojourn = u, age_EAC = u+tau))
 }

# sims = 50000
# indiv_test= rep(0,sims)
# total_test = rep(0,sims)
# total_M_cells = 100000

# for ( k in 1:sims){
# 	temp = malig_sojourn(params,tau=0)
# 	M_to_EAC_time = temp$age_EAC
# 	for (m in 1:(total_M_cells-1)){
# 		temp = malig_sojourn(params,tau=0)
# 		M_to_EAC_time = min(M_to_EAC_time,temp$age_EAC)
# 	}
# 	indiv_test[k] <- M_to_EAC_time
# 	total_test[k] <- malig_sojourn_m(params,tau=0,total_M_cells)$age_EAC
# }

# quartz()
# hist(indiv_test,col='red')
# quartz()
# hist(total_test,col='blue')
#  