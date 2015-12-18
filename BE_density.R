######################################################################################################
##  Function to generate a GERD-dependent BE onset rate nu, the cumulative distribution up until   ###
##         		 index screening age, the density funtion up until index screening age,				##
##  							and the corresponding prevalence of GERD							##
##	This example for the Screening paper 2014 uses an exponential distribution for BE onset			##
######################################################################################################

BEdensity <- function(screen_age,gerd1,gerd2,gerd3,nu0){
	ages=1:screen_age
	RR <- 5
	nu <- rep(0,length(ages))
	nu_cum <- rep(0,length(ages))
	GERDpr1 <- rep(0,length(ages))
	for ( i in 1:length(ages)){
		GERDpr1[i]=1-exp(-gerd1*min(gerd3,ages[i])-gerd2*max(0,ages[i]-gerd3))
		nu[i] <- nu0*((1-GERDpr1[i])+RR*GERDpr1[i])
		GLquad1<-legauss(0 ,ages[i],20)
		wg <- GLquad1$weights
		xg <- GLquad1$mesh
		int1 <- rep(0,20)
		for ( j in 1:20){
			GERDpr=1-exp(-gerd1*min(gerd3,xg[j])-gerd2*max(0,xg[j]-gerd3))
			nu_temp <-nu0*((1-GERDpr)+RR*GERDpr)
			int1[j] <- nu_temp
		}
		nu_cum[i]<-1-exp(-sum(wg*int1))
	}
	nu_dens <- nu*(1-nu_cum)
	return(list(nu=nu,nu_cum=nu_cum,nu_dens=nu_dens,GERDpr=GERDpr1))
}


BEdensity_exactage <- function(screen_age,gerd1,gerd2,gerd3,nu0){
	ages=screen_age
	RR <- 5
	nu <- rep(0,length(ages))
	nu_cum <- rep(0,length(ages))
	GERDpr1 <- rep(0,length(ages))
	for ( i in 1:length(ages)){
		GERDpr1[i]=1-exp(-gerd1*min(gerd3,ages[i])-gerd2*max(0,ages[i]-gerd3))
		nu[i] <- nu0*((1-GERDpr1[i])+RR*GERDpr1[i])
		GLquad1<-legauss(0 ,ages[i],20)
		wg <- GLquad1$weights
		xg <- GLquad1$mesh
		int1 <- rep(0,20)
		for ( j in 1:20){
			GERDpr=1-exp(-gerd1*min(gerd3,xg[j])-gerd2*max(0,xg[j]-gerd3))
			nu_temp <-nu0*((1-GERDpr)+RR*GERDpr)
			int1[j] <- nu_temp
		}
		nu_cum[i]<-1-exp(-sum(wg*int1))
	}
	nu_dens <- nu*(1-nu_cum)
	return(list(nu=nu,nu_cum=nu_cum,nu_dens=nu_dens,GERDpr=GERDpr1))
}
