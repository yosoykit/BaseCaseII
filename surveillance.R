
surv_period <- function(currentage, time_interval, params, preinit, Pnumbers, Psizes, Msizes,extramalig,BEdist, EAdetect){
	BEdistmm = 10*BEdist
	screen_age1 <-current_age
	screen_age2 <-current_age+time_interval
	## 65-60, 70-65 etc.
	Tfinal = screen_age2-screen_age1
	### KEEP TRACK OF FOLLOWING OUTPUT: ####
	### Premalignant clone number and sizes at start time AFTER ABLATION
	Pnumbers2 = 0
	Psizes2 = rep(0,500)
	### Malignant clones : Matrix with size of malignant clone, which individual k it is in, ###
	####  and which premalignant clone (characterized by its size) it is in person k       ####
	Msizes_new <- matrix(rep(0,500*3), nrow=3, ncol=500)
	## detection of EA check
	EAdetect2 <- EAdetect 
	EAC_age <-NULL
	extramalig_new=NULL

	if (EAdetect2==0 && BEdistmm>0){
		### WHOLE PROCESS SEQUENTIALLY retaining heterogeneity in clone number and size
		## 1) Poisson events for progenitor, preinitiated cells N1
		BE = BEdistmm*kstem*5000
		X=BE
		lambda1 = params[1]*X*Tfinal
		N1 <- rpois(1,lambda1)
		preinit2 <- N1
		## times of preinitiations are uniformly distributed:
		tau1=0
		tau1 <- runif(N1,screen_age1, screen_age2)
		## Also have the number of preintitiated immortal stem cells at time of screen, so add # of preinitiations at screen_age1 to be possibly initiated rather than all pre-initiations de-novo
		#print(preinit)

		tau1<- c(rep(screen_age1,preinit),tau1)
		Ptimes=0
		Ptimes<- initiations(tau1,params[2],screen_age2)
		if (Ptimes[1]>0){
			#totalinit2 <- length(Ptimes)
			tmax = rep(0,length(Ptimes))
			tmax = screen_age2-Ptimes
			pclonesizesnew <- rep(0,length(Ptimes))
			mclonesizesnew <-rep(0,length(Ptimes))
			for (i in 1:length(Ptimes)){
				out=0
				out=pclone(params,n0=1,M=1000,tmax[i])
				pclonesizesnew[i] <- out$finalsize
				mclonesizesnew[i] <- out$maligclone
				if (out$EAdetect==1){ 
					#cat('EA detected in person', '\n')
					## keep track of EA detections to condition
					EAdetect2<- 1
					EAC_age <- min(EAC_age,(out$p_malig_soj+Ptimes[i]+screen_age1))
					break
				}	    
			}
		}
		else{
			pclonesizesnew <- NULL
			mclonesizesnew <- NULL
		}
		#### CONTINUE TAU LEAP GROWTH OF REMAINING P AND M CLONES IN HGD/EAC negative screened individuals
		m2<-Msizes[2,]
		if (Pnumbers>0){
			clonesize1 <- Psizes
			pclonesizes <- c(pclonesizesnew,rep(0,length(clonesize1)))
			mclonesizes_1<- c(mclonesizesnew, rep(0,length(clonesize1)))
			for (i in 1:length(clonesize1)){
				out=0
				out=pclone(params,n0=clonesize1[i],M=1000,(screen_age2-screen_age1))
				pclonesizes[(length(pclonesizesnew)+i)] <- out$finalsize
				mclonesizes_1[(length(mclonesizesnew)+i)] <- out$maligclone
				if (out$EAdetect==1){ 
					#cat('EA detected in person', k, '\n')
					# keep track of EA detections to condition 
					EAdetect2 <- 1
					EAC_age <- min(EAC_age,(out$p_malig_soj+screen_age1))
					break
				}	    
			}
		}
		else {
			mclonesizes_1=NULL
		}
		if (length(which((m2==1)))>0 && EAdetect2==0){
			### TAU LEAP the small malignant clones in non HGD/EAC screened population (part of missed malig population)
			epsilon <- .005
			tau <- epsilon/(params[5]-params[3])
			leaprates<- params[6:8]
			v<- c(-1,0,1) 
			mcloneind <- which(m2==1)
			mclonesizes <- c(mclonesizes_1,rep(0,length(mcloneind)))
			m_sizes_old <- Msizes[1,mcloneind]
			past_p_sizes <- Msizes[3,mcloneind]
			pclone_reference <- c(pclonesizes,past_p_sizes)
			for (p in 1:length(past_p_sizes)){
				if (past_p_sizes[p]>0){
					for (r in 1: length(clonesize1)){
						if (past_p_sizes[p]==clonesize1[r]){
							pclone_reference[length(pclonesizes)+p]<-pclonesizes[(length(pclonesizesnew)+r)]

						}
					}
				}
			}
			for (j in 1:length(mcloneind)){
				tfinal <- screen_age2 - screen_age1
				time0 <- 0   #running time since start of tau leap algorithm
				## create new variable x to track M clone size
				x <- m_sizes_old[j]     
				while (time0 < tfinal){
					lambda <- x*leaprates
					N<- rep(0,3)
					N <- rpois(3,lambda*tau)
					if (N[2]>0){
						EAdetect2 <- 1
						EAC_age <- min(EAC_age,(time0+screen_age1))
					}
					if (!(is.na(N[3])) && !(is.na(N[1]))){
						add<-0
						add <- N*v
						xtemp <- x + sum(add)
						## If count goes negative, discard this change and try smaller tau
						if (xtemp <0){
							tau=tau/2
						}
						else{	
							x <-xtemp
							time0 <- time0 + tau
						}
					}
					else {time0=time0+tau}
					if (x==0){ break}
				}
				mclonesizes[(length(mclonesizes_1)+j)] <- x
			}		
		}
		else {
			pclone_reference <- pclonesizes
			mclonesizes <- mclonesizes_1
		}

		## Tau leap extra malignant clones with no surviving premalignant progenitor clone, if any
		if (length(extramalig)>0 && EAdetect2==0){
			print('extra malig!')
			### TAU LEAP the small malignant clones in non HGD/EAC screened population (part of missed malig population)
			epsilon <- .005
			tau <- epsilon/(params[5]-params[3])
			leaprates<- params[6:8]
			v<- c(-1,0,1) 
			extramalig_new1 <- rep(0,length(extramalig))
			for (j in 1:length(extramalig_new1)){
				tfinal <- screen_age2 - screen_age1
				time0 <- 0   #running time since start of tau leap algorithm
				## create new variable x to track M clone size
				x <- extramalig[j]     
				while (time0 < tfinal){
					lambda <- x*leaprates
					N<- rep(0,3)
					N <- rpois(3,lambda*tau)
					if (N[2]>0){
						EAdetect2 <- 1
						EAC_age <- min(EAC_age,(time0+screen_age1))
					}
					if (!(is.na(N[3])) && !(is.na(N[1]))){
						add<-0
						add <- N*v
						xtemp <- x + sum(add)
						## If count goes negative, discard this change and try smaller tau
						if (xtemp <0){
							tau=tau/2
						}
						else{	
							x <-xtemp
							time0 <- time0 + tau
						}
					}
					else {time0=time0+tau}
					if (x==0){ break}
				}
				extramalig_new1[j] <- x
			}		
		}
		else {
			extramalig_new1 <- NULL
		}

		## keep track of malig clone sizes from which individual in which clone
		mcloneind <-0
		mtotal <- sum(mclonesizes)
		if (mtotal>0){
			mcloneind <- which(mclonesizes>0)
			Msizes_new[1,1:length(mcloneind)] <- mclonesizes[mcloneind]
			Msizes_new[2,1:length(mcloneind)] <- 1
			Msizes_new[3,1:length(mcloneind)] <- pclone_reference[mcloneind]
		}

		extramalig_new<- NULL
		mtemp = Msizes_new[2,]
		mlast = min(which(mtemp==0))
		if (mlast>1){
			Msizes_temp=matrix(Msizes_new[1:3,1:(mlast-1)],nrow=3)
			mtemp2 <- Msizes_temp[3,]
			if(length(which(mtemp2==0))>0){
				nopremaligind <- which(mtemp2==0)
				extramalig_new <- Msizes_temp[1,nopremaligind]
			}
		}

		extramalig_new =c(extramalig_new1,extramalig_new)
		## keep track of premalig clone numbers and sizes
		pcloneind <- 0
		if (sum(pclonesizes)>0){
			pcloneind <- which(pclonesizes>0)
			Pnumbers2 <- length(pcloneind)
			Psizes2[1:Pnumbers2] <-pclonesizes[pcloneind]
		}
	}
	#print(Psizes2[Psizes2>0])
	#print(Msizes_new[3,])
	return(list(pnumbers= Pnumbers2, preinit = length(tau1), psizesnew= Psizes2,msizesnew=Msizes_new, extramalignew=extramalig_new, clinicalEA=EAdetect2, current_age = screen_age2, EAC_age=EAC_age))
}
