
surv_period_post80 <- function(currentage, time_interval, params, preinit, Pnumbers, Psizes, Msizes,extramalig,BEdist, EAdetect){
	BEdistmm = 10*BEdist
	screen_age1 <-current_age
	screen_age2 <-current_age+time_interval
	## 65-60, 70-65 etc.
	Tfinal = screen_age2-screen_age1
	### KEEP TRACK OF FOLLOWING OUTPUT: ####
	## detection of EA check
	EAdetect2 <- EAdetect 
	EAC_age <-NULL
	EAC_time=NULL
	EAC_time2=NULL
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
		tau1<- c(rep(screen_age1,preinit),tau1)
		Ptimes=0
		Ptimes<- initiations(tau1,params[2],screen_age2)
		if (Ptimes[1]>0){
			#totalinit2 <- length(Ptimes)
			tmax = rep(0,length(Ptimes))
			tmax = screen_age2-Ptimes
			for (i in 1:length(Ptimes)){
				out=0
				out=pclone(params,n0=1,M=1000,tmax[i])
				if (out$EAdetect==1){ 
					#cat('EA detected in person', '\n')
					## keep track of EA detections to condition
					EAdetect2<- 1
					EAC_age <- min(EAC_age,(out$p_malig_soj+Ptimes[i]))
					break
				}	    
			}
		}
		#### CONTINUE TAU LEAP GROWTH OF REMAINING P AND M CLONES IN HGD/EAC negative screened individuals
		if (Pnumbers>0){
			clonesize1 <- Psizes
			for (i in 1:length(clonesize1)){
				out=0
				out=pclone(params,n0=clonesize1[i],M=1000,(screen_age2-screen_age1))
				if (out$EAdetect==1){ 
					# keep track of EA detections to condition 
					EAdetect2 <- 1
					EAC_age <- min(EAC_age,(out$p_malig_soj+screen_age1))
					break
				}	    
			}
		}
		m2<-Msizes[2,]
		if (length(which((m2==1)))>0 && EAdetect2==0){
			### check to see if M close becomes detectable cancer by age 100
			total_M_cells= sum(Msizes[1,])
			EAC_time=malig_sojourn_m(params,tau=current_age,total_M_cells)$age_EAC
			EAC_age = min(EAC_age, EAC_time)
		}
		## check extra malignant clones with no surviving premalignant progenitor clone, if any, if they develop EAC
		if (length(extramalig)>0 && EAdetect2==0){
			total_M_cells= sum(extramalig)
			EAC_time2=malig_sojourn_m(params,tau=current_age,total_M_cells)$age_EAC
			EAC_age = min(EAC_age,EAC_time2,EAC_time)			
		}
	}
	return(list(clinicalEA=EAdetect2, current_age = screen_age2, EAC_age=EAC_age))
}




# total clinical, total screen, EAC deaths, endoscopies, RFA treatments, 

# start0_old=c(137.50978, 0,81.33647,0,0)
# start3_old = c(0.4177764,224.97258,72.64086,9699.828,0 )
# start4_old=c(23.422765,69.691888,35.422339,8158.457,253.82523 )
# strat6_old=c(27.55129,3.943512,18.27871,8742.02,1000)

# #dev.new()

# compare = matrix( c(start0_old,start3_old,start4_old,strat6_old, totals_0, totals_3,totals_4,totals_6),nrow=2)


