## 1 entry: Strategies- input 4 for Strat 4, and 5 for Strat 6 for sensitivity analyses
## 2 entry: Age at time of index endoscopsy
temp = c(4,60)
sex="M"
#stop surveilling CE-IM patients after certain years of surveillance for sensitivity analysis
sens_surv = 5
totalpop <- 1000
## number of stem cells in one crypt
kstem_biop <- 100
# fraction leftover of each cell-type post ablation 
omega=.6
## sensitivity of biopsy   4=30% HGD, 90% malig threshold for Base Case II
sens = 4
## BASE CASE 2: WHITE MALES, 1950 cohort
# surv_input<- function() {
#   ANSWER <- readline("Which surveillance strategy to run? ")
#   ## a better version would check the answer less cursorily, and
#   ## perhaps re-prompt
#   if (ANSWER==0)
#     temp=1
#   if (ANSWER==1)
#     temp=2
#   if (ANSWER==3)
#   	temp = 3
#   if (ANSWER==4)
#   	temp=4
#   if (ANSWER==6)
#   	temp=5
#   	ANSWER2 <- readline("What is the initial screen age?")
#    return(as.numeric(c(temp,ANSWER2)))
# }
#if(interactive()){ 
#	temp = surv_input()
#}
#else {


#}
print(temp[2])
print(temp[1])
strategy =temp[1]
screen_age=temp[2]

source('AllMortalityGenUS.R')
source('biopsy_functions.R')
source('am_parameter_list.R')
source('RFA_treatment.R')
source('surveillance.R')
source('surveillance_post80.R')
source('BE_density.R')
source('EACdeathGen.R')
source('EACdeathGen2.R')
source('EACdeathGen2_Wani.R')
source("EACstageGen2.R")
params = amparams	

## (initial Maligs, initial EAC, screen detected EAC, clinical EAC, deaths by EAC, deaths from surgery, initial ND (totalpop-initial HGD),initial HGD, screen EAC from ND, screen EAC from HGD, clin EAC from ND, clin EAC from HGD, # of endoscopies, death of other causes, # of EAC surgeries)
totals <- rep(0,28)
#names(totals)=c('initial maligs', 'initial EAC','total EAC screened','total clinical EAC', 'EAC deaths', 'surgery deaths', 'initial ND', 'initial HGD','EAC screened: ND','EAC screened: HGD', 'clinical EAC:ND', 'clinical EAC:HGD', 'num. endoscopies','OC deaths', 'num EAC surgeries', 'EAC death:ND', 'EAC death:HGD', 'OC death:ND', 'OC death:HGD', 'no. endoscopies:ND', 'no.endscopies:HGD')
names(totals)=c('initial maligs', 'initial EAC','total EAC screened','total clinical EAC', 'EAC deaths', 'surgery deaths', 'initial ND', 'initial HGD','EAC screened: ND','EAC screened: HGD', 'clinical EAC:ND', 'clinical EAC:HGD', 'num. endoscopies','OC deaths', 'num EAC surgeries', 'EAC death:ND', 'EAC death:HGD', 'OC death:ND', 'OC death:HGD', 'no. endoscopies:ND', 'no.endscopies:HGD','RFA touchups', 'RFA:ND', 'RFA:HGD', 'RFA touchups:ND', 'RFA touchups:HGD', 'Lifeyears:ND', 'Lifeyears:HGD')
end_age = 80
final_surv = 100
for (run in 1:totalpop){
	## time of first screen of BE population
	survey_int = matrix(c(end_age-screen_age,end_age-screen_age,3,3,3, end_age-screen_age,1/4,1/4,1/4,1/4),nrow=2, byrow=T)
	## for ablation strategies
	## make sure person with EAC before screen_age does not die before screen_age to be counted as a prevalent case
	EAC_death_check=1
	## assuming initial treatment efficacy ~ efficacy of 3.55 touchups
	omega_TU = omega^(1/3.55) 
	while (EAC_death_check==1){
		death_by_EAC = 0
		death_by_OC=0
		EAC_age = NULL
		death_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=screen_age)$age
		## will change as individual steps through surveillance strategy
		current_age = screen_age
		BEdist <- rbeta(1,.5,4)*15+1 #in cm
		BEdistmm <- BEdist*10   #in mm
		## male individual, Tfinal = screen age - BE onset age : total time for process
		## must develop BE before screening age
		BEparams=BEdensity(screen_age,gerd1,gerd2,gerd3,nu0)
		nu_cum = BEparams$nu_cum
		ages = 1:screen_age
		BEstart <- 100

		while (BEstart >= screen_age){
			x <- runif(1,0,nu_cum[length(nu_cum)])
			if (x<nu_cum[1]){
				x <- nu_cum[1]
			}
			temp_BE <- approx(nu_cum, ages, xout=x)$y
			BEstart <-round(temp_BE)
		}
		Tfinal = screen_age-BEstart
		### KEEP TRACK OF FOLLOWING OUTPUT: ####
		### Premalignant clone number and sizes
		preinit =0
		Pnumbers = 0
		Psizes = rep(0,500)
		### Malignant clones : Matrix with size of malignant clone, which individual k it is in, ###
		####  and which premalignant clone (characterized by its size) it is in person k       ####
		Msizes <- matrix(rep(0,50*3), nrow=3, ncol=50)
		## detection of EA check
		EAdetect <- 0
		### WHOLE PROCESS SEQUENTIALLY retaining heterogeneity in clone number and size
		## 1) Poisson events for progenitor, preinitiated cells N1
		BE = BEdistmm*kstem*5000
		X=BE
		lambda1 = params[1]*X*Tfinal
		N1 <- rpois(1,lambda1)
		preinit <- N1
		## times of preinitiations are uniformly distributed:
		tau1=0
		tau1 <- runif(N1,BEstart,screen_age)
		Ptimes=0
		Ptimes<- initiations(tau1,params[2], screen_age)
		if (Ptimes[1]>0){
			tmax = rep(0,length(Ptimes))
			tmax = screen_age-Ptimes
			pclonesizes <- rep(0,length(Ptimes))
			mclonesizes <-rep(0,length(Ptimes))
			for (i in 1:length(Ptimes)){
				out=0
				out=pclone(params,n0=1,M=10000,tmax[i])
				pclonesizes[i] <- out$finalsize
				mclonesizes[i] <- out$maligclone
				if (out$EAdetect==1){ 
					EAdetect <- 1
					EAC_age <- min(EAC_age,(out$p_malig_soj+Ptimes[i]))
					EACstage <- EACstageGen(n=1,age=EAC_age,sex=sex,race="A",screen=F)
					death_by_EAC = EACdeathGen_stage(n=1, age = EAC_age, sex = sex, yrdx = (birth_cohort+EAC_age), stage = EACstage)$time
					EACdeath = death_by_EAC/12+EAC_age 
					if (EACdeath > screen_age){
						EAC_death_check=0
					}
					else if (EACdeath <= screen_age){
						EAC_death_check=1
						break
					}
				}	    
			}
			## keep track of malig clone sizes from which individual in which clone
			mcloneind <-0
			mtotal <- sum(mclonesizes)
			if (mtotal>0){
				mcloneind <- which(mclonesizes>0)
				Msizes[1,1:length(mcloneind)] <- mclonesizes[mcloneind]
				Msizes[2,1:length(mcloneind)] <- 1
				Msizes[3,1:length(mcloneind)] <- pclonesizes[mcloneind]
			}	
			## keep track of premalig clone numbers and sizes
			pcloneind <- 0
			if (sum(pclonesizes)>0){
				pcloneind <- which(pclonesizes>0)
				Pnumbers <- length(pcloneind)
				Psizes[1:Pnumbers] <-pclonesizes[pcloneind]
			}
		}
		if (EAdetect==0){
			EAC_death_check=0
		}	
	}
	if (EAdetect ==0){
		plast = min(which(Psizes==0))
		Psizes=Psizes[1:(plast-1)]
		mtemp = Msizes[2,]
		mlast = min(which(mtemp==0))
		Msizes2 <- rep(0,length(Psizes))
		extramalig <- NULL
		if (mlast>1){
			Msizes=matrix(Msizes[1:3,1:(mlast-1)],nrow=3)
			mtemp2 <- Msizes[3,]
			if(length(which(mtemp2==0))>0){
				nopremaligind <- which(mtemp2==0)
				extramalig <- Msizes[1,nopremaligind]
			}
			clonesize1 <- Psizes
			for (j in 1:length(Msizes[1,])){
				for (k in 1:length(clonesize1)){
					if (Msizes[3,j]==clonesize1[k]){
						kspot=k
						if (Msizes2[kspot]==0){
						 	Msizes2[kspot]<- Msizes[1,j]
						 	break
						}		 			 	
					}
				}
			}
		}
		else {
			Msizes = matrix(rep(0,3),nrow=3)
		}
		pos.biops <- biop.sample_ind(Pnumbers, Psizes, Msizes2, BEdist,BEdistmm, extramalig,kstem_biop)
		premaligbiops <- pos.biops$pbiop
		maligbiops <- pos.biops$mbiop
		number.posbiops<-rep(0,10)
		mnumber.posbiops<-rep(0,10)
		Psizes=pos.biops$Psizes
		Msizes2 = pos.biops$Msizes
		if (sum(Msizes2)>0){
			msizes_ind = which(Msizes2>0)
			Msizes=matrix(c(Msizes2[msizes_ind],rep(1,length(msizes_ind)),Psizes[msizes_ind]),nrow=3,byrow=T)
		}
		Pnumbers = pos.biops$Pnumbers
		extramalig = pos.biops$extramalig
		for (k in 1:10){
		 	kbiops <- premaligbiops[k,]
		 	number.posbiops[k] <- length(kbiops[kbiops>0])
		 	mkbiops <- maligbiops[k,]
		 	mnumber.posbiops[k] <- length(mkbiops[mkbiops>0])
		}
	}
	# Adjust future RFA failure rates as higher for longer BE segments
	BEdist_RFA_effect = BEdist/(15*.5/4.5+1)
	current_age=screen_age
	## DECIDE IF THIS PERSON IS ND,HGD, Malignant. 2 = EAC or malignancy at time of first screen (not included in eventual denominator)
	grade_status=c(0,0,0,EAdetect)
	# touchup counter
	touchup=0
	BE_recur =0
	## RFA_success set to 3, will only be 2,1 or 0 if RFA occurs
	RFA_success=3
	## check if its first surveillance after RFA to add up extra screens not explicitly modeled
	postRFA_screens=0
	## check if previous screen was CE-IM for surveillance interval sensitivity	
	last_screen_CEIM = 0
	## count time (years) of CE-IM
	CEIM = 0
	## sensitivity of biopsy   5=40% HGD, 20% malig
	sens = 4
	if (EAdetect==1){
		grade_status[4]=2
	}
	if (mnumber.posbiops[sens]>0 && EAdetect==0){
		grade_status[3]=2
	} 
	else if(number.posbiops[sens]>0 && EAdetect==0){
		totals[8]=totals[8]+1
		grade_status[2]=1
		surv_int = survey_int[2,strategy]
		## HGD ablation
		if (strategy >=4){
			totals[24]=totals[24]+1
			RFA_outcome<- RFA_treatment(omega, Psizes, Msizes, Pnumbers, preinit, BEdist) 
			Psizes = RFA_outcome$Psizes
			Pnumbers=RFA_outcome$Pnumbers
			Msizes=RFA_outcome$Msizes
			preinit = RFA_outcome$preinit
			BEdist=RFA_outcome$BEdist
			# 2 years of RFA treatment
			current_age=current_age+2
			# add 4 endoscopies for each initial treatment during those two years
			totals[21] = totals[21]+4
			# decide if patient is success or failure: base values
			successes = c((.15*BEdist_RFA_effect),.68,.17)/sum(c((.15*BEdist_RFA_effect),.68,.17))
			RFA_success = sample(c(0,1,2),1,prob=successes)
			## Sensitivity analysis for effectiveness of RFA
			#successes = c((.178*BEdist_RFA_effect),.644,.178)/sum(c((.178*BEdist_RFA_effect),.644,.178))
			#RFA_success = sample(c(0,1,2),1,prob=successes)
			#successes = c((.074*BEdist_RFA_effect),.889,.037)/sum(c((.074*BEdist_RFA_effect),.889,.037))
			#RFA_success = sample(c(0,1,2),1,prob=successes)
			if (RFA_success==0){
				# put to max value, non CE-IM and non CE-D never get touchups
				touchup = 3
			}
			else{
				# yearly screens for HGD ablated, CE-D
				surv_int=1
				## Add extra endoscopies for RFA_success cases.  Failures are modeled and counted explicitly 
				totals[13]=totals[13]+4
				totals[21]=totals[21]+4
				if (RFA_success==1){
					last_screen_CEIM =1
				}
			}
		}
	} 
	else if(EAdetect==0){
		grade_status[1]=1
		surv_int = survey_int[1,strategy]
		## ND ablation
		if (strategy ==5){
			totals[23]=totals[23]+1
			RFA_outcome<- RFA_treatment(omega, Psizes, Msizes, Pnumbers, preinit, BEdist) 
			Psizes = RFA_outcome$Psizes
			Pnumbers=RFA_outcome$Pnumbers
			Msizes=RFA_outcome$Msizes
			preinit = RFA_outcome$preinit
			BEdist=RFA_outcome$BEdist
			# 2 years of RFA treatment
			current_age=current_age+2
			# add 4 endoscopies for each initial treatment during those two years
			totals[20] = totals[20]+4
			# decide if patient is success or failure: base case values
			successes=c((.189*BEdist_RFA_effect),.811)/sum(c((.189*BEdist_RFA_effect),.811))
			RFA_success = sample(c(0,1),1,prob=successes)
			### Sensitivity analysis for effectiveness of RFA
			#successes=c((.315*BEdist_RFA_effect),.685)/sum(c((.315*BEdist_RFA_effect),.685))
			#RFA_success = sample(c(0,1),1,prob=successes)
			#successes=c((.016*BEdist_RFA_effect),.984)/sum(c((.016*BEdist_RFA_effect),.984))
			#RFA_success = sample(c(0,1),1,prob=successes)
			if (RFA_success==0){
				# put to max value, non CE-IM never get touchups
				touchup = 3
			}
			else{
				surv_int = 3
				last_screen_CEIM =1
			}
		}
	}
	EAC_age=0
	while (EAdetect ==0 && current_age < end_age && current_age < death_by_OC){
		if (grade_status[3]>0){
			break
		}
		else{
			next_screen <- surv_period(current_age, surv_int, params, preinit, Pnumbers, Psizes, Msizes,extramalig,BEdist, EAdetect)
			Pnumbers<-next_screen$pnumbers
			preinit <- next_screen$preinit
			Psizes <- next_screen$psizesnew
			Msizes<- next_screen$msizesnew
			EAdetect <- next_screen$clinicalEA
			extramalig <- next_screen$extramalignew
			if (EAdetect==1){
				EAC_age=next_screen$EAC_age
			}
			time_lapsed <- next_screen$current_age-current_age
			current_age <- next_screen$current_age

		}
		if (EAdetect ==0 && current_age < end_age){
			# perform endoscopy at current_age
			totals[13]=totals[13]+1
			if (grade_status[1]==1){
				totals[20]=totals[20]+1
			}
			else if(grade_status[2]==1){
				totals[21]=totals[21]+1
			}
			plast = min(which(Psizes==0))
			Psizes=Psizes[1:(plast-1)]
			mtemp = Msizes[2,]
			mlast = min(which(mtemp==0))
			Msizes2 <- rep(0,length(Psizes))
			if (mlast>1){
				Msizes=matrix(Msizes[1:3,1:(mlast-1)],nrow=3)
				mtemp2 <- Msizes[3,]
				clonesize1 <- Psizes
				for (j in 1:length(Msizes[1,])){
					for (k in 1:length(clonesize1)){
						if (Msizes[3,j]==clonesize1[k]){
							kspot=k
							Msizes2[kspot]=Msizes2[kspot]+Msizes[1,j]
							break 			 	
						}
					}
				}
			}
			if (RFA_success==1 && touchup <3){
				if (grade_status[2]==1){
					BE_recur = sample(c(0,1),1,prob=c(.9,.1))
					## sensitivity analysis
					#BE_recur = sample(c(0,1),1,prob=c(.95,.5))
					#BE_recur = sample(c(0,1),1,prob=c(.80,.20))
				}
				else if (grade_status[1]==1){
					## BE RECUR from initially ND at year 3 = 21% for Strat 6
					BE_recur = sample(c(0,1),1,prob=c(.79,.21))
					## sensitivity analysis
					#BE_recur = sample(c(0,1),1,prob=c(.895,.105))
					#BE_recur = sample(c(0,1),1,prob=c(.58,.42))
				}
				# ablate recurrent BE cases up to 3 times from successful HGD or ND RFA treatment
				if (BE_recur==1){
					RFA_outcome<- RFA_treatment(omega_TU, Psizes, Msizes, Pnumbers, preinit, BEdist) 
					Psizes = RFA_outcome$Psizes
					Pnumbers=RFA_outcome$Pnumbers
					Msizes=RFA_outcome$Msizes
					preinit = RFA_outcome$preinit
					BEdist=RFA_outcome$BEdist
					touchup=touchup+1
					totals[22]=totals[22]+1
					if (grade_status[1]==1){
						totals[25]=totals[25]+1
					}
					else if(grade_status[2]==1){
						totals[26]=totals[26]+1
					}
					if (grade_status[2]==1){
						# decide if patient is success or failure:base case values
						successes = c((.15*BEdist_RFA_effect),.68,.17)/sum(c((.15*BEdist_RFA_effect),.68,.17))
						RFA_success = sample(c(0,1,2),1,prob=successes)
						## Sensitivity analysis for effectiveness of RFA
						#successes = c((.178*BEdist_RFA_effect),.644,.178)/sum(c((.178*BEdist_RFA_effect),.644,.178))
						#RFA_success = sample(c(0,1,2),1,prob=successes)
						#successes = c((.074*BEdist_RFA_effect),.889,.037)/sum(c((.074*BEdist_RFA_effect),.889,.037))
						#RFA_success = sample(c(0,1,2),1,prob=successes)
						if (RFA_success==0){
							# put to max value, non CE-IM and non CE-D never get touchups
							touchup = 3
						}
					}
					else if (grade_status[1]==1){
						# decide if patient is success or failure : base case values
						successes=c((.189*BEdist_RFA_effect),.811)/sum(c((.189*BEdist_RFA_effect),.811))
						RFA_success = sample(c(0,1),1,prob=successes)
						### Sensitivity analysis for effectiveness of RFA
						#successes=c((.315*BEdist_RFA_effect),.685)/sum(c((.315*BEdist_RFA_effect),.685))
						#RFA_success = sample(c(0,1),1,prob=successes)
						#successes=c((.016*BEdist_RFA_effect),.984)/sum(c((.016*BEdist_RFA_effect),.984))
						#RFA_success = sample(c(0,1),1,prob=successes)
						if (RFA_success==0){
							# put to max value, non CE-IM never get touchups
							touchup = 3
						}
					}
					## Add extra endoscopies for RFA_success cases (looking AHEAD).  Failures are modeled and counted explicitly 
					if (grade_status[2]==1 && grade_status[1]==1 && RFA_success >=1){
						## ablate HGD, started as ND
						totals[13]=totals[13]+4
						totals[20]=totals[20]+4
					}
					else if (grade_status[2]==1 && grade_status[1]==0 && RFA_success >=1){
						## ablate HGD, started as HGD
						totals[13]=totals[13]+4
						totals[21]=totals[21]+4
					}
				}
			}
			else {
				pos.biops <- biop.sample_ind(Pnumbers, Psizes, Msizes2, BEdist,BEdistmm, extramalig,kstem_biop)
				premaligbiops <- pos.biops$pbiop
				maligbiops <- pos.biops$mbiop
				Psizes=pos.biops$Psizes
				Msizes2 = pos.biops$Msizes
				if (sum(Msizes2)>0){
					msizes_ind = which(Msizes2>0)
					Msizes=matrix(c(Msizes2[msizes_ind],rep(1,length(msizes_ind)),Psizes[msizes_ind]),nrow=3,byrow=T)
				}
				Pnumbers = pos.biops$Pnumbers
				extramalig = pos.biops$extramalig
				number.posbiops<-rep(0,10)
				mnumber.posbiops<-rep(0,10)
				for (k in 1:10){
				 	kbiops <- premaligbiops[k,]
				 	number.posbiops[k] <- length(kbiops[kbiops>0])
				 	mkbiops <- maligbiops[k,]
				 	mnumber.posbiops[k] <- length(mkbiops[mkbiops>0])
				}

				if (mnumber.posbiops[sens]>0 && EAdetect==0){
					grade_status[3]=1
				}
				else if (number.posbiops[sens]>0 && EAdetect==0){
					grade_status[2]=1
					## HGD caught during surveillance
					if (RFA_success==1 && touchup>=3){ 
						surv_int = 1/4
					}
					else if( RFA_success==2&& touchup>=3){ 
						surv_int = 1/4
					}
					if (strategy >=4 && RFA_success >=2 && touchup<3){
						# this should be only originally ND patients who are found with HGD on subsequent screen or RFA_success==2 is HGD -> non CE-IM, CE-D or RFA_success==1 and you've had a BE_recur==1
						if (RFA_success==3){
							## originally ND patients who are found with HGD on subsequent screen receive RFA for the first time now
							totals[23]=totals[23]+1
							## Ablate newly found HGD case
							RFA_outcome<- RFA_treatment(omega, Psizes, Msizes, Pnumbers, preinit, BEdist) 
								successes = c((.15*BEdist_RFA_effect),.68,.17)/sum(c((.15*BEdist_RFA_effect),.68,.17))
								RFA_success = sample(c(0,1,2),1,prob=successes)
								## Sensitivity analysis for effectiveness of RFA
								#successes = c((.178*BEdist_RFA_effect),.644,.178)/sum(c((.178*BEdist_RFA_effect),.644,.178))
								#RFA_success = sample(c(0,1,2),1,prob=successes)
								#successes = c((.074*BEdist_RFA_effect),.889,.037)/sum(c((.074*BEdist_RFA_effect),.889,.037))
								#RFA_success = sample(c(0,1,2),1,prob=successes)
							if (RFA_success==0){
								# put to max value, non CE-IM and non CE-D never get touchups
								touchup = 3
							}
							else{
								# yearly screens for HGD ablated, CE-D
								surv_int=1
								## Add extra endoscopies for RFA_success cases (looking AHEAD).  Failures are modeled and counted explicitly 
								## ablate HGD, started as ND
								totals[13]=totals[13]+4
								totals[20]=totals[20]+4
							}
							current_age=current_age+2
							# add 4 endoscopies for each initial treatment during those two years
							totals[20] = totals[20]+4
						}
						else if (RFA_success==2){
							totals[22]=totals[22]+1
							touchup = touchup +1
							if (grade_status[1]==1){
								totals[25]=totals[25]+1
							}		
							else if(grade_status[2]==1){
								totals[26]=totals[26]+1
							}
							## Ablate newly found HGD case
							RFA_outcome<- RFA_treatment(omega_TU, Psizes, Msizes, Pnumbers, preinit, BEdist) 
							successes = c((.15*BEdist_RFA_effect),.68,.17)/sum(c((.15*BEdist_RFA_effect),.68,.17))
							RFA_success = sample(c(0,1,2),1,prob=successes)
							## Sensitivity analysis for effectiveness of RFA
							#successes = c((.178*BEdist_RFA_effect),.644,.178)/sum(c((.178*BEdist_RFA_effect),.644,.178))
							#RFA_success = sample(c(0,1,2),1,prob=successes)
							#successes = c((.074*BEdist_RFA_effect),.889,.037)/sum(c((.074*BEdist_RFA_effect),.889,.037))
							#RFA_success = sample(c(0,1,2),1,prob=successes)
							if (RFA_success==0){
								# put to max value, non CE-IM and non CE-D never get touchups
								touchup = 3
							}
							else{
								if (grade_status[2]==1 && grade_status[1]==1){
									## ablate HGD, started as ND
									totals[13]=totals[13]+4
									totals[20]=totals[20]+4
								}
								else if (grade_status[2]==1 && grade_status[1]==0){
									## ablate HGD, started as HGD
									totals[13]=totals[13]+4
									totals[21]=totals[21]+4
								}
								# yearly screens for HGD ablated, CE-D
								surv_int=1
							}
						}
						Psizes = RFA_outcome$Psizes
						Pnumbers=RFA_outcome$Pnumbers
						Msizes=RFA_outcome$Msizes
						preinit = RFA_outcome$preinit
						BEdist=RFA_outcome$BEdist
					}
				}
			}
		}
		grade_status[4]=EAdetect
		if (RFA_success==1 && last_screen_CEIM==1){
			CEIM = CEIM + time_lapsed
		}
		if (RFA_success==1){
			last_screen_CEIM=1
		}
		else {
			CEIM=0
		}
		if (CEIM >=sens_surv){
			break
		}
	}
	if (grade_status[4]==0 && grade_status[3]==0 && current_age >=end_age && BEdist>0){
		if (current_age >=end_age || CEIM >=sens_surv ){
			if (!all(Psizes>0)){
				plast = min(which(Psizes==0))
				Psizes=Psizes[1:(plast-1)]
				mtemp = Msizes[2,]
				mlast = min(which(mtemp==0))
				Msizes2 <- rep(0,length(Psizes))
				if (mlast>1){
					Msizes=matrix(Msizes[1:3,1:(mlast-1)],nrow=3)
					mtemp2 <- Msizes[3,]
					clonesize1 <- Psizes
					for (j in 1:length(Msizes[1,])){
						for (k in 1:length(clonesize1)){
							if (Msizes[3,j]==clonesize1[k]){
								kspot=k
								Msizes2[kspot]=Msizes2[kspot]+Msizes[1,j]
								break 			 	
							}
						}
					}
				}
			}
			surv_int = final_surv - current_age
			post80_outcome <- surv_period_post80(current_age, surv_int, params, preinit, Pnumbers, Psizes, Msizes,extramalig,BEdist, EAdetect)
			EAdetect = post80_outcome$clinicalEA
			grade_status[4]=EAdetect
			if (EAdetect==1){
				EAC_age = post80_outcome$EAC_age
			}
		}
	}
	EACdeath=NULL
	### Explain grade_status
	if (grade_status[3]==2){
		totals[1]=totals[1]+1
	}
	else if (grade_status[4]==2){
		totals[2]=totals[2]+1
	}
	else if (grade_status[4]==1){
		# rho detection event cisnet occur sometime interval after death by OC 
		if (EAC_age < death_by_OC){
			#cat('clinical EAC in person', run, '\n')
			totals[4]=totals[4]+1
			if (grade_status[1]==1){
				totals[11]=totals[11]+1
			}
			else if(grade_status[2]==1){
				totals[12]=totals[12]+1
			}
			EACstage <- EACstageGen(n=1,age=EAC_age,sex=sex,race="A",screen=F)
			death_by_EAC = EACdeathGen_stage(n=1, age = EAC_age, sex = sex, yrdx = (birth_cohort+EAC_age), stage = EACstage)$time
			EACdeath = death_by_EAC/12+EAC_age 
			if (EACdeath < final_surv && EACdeath < death_by_OC){
				totals[5]=totals[5]+1
				if (grade_status[1]==1){
					totals[16]=totals[16]+1
				}
				else if(grade_status[2]==1){
					totals[17]=totals[17]+1
				}
			}
		}
		else if (EAC_age >= death_by_OC && death_by_OC<final_surv){
			totals[14]=totals[14]+1
			if (grade_status[1]==1){
					totals[18]=totals[18]+1
			}
			else if(grade_status[2]==1){
				totals[19]=totals[19]+1
			}
		}	
	}
	else if (grade_status[3]==1){
		totals[3]=totals[3]+1
		if (grade_status[1]==1){
			totals[9]=totals[9]+1
		}
		else if (grade_status[2]==1){
			totals[10]=totals[10]+1
		}
		EACstage <- EACstageGen(n=1,age=current_age,sex=sex,race="A",screen=T)
		total_M_cells= sum(Msizes[1,])
		EAC_hypoth_time=malig_sojourn_m(params,tau=current_age,total_M_cells)$age_EAC
		if (EAC_hypoth_time < Inf){
			death_by_EAC = EACdeathGen_stage(n=1, age = EAC_hypoth_time, sex = sex, yrdx = (birth_cohort+EAC_hypoth_time), stage = EACstage)$time
			EACdeath = death_by_EAC/12+EAC_hypoth_time
		}
		else{ EACdeath=Inf}
		if (EACdeath < final_surv && EACdeath < death_by_OC){
				totals[5]=totals[5]+1
				if (grade_status[1]==1){
					totals[16]=totals[16]+1
				}
				else if(grade_status[2]==1){
					totals[17]=totals[17]+1
				}
			}
		else if (EAC_age >= death_by_OC && death_by_OC<final_surv){
			totals[14]=totals[14]+1
			if (grade_status[1]==1){
					totals[18]=totals[18]+1
			}
			else if(grade_status[2]==1){
				totals[19]=totals[19]+1
			}
		}
	}
	else if (grade_status[2]==1){
		if (final_surv > death_by_OC){
			totals[14]=totals[14]+1
			if (grade_status[1]==1){
					totals[18]=totals[18]+1
			}
			else if(grade_status[2]==1){
				totals[19]=totals[19]+1
			}
		}
	}
	else if (grade_status[1]==1){
		if (final_surv > death_by_OC){
			totals[14]=totals[14]+1
			if (grade_status[1]==1){
					totals[18]=totals[18]+1
			}
			else if(grade_status[2]==1){
				totals[19]=totals[19]+1
			}
		}
	}
	if (grade_status[3]<2 && grade_status[4]<2){
		life_years= min(EACdeath, final_surv, death_by_OC)-screen_age
		if (grade_status[1]==1){
			totals[27]=totals[27]+life_years
		}
		else if (grade_status[2]==1){
			totals[28]=totals[28]+life_years
		}
	}
	if (run %% 100==0){print(run)}
}


nonEAC_pop = totalpop-sum(totals[1:2])

totals[7]=nonEAC_pop-totals[8]
## first 6 columns in template include screen detected cancers on initial screen... rest do not
totals[1:2] = totals[1:2]/totalpop*1000
totals[7:8] = totals[7:8]/totalpop*1000

totals[3:6]=totals[3:6]/nonEAC_pop*1000
totals[9:28]=totals[9:28]/nonEAC_pop*1000

if (strategy==1){
	write(totals,file='am_basecase2_strat0_60_020516.txt')
}
if (strategy==3){
	write(totals,file='am_basecase2_strat3_60_020516.txt')
}
if (strategy==4){
	write(totals,file='am_basecase2_strat4_60_020516.txt')
}
if (strategy==5){
	write(totals,file='am_basecase2_strat6_60_020516.txt')
}
