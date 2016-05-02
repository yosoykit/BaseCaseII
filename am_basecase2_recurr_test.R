## STRATEGY 1 WHITE MALES Recurrence of HGD 

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
sex="M"
totalpop <- 100000
totals <- c(0,0,0,0,0)
# ages for recurrences: 
recurr_starts <- 50
# omega efficiencies : 40%, 50%, 70%,90%
#recurr_omegas <- c(.6,.5,.3,.1)
recurr_omegas <- c(.6,.3,.5)
recurr_m <- rep(0,length(recurr_omegas))
surv_int = 1
screen_age=50
count_o <-1
for (omega in recurr_omegas){
	hgdpop = 0
	hgdpop2=0
	end_age = screen_age+3
	death_by_EAC = 0
	death_by_OC=0
	for (run in 1:totalpop){
		## time of first screen of BE population
		#screen_age = 55
		## Draw time of death by other causes to not include EAC cases for those who have died
		EAC_death_check=1
		## assuming initial treatment efficacy ~ efficacy of 3.55 touchups
		omega_TU = omega^(1/3.55) 
		while (EAC_death_check==1){
			death_by_EAC = 0
			death_by_OC=0
			EAC_age = NULL
			## Draw time of death by other causes to not include EAC cases for those who have died
			#while (death_by_OC < screen_age){
			#	death_by_OC <-OCdeathGenUS2K(race='white')
			#}
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
			kstem_biop <- 100
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
						##cat('EA detected in person', '\n')
						## keep track of EA detections to condition
						#cat('clinical EAC in person', run, '\n')
						EAdetect <- 1
						EAC_age <- min(EAC_age,(out$p_malig_soj+Ptimes[i]))
						# cat('EAC age before screening', EAC_age, '\n')
						# surg = sample(c(0,1),1,prob=c(.2,.8))
						# resect=0
						# if (surg==1){
						# 	resect = sample(c(0,1),1,prob=c(.66,.33))
						# }
						# ## draw time of death by EAC for this individual
						# if (resect==1){
						# 	esoph_mort = sample(c(0,1),1,prob=c(.95,.05))
						# 	if (esoph_mort==1){
						# 		EACdeath= Inf
						# 	}
						# 	else{
						# 		death_by_EAC = EACdeathGen_resect(age=EAC_age)
						# 		EACdeath = death_by_EAC$age/12+EAC_age
						# 	}
						# }
						#  else{
							EACstage <- EACstageGen(n=1,age=EAC_age,sex=sex,race="A",screen=F)
							death_by_EAC = EACdeathGen_stage(n=1, age = EAC_age, sex = sex, yrdx = (birth_cohort+EAC_age), stage = EACstage)$time
							EACdeath = death_by_EAC/12+EAC_age 
						#} 
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
					#Msizes[1,(msizecount+1):(msizecount+length(mcloneind))] <- mclonesizes[mcloneind]
					Msizes[1,1:length(mcloneind)] <- mclonesizes[mcloneind]
					Msizes[2,1:length(mcloneind)] <- 1
					Msizes[3,1:length(mcloneind)] <- pclonesizes[mcloneind]
					#msizecount <- msizecount + length(mcloneind)
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
			missedmalig <- pos.biops$missedmalig
			number.posbiops<-rep(0,10)
			mnumber.posbiops<-rep(0,10)
			missednumber<- rep(0,10)
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
			 	missed<- missedmalig[k,]
			 	## to get missed malignancy percentage, don't count individuals with one missed malignancy but another detected
			 	if (missed>0 && mkbiops>0){
			 		missed<-0
			 	}
			 	missednumber[k] <- length(missed[missed>0])
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
			## HGD ablation
				# totals 22 is RFA touchups only
				#totals[22]=totals[22]+1
				RFA_outcome<- RFA_treatment(omega, Psizes, Msizes, Pnumbers, preinit, BEdist) 
				Psizes = RFA_outcome$Psizes
				Pnumbers=RFA_outcome$Pnumbers
				Msizes=RFA_outcome$Msizes
				preinit = RFA_outcome$preinit
				BEdist=RFA_outcome$BEdist
				# 2 years of RFA treatment
				# add 4 endoscopies for each initial treatment during those two years
				# decide if patient is success or failure: base values
				successes = c((.15*BEdist_RFA_effect),.68,.17)/sum(c((.15*BEdist_RFA_effect),.68,.17))
				RFA_success = sample(c(0,1,2),1,prob=successes)
				## Sensitivity analysis for effectiveness of RFA
				#RFA_success = sample(c(0,1,2),1,prob=c(.178,.644,.178))
				#RFA_success = sample(c(0,1,2),1,prob=c(.074,.889,.037))
				if (RFA_success==2){
					hgdpop=hgdpop+1
				}
				##cat('M clones post RFA', Msizes[3,], '\n')
				##cat('P clones post RFA', Psizes,'\n')
			#Strategy 1 & 3
			#surv_int = 1/4
		} 
		EAC_age=0
		if (EAdetect ==0  && RFA_success==2){
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
				if (EAdetect==1){
					EAC_age=next_screen$EAC_age
					cat('EAC age onset',next_screen$EAC_age, '\n')
				}
				current_age <- next_screen$current_age
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
					nopremaligind <- which(mtemp2==0)
					extramalig <- matrix(rep(0,3*length(nopremaligind)), nrow=3,ncol=length(nopremaligind))
					extramalig <- matrix(Msizes[,nopremaligind], nrow=3,ncol=length(nopremaligind),byrow=FALSE)
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
				kstem_biop <- 100
				pos.biops <- biop.sample_ind(Pnumbers, Psizes, Msizes2, BEdist,BEdistmm, extramalig,kstem_biop)
				premaligbiops <- pos.biops$pbiop
				maligbiops <- pos.biops$mbiop
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
					hgdpop2 = hgdpop2 +1
					## HGD caught during surveillance
					#surv_int = 1/4
				}
			}
			grade_status[4]=EAdetect
		}
		if (run %% 100==0){print(run); print(omega)}
	}
	print(hgdpop2)
	print(hgdpop)
	recurr_m[count_o] <- hgdpop2/hgdpop
	count_o = count_o+1
}

# dev.new()
# plot(seq(30,70,by=10),recurr_m[1,],type='o',pch=19,col="blue",lwd=2,xlab="Age", ylab="HGD Recurrence Rate among HGD", main='HGD Recurrence by RFA Efficacy',cex=2,cex.lab=1.5,cex.axis=1.5)
# lines(seq(30,70,by=10),recurr_m[2,],type='o',pch=19,col="red",lwd=2)
# lines(seq(30,70,by=10),recurr_m[3,],type='o',pch=19,col="forestgreen",lwd=2)
# lines(seq(30,70,by=10),recurr_m[4,],type='o',pch=19,col="slateblue",lwd=2)
# legend('topleft',c("40% Cells Removed","50% Cells Removed","70% Cells Removed","90% Cells Removed"), col=c('blue','red','forestgreen','slateblue'),pch=19,lwd=2,cex=2)

write(recurr_m,file="am_recurr_HGD_300316_sens4_100K.txt")
# nonEACpop0=totalpop-totals_0[1]
# nonEACpop1=totalpop-totals_1[1]
# nonEACpop3=totalpop-totals_3[1]

# # per 1000 BE patients
# eac_values = c((totals_0[2:3])/nonEACpop0,0,(totals_1[2:3])/nonEACpop1,0,(totals_3[2:3])/nonEACpop3,0)
# strategies <- matrix(heights1, ncol = 3, byrow = FALSE, dimnames = list(c("Detected EAC", "Clinical EAC", "Deaths by EAC"), c("Strategy 0", "Strategy 1", "Strategy 3")))
# colors <- c("darkblue", "red", "darkgreen")
# barplot(t(strategies), beside = TRUE, col = colors, ylim = c(0, 6), axes = FALSE, xlab = "Outcome by 85", main = "All Males: intervention begins at age 55")
# axis(2, at = 0:5, labels = 0:5)
# legend("topright", colnames(strategies), fill = colors, bty = "n")


