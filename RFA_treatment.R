RFA_treatment <- function(omega, Psizes, Msizes, Pnumbers, preinit, BEdist){
	##  OPTION 1: deplete all cell types by same fraction
	Psizes = floor(Psizes*omega)
	Psizes=Psizes[Psizes>0]
	Pnumbers=length(Psizes)
	if (length(Psizes)==0){
		Psizes=0
		Pnumbers=0
	}
	Mtemp = floor(Msizes[1,]*omega)
	if (length(which(Mtemp>0))>0){
		zero_ind<- which(Mtemp>0)
		Msizes = matrix(Msizes[,zero_ind],nrow=3)
		Msizes[1,]=floor(Msizes[1,]*omega)
		Msizes[3,]=floor(Msizes[3,]*omega)
		Mtemp3 <- Msizes[3,]
		if (!all(Msizes[3,]>0)){
			zero_ind <- which(Mtemp3==0)
			if (length(extramalig)>0){
				extramalig <- floor(omega*extramalig)
			}
			extramalig <- c(extramalig,Msizes[1,zero_ind])
			Msizes=matrix(Msizes[,-zero_ind],nrow=3)
		}
	}
	else{
		Msizes=matrix(rep(0,3),nrow=3)
	}
	preinit = floor(preinit*omega)
	BEdist=floor(BEdist*omega)
	return(list(Psizes=Psizes, Pnumbers=Pnumbers,Msizes = Msizes, preinit=preinit, BEdist=BEdist))
}