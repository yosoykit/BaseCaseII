EACstageGen = function(n=1,age=60,sex="A",race="A",screen=FALSE) {

  # stage distribution based on SEER18 EAC 2010-2012 incident cases (normalized to exclude undefined stage) 
    age.ind = 1; if(age >=60) age.ind = 2; if(age >=70) age.ind = 3

    ## All races (files SEER18AllRaceEACbyStage-2010-2012-Males.csv; etc)
    if(toupper(sex)=="M" | toupper(sex)=="MALE") {
      if(age.ind==1) {prob.stage = c(0.1682722, 0.2924559, 0.4942483, 0.0450236)}
      if(age.ind==2) {prob.stage = c(0.1952889, 0.3319850, 0.4104304, 0.0622957)}
      if(age.ind==3) {prob.stage = c(0.2432828, 0.3101813, 0.3408756, 0.1056603)}
    } else if (toupper(sex)=="F" | toupper(sex)=="FEMALE") {
      if(age.ind==1) {prob.stage = c(0.1862843, 0.2753048, 0.4996771, 0.0387338)}
      if(age.ind==2) {prob.stage = c(0.2527574, 0.2714176, 0.4247218, 0.0511032)}
      if(age.ind==3) {prob.stage = c(0.2123367, 0.2929681, 0.3253388, 0.1693563)}
    } else {
      if(age.ind==1) {prob.stage = c(0.1701451, 0.2926728, 0.4927897, 0.0443924)}
      if(age.ind==2) {prob.stage = c(0.2010207, 0.3261406, 0.4116283, 0.0612104)}
      if(age.ind==3) {prob.stage = c(0.2383783, 0.3082699, 0.3385619, 0.1147899)}
    }
    ## stage-shifting based on study by Cooper, G.S. et al., 2002.
    ## Association of prediagnosis endoscopy with stage and survival in adenocarcinoma of the esophagus
    ## and gastric cardia. Cancer, 95(1), pp.32–38.Cooper et al. 2002 (Fig.1B)
    if(screen==TRUE) {
      prob.stage[1] = prob.stage[1]*75/35
      prob.stage[2] = prob.stage[2]*13/32
      prob.stage[3] = prob.stage[3]*12/33
      prob.stage[4] = 1-sum(prob.stage[1:3])
    }
     return(sample(c("L","R","D","U"),size=n,prob=prob.stage,replace=T)) 
  }
