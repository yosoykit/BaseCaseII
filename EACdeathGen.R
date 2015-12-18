EACdeathGen = function(n=1, age = 60, sex = 1, yrdx = 2000, stage = "L") {

  # This routine uses output from CanSurv fits for EAC survival.
  # Stratified variables: age, sex, stage. Note, we have not found significant differences between the sexes. This variable is therefore
  # ignored. Range input is in months.

  ## all stages (A)

  # Cure (for 3 age groups: 50-59, 60-69, 70-84)
  b0 = -c(1.7763, 1.7748, 2.1619); b1 = c(0.0086, 0.0100, 0.0089)
  sd.b0 = c(0.0782, 0.0687, 0.0822); sd.b1 = c(0.0042, 0.0039, 0.0045)

  Cure = 1/(1+exp(-(b0+b1*(yrdx-1975)))) # uses logit
  sd.Cure = Cure*(1-Cure)*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # mu
  b0 = c(2.4053, 2.2644, 2.1014); b1 = c(0.0013, 0.0028, 0.0021)
  sd.b0 = c(0.0283, 0.0257, 0.0250); sd.b1 = c(0.0016, 0.0015, 0.0014)

  mu = (b0+b1*(yrdx-1975)) # log transformation
  msurv = exp(mu); sd.msurv = msurv*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # sigma
  b0 = c(0.0693, 0.1156, 0.1652); b1 = c(0.0003, 0.0003, 0.0004)
  sigma =  exp(b0+b1*(yrdx-1975))

  age.ind = 1; if(age >=60) age.ind = 2; if(age >=70) age.ind = 3

  # n random uniform deviates
  r = runif(n)
  non.cured = length(r[r>Cure[age.ind]])
  age.candth = rep(Inf,n)
  age.candth[r >  Cure[age.ind]] =  rlnorm(non.cured, meanlog = mu[age.ind], sdlog = sigma[age.ind])

  return(list(age = age.candth, Cure = Cure[age.ind], sd.Cure = sd.Cure[age.ind], msurv = msurv[age.ind], sd.msurv = sd.msurv[age.ind]))
}


EACdeathGen_resect = function(n=1, age = 60, sex = 1, yrdx = 2000, stage = "L",surgcurefactor=0.26042) {

  # This routine uses output from CanSurv fits for EAC survival.
  # Stratified variables: age, sex, stage. Note, we have not found significant differences between the sexes. This variable is therefore
  # ignored. Range input is in months.

  ## all stages (A)

  # Cure (for 3 age groups: 50-59, 60-69, 70-84)
  #b0 = -surgcurefactor*c(1.7763, 1.7748, 2.1619); b1 = c(0.0086, 0.0100, 0.0089)
  b0 = -c(1.7763, 1.7748, 2.1619); b1 = c(0.0086, 0.0100, 0.0089)
  sd.b0 = c(0.0782, 0.0687, 0.0822); sd.b1 = c(0.0042, 0.0039, 0.0045)

  Cure1 = 1/(1+exp(-(b0+b1*(yrdx-1975)))) # uses logit
  Cure = (1-surgcurefactor*(1-Cure1))
  sd.Cure = Cure*(1-Cure)*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # mu
  b0 = c(2.4053, 2.2644, 2.1014); b1 = c(0.0013, 0.0028, 0.0021)
  sd.b0 = c(0.0283, 0.0257, 0.0250); sd.b1 = c(0.0016, 0.0015, 0.0014)

  mu = (b0+b1*(yrdx-1975)) # log transformation
  msurv = exp(mu); sd.msurv = msurv*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # sigma
  b0 = c(0.0693, 0.1156, 0.1652); b1 = c(0.0003, 0.0003, 0.0004)
  sigma =  exp(b0+b1*(yrdx-1975))

  age.ind = 1; if(age >=60) age.ind = 2; if(age >=70) age.ind = 3

  # n random uniform deviates
  r = runif(n)
  non.cured = length(r[r>Cure[age.ind]])
  age.candth = rep(Inf,n)
  age.candth[r >  Cure[age.ind]] =  rlnorm(non.cured, meanlog = mu[age.ind], sdlog = sigma[age.ind])

  return(list(age = age.candth, Cure = Cure[age.ind], sd.Cure = sd.Cure[age.ind], msurv = msurv[age.ind], sd.msurv = sd.msurv[age.ind]))
}

# age_test=60
# runs=100000
# EACdeaths_nores = rep(0,runs)
# EACdeaths_withres = rep(0,runs)
# OCdeath = rep(0,runs)
# for (i in 1:runs){
#   EACdeaths_no = EACdeathGen(age=age_test)
#   EACdeaths_nores[i] = EACdeaths_no$age/12+age_test
#   EACdeaths_with = EACdeathGen_resect(age=age_test)
#   EACdeaths_withres[i]=EACdeaths_with$age/12+age_test
#   while (OCdeath[i] < age_test){
#     OCdeath[i] = OCdeathGenUS2K(race='white')
#   }
# }

