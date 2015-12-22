EACdeathGen = function(n=1, age = 60, sex = 1, yrdx = 2000, stage = "U") {

  # This routine uses outputs from CanSurv fits for EAC survival.
  # Stratified variables: age, sex, stage. Note, no significant differences in EAC-specific survival
  # were found between the sexes. This variable is therefore ignored.
  # Range input is in months.

  # need stage info
  if(length(intersect(toupper(stage),c("A","L","R","D","U")))==0) {warning("stage variable needs to be one of L,R,D,A,U"); return}

  if(toupper(stage)=="A" | toupper(stage)=="U") {
  # Cure ALL (for 3 age groups: 50-59, 60-69, 70-84)
  b0 = -c(1.7763, 1.7748, 2.1619); b1 = c(0.0086, 0.0100, 0.0089)
  sd.b0 = c(0.0782, 0.0687, 0.0822); sd.b1 = c(0.0042, 0.0039, 0.0045)

  Cure = 1/(1+exp(-(b0+b1*(yrdx-1975)))) # uses logit
  sd.Cure = Cure*(1-Cure)*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # mu
  b0 = c(2.4053, 2.2644, 2.1014); b1 = c(0.0013, 0.0028, 0.0021)
  sd.b0 = c(0.0283, 0.0257, 0.0250); sd.b1 = c(0.0016, 0.0015, 0.0014)

  mu = (b0+b1*(yrdx-1975)) # log transformation
  msurv = exp(mu); sd.msurv = msurv*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)
}

  if(toupper(stage)=="L") {
  # Cure Localized (for 3 age groups: 50-59, 60-69, 70-84)
  b0 = -c(3.204e-01,5.902e-01,1.405e+00); b1 = c(1.119e-02,8.698e-03,1.550e-02) ## using logit Cure regressions
  sd.b0 = c(0.0782, 0.0687, 0.0822); sd.b1 = c(0.0042, 0.0039, 0.0045) ## using All

  Cure = 1/(1+exp(-(b0+b1*(yrdx-1975)))) # uses logit
  sd.Cure = Cure*(1-Cure)*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # mu
  b0 = c(3.079e+00,3.020e+00,2.767e+00); b1 = c(4.692e-03,4.604e-03,-6.089e-04)
  sd.b0 = c(0.0283, 0.0257, 0.0250); sd.b1 = c(0.0016, 0.0015, 0.0014) ## using All

  mu = (b0+b1*(yrdx-1975)) # log transformation
  msurv = exp(mu); sd.msurv = msurv*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)
}

  if(toupper(stage)=="R") {
  # Cure Regional (for 3 age groups: 50-59, 60-69, 70-84)
  b0 = -c(2.168e+00,2.018e+00,2.547e+00); b1 = c(1.905e-02,1.717e-02,6.802e-03) ## using logit Cure regressions
  sd.b0 = c(0.0782, 0.0687, 0.0822); sd.b1 = c(0.0042, 0.0039, 0.0045) ## using All

  Cure = 1/(1+exp(-(b0+b1*(yrdx-1975)))) # uses logit
  sd.Cure = Cure*(1-Cure)*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # mu
  b0 = c(2.858e+00,2.638e+00,2.420e+00); b1 = c(1.040e-03,1.834e-03,2.895e-03)
  sd.b0 = c(0.0283, 0.0257, 0.0250); sd.b1 = c(0.0016, 0.0015, 0.0014) ## using All

  mu = (b0+b1*(yrdx-1975)) # log transformation
  msurv = exp(mu); sd.msurv = msurv*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)
}
  
  if(toupper(stage)=="D") {
  # Cure Distant (for 3 age groups: 50-59, 60-69, 70-84)
  b0 = -c(3.3500,3.77494,4.265); b1 = c(0.00598,0.0130,1.021e-03) ## using logit Cure regressions
  sd.b0 = c(0.0782, 0.0687, 0.0822); sd.b1 = c(0.0042, 0.0039, 0.0045) ## using All

  Cure = 1/(1+exp(-(b0+b1*(yrdx-1975)))) # uses logit
  sd.Cure = Cure*(1-Cure)*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)

  # mu
  b0 = c(1.901e+00,1.743e+00,1.546e+00); b1 = c(1.741e-03,3.818e-03,3.072e-03)
  sd.b0 = c(0.0283, 0.0257, 0.0250); sd.b1 = c(0.0016, 0.0015, 0.0014) ## using All

  mu = (b0+b1*(yrdx-1975)) # log transformation
  msurv = exp(mu); sd.msurv = msurv*sqrt(sd.b0^2+(sd.b1*(yrdx-1975))^2)
}

  # sigma from All
  b0 = c(0.0693, 0.1156, 0.1652); b1 = c(0.0003, 0.0003, 0.0004)
  sigma =  exp(b0+b1*(yrdx-1975))

  age.ind = 1; if(age >=60) age.ind = 2; if(age >=70) age.ind = 3

  # n random uniform deviates
  r = runif(n)
  non.cured = length(r[r>Cure[age.ind]])
  time.candth = rep(Inf,n)
  time.candth[r >  Cure[age.ind]] =  rlnorm(non.cured, meanlog = mu[age.ind], sdlog = sigma[age.ind])

  return(list(time = time.candth, Cure = Cure[age.ind], sd.Cure = sd.Cure[age.ind], msurv = msurv[age.ind], sd.msurv = sd.msurv[age.ind]))
}

