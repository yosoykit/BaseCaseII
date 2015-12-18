AllMortalityGenUS = function(n=1, byear=2000, start=60, sex="male", race="all") {

  if(byear > 2008) {warning("birthyear < 2009!"); return}
  check = 0
  range = 0:119
  birthyears = 1860:2008
  ages = floor(start):max(range)
  if(toupper(race) == "WHITE") {
    if(toupper(sex) == "MALE") {
      check=1
      dat = read.csv('AllCauseMortRatesFromUSLifeTablesByBirthYear-UpdateTo2008Table-Age120-whitemales.csv',sep=",",header=F)
      # surv = c(100000, 99243, 99163, 99052, 98615, 98002, 97434, 96772, 95855, 94522, 92573, 89854, 85710, 79515, 70912, 59139, 44043, 27376, 12817, 3892, 624, 41)
    } else {
      check=1
      dat = read.csv('AllCauseMortRatesFromUSLifeTablesByBirthYear-UpdateTo2008Table-Age120-whitefemales.csv',sep=",",header=F)
      # surv = c(100000, 99385, 99319, 99245, 99049, 98835, 98601, 98282, 97789, 97047, 95958, 94284, 91591, 87391, 81346, 72546, 59681, 42820, 24475, 9673, 2145, 201)
    } 
  }
  if(toupper(race) == "ALL") {
    if(toupper(sex) == "MALE") {
      check=1
      dat = read.csv('AllCauseMortRatesFromUSLifeTablesByBirthYear-UpdateTo2008Table-Age120-males.csv',sep=",",header=F)
      # surv = c(100000, 99095, 99008, 98890, 98426, 97746, 97112, 96382, 95384, 93931, 91800, 88862, 84478, 78083, 69350, 57572, 42683, 26473, 12447, 3847, 643, 46)
    } else {
      check=1
      dat = read.csv('AllCauseMortRatesFromUSLifeTablesByBirthYear-UpdateTo2008Table-Age120-females.csv',sep=",",header=F)
      # surv = c(100000, 99261, 99190, 99111, 98915, 98682, 98418, 98052, 97492, 96645, 95420, 93597, 90739, 86367, 80158, 71257, 58411, 41798, 23918, 9553, 2181, 217)
    }
  }

  if(check == 0) {warning("input not specified correctly"); return}

  dimnames(dat)[[2]]=birthyears
  dimnames(dat)[[1]]=range
  # integrate rates:
  cumhaz = cumsum(dat[(ages+1),(byear-1859)])
  surv=exp(-cumhaz)
  prob=1-surv
  return(list(age = approx(x=prob, y=ages, xout=runif(n), yleft=start, yright=(max(range)+1))$y, byear=byear, start=start, surv=surv))
}
  


