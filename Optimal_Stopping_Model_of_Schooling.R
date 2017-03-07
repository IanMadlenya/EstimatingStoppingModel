##########################################
## Preliminaries
##########################################
## 1- Work on data

# Do not forget to change the current directory
setwd("...")

## read data as csv (comma separated values)
x = read.csv(file = "g98_5years2.csv",sep=",")
x <- na.omit(x) #incomplete cases removed 

## open the dataset, dimensions, list of variables, summarize a variable, plots...
# fix(x) # open data
dim(x) #dimensions
names(x) #list of variables
table(x$sch_lev) # frequencies of level of schooling
#hist(x$wage_60) # histogram 
plot(x$sch_lev,x$wage_first) 

## reading the variables keeping identifier 
id	<- x$numenq  #identifier
S	<- cbind(id, x$sch_lev)
# father, reference x$blc_f
Xf	<- cbind(x$farm_f,x$trad_f,x$exec_f,x$tech_f,x$whc_f,x$noocc_f) 
# mother, reference x$blc_m
Xm	<- cbind(x$farm_m,x$trad_m,x$exec_m,x$tech_m,x$whc_m,x$noocc_m) 

# origin, reference x$orfr
Xo	<- cbind(x$oroecd,x$orother,x$foreign) 
# other individual variable
Xi	<- cbind(x$late6,x$urban,x$male)
# full time and transformed wages - now yearly
wf 	<- cbind(x$full_first, x$wage_first)
wag36	<- x$wage_36*12
w36	<- cbind(x$full_36,wag36)
wag60	<- x$wage_60*12
w60	<- cbind(x$full_60,wag60)
# employment and experience
datef	<- x$date_first
exp36	<- x$experience_36/12
e36	<- cbind(x$employed_36,exp36)
exp60	<- x$experience_60/12
e60	<- cbind(x$employed_60,exp60)

## 2- Descriptive work
## Independent linear estimations of Mincer equations on full employed

T <- 50 # number of periods

s	<- S[,2] #schooling
s1	<- (S[,2]==1)*1
s2	<- (S[,2]==2)*1
s3	<- (S[,2]==3)*1
s4	<- (S[,2]==4)*1
s5	<- (S[,2]==5)*1
s6	<- (S[,2]==6)*1
s7	<- (S[,2]==7)*1
s8	<- (S[,2]==8)*1
late	<- Xi[,1] # repeating variable
# experience
ex36	<- e36[,2] 
ex36sq	<- ex36*ex36
ex60	<- e60[,2]
ex60sq	<- ex60*ex60

# linear models

y36	<-log(w36[,2]/1000) # wage at month 36
U36	<- lm(y36 ~ s2+s3+s4+s5+s6+s7+s8 + late + ex36 + ex36sq, subset= w36[,1]==1)
c36	<- coefficients(U36)
summary(U36)

y60	<-log(w60[,2]/1000) # wage at month 36
U60	<- lm(y60 ~ s2+s3+s4+s5+s6+s7+s8 + late + ex60 + ex60sq, subset= w60[,1]==1)
c60	<- coefficients(U60)
summary(U60)


## 3- Value functions of labor

# regular age for schooling s
A 	<- function(s) {
  As <- 16*(s==1)+17*(s==2)+18*(s==3)+19*(s==4)+20*(s==5)+21*(s==6)+22*(s==7)+23*(s==8)
  return(As)
}


# simplified version of VW(s, As) :
# with 0.01 for experience and no squared experience

# VW	<- function(s,a,r,b){
#   VS36	<- c36[1]+c36[2]*(s==2)+c36[3]*(s==3)+c36[4]*(s==4)+c36[5]*(s==5)+c36[6]*(s==6)+c36[7]*(s==7)+c36[8]*(s==8)
#   VS	<- VS36*(b^(T-a)-1)/(b-1) #schooling
#   R36	<- c36[9]*(r==1)
#   R	<- R36*(b^(T-a)-1)/(b-1) #repeat
#   i	<- seq(0,(T-a-1),1)
#   E	<- sum(b^i*(c36[10]*i)) #experience
#   VW	<- VS + R + E
#   return(VW)
# }

VW	<- function(s,a,r,b){
  VS60	<- c60[1]+c60[2]*(s==2)+c60[3]*(s==3)+c60[4]*(s==4)+c60[5]*(s==5)+c60[6]*(s==6)+c60[7]*(s==7)+c60[8]*(s==8)
  VS	<- VS60*(b^(T-a)-1)/(b-1) #schooling
  R60	<- c60[9]*(r==1)
  R	<- R60*(b^(T-a)-1)/(b-1) #repeat
  i	<- seq(0,(T-a-1),1)
  E	<- sum(b^i*(c60[10]*i)) #experience
  VW	<- VS + R + E
  return(VW)
}


##########################################
## 1 - Maximum likelihood for Vs
##########################################

### value functions for each steps
## Px for probability of schooling in step x
## EVx for Emax (EV(x) in the problem set)
#########

P1 <- function(X,r,dx,ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta) {
  VS1 	<- X %*% dx+ds1+dl*(r==1) + beta*EV2(X,r,dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta)

  VW0 	<- 0
  c		<- VW0 - VS1
  return(1-pnorm(c))		
}

EV1 <- function(X,r,dx,ds1,ds2,ds3, ds4,ds5,ds6,ds7,ds8,dl,beta) {
  VS1		<- X %*% dx+ds1+dl*(r==1) + beta*EV2(X,r,dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta)
  VW0 	<- 0
  c		<- VW0 - VS1
  P1		<- 1-pnorm(c)
  V1T 	<- VS1 + dnorm(c)/(1-pnorm(c))
  EV		<- P1*V1T + (1-P1)*VW0
  return(EV)		
}




P2 <- function(X,r,dx,ds2,ds3, ds4,ds5,ds6,ds7,ds8,dl,beta) {
  VS2 	<- X %*% dx+ds2+dl*(r==1) + beta*EV3(X,r,dx,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta)
  VW1 	<- VW(s=1,a=A(1),r=r,b=beta)
  c		<- VW1 - VS2
  return(1-pnorm(c))		
}


EV2 <- function(X,r,dx,ds2,ds3, ds4,ds5,ds6,ds7,ds8,dl,beta) {
  VS2		<- X %*% dx+ds2+dl*(r==1) + beta*EV3(X,r,dx,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta)
  VW1 	<- VW(s=1,a=A(1),r=r,b=beta)
  c		<- VW1 - VS2
  P2		<- 1-pnorm(c)
  V2T 	<- VS2 + dnorm(c)/(1-pnorm(c))
  EV		<- P2*V2T + (1-P2)*VW1
  return(EV)		
}


P3 <- function(X,r,dx,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta) {
  VS3 	<- X %*% dx+ds3+dl*(r==1) + beta*EV4(X,r,dx,ds4,ds5,ds6,ds7,ds8,dl,beta)
  VW2 	<- VW(s=2,a=A(2),r=r,b=beta)
  c		<- VW2 - VS3
  return(1-pnorm(c))		
}

EV3 <- function(X,r,dx,ds3, ds4,ds5,ds6,ds7,ds8,dl,beta) {
  VS3		<- X %*% dx+ds3+dl*(r==1) + beta*EV4(X,r,dx,ds4,ds5,ds6,ds7,ds8,dl,beta)
  VW2 	<- VW(s=2,a=A(2),r=r,b=beta)
  c		<- VW2 - VS3
  P3		<- 1-pnorm(c)
  V3T 	<- VS3 + dnorm(c)/(1-pnorm(c))
  EV		<- P3*V3T + (1-P3)*VW2
  return(EV)		
}

P4 <- function(X,r,dx,ds4,ds5,ds6,ds7,ds8,dl,beta) {
  VS4 	<- X %*% dx+ds4+dl*(r==1) + beta*EV5(X,r,dx,ds5,ds6,ds7,ds8,dl,beta)
  VW3 	<- VW(s=3,a=A(3),r=r,b=beta)
  c		<- VW3 - VS4
  return(1-pnorm(c))		
}

EV4 <- function(X,r,dx,ds4, ds5,ds6,ds7,ds8,dl,beta) {
  VS4		<- X %*% dx+ds4+dl*(r==1) + beta*EV5(X,r,dx,ds5,ds6,ds7,d8,dl,beta)
  VW3 	<- VW(s=3,a=A(3),r=r,b=beta)
  c		<- VW3 - VS4
  P4		<- 1-pnorm(c)
  V4T 	<- VS4 + dnorm(c)/(1-pnorm(c))
  EV		<- P4*V4T + (1-P4)*VW3
  return(EV)		
}



P5 <- function(X,r,dx,ds5,ds6,ds7,d8,dl,beta) {
  VS5 	<- X %*% dx+ds5+dl*(r==1) + beta*EV6(X,r,dx,ds6,ds7,d8,dl,beta)
  VW4 	<- VW(s=4,a=A(4),r=r,b=beta)
  c		<- VW4 - VS5
  return(1-pnorm(c))		
}

EV5 <- function(X,r,dx,ds5,ds6,ds7,d8,dl,beta) {
  VS5		<- X %*% dx+ds5+dl*(r==1) + beta*EV6(X,r,dx,ds6,ds7,ds8,dl,beta)
  VW4 	<- VW(s=4,a=A(4),r=r,b=beta)
  c		<- VW4 - VS5
  P5		<- 1-pnorm(c)
  V5T 	<- VS5 + dnorm(c)/(1-pnorm(c))
  EV		<- P5*V5T + (1-P5)*VW4
  return(EV)		
}



P6 <- function(X,r,dx,ds6,ds7,ds8,dl,beta) {
  VS6 	<- X %*% dx+ds6+dl*(r==1) + beta*EV7(X,r,dx,ds7,ds8,dl,beta)
  VW5 	<- VW(s=5,a=A(5),r=r,b=beta)
  c		<- VW5 - VS6
  return(1-pnorm(c))		
}

EV6 <- function(X,r,dx,ds6,ds7,ds8,dl,beta) {
  VS6		<- X %*% dx+ds6+dl*(r==1) + beta*EV7(X,r,dx,ds7,ds8,dl,beta)
  VW5 	<- VW(s=5,a=A(5),r=r,b=beta)
  c		<- VW5 - VS6
  P6		<- 1-pnorm(c)
  V6T 	<- VS6 + dnorm(c)/(1-pnorm(c))
  EV		<- P6*V6T + (1-P6)*VW5
  return(EV)		
}



P7 <- function(X,r,dx,ds7,ds8,dl,beta) {
  VS7 	<- X %*% dx+ds7+dl*(r==1) + beta*EV8(X,r,dx,ds8,dl,beta)
  VW6 	<- VW(s=6,a=A(6),r=r,b=beta)
  c		<- VW6 - VS7
  return(1-pnorm(c))		
}

EV7 <- function(X,r,dx,ds7, ds8,dl,beta) {
  VS7		<- X %*% dx+ds7+dl*(r==1) + beta*EV8(X,r,dx,ds8,dl,beta)
  VW6 	<- VW(s=6,a=A(6),r=r,b=beta)
  c		<- VW6 - VS7
  P7		<- 1-pnorm(c)
  V7T 	<- VS7 + dnorm(c)/(1-pnorm(c))
  EV		<- P7*V7T + (1-P7)*VW6
  return(EV)		
}


P8 <- function(X,r,dx,ds8,dl,beta) {
  VS8 	<- X %*% dx+ds8+dl*(r==1) + beta*VW(s=8,a=A(8),r=r,b=beta)
  VW7 	<- VW(s=7,a=A(7),r=r,b=beta)
  c		<- VW7 - VS8
  return(1-pnorm(c))		
}

EV8 <- function(X,r,dx,ds8,dl,beta) {
  VS8		<- X %*% dx+ds8+dl*(r==1) + beta*VW(s=8,a=A(8),r=r,b=beta)
  VW7 	<- VW(s=7,a=A(7),r=r,b=beta)
  c		<- VW7 - VS8
  P8 		<- 1-pnorm(c)
  V8T 	<- VS8 + dnorm(c)/(1-pnorm(c))
  EV		<- P8*V8T + (1-P8)*VW7
  return(EV)		
}

### TO CONTINUE  : P7, EV7, P6... + logLik, etc.

######## full log L
fll <- function(dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta1) {
  X	<- cbind(Xf, Xm, Xo, Xi)
  r 	<- late
 # pi1   	<- P1(X,r,dx,ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta)=1
  pi1 <- 1
  pi2 <- P2(X,r,dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta1)
  pi3 <- P3(X,r,dx,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta1)
  pi4 <- P4(X,r,dx,ds4,ds5,ds6,ds7,ds8,dl,beta1)
  pi5 <- P5(X,r,dx,ds5,ds6,ds7,ds8,dl,beta1)
  pi6 <- P6(X,r,dx,ds6,ds7,ds8,dl,beta1)
  pi7 <- P7(X,r,dx,ds7,ds8,dl,beta1)
  pi8 <- P8(X,r,dx,ds8,dl,beta1)
  
  proba<-c(1,P2,P3,P4,P5,P6,P7,P8)
  
  # (...)
  li1 	<- pi1*(1-pi2)
  li2 <- pi1*pi2*(1-pi3)
  li3 <- pi1*pi2*pi3*(1-pi4)
  li4 <- pi1*pi2*pi3*pi4*(1-pi5)
  li5 <- pi1*pi2*pi3*pi4*pi5*(1-pi6)
  li6 <- pi1*pi2*pi3*pi4*pi5*pi6*(1-pi7)
  li7 <- pi1*pi2*pi3*pi4*pi5*pi6*pi7*(1-pi8)
  li8 <- pi1*pi2*pi3*pi4*pi5*pi6*pi7*pi8
  # (...)
  li	<- li1^(s==1)*li2^(s==2)*li3^(s==3)*li4^(s==4)*li5^(s==5)*li6^(s==6)*li7^(s==7)*li8^(s==8) # (...)
    logL  <- log(li)
  return(logL)
}

#Starting values of parameters

dx<-rep(0.5, 18)

ds2<-0.5
ds3<-0.5
ds4<-0.5
ds5<-0.5
ds6<-0.5
ds7<-0.5
ds8<-0.5
dl<-0.5
beta1<-0.6




sum(fll(dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta1))

# parameters as vector, beta estimated remove ds1
fllv <- function(p){
  dx	<- p[1:18]
  
  ds2	<- p[19]
  ds3	<- p[20]
  ds4	<- p[21]
  ds5	<- p[22]
  ds6	<- p[23]
  ds7	<- p[24]
  ds8	<- p[25]
  dl	<- p[26]
  #beta	<- p[27] beta fixed below
  beta	<- 0.5
  return(fll(dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta))
}

param_start1 <- c(dx, ds2, ds3, ds4, ds5, ds6, ds7, ds8, dl, beta1)
#sum(fllv(param_start1))

#parameter vector for q3 uniquement
param_startq3 <- c(dx, ds2, ds3, ds4, ds5, ds6, ds7, ds8, dl)
##### optimization routines ####

#First method with optim
# be careful : optim runs minimization 
flls <- function(par){
  - sum(fllv(par))
}

time 	<- proc.time()
est1	<- optim(param_startq3,flls, method = c("BFGS"),control=list(trace=2,maxit=2000), hessian=FALSE)
proc.time()- time
est1$convergence

time 	<- proc.time()
est2	<- optim(est1$par,flls, method = c("BFGS"),control=list(trace=2,maxit=2000), hessian=TRUE)
proc.time()- time
est2$convergence

#reltol=sqrt(.Machine$double.eps)
#reltol=10^(-4)

# computation of standard errors
hatp	<- est2$par
n	<- length(id)
HINV	<- (est2$hessian)^-1
S2 	<- sqrt(diag(HINV)/n)

#################
###second method with maxLik (works bad, inverse matrix for SE is not computed)
library(maxLik)

ttt <- proc.time()
opti1 <- maxBFGS(fllv, start = param_start1, print.level = 4, finalHessian = F)
ttt2 <- proc.time() - ttt
ttt2

param_opti1= opti1$estimate

#######Calculate the probability of achieving a specific schooling level

plp1 <- function(dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta1) {
  X	<- cbind(Xf, Xm, Xo, Xi)
  r 	<- late
  # pi1   	<- P1(X,r,dx,ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta)=1
  pi1 <- 1
  pi2 <- P2(X,r,dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta1)
  pi3 <- P3(X,r,dx,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta1)
  pi4 <- P4(X,r,dx,ds4,ds5,ds6,ds7,ds8,dl,beta1)
  pi5 <- P5(X,r,dx,ds5,ds6,ds7,ds8,dl,beta1)
  pi6 <- P6(X,r,dx,ds6,ds7,ds8,dl,beta1)
  pi7 <- P7(X,r,dx,ds7,ds8,dl,beta1)
  pi8 <- P8(X,r,dx,ds8,dl,beta1)
  
  li1 	<- pi1*(1-pi2)
  li2 <- pi1*pi2*(1-pi3)
  li3 <- pi1*pi2*pi3*(1-pi4)
  li4 <- pi1*pi2*pi3*pi4*(1-pi5)
  li5 <- pi1*pi2*pi3*pi4*pi5*(1-pi6)
  li6 <- pi1*pi2*pi3*pi4*pi5*pi6*(1-pi7)
  li7 <- pi1*pi2*pi3*pi4*pi5*pi6*pi7*(1-pi8)
  li8 <- pi1*pi2*pi3*pi4*pi5*pi6*pi7*pi8
  
  
  #proba1<-c(mean(pi1),mean(pi2),mean(pi3),mean(pi4),mean(pi5),mean(pi6),mean(pi7),mean(pi8))
  proba2<-c(mean(li1),mean(li2),mean(li3),mean(li4),mean(li5),mean(li6),mean(li7),mean(li8))
  return(proba2)
  }

plp2 <- function(p){
  dx	<- p[1:18]
  
  ds2	<- p[19]
  ds3	<- p[20]
  ds4	<- p[21]
  ds5	<- p[22]
  ds6	<- p[23]
  ds7	<- p[24]
  ds8	<- p[25]
  dl	<- p[26]
  #beta	<- p[27]
  beta<-0.5
  return(plp1(dx,ds2,ds3,ds4,ds5,ds6,ds7,ds8,dl,beta))
}
prtoplot<-plp2(est2$par)
more4<-prtoplot[4]+prtoplot[5]+prtoplot[6]+prtoplot[7]+prtoplot[8]

prtoplot

more4

#cbind(c(0.05159809,0.20224289,0.18006186,0.09424392,0.25223088,0.05534522,0.06889206,0.09538508),c(0.05173560,0.20275089,0.17895821,0.09317897,0.25279021,0.05581789,0.06972586,0.09504237),c(0.05173560,0.20275089,0.17895821,0.09317897,0.25279021,0.05581789,0.06972586,0.09504237),c(0.05219953,0.20397420,0.17716853,0.09147798,0.25273482,0.05683339,0.07189813,0.09371341),c(0.05244858,0.20452779,0.17675385,0.09101501,0.25242339,0.05716355,0.07276864,0.09289919))

#rm(list=ls())
##########################
save(file = "opti.txt", opti1)

# standard errors
sigma = sqrt(diag(solve(t(opti1$gradientObs) %*% opti1$gradientObs)))

#to load the estimated objects
load("opti.txt")

param_start2 <- param_opti1

# # R net search : http://rseek.org
# # http://cran.r-project.org/web/views/Optimization.html
# 
# 
# ######## STARTING VALUES 
# # population at risk and hazard rates
# pr <- function(tt){
#   p1 <- n
#   p2 <- n-tt[1]
#   p3 <- n-tt[1]-tt[2]
#   p4 <- n-tt[1]-tt[2]-tt[3]
#   p5 <- n-tt[1]-tt[2]-tt[3]-tt[4]
#   p6 <- tt[6]+tt[7]+tt[8]
#   p7 <- tt[7]+tt[8]
#   p8 <- tt[8]	
#   c(1-p2/p1, p2/p1, p3/p2, p4/p3, p5/p4, p6/p5, p7/p6, p8/p7)
# }
# # descriptives
# ts		<- table(s)
# # hazard rates 
# haz	<- pr(ts)
# 
# # for utilities of schooling
# n		<- length(id)
# XXX 	<- cbind(Xf, Xm, Xo, Xi)
# rrr 	<- late
# ddd 	<- array(0, c(18, 1)) 
# dl		<- -0.02
# beta 	<- 0.95
# 
# # TEST FOR VALUE ds8
# ds8	<- -0.32
# 
# tVX8	<- XXX %*% ddd +ds8+dl*(rrr==1)
# tV8 	<- beta*VW(s=8,a=A(8),r=rrr,b=beta)
# tV7 	<- VW(s=7,a=A(7),r=rrr,b=beta)
# tc8	<- tV7 - tV8-tVX8
# tp8	<- 1-pnorm(tc8) 
# #mean(tc8)
# ## to compare to hazard rate 1301/(749+1301)=0.63
# haz[8]
# mean(tp8)
# 
# # and so on for ds7, ds6 etc.
# # caution with ds1 which is similar to a constant in the model
# ### to be compared with the sum of phi_1*beta^i on the whole period
# ind	<- seq(0,(T-1),1)
# sum(beta^ind)*c36[1] 