# This code was run during winter break 2009-2010.


require(mstate)
require(diagram)
require(msm)
##------------          

# SOLID Example: Rates triple for first 30 days following MACE, double for remainder of year, then return back to 0.075
#p.Z1              = c(.075*3*.8, .075*2.8, .075*.8)    # Annual rates for non-fatal MI/Stroke
#p.Z2              = c(.075*3*.2, .075*2.2, .075*.2)    # Anndual rates for CVD
#tp1                = 30                                # Timepoint for rate switch 1
#tp2                = 365                               # Timepoint for rate switch 2
#narm              = 5750                               # patients per arm
#Exp.No.ACM        = 371                                # Following determines the NCVD rate (500 assumed in Stability 371 = 500*11500/15500)
#prop.NCVD         = 1/3                                # Proportion of MACE that are no fatal
#Exp.med.exposure  = 3.00                               # Median Exposure (Approximated from East)
#EnrollYears       = 24/12                              # enrollment period 
#eshape            = 1/3
#HR.adj            = .845

#---------- # LOCK  OPEN ---------------------------
#--------- 3-piece exponential cdf

#F0      <- function(t, lambda1=-log(1-.15/365), lambda2=-log(1-.10/365), lambda3=-log(1-.20/365), t1=365, t2=2*365){
#        1 - (exp(-1*lambda1*t)*(t < t1) +
#             exp(-1*lambda1*t1)*exp(-1*lambda2*(t - t1))*(t >= t1 & t < t2) +
#             exp(-1*lambda1*t1)*exp(-1*lambda2*(t2 - t1))*exp(-lambda3*(t- t2))*(t >= t2))}

#----------Inverse of 3-piece exponential cdf
#F0_inv  <- function(t, lambda1=-log(1-.15/365), lambda2=-log(1-.10/365), lambda3=-log(1-.20/365),  t1=365, t2=2*365){
#        log(1-t)              /-lambda1                         *(t <  F0(t1, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, t1=t1,t2=t2)) +
#      ((log(1-t) + lambda1*t1)/-lambda2 + t1)                   *(t >= F0(t1, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, t1=t1,t2=t2) & t < F0(t2, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, t1=t1,t2=t2)) +
#      ((log(1-t) + lambda1*t1 + lambda2*(t2-t1))/-lambda3 + t2) *(t >= F0(t2, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, t1=t1,t2=t2)) }
# LOCK CLOSE---------------------------

get.event.hist <- function(
    p.Z1              = c(.075*3*.8, .075*2*.8, .075*.8),  # Annual rates for non-fatal MI/Stroke
    p.Z2              = c(.075*3*.2, .075*2*.2, .075*.2),  # Anndual rates for CVD
    tp1               = 30,                          # Timepoint for rate switch 1
    tp2               = 365,                         # Timepoint for rate switch 2
    narm              = 5750,                        # patients per arm
    Exp.No.ACM        = 371,                         # Following determines the NCVD rate
    prop.NCVD         = 1/3,                         
    Exp.med.exposure  = 3,                           
    EnrollYears       = 24,                          #  enrollment period
    eshape            = 1/3,                         # enrollment shape
    HR.adj            = 1) {                         # Take HR.adj = 1 for PBO and .845 for SB

EnrollDay     <- ceiling((runif(narm)^(eshape))*365*EnrollYears)
p.Z12         <- rep((Exp.No.ACM*prop.NCVD)/(2*narm*Exp.med.exposure),3)
# Determines the exponential parameters
Z12.rate      <- -log(1-p.Z12)             # Drug has no effect on NON-CV Death
Z1.rate       <- -log(1-p.Z1)*HR.adj       # Take HR.adj = 1 for PBO and .845 for SB
Z2.rate       <- -log(1-p.Z2)*HR.adj       # Take HR.adj = 1 for PBO and .845 for SB

#------------------ Transitions into Non-CV ACM Absorbing States------------#
E1E7 <- rpexp(n=narm, rate=Z12.rate, t=c(0, tp1, tp2))
E2E7 <- rpexp(n=narm, rate=Z12.rate, t=c(0, tp1, tp2))
E3E7 <- rpexp(n=narm, rate=Z12.rate, t=c(0, tp1, tp2))
E4E7 <- rpexp(n=narm, rate=Z12.rate, t=c(0, tp1, tp2))
E5E7 <- rpexp(n=narm, rate=Z12.rate, t=c(0, tp1, tp2))


#------------------ Transitions into CVDeath Absorbing States------------#
# Hardcoding E1E6 to install the Protocol assumptions on time to first FATAL MACE
E1E6 <- rpexp(n=narm, rate=c(-log(1-.075*.2), -log(1-.075*.2), -log(1-.035*.2))*HR.adj, t=c(0, tp1, tp2))
E2E6 <- rpexp(n=narm, rate=Z2.rate, t=c(0, tp1, tp2))
E3E6 <- rpexp(n=narm, rate=Z2.rate, t=c(0, tp1, tp2))
E4E6 <- rpexp(n=narm, rate=Z2.rate, t=c(0, tp1, tp2))
E5E6 <- rpexp(n=narm, rate=Z2.rate, t=c(0, tp1, tp2))


#--------------Walk Through Transient States
# Round 1: Non-Fatal Stroke/MI rate is the same regardless of transitional state.
# Hardcoding E1E2 to install the Protocol assumptions on time to first NON-FATAL MACE
E1E2 <- rpexp(n=narm, rate=c(-log(1-.075*.8), -log(1-.075*.8), -log(1-.075*.8))*HR.adj, t=c(0, tp1, tp2))
E2E3 <- rpexp(n=narm, rate=Z1.rate, t=c(0, tp1, tp2))
E3E4 <- rpexp(n=narm, rate=Z1.rate, t=c(0, tp1, tp2))
E4E5 <- rpexp(n=narm, rate=Z1.rate, t=c(0, tp1, tp2))
# This just gives a quick check on transition times.
#par(mfrow=c(3,5))
#hist(E1E7, xlim=c(0,2000), ylim=c(0,.005),breaks=10, freq=F)
#hist(E2E7, xlim=c(0,2000), ylim=c(0,.005),breaks=10, freq=F)
#hist(E3E7, xlim=c(0,2000), ylim=c(0,.005),breaks=10, freq=F)
#hist(E4E7, xlim=c(0,2000), ylim=c(0,.005),breaks=10, freq=F)
#hist(E5E7, xlim=c(0,2000), ylim=c(0,.005),breaks=10, freq=F)
#hist(E1E6, xlim=c(0,1500), ylim=c(0,.015),breaks=10, freq=F)
#hist(E2E6, xlim=c(0,1500), ylim=c(0,.015),breaks=10, freq=F)
#hist(E3E6, xlim=c(0,1500), ylim=c(0,.015),breaks=10, freq=F)
#hist(E4E6, xlim=c(0,1500), ylim=c(0,.015),breaks=10, freq=F)
#hist(E5E6, xlim=c(0,1500), ylim=c(0,.015),breaks=10, freq=F)
#hist(E1E2, xlim=c(0,500), ylim=c(0,.1),breaks=10, freq=F)
#hist(E2E3, xlim=c(0,500), ylim=c(0,.1),breaks=10, freq=F)
#hist(E3E4, xlim=c(0,500), ylim=c(0,.1),breaks=10, freq=F)
#hist(E4E5, xlim=c(0,500), ylim=c(0,.1),breaks=10, freq=F)
#

# So each row corresponds to times generated for a single patient
trans.time <- cbind(E1E2, E1E6, E1E7, E2E3, E2E6, E2E7, E3E4, E3E6, E3E7, E4E5, E4E6, E4E7, E5E6, E5E7)
colnames(trans.time) <- c("E1E2", "E1E6", "E1E7", "E2E3", "E2E6", "E2E7", "E3E4", "E3E6", "E3E7", "E4E5", "E4E6", "E4E7", "E5E6", "E5E7")

E1trans <- c()
E2trans <- c()
E3trans <- c()
E4trans <- c()
E5trans <- c()

for(i in 1:narm){
E1trans <- c(E1trans,which.min(trans.time[i,1:3]))                 # Starting from 1 which state tranisition time occurs first
E2trans <- c(E2trans,which.min(trans.time[i,4:6]))                 # Starting from 2 which state transition time occurs first
E3trans <- c(E3trans,which.min(trans.time[i,7:9]))                 # Starting from 3 ....
E4trans <- c(E4trans,which.min(trans.time[i,10:12]))               # Starting from 4 ....
E5trans <- c(E5trans,which.min(trans.time[i,13:14]))               # Starting from 5 ....
}

E1trans <- (E1trans==1)*2 +(E1trans==2)*6+(E1trans==3)*7           # E1trans now holds which state we go to next
E2trans <- (E2trans==1)*3 +(E2trans==2)*6+(E2trans==3)*7           # Assuming we start at 2, which state do we go to next
E3trans <- (E3trans==1)*4 +(E3trans==2)*6+(E3trans==3)*7           # Assuming we start at 3, ...
E4trans <- (E4trans==1)*5 +(E4trans==2)*6+(E4trans==3)*7           # Assuming we start at 4, ....
E5trans <-                 (E5trans==1)*6+(E5trans==2)*7           # Assuming we start at 5, ....

trans <- (cbind(E1trans, E2trans, E3trans, E4trans, E5trans))
rownames(trans) <- c()
colnames(trans) <- paste("Next State from", 1:5)

# Following code determines subject paths to absorbing states
c1 <- trans[,1]
c2 <- c()
for(i in 1:narm)  {
  if (c1[i] %in% c(6,7))
    c2[i] <- NA
  if (c1[i] == 2)
    c2[i] <- trans[i,2] }

c3 <- c()
for(i in 1:narm)  {
  if (c2[i] %in% c(6,7, NA))
    c3[i] <- NA
  if (is.na(c2[i]) == F & c2[i] == 3)
    c3[i] <- trans[i,3] }

c4 <- c()
for(i in 1:narm)  {
  if (c3[i] %in% c(6,7, NA))
    c4[i] <- NA
  if (is.na(c3[i])==F & c3[i] == 4)
    c4[i] <- trans[i,4]  }

c5 <- c()
for(i in 1:narm)  {
  if (c4[i] %in% c(6,7, NA))
    c5[i] <- NA
  if (is.na(c4[i])==F & c4[i] == 5)
    c5[i] <- trans[i,5]   }

subj.paths            <- cbind(c1, c2, c3, c4, c5)
colnames(subj.paths)  <- c(paste("Step", 1:5,sep=""))
# Subj.paths holds the actual paths of the subjects, assuming we have the ability to follow til infinity

no.trans <- c()
for(i in 1:narm){
  no.trans <- c(no.trans, 5 - sum(is.na(subj.paths[i,]==T)))}
# no.trans holds the number of transitions to get to an absorbing state, assuming we have the ability to follow til infinity

pull.times2 <- rep(NA, narm)
pull.times3 <- rep(NA, narm)
pull.times4 <- rep(NA, narm)
pull.times5 <- rep(NA, narm)

# Grab the tranistion time from State 1 to state (2, 6, 7) according to subj.path (i.e., based on which transit time is minimum)
pull.times  <- trans.time[,1]*(subj.paths[,1]  == 2) + trans.time[,2]*(subj.paths[,1]  == 6) + trans.time[,3]*(subj.paths[,1]  == 7) + .0001
pull.times2 <- trans.time[,4]*(subj.paths[,2]  == 3) + trans.time[,5]*(subj.paths[,2]  == 6) + trans.time[,6]*(subj.paths[,2]  == 7) + .0001
pull.times3 <- trans.time[,7]*(subj.paths[,3]  == 4) + trans.time[,8]*(subj.paths[,3]  == 6) + trans.time[,9]*(subj.paths[,3]  == 7) + .0001
pull.times4 <- trans.time[,10]*(subj.paths[,4] == 5) + trans.time[,11]*(subj.paths[,4] == 6) + trans.time[,12]*(subj.paths[,4] == 7) + .0001
pull.times5 <-                                         trans.time[,13]*(subj.paths[,5] == 6) + trans.time[,14]*(subj.paths[,5] == 7) + .0001

times           <- cbind(pull.times, pull.times2, pull.times3, pull.times4, pull.times5)
colnames(times) <- c(paste("SubjTransTime",1:5,sep=""))
# So times holds the transit time spent in a state before moving to next state.  
# Each row corresponds to a subject
# The subject's path is given by subj.path
# NA's indicate the subject had entered an absorbing state at last transit.

study.time1 <- EnrollDay/365
study.time2 <- EnrollDay/365 + times[,1]
study.time3 <- EnrollDay/365 + times[,1] + times[,2]
study.time4 <- EnrollDay/365 + times[,1] + times[,2] + times[,3]
study.time5 <- EnrollDay/365 + times[,1] + times[,2] + times[,3] + times[,4]
study.time6 <- EnrollDay/365 + times[,1] + times[,2] + times[,3] + times[,4] + times[,5]

# study.time assumes we observed full data
study.time            <- cbind(study.time1, study.time2, study.time3,  study.time4,  study.time5, study.time6)
colnames(study.time)  <- c("EnrollDay", "Trans1", "Trans2", "Trans3", "Trans4", "Trans5")
# So study.times holds the clock time relative to study start. 
# First column is enrollment time.
# Next columns are cumulative times, relative to the study start 

# The subject's path is given by subj.path
# NA's indicate the subject had entered an absorbing state at last transit.

list(subj.paths = subj.paths, times = times, study.time = study.time) 
# Subj.paths holds the actual paths of the subjects, assuming we  follow til infinity
# So times holds the transit time spent in a state before moving to next state, assuming follow til infinity.  
# So study.times includes the enrollment time and shifts gap times from "times" relative to study start. 
}

#----------------------------------
#----------------------------------

# Event history calls get.event.hist twice to get pbo and sb data
# It determines the study duration (i.e., time to capture 1500 events
# 
event.hist <- function(
    p.Z1              = c(.075*3*.8, .075*2*.8, .075*.8),  # Annual rates for non-fatal MI/Stroke
    p.Z2              = c(.075*3*.2, .075*2*.2, .075*.2),  # Anndual rates for CVD
    tp1               = 30,                         # Timepoint for rate switch 1
    tp2               = 365,                        # Timepoint for rate switch 2
    narm              = 5750,                       # patients per arm
    Exp.No.ACM        = 371,                        # Following determines the NCVD rate
    prop.NCVD         = 1/3,
    Exp.med.exposure  = 3.00,
    EnrollYears       = 24/12,                      # enrollment period
    eshape            = 1/3,
    HR.adjustment     = 0.845
    ) {        

pbo <- get.event.hist(
    p.Z1              = p.Z1,  # Annual rates for non-fatal MI/Stroke
    p.Z2              = p.Z2,  # Anndual rates for CVD
    tp1               = tp1,                         # Timepoint for rate switch 1
    tp2               = tp2,                        # Timepoint for rate switch 2
    narm              = narm,                       # patients per arm
    Exp.No.ACM        = Exp.No.ACM,                        # Following determines the NCVD rate
    prop.NCVD         = prop.NCVD ,
    Exp.med.exposure  = Exp.med.exposure,
    EnrollYears       = EnrollYears,                      # enrollment period
    eshape            = eshape,
    HR.adj            = 1)

sb <- get.event.hist(
    p.Z1              = p.Z1,  # Annual rates for non-fatal MI/Stroke
    p.Z2              = p.Z2,  # Anndual rates for CVD
    tp1               = tp1,                         # Timepoint for rate switch 1
    tp2               = tp2,                        # Timepoint for rate switch 2
    narm              = narm,                       # patients per arm
    Exp.No.ACM        = Exp.No.ACM,                        # Following determines the NCVD rate
    prop.NCVD         = prop.NCVD ,
    Exp.med.exposure  = Exp.med.exposure,
    EnrollYears       = EnrollYears,                      # enrollment period
    eshape            = eshape,
    HR.adj            = HR.adjustment)


# Determine the study duration required to capture 1500 events
for(i in seq(2.5,7,.025)){

# Check on first instances of MACE (non-fatal = 2, fatal = 6) and count when these occur at times < i*365
  StudyDur <- i*365
    if(nrow(pbo$subj.paths[pbo$subj.paths[,1] %in% c(2,6) & pbo$study.time[,2] < i,])  + 
       nrow( sb$subj.paths[sb$subj.paths[,1]  %in% c(2,6) & sb$study.time[,2] < i,]) > 1500) break}


  pbo$study.time2 <- pbo$study.time3 <- (pmin(pbo$study.time,StudyDur/365))
  sb$study.time2 <- sb$study.time3 <- (pmin(sb$study.time,StudyDur/365))

for(i in 1:narm){
if (is.na(pbo$study.time2[i,2]) == F & pbo$study.time2[i,2] == StudyDur/365)   pbo$study.time3[i,3] <- NA
if (is.na(sb$study.time2[i,2])  == F &  sb$study.time2[i,2] == StudyDur/365)   sb$study.time3[i,3]  <- NA  }

for(i in 1:narm){
if (is.na(pbo$study.time2[i,3]) == F & pbo$study.time2[i,3] == StudyDur/365)  pbo$study.time3[i,4] <- NA
if (is.na(pbo$study.time2[i,3]) == T)                                         pbo$study.time3[i,4] <- NA

if (is.na(sb$study.time2[i,3]) == F & sb$study.time2[i,3] == StudyDur/365)  sb$study.time3[i,4] <- NA
if (is.na(sb$study.time2[i,3]) == T)                                        sb$study.time3[i,4] <- NA  }

for(i in 1:narm){
if (is.na(pbo$study.time2[i,4]) == F & pbo$study.time2[i,4] == StudyDur/365)  pbo$study.time3[i,5] <- NA
if (is.na(pbo$study.time2[i,4]) == T)                                         pbo$study.time3[i,5] <- NA

if (is.na(sb$study.time2[i,4]) == F & sb$study.time2[i,4] == StudyDur/365)  sb$study.time3[i,5] <- NA
if (is.na(sb$study.time2[i,4]) == T)                                        sb$study.time3[i,5] <- NA }

for(i in 1:narm){
if (is.na(pbo$study.time2[i,5]) == F & pbo$study.time2[i,5] == StudyDur/365)  pbo$study.time3[i,6] <- NA
if (is.na(pbo$study.time2[i,5]) == T)                                         pbo$study.time3[i,6] <- NA

if (is.na(sb$study.time2[i,5]) == F & sb$study.time2[i,5] == StudyDur/365)  sb$study.time3[i,6] <- NA
if (is.na(sb$study.time2[i,5]) == T)                                        sb$study.time3[i,6] <- NA  }

# So study.time2 contains superfluous censored values, whereas study.time3 cleans these up to that they match
# better with observed subj.paths AND censoring due to study duration (Which was based on capturing 1500 MACE)



#  Observed Transition Matrix
pbo$mat.obs <- matrix(NA,7,7)
  pbo$mat.obs[1,2] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,2] != StudyDur/365 & is.na(pbo$study.time3[,2]) ==F ,]),0)
  pbo$mat.obs[1,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 6 & pbo$study.time3[,2] != StudyDur/365 & is.na(pbo$study.time3[,2]) ==F ,]),0)
  pbo$mat.obs[1,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 7 & pbo$study.time3[,2] != StudyDur/365 & is.na(pbo$study.time3[,2]) ==F ,]),0)
  pbo$mat.obs[2,3] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,3] != StudyDur/365 & is.na(pbo$study.time3[,3]) ==F & pbo$subj.paths[,2]==3,]),0)
  pbo$mat.obs[2,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,3] != StudyDur/365 & is.na(pbo$study.time3[,3]) ==F & pbo$subj.paths[,2]==6,]),0)
  pbo$mat.obs[2,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,3] != StudyDur/365 & is.na(pbo$study.time3[,3]) ==F & pbo$subj.paths[,2]==7,]),0)
  pbo$mat.obs[3,4] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,4] != StudyDur/365 & is.na(pbo$study.time3[,4]) ==F & pbo$subj.paths[,3]==4,]),0)
  pbo$mat.obs[3,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,4] != StudyDur/365 & is.na(pbo$study.time3[,4]) ==F & pbo$subj.paths[,3]==6,]),0)
  pbo$mat.obs[3,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,4] != StudyDur/365 & is.na(pbo$study.time3[,4]) ==F & pbo$subj.paths[,3]==7,]),0)
  pbo$mat.obs[4,5] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,5] != StudyDur/365 & is.na(pbo$study.time3[,5]) ==F & pbo$subj.paths[,4]==5,]),0)
  pbo$mat.obs[4,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,5] != StudyDur/365 & is.na(pbo$study.time3[,5]) ==F & pbo$subj.paths[,4]==6,]),0)
  pbo$mat.obs[4,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,5] != StudyDur/365 & is.na(pbo$study.time3[,5]) ==F & pbo$subj.paths[,4]==7,]),0)
  pbo$mat.obs[5,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,6] != StudyDur/365 & is.na(pbo$study.time3[,6]) ==F & pbo$subj.paths[,5]==6,]),0)
  pbo$mat.obs[5,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,6] != StudyDur/365 & is.na(pbo$study.time3[,6]) ==F & pbo$subj.paths[,5]==7,]),0)
  rownames(pbo$mat.obs) <- paste("From E", 1:7,sep="")
  colnames(pbo$mat.obs) <- paste("To E", 1:7,sep="")

sb$mat.obs <- matrix(NA,7,7)
  sb$mat.obs[1,2] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,2] != StudyDur/365 & is.na(sb$study.time3[,2]) ==F ,]),0)
  sb$mat.obs[1,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 6 & sb$study.time3[,2] != StudyDur/365 & is.na(sb$study.time3[,2]) ==F ,]),0)
  sb$mat.obs[1,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 7 & sb$study.time3[,2] != StudyDur/365 & is.na(sb$study.time3[,2]) ==F ,]),0)
  sb$mat.obs[2,3] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,3] != StudyDur/365 & is.na(sb$study.time3[,3]) ==F & sb$subj.paths[,2]==3,]),0)
  sb$mat.obs[2,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,3] != StudyDur/365 & is.na(sb$study.time3[,3]) ==F & sb$subj.paths[,2]==6,]),0)
  sb$mat.obs[2,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,3] != StudyDur/365 & is.na(sb$study.time3[,3]) ==F & sb$subj.paths[,2]==7,]),0)
  sb$mat.obs[3,4] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,4] != StudyDur/365 & is.na(sb$study.time3[,4]) ==F & sb$subj.paths[,3]==4,]),0)
  sb$mat.obs[3,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,4] != StudyDur/365 & is.na(sb$study.time3[,4]) ==F & sb$subj.paths[,3]==6,]),0)
  sb$mat.obs[3,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,4] != StudyDur/365 & is.na(sb$study.time3[,4]) ==F & sb$subj.paths[,3]==7,]),0)
  sb$mat.obs[4,5] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,5] != StudyDur/365 & is.na(sb$study.time3[,5]) ==F & sb$subj.paths[,4]==5,]),0)
  sb$mat.obs[4,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,5] != StudyDur/365 & is.na(sb$study.time3[,5]) ==F & sb$subj.paths[,4]==6,]),0)
  sb$mat.obs[4,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,5] != StudyDur/365 & is.na(sb$study.time3[,5]) ==F & sb$subj.paths[,4]==7,]),0)
  sb$mat.obs[5,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,6] != StudyDur/365 & is.na(sb$study.time3[,6]) ==F & sb$subj.paths[,5]==6,]),0)
  sb$mat.obs[5,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,6] != StudyDur/365 & is.na(sb$study.time3[,6]) ==F & sb$subj.paths[,5]==7,]),0)
  rownames(sb$mat.obs) <- paste("From E", 1:7,sep="")
  colnames(sb$mat.obs) <- paste("To E", 1:7,sep="")

# Full transition matrix assuming no censoring
pbo$mat.Full <- matrix(NA,7,7)
  pbo$mat.Full[1,2] <- sum(pbo$subj.paths[,1] == 2)
  pbo$mat.Full[1,6] <- sum(pbo$subj.paths[,1] == 6)
  pbo$mat.Full[1,7] <- sum(pbo$subj.paths[,1] == 7)
  pbo$mat.Full[2,3] <- sum(pbo$subj.paths[,2] == 3,na.rm=T)
  pbo$mat.Full[2,6] <- sum(pbo$subj.paths[,2] == 6,na.rm=T)
  pbo$mat.Full[2,7] <- sum(pbo$subj.paths[,2] == 7,na.rm=T)
  pbo$mat.Full[3,4] <- sum(pbo$subj.paths[,3] == 4,na.rm=T)
  pbo$mat.Full[3,6] <- sum(pbo$subj.paths[,3] == 6,na.rm=T)
  pbo$mat.Full[3,7] <- sum(pbo$subj.paths[,3] == 7,na.rm=T)
  pbo$mat.Full[4,5] <- sum(pbo$subj.paths[,4] == 5,na.rm=T)
  pbo$mat.Full[4,6] <- sum(pbo$subj.paths[,4] == 6,na.rm=T)
  pbo$mat.Full[4,7] <- sum(pbo$subj.paths[,4] == 7,na.rm=T)
  pbo$mat.Full[5,6] <- sum(pbo$subj.paths[,5] == 6,na.rm=T)
  pbo$mat.Full[5,7] <- sum(pbo$subj.paths[,5] == 7,na.rm=T)
  rownames(pbo$mat.Full) <- paste("From E", 1:7,sep="")
  colnames(pbo$mat.Full) <- paste("To E", 1:7,sep="")

sb$mat.Full <- matrix(NA,7,7)
  sb$mat.Full[1,2] <- sum(sb$subj.paths[,1] == 2)
  sb$mat.Full[1,6] <- sum(sb$subj.paths[,1] == 6)
  sb$mat.Full[1,7] <- sum(sb$subj.paths[,1] == 7)
  sb$mat.Full[2,3] <- sum(sb$subj.paths[,2] == 3,na.rm=T)
  sb$mat.Full[2,6] <- sum(sb$subj.paths[,2] == 6,na.rm=T)
  sb$mat.Full[2,7] <- sum(sb$subj.paths[,2] == 7,na.rm=T)
  sb$mat.Full[3,4] <- sum(sb$subj.paths[,3] == 4,na.rm=T)
  sb$mat.Full[3,6] <- sum(sb$subj.paths[,3] == 6,na.rm=T)
  sb$mat.Full[3,7] <- sum(sb$subj.paths[,3] == 7,na.rm=T)
  sb$mat.Full[4,5] <- sum(sb$subj.paths[,4] == 5,na.rm=T)
  sb$mat.Full[4,6] <- sum(sb$subj.paths[,4] == 6,na.rm=T)
  sb$mat.Full[4,7] <- sum(sb$subj.paths[,4] == 7,na.rm=T)
  sb$mat.Full[5,6] <- sum(sb$subj.paths[,5] == 6,na.rm=T)
  sb$mat.Full[5,7] <- sum(sb$subj.paths[,5] == 7,na.rm=T)
  rownames(sb$mat.Full) <- paste("From E", 1:7,sep="")
  colnames(sb$mat.Full) <- paste("To E", 1:7,sep="")

# This counts the number censored.  Note the difference here is use of '!=' 
pbo$mat.cen <- matrix(NA,7,7)
  pbo$mat.cen[1,2] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,3] == StudyDur/365 & is.na(pbo$study.time3[,3]) ==F ,]),0)
  pbo$mat.cen[1,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 6 & pbo$study.time3[,3] == StudyDur/365 & is.na(pbo$study.time3[,3]) ==F ,]),0)
  pbo$mat.cen[1,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 7 & pbo$study.time3[,3] == StudyDur/365 & is.na(pbo$study.time3[,3]) ==F ,]),0)
  pbo$mat.cen[2,3] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,3] == StudyDur/365 & is.na(pbo$study.time3[,3]) ==F & pbo$subj.paths[,2]==3,]),0)
  pbo$mat.cen[2,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,3] == StudyDur/365 & is.na(pbo$study.time3[,3]) ==F & pbo$subj.paths[,2]==6,]),0)
  pbo$mat.cen[2,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,3] == StudyDur/365 & is.na(pbo$study.time3[,3]) ==F & pbo$subj.paths[,2]==7,]),0)
  pbo$mat.cen[3,4] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,4] == StudyDur/365 & is.na(pbo$study.time3[,4]) ==F & pbo$subj.paths[,3]==4,]),0)
  pbo$mat.cen[3,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,4] == StudyDur/365 & is.na(pbo$study.time3[,4]) ==F & pbo$subj.paths[,3]==6,]),0)
  pbo$mat.cen[3,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,4] == StudyDur/365 & is.na(pbo$study.time3[,4]) ==F & pbo$subj.paths[,3]==7,]),0)
  pbo$mat.cen[4,5] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,5] == StudyDur/365 & is.na(pbo$study.time3[,5]) ==F & pbo$subj.paths[,4]==5,]),0)
  pbo$mat.cen[4,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,5] == StudyDur/365 & is.na(pbo$study.time3[,5]) ==F & pbo$subj.paths[,4]==6,]),0)
  pbo$mat.cen[4,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,5] == StudyDur/365 & is.na(pbo$study.time3[,5]) ==F & pbo$subj.paths[,5]==7,]),0)
  pbo$mat.cen[5,6] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,6] == StudyDur/365 & is.na(pbo$study.time3[,6]) ==F & pbo$subj.paths[,5]==6,]),0)
  pbo$mat.cen[5,7] <- max(nrow(pbo$study.time3[pbo$subj.paths[,1] == 2 & pbo$study.time3[,6] == StudyDur/365 & is.na(pbo$study.time3[,6]) ==F & pbo$subj.paths[,5]==7,]),0)
  rownames(pbo$mat.cen) <- paste("From E", 1:7,sep="")
  colnames(pbo$mat.cen) <- paste("To E", 1:7,sep="")

sb$mat.cen <- matrix(NA,7,7)
  sb$mat.cen[1,2] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,3] == StudyDur/365 & is.na(sb$study.time3[,3]) ==F ,]),0)
  sb$mat.cen[1,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 6 & sb$study.time3[,3] == StudyDur/365 & is.na(sb$study.time3[,3]) ==F ,]),0)
  sb$mat.cen[1,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 7 & sb$study.time3[,3] == StudyDur/365 & is.na(sb$study.time3[,3]) ==F ,]),0)
  sb$mat.cen[2,3] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,3] == StudyDur/365 & is.na(sb$study.time3[,3]) ==F & sb$subj.paths[,2]==3,]),0)
  sb$mat.cen[2,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,3] == StudyDur/365 & is.na(sb$study.time3[,3]) ==F & sb$subj.paths[,2]==6,]),0)
  sb$mat.cen[2,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,3] == StudyDur/365 & is.na(sb$study.time3[,3]) ==F & sb$subj.paths[,2]==7,]),0)
  sb$mat.cen[3,4] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,4] == StudyDur/365 & is.na(sb$study.time3[,4]) ==F & sb$subj.paths[,3]==4,]),0)
  sb$mat.cen[3,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,4] == StudyDur/365 & is.na(sb$study.time3[,4]) ==F & sb$subj.paths[,3]==6,]),0)
  sb$mat.cen[3,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,4] == StudyDur/365 & is.na(sb$study.time3[,4]) ==F & sb$subj.paths[,3]==7,]),0)
  sb$mat.cen[4,5] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,5] == StudyDur/365 & is.na(sb$study.time3[,5]) ==F & sb$subj.paths[,4]==5,]),0)
  sb$mat.cen[4,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,5] == StudyDur/365 & is.na(sb$study.time3[,5]) ==F & sb$subj.paths[,4]==6,]),0)
  sb$mat.cen[4,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,5] == StudyDur/365 & is.na(sb$study.time3[,5]) ==F & sb$subj.paths[,5]==7,]),0)
  sb$mat.cen[5,6] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,6] == StudyDur/365 & is.na(sb$study.time3[,6]) ==F & sb$subj.paths[,5]==6,]),0)
  sb$mat.cen[5,7] <- max(nrow(sb$study.time3[sb$subj.paths[,1] == 2 & sb$study.time3[,6] == StudyDur/365 & is.na(sb$study.time3[,6]) ==F & sb$subj.paths[,5]==7,]),0)
  rownames(sb$mat.cen) <- paste("From E", 1:7,sep="")
  colnames(sb$mat.cen) <- paste("To E", 1:7,sep="")

pbo.counts <- c(11500 - sum(c(pbo$mat.obs[1,2] + pbo$mat.obs[1,6],
pbo$mat.obs[2,3] + pbo$mat.obs[2,6],
pbo$mat.obs[3,4] + pbo$mat.obs[3,6],
pbo$mat.obs[4,5] + pbo$mat.obs[4,6], 
pbo$mat.obs[5,6] )), c(pbo$mat.obs[1,2] + pbo$mat.obs[1,6],
pbo$mat.obs[2,3] + pbo$mat.obs[2,6],
pbo$mat.obs[3,4] + pbo$mat.obs[3,6],
pbo$mat.obs[4,5] + pbo$mat.obs[4,6],
pbo$mat.obs[5,6] ))
sb.counts <- c(11500 - sum(c(sb$mat.obs[1,2] + sb$mat.obs[1,6],
sb$mat.obs[2,3] + sb$mat.obs[2,6],
sb$mat.obs[3,4] + sb$mat.obs[3,6],
sb$mat.obs[4,5] + sb$mat.obs[4,6],
sb$mat.obs[5,6] )), c(sb$mat.obs[1,2] + sb$mat.obs[1,6],
sb$mat.obs[2,3] + sb$mat.obs[2,6],
sb$mat.obs[3,4] + sb$mat.obs[3,6],
sb$mat.obs[4,5] + sb$mat.obs[4,6],
sb$mat.obs[5,6] ))

list(studydur = StudyDur/365, PBO=pbo, SB = sb, counts = rbind(pbo.counts, sb.counts)) }


# Simulation Testing
# Single Test with seed = 123 worked: time =  1.182217 mins  
# Test with 10 under V1 assumptions
#------------------------------

p.Z1.v <- list()
p.Z2.v <- list()
p.Z1.v[[1]] <- c(.075*1*.8, .075*1*.8, .075*1*.8)     #1x, 1x, 1x        
p.Z1.v[[2]] <- c(.075*2*.8, .075*1*.8, .075*1*.8)     #2x, 1x, 1x          
p.Z1.v[[3]] <- c(.075*2*.8, .075*2*.8, .075*1*.8)     #2x, 2x, 1x       
p.Z1.v[[4]] <- c(.075*3*.8, .075*2*.8, .075*1*.8)     #3x, 2x, 1x          
p.Z1.v[[5]] <- c(.075*3*.8, .075*2*.8, .075*2*.8)     #3x, 2x, 2x           
p.Z1.v[[6]] <- c(.075*3*.8, .075*3*.8, .075*2*.8)     #3x, 3x, 2x           

p.Z2.v[[1]] <- c(.075*1*.2, .075*1*.2, .075*1*.2)
p.Z2.v[[2]] <- c(.075*2*.2, .075*1*.2, .075*1*.2)
p.Z2.v[[3]] <- c(.075*2*.2, .075*2*.2, .075*1*.2)
p.Z2.v[[4]] <- c(.075*3*.2, .075*2*.2, .075*1*.2)
p.Z2.v[[5]] <- c(.075*3*.2, .075*2*.2, .075*2*.2)
p.Z2.v[[6]] <- c(.075*3*.2, .075*3*.2, .075*2*.2)


big.sim <- function(version.num, simsize, HRsend=.845,startme=0)    {

tmat <- (matrix(NA, 7, 7))
tmat[1,c(2,6,7)] <- c(2,6,7)
tmat[2,c(3,6,7)] <- c(3,6,7)
tmat[3,c(4,6,7)] <- c(4,6,7)
tmat[4,c(5,6,7)] <- c(5,6,7)
tmat[5,c(6,7)] <- c(6,7)
narm    <-5750

start <- Sys.time()
par(las=1)
study.dur.report              <- c()
ag.common.report              <- c()
pwp.total.common.report       <- c()
pwp.gap.common.report         <- c()
pwp.total.uncommon.report     <- c()
pwp.gap.uncommon.report       <- c()
wlw.common.report             <- c()
wlw.uncommon.report           <- c()
poisreg.report                <- c()
ttest1.report                 <- c()
ttest2.report                 <- c()
chisq.report                  <- c()
landmark.report               <- c()
pois.reg.report               <- c()
simple.report                 <- c()
for(iteration in 1:simsize){
  #-------------------- Displays graph to alert about remaining time in sim
  plot(iteration, (difftime(Sys.time(),start,units="sec"))/(iteration)*(simsize-(iteration))/60,
       xlim=c(0,simsize), ylim=c(0,simsize*1.25),cex=.25, ylab="Minutes Left", xlab="Simulation", main="Time Left in Sim",axes=F)
       axis(1)
       axis(2,at=seq(0,simsize*1.5, by=30),cex.axis=.4)
       abline(lty=3,h=seq(0,simsize*1.5, by=60)) 
  par(new=T) 
  #--------------------------
  
  temp <- event.hist(
      p.Z1              = p.Z1.v[[version.num-6]],                    # Annual rates for non-fatal MI/Stroke
      p.Z2              = p.Z2.v[[version.num-6]],                    # Anndual rates for CVD
      tp1               = 30,                         # Timepoint for rate switch 1
      tp2               = 365,                        # Timepoint for rate switch 2
      narm              = 5750,                       # patients per arm
      Exp.No.ACM        = 371,                        # Following determines the NCVD rate
      prop.NCVD         = 1/3,
      Exp.med.exposure  = 3,
      EnrollYears       = 24/12,                      # enrollment period
      eshape            = 1/3,
      HR.adjustment     = HRsend
      ) 
  
  
  # Write the recurrent event data to a csv file - replace 999999 with iteration number in the full sim
  write.csv(rbind(cbind(rep(1,nrow(temp[[2]]$study.time3)),temp[[2]]$study.time3),
                  cbind(rep(2,nrow(temp[[3]]$study.time3)),temp[[3]]$study.time3)), 
                  paste("c:/recurrenteventdummies/version ",version.num,"/clintrial_", iteration+startme,".csv", sep=""))
  
                   
  
  study.dur.report    <- c(study.dur.report, temp[[1]])
  
  #-----------------------RECURRENT EVENTS ANALYSIS--------------------------#
  # Reorganizes the data in the style of bladder1
  # ID, Tstart, Tstop, status, trt, visit
  # where there's only rows avail when the subject is at risk for next event.
  # E.g., sub with 0,1,2,3,4 events has 1,2,3,4,5 rows.
  
  temp2 <- temp[[2]]$study.time3
  reorg.v1.pbo <- c()
  for(i in 1:nrow(temp2)){
  for(j in 0:(4 - sum(is.na(temp2[i,])))){
  reorg.v1.pbo <- rbind(reorg.v1.pbo,
  c(i, temp2[i,j+1], temp2[i,j+2], temp2[i,j+2] < temp[[1]], 0, j+1))   }}
  
  reorg.v2.pbo <- c()
  for(i in 1:nrow(temp2)){
  for(j in 1:4){
  reorg.v2.pbo <- rbind(reorg.v2.pbo,
  c(i, min(na.rm=T, temp2[i,j+1], temp[[1]]) - temp2[i,1],  min(na.rm=T, temp2[i,j+1], temp[[1]]) < temp[[1]], 0, j))}}
  
  temp2 <- temp[[3]]$study.time3
  reorg.v1.sb <- c()
  for(i in 1:nrow(temp2)){
  for(j in 0:(4 - sum(is.na(temp2[i,])))){
  reorg.v1.sb <- rbind(reorg.v1.sb,
  c(i+nrow(reorg.v1.pbo), temp2[i,j+1], temp2[i,j+2], temp2[i,j+2] < temp[[1]], 1, j+1))   }}
  
  reorg.v2.sb <- c()
  for(i in 1:nrow(temp2)){
  for(j in 1:4){
  reorg.v2.sb <- rbind(reorg.v2.sb,
  c(i, min(na.rm=T, temp2[i,j+1], temp[[1]]) - temp2[i,1],  min(na.rm=T, temp2[i,j+1], temp[[1]]) < temp[[1]], 1, j))}}
  
  reorg.v1 <- as.data.frame(rbind(reorg.v1.pbo, reorg.v1.sb))
  reorg.v2 <- as.data.frame(rbind(reorg.v2.pbo, reorg.v2.sb))
  colnames(reorg.v1) <- c("ID", "Tstart", "Tstop", "status", "trt", "visit")
  colnames(reorg.v2) <- c("ID", "Followup", "status", "trt", "visit")

if(nrow(reorg.v2[reorg.v2$visit==4 & reorg.v2$status==1,]) < 10) reorg.v2 <- reorg.v2[reorg.v2$visit < 4,]
if(nrow(reorg.v2[reorg.v2$visit==3 & reorg.v2$status==1,]) < 10) reorg.v2 <- reorg.v2[reorg.v2$visit < 3,]
  
  
  #COMMON TREATMENT EFFECT
  ag.common <- coxph(Surv(Tstart,Tstop,status)~as.factor(trt)+cluster(ID),data=reorg.v1, method="breslow")
  ag.common.report <- rbind(ag.common.report, summary(ag.common)$coef)
  
  #Conditional model: uses time intervals; only at risk for event 2 after experience event 1; tends to be much more efficient than WLW; stratify on event type like WLW
  # NOTE from R: If there are no tied death times all the methods are equivalent. Nearly all Cox regression programs 
  # use the Breslow method by default, but not this one. The Efron approximation is used as the default here, as it is 
  # much more accurate when dealing with tied death times, and is as efficient computationally. The exact method 
  # computes the exact partial likelihood, which is equivalent to a conditional logistic model. If there are a large 
  # number of ties the computational time will be excessive. 
  
  pwp.total.common          <- coxph(Surv(Tstart,Tstop,status)~as.factor(trt)+strata(visit)+cluster(ID),data=reorg.v1, method="breslow")
  pwp.total.common.report   <- rbind(pwp.total.common.report , summary(pwp.total.common)$coef)
  pwp.gap.common            <- coxph(Surv(Tstop-Tstart,status)~as.factor(trt)+strata(visit)+cluster(ID),data=reorg.v1,method="breslow")
  pwp.gap.common.report     <- rbind(pwp.gap.common.report, summary(pwp.gap.common)$coef)

  # UNCOMMON TREATMENT EFFECT
  pwp.total.uncommon        <- coxph(Surv(Tstart,Tstop,status)~as.factor(trt*visit)+strata(visit)+cluster(ID),data=reorg.v1, method="breslow")
  pwp.total.uncommon.report <- rbind(pwp.total.uncommon.report, summary(pwp.total.uncommon)$coef)
  pwp.gap.uncommon          <- coxph(Surv(Tstop-Tstart,status)~as.factor(trt*visit)+strata(visit)+cluster(ID),data=reorg.v1,method="breslow")
  pwp.gap.uncommon.report   <- rbind(pwp.gap.uncommon.report, summary(pwp.gap.uncommon)$coef)
  
  # WLW
  wlw.common                <- coxph(Surv(Followup,status)~as.factor(trt)+strata(visit)+cluster(ID),method="breslow", data=reorg.v2)
  wlw.common.report         <- rbind(wlw.common.report, summary(wlw.common)$coef)
  
  #WLW uncommon
  wlw.uncommon <- coxph(Surv(Followup,status)~as.factor(trt*visit)+strata(visit)+cluster(ID),method="breslow", data=reorg.v2)
  wlw.uncommon.report <- rbind(wlw.uncommon.report, summary(wlw.uncommon)$coef)

  
  #-----------------------LANDMARK ANALYSIS--------------------------#
  # TIME TO FIRST EVENT
  
  temp2 <- temp[[2]]$study.time3
  time.to.first.event.pbo <- cbind(
  Tstart=temp2[,1], Tstop=temp2[,2], status = temp2[,2] < temp[[1]])
  time.to.first.event.pbo <- cbind(time.to.first.event.pbo, treatment=(rep(0, nrow(time.to.first.event.pbo ))))
  
  temp2 <- temp[[3]]$study.time3
  time.to.first.event.sb <- cbind(
  Tstart=temp2[,1], Tstop=temp2[,2], status = temp2[,2] < temp[[1]])
  time.to.first.event.sb <- cbind(time.to.first.event.sb, treatment=(rep(1, nrow(time.to.first.event.sb ))))
  
  time.to.first.event <- as.data.frame(rbind(time.to.first.event.pbo, time.to.first.event.sb))
  #head(time.to.first.event.pbo)
  #head(time.to.first.event.sb)
  c1 <- coxph(Surv(time.to.first.event$Tstart, time.to.first.event$Tstop, time.to.first.event$status) ~ as.factor(time.to.first.event$treatment))
  
  # TIME TO SECOND EVENT
  temp2 <- temp[[2]]$study.time3[is.na(temp[[2]]$study.time3[,3]) == F,]
  time.to.second.event.pbo <- cbind(
  Tstart=temp2[,2], Tstop=temp2[,3], status = temp2[,3] < temp[[1]])
  time.to.second.event.pbo <- cbind(time.to.second.event.pbo, treatment=(rep(0, nrow(time.to.second.event.pbo ))))
  
  temp2 <- temp[[3]]$study.time3[is.na(temp[[3]]$study.time3[,3]) == F,]
  time.to.second.event.sb <- cbind(
  Tstart=temp2[,2], Tstop=temp2[,3], status = temp2[,3] < temp[[1]])
  time.to.second.event.sb <- cbind(time.to.second.event.sb, treatment=(rep(1, nrow(time.to.second.event.sb ))))
  
  time.to.second.event <- as.data.frame(rbind(time.to.second.event.pbo, time.to.second.event.sb))
  #head(time.to.second.event.pbo)
  #head(time.to.second.event.sb)
  c2 <- coxph(Surv(time.to.second.event$Tstart, time.to.second.event$Tstop, time.to.second.event$status) ~as.factor(time.to.second.event$treatment))
  
  # TIME TO Third EVENT


if (nrow(reorg.v2[reorg.v2$visit==3 & reorg.v2$status==1,]) >= 10)   {
  temp2 <- temp[[2]]$study.time3[is.na(temp[[2]]$study.time3[,4]) == F,]
  time.to.third.event.pbo <- cbind(
  Tstart=temp2[,3], Tstop=temp2[,4], status = temp2[,4] < temp[[1]])
  time.to.third.event.pbo <- cbind(time.to.third.event.pbo, treatment=rep(0, nrow(time.to.third.event.pbo )))
  
  temp2 <- temp[[3]]$study.time3[is.na(temp[[3]]$study.time3[,4]) == F,]
  time.to.third.event.sb <- cbind(
  Tstart=temp2[,3], Tstop=temp2[,4], status = temp2[,4] < temp[[1]])
  time.to.third.event.sb <- cbind(time.to.third.event.sb, treatment=(rep(1, nrow(time.to.third.event.sb ))))
  
  time.to.third.event <- as.data.frame(rbind(time.to.third.event.pbo, time.to.third.event.sb))
  #head(time.to.third.event.pbo)
  #head(time.to.third.event.sb)
  c3 <- coxph(Surv(time.to.third.event$Tstart, time.to.third.event$Tstop, time.to.third.event$status) ~ as.factor(time.to.third.event$treatment))
  # pick off summary(c1)$coef
  landmark <-  cbind(1:3, rbind(summary(c1)$coef,summary(c2)$coef,summary(c3)$coef))
           }
  
  
if (nrow(reorg.v2[reorg.v2$visit==3 & reorg.v2$status==1,]) < 10) {  
  landmark <-  cbind(1:3, rbind(summary(c1)$coef,summary(c2)$coef,rep(NA, 5)))}
    
  landmark.report <- rbind(landmark.report, landmark)
  
  #-----------------------POISSON REGRESSION ANALYSIS--------------------------#
  total.events.pbo <- c(rep(0, temp[[4]][1,1]), rep(1, temp[[4]][1,2]), rep(2, temp[[4]][1,3]), rep(3, temp[[4]][1,4]), rep(4, temp[[4]][1,5]), rep(5, temp[[4]][1,6]))
  total.events.sb  <- c(rep(0, temp[[4]][2,1]), rep(1, temp[[4]][2,2]), rep(2, temp[[4]][2,3]), rep(3, temp[[4]][2,4]), rep(4, temp[[4]][2,5]), rep(5, temp[[4]][2,6]))
  totals           <- as.data.frame(cbind(events = c(total.events.pbo,total.events.sb), treatment = factor(c(rep(0, length(total.events.pbo)), rep(1, length(total.events.sb))))))
  
  
  
  pois.reg <-  summary(glm(events ~ treatment, family = "poisson", data =totals))$coefficients[2,]
  pois.reg.report <- rbind(pois.reg.report, pois.reg)
  
  t.test2 <- NA
  if(sd(totals[totals$events>1 & totals$treatment == 1,]$events) != 0 ||
     sd(totals[totals$events>1 & totals$treatment == 2,]$events) != 0 ){
  t.test2 <- t.test(events~treatment, data=totals[totals$events>1,])$p.value}
  simple.report <- rbind(simple.report, c(t.test(events~treatment, data=totals)$p.value, t.test2, chisq.test(with(totals, table(treatment, events)))$p.value))
  

write.csv(append=F,  study.dur.report            ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/studydur",startme,".csv",sep=""))
write.csv(append=F,  ag.common.report           ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/ag.common",startme,".csv",sep=""))
write.csv(append=F,  pwp.total.common.report    ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/pwp.total.common",startme,".csv",sep=""))
write.csv(append=F,  pwp.gap.common.report      ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/pwp.gap.common",startme,".csv",sep=""))
write.csv(append=F,  pwp.total.uncommon.report  ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/pwp.total.uncommon",startme,".csv",sep=""))
write.csv(append=F,  pwp.gap.uncommon.report    ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/pwp.gap.uncommon",startme,".csv",sep=""))
write.csv(append=F,  wlw.common.report          ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/wlw.common",startme,".csv",sep=""))
write.csv(append=F, wlw.uncommon.report        ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/wlw.uncommon",startme,".csv",sep=""))
write.csv(append=F,  pois.reg.report            ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/poisreg",startme,".csv",sep=""))
write.csv(append=F, simple.report              ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/simple",startme,".csv", sep=""))
write.csv(append=F, landmark.report            ,   paste("c:/recurrenteventdummies/version ",version.num,"/results/landmark",startme,".csv",sep=""))
}

stop <- Sys.time()
list(stop - start, temp[[2]], temp[[3]])    
}


 

  

# Initial code caused occasional errors: Each time simulation kicks its due to:
#     Error in Surv(Tstart, Tstop, status) : Stop time must be > start time
# Added an offset of .0001 in the pull.times in get.event.hist function ; jumping to test.v3 to see if it works
# This did the trick.  Rerunning all simulations


set.seed(111222)                                            
test.v1 <- big.sim(version.num=1, simsize=1000, HRsend=.845,startme=0)
set.seed(222333)                       
test.v2 <- big.sim(version.num=2, simsize=1000, HRsend=.845)
set.seed(333444)
test.v3 <- big.sim(version.num=3, simsize=1000, HRsend=.845)
set.seed(444555)
test.v4 <- big.sim(version.num=4, simsize=1000, HRsend=.845)
set.seed(555666)
test.v5 <- big.sim(version.num=5, simsize=1000, HRsend=.845)
set.seed(666777)
test.v6 <- big.sim(version.num=6, simsize=1000, HRsend=.845)


# checking under the null hypothesis... added 6 to version.num in code above

set.seed(111223)                                            
test.v1 <- big.sim(version.num=7, simsize=1000, HRsend=1,startme=0)
par(new=F)
set.seed(222334)                       
test.v2 <- big.sim(version.num=8, simsize=1000, HRsend=1)
par(new=F)
set.seed(333445)
par(new=F)
test.v3 <- big.sim(version.num=9, simsize=1000, HRsend=1)
set.seed(444556)
par(new=F)
test.v4 <- big.sim(version.num=10, simsize=1000, HRsend=1)
set.seed(555667)
par(new=F)
test.v5 <- big.sim(version.num=11, simsize=1000, HRsend=1)
set.seed(666778)
par(new=F)
test.v6 <- big.sim(version.num=12, simsize=1000, HRsend=1)






ag.check <- list()
pwp.total.common.check <- list()
pwp.gap.common.check <- list()
pwp.total.uncommon.check <- list()
pwp.gap.uncommon.check <- list()
wlw.common.check <- list()
wlw.uncommon.check <- list()
pois.reg.check <- list()
simple.check <- list()
landmark.check <- list()
  for(i in 1:12){
  ag.check[[i]]                     <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/ag.common0.csv",sep=""))
  pwp.total.common.check[[i]]       <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/pwp.total.common0.csv",sep=""))
  pwp.gap.common.check[[i]]         <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/pwp.gap.common0.csv",sep=""))  
  pwp.total.uncommon.check[[i]]   <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/pwp.total.uncommon0.csv",sep="")) 
  pwp.gap.uncommon.check[[i]]   <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/pwp.gap.uncommon0.csv",sep="")) 
  wlw.common.check[[i]]           <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/wlw.common0.csv",sep="")) 
  wlw.uncommon.check[[i]]         <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/wlw.uncommon0.csv",sep="")) 
  pois.reg.check[[i]]             <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/poisreg0.csv",sep="")) 
  simple.check[[i]]               <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/simple0.csv",sep="")) 
  landmark.check[[i]]             <- read.csv(paste("c:/recurrenteventdummies/version ",i,"/results/landmark0.csv",sep="")) 
  }


versions <- rbind(
c(.075*1, .075*1, .075*1),     #1x, 1x, 1x        
c(.075*2, .075*1, .075*1),     #2x, 1x, 1x          
c(.075*2, .075*2, .075*1),     #2x, 2x, 1x       
c(.075*3, .075*2, .075*1),     #3x, 2x, 1x          
c(.075*3, .075*2, .075*2),     #3x, 2x, 2x           
c(.075*3, .075*3, .075*2))     #3x, 3x, 2x

colnames(versions) <- c("1-30","31-365", "365+")
rownames(versions) <- paste("Version", 1:6)

#-------------------COLLECTING THE RESULTS ------------------------------#




pdf(file = "C:/Superduper Sim Results.pdf",
    width = 12, height = 8.5,
    onefile = TRUE, family = "Helvetica",
    title = "R Graphics Output", fonts = NULL, version = "1.1",
    paper = "special")

storeit <- list()


#---------------------------- Anderson Gill
holdme <- ag.check 
results <- c()
par(mfrow=c(2,3), ask=F)
for(i in 1:12){
hist(holdme[[i]]$p, main="Anderson-Gill P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = paste("Prop P-values < .05:", 
round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4)))
results <- c(results,round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4))
}
storeit <- list(ag.check=results)
#---------------------------- 

#----------------------------  PWP - Total Time - Common Effect
holdme <- pwp.total.common.check 
results <- c()
par(mfrow=c(2,3), ask=F)
for(i in 1:12){
hist(holdme[[i]]$p, main="PWP - Total Time, Common Effect P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = paste("Prop P-values < .05:", 
round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4)))
results <- c(results,round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4))
}
storeit$pwp.total.common.check = results
#---------------------------- 



# PWP - Total Time - Uncommon Effect

results <- c()
for(j in 1:12){
if (j <= 6) label <- c("Alternative")
if (j > 6) label <- c("Null")
par(mfrow=c(2,3), ask=F)
for(i in 1:5){
holdme <- pwp.total.uncommon.check[[j]][pwp.total.uncommon.check[[j]][,1] == paste("as.factor(trt * visit)",i,sep=""),]$p
if(sum(is.na(holdme)) < length(holdme)){
hist(holdme, main = "PWP - Total Time - Uncommon Effect P-value", xlab="", freq=F, sub=
paste("Version", j*(j <= 6) + (j-6)*(j > 6), "Treatment Effect", i, "Under", label))
legend("topright", bty="n", legend = c(
paste("Missing:", 1000 - sum(is.na(holdme) == F)), 
paste("Prop P-values < .05:", round(sum(holdme < .05, na.rm=T)/length(holdme),4))))     
results <- rbind(results, c(j, i, 1000 - sum(is.na(holdme) == F), round(sum(holdme < .05, na.rm=T)/length(holdme),4)))
}}}

storeit$pwp.total.uncommon.check=results


# PWP - Total Time - Common Effect
holdme <- pwp.gap.common.check 
results <- c()
par(mfrow=c(2,3), ask=F)
for(i in 1:12){
hist(holdme[[i]]$p, main="PWP - Gap Time, Common Effect P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = paste("Prop P-values < .05:", 
round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4)))
results <- c(results,round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4))
}

storeit$pwp.gap.common.check = results


# PWP - Total Time - Uncommon Effect

results <- c()
for(j in 1:12){
if (j <= 6) label <- c("Alternative")
if (j > 6) label <- c("Null")
par(mfrow=c(2,3), ask=F)
for(i in 1:5){
holdme <- pwp.gap.uncommon.check[[j]][pwp.gap.uncommon.check[[j]][,1] == paste("as.factor(trt * visit)",i,sep=""),]$p
if(sum(is.na(holdme)) < length(holdme)){
hist(holdme, main = "PWP - Gap Time - Uncommon Effect P-value", xlab="", freq=F, sub=
paste("Version", j*(j <= 6) + (j-6)*(j > 6), "Treatment Effect", i, "Under", label))
legend("topright", bty="n", legend = c(
paste("Missing:", 1000 - sum(is.na(holdme) == F)), 
paste("Prop P-values < .05:", round(sum(holdme < .05, na.rm=T)/length(holdme),4))))     
results <- rbind(results, c(j, i, 1000 - sum(is.na(holdme) == F), round(sum(holdme < .05, na.rm=T)/length(holdme),4)))
}}}

storeit$pwp.gap.uncommon.check=results






# WLW - Total Time - Common Effect
holdme <- wlw.common.check 
results <- c()
par(mfrow=c(2,3), ask=F)
for(i in 1:12){
hist(holdme[[i]]$p, main="WLW - Common Effect P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = paste("Prop P-values < .05:", 
round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4)))
results <- c(results,round(sum(holdme[[i]]$p < .05)/length(holdme[[i]]$p),4))
}

storeit$wlw.common.check =results



# WLW - Total Time - Uncommon Effect
results <- c()
for(j in 1:12){
if (j <= 6) label <- c("Alternative")
if (j > 6) label <- c("Null")
par(mfrow=c(2,3), ask=F)
for(i in 1:5){
holdme <- wlw.uncommon.check[[j]][wlw.uncommon.check[[j]][,1] == paste("as.factor(trt * visit)",i,sep=""),]$p
if(sum(is.na(holdme)) < length(holdme)){
hist(holdme, main = "WLW - Uncommon Effect P-value", xlab="", freq=F, sub=
paste("Version", j*(j <= 6) + (j-6)*(j > 6), "Treatment Effect", i, "Under", label))
legend("topright", bty="n", legend = c(
paste("Missing:", 1000 - sum(is.na(holdme) == F)), 
paste("Prop P-values < .05:", round(sum(holdme < .05, na.rm=T)/length(holdme),4))))
results <- rbind(results, c(j, i, 1000 - sum(is.na(holdme) == F), round(sum(holdme < .05, na.rm=T)/length(holdme),4)))
}}}

storeit$wlw.uncommon.check =results


# Poisson Regression
results <- c()
holdme <- pois.reg.check
par(mfrow=c(2,3), ask=F)
for(i in 1:12){
hist(holdme[[i]]$Pr...z.., main="Poisson Regression P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = paste("Prop P-values < .05:", 
round(sum(holdme[[i]]$Pr...z.. < .05)/length(holdme[[i]]$Pr...z..),4)))
results <- c(results,round(sum(holdme[[i]]$Pr...z..< .05)/length(holdme[[i]]$Pr...z..),4))
}

storeit$pois.reg.check =results


holdme <- simple.check
results <- c()
par(mfrow=c(2,3), ask=F)
for(i in 1:12){
hist(holdme[[i]][,2], main="T-test on MACE Incidence - All Patients - P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = paste("Prop P-values < .05:", 
round(sum(holdme[[i]][,2] < .05)/length(holdme[[i]][,2]),4)))
results <- c(results,round(sum(holdme[[i]][,2] < .05)/length(holdme[[i]][,2]),4))
}

storeit$simple.check1 =results

results <- c()
par(mfrow=c(2,3), ask=F)
for(i in 1:12){
hist(holdme[[i]][,3], main="T-test on MACE Incidence - Patients with 2+ MACE - P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = c(
paste("Missing:", 1000 - sum(is.na(holdme[[i]][,3]) == F)), 
paste("Prop P-values < .05:", round(sum(holdme[[i]][,3] < .05, na.rm=T)/length(holdme[[i]][,3]),4))))
results <- rbind(results, c(1000 - sum(is.na(holdme[[i]][,3]) == F), 
round(sum(holdme[[i]][,3] < .05, na.rm=T)/length(holdme[[i]][,3]),4))) 
}

storeit$simple.check2 =results


par(mfrow=c(2,3), ask=F)
results <- c()
for(i in 1:12){
hist(holdme[[i]][,4], main="Chi-Sq test on MACE Incidence - P-values", xlab="", freq = F, sub=
paste("Version", c(1:6,1:6)[i],"Under the", c(rep("Alternative",6), rep("Null", 6))[i])) 
legend("topright", bty="n", legend = c(
paste("Missing:", 1000 - sum(is.na(holdme[[i]][,4]) == F)), 
paste("Prop P-values < .05:", round(sum(holdme[[i]][,4] < .05, na.rm=T)/length(holdme[[i]][,4]),4))))
results <- c(results, round(sum(holdme[[i]][,4] < .05, na.rm=T)/length(holdme[[i]][,4]),4)) 
}
storeit$simple.check3 =results

head(landmark.check[[1]])
results <- c()
# Landmark - Total Time - Uncommon Effect
for(j in 1:12){
if (j <= 6) label <- c("Alternative")
if (j > 6) label <- c("Null")
par(mfrow=c(2,3), ask=F)
for(i in 1:3){
holdme <- landmark.check[[j]][landmark.check[[j]][,2] == i,]$p
if(sum(is.na(holdme)) < length(holdme)){
hist(holdme, main = "Landmark Analysis P-value", xlab="", freq=F, sub=
paste("Version", j*(j <= 6) + (j-6)*(j > 6), "Treatment Effect", i, "Under", label))
legend("topright", bty="n", legend = c(
paste("Missing:", 1000 - sum(is.na(holdme) == F)), 
paste("Prop P-values < .05:", round(sum(holdme < .05, na.rm=T)/length(holdme),4))))
results <- rbind(results, c(j, i, 1000 - sum(is.na(holdme) == F), round(sum(holdme < .05, na.rm=T)/length(holdme),4)))
}}}
storeit$landmark.check =results

names(storeit)


creport <- rbind(
storeit[[1]], 
storeit[[2]], 
storeit[[3]][storeit[[3]][,2] == 1,][,4],
storeit[[3]][storeit[[3]][,2] == 2,][,4],
storeit[[3]][storeit[[3]][,2] == 3,][,4],
storeit[[3]][storeit[[3]][,2] == 3,][,3],
storeit[[3]][storeit[[3]][,2] == 4,][,4],
storeit[[3]][storeit[[3]][,2] == 4,][,3],
storeit[[4]], 
storeit[[5]][storeit[[5]][,2] == 1,][,4],
storeit[[5]][storeit[[5]][,2] == 2,][,4],
storeit[[5]][storeit[[5]][,2] == 3,][,4],
storeit[[5]][storeit[[5]][,2] == 3,][,3],
storeit[[5]][storeit[[5]][,2] == 4,][,4],
storeit[[5]][storeit[[5]][,2] == 4,][,3],
storeit[[6]],                    
storeit[[7]][storeit[[7]][,2] == 1,][,4],
storeit[[7]][storeit[[7]][,2] == 2,][,4],
storeit[[7]][storeit[[7]][,2] == 3,][,4],
storeit[[7]][storeit[[7]][,2] == 3,][,3],
storeit[[7]][storeit[[7]][,2] == 4,][,4],
storeit[[7]][storeit[[7]][,2] == 4,][,3],
storeit[[8]],
storeit[[9]],
storeit[[11]],
storeit[[12]][storeit[[12]][,2] == 1,][,4],
storeit[[12]][storeit[[12]][,2] == 2,][,4],
storeit[[12]][storeit[[12]][,2] == 3,][,4],
storeit[[12]][storeit[[12]][,2] == 3,][,3])


rownames(creport) <- 
c(
"Anderson-Gill", 
"PWP Total Common", 
"PWP Total Uncommon 1",
"PWP Total Uncommon 2",
"PWP Total Uncommon 3",
"   Missing",
"PWP Total Uncommon 4",
"   Missing",
"PWP Gap Common",
"PWP Gap Uncommon 1",
"PWP Gap Uncommon 2",
"PWP Gap Uncommon 3",
"   Missing",
"PWP Gap Uncommon 4",
"   Missing",
"WLW Common", 
"WLW Uncommon 1",
"WLW Uncommon 2",
"WLW Uncommon 3",
"   Missing",
"WLW Uncommon 4",
"   Missing",
"Poisson Regression",
"T-test",
"Chi-Square",
"Landmark 1",
"Landmark 2",
"Landmark 3",
"   Missing")

write.csv(creport, "c:/creport.csv")



dev.off()




## Digraphs
#make.graphics <- function(arminfo){
#tmat <- (matrix(NA, 7, 7))
#tmat[1,c(2,6,7)] <- c(2,6,7)
#tmat[2,c(3,6,7)] <- c(3,6,7)
#tmat[3,c(4,6,7)] <- c(4,6,7)
#tmat[4,c(5,6,7)] <- c(5,6,7)
#tmat[5,c(6,7)] <- c(6,7)
#trans.n <- rbind(
#nrow(arminfo[arminfo$trans %in% c(1) & arminfo$status==1,]),             #folks with observable 1st mace
#nrow(arminfo[arminfo$trans %in% c(2) & arminfo$status==1,]),             #folks with observable fatal 1st mace
#nrow(arminfo[arminfo$trans %in% c(3) & arminfo$status==1,]),              # folks who then die from non-CV death
#nrow(arminfo[arminfo$trans %in% c(4) & arminfo$status==1,]),                # folks with observable 2nd non-fatal Mace
#nrow(arminfo[arminfo$trans %in% c(5) & arminfo$status==1,]),                # folks with observable 2nd fatal MACE
#nrow(arminfo[arminfo$trans %in% c(6) & arminfo$status==1,]),                # folks with censored time to 2nd mace due to observed nonCV death
#nrow(arminfo[arminfo$trans %in% c(7) & arminfo$status==1,]),                # folks 3rd non-fatal MACE
#nrow(arminfo[arminfo$trans %in% c(8) & arminfo$status==1,]),                # folks 3rd fatal MACE
#nrow(arminfo[arminfo$trans %in% c(9) & arminfo$status==1,]),                # folks with censored time to 3rd mace due to observed nonCV death
#nrow(arminfo[arminfo$trans %in% c(10) & arminfo$status==1,]),                # folks 4th non-fatal MACE
#nrow(arminfo[arminfo$trans %in% c(11) & arminfo$status==1,]),                # folks 4th fatal MACE
#nrow(arminfo[arminfo$trans %in% c(12) & arminfo$status==1,]),                # folks with censored time to 4th mace due to observed nonCV death
#nrow(arminfo[arminfo$trans %in% c(13) & arminfo$status==1,]),                # folks 5th fatal MACE
#nrow(arminfo[arminfo$trans %in% c(14) & arminfo$status==1,]))                # folks with censored time to 4th mace due to observed nonCV death
#rownames(trans.n) <- c(
#"Folks with observable 1st MACE",
#"Folks with observable fatal 1st MACE",
#"Folks who then die from non-CV death",
#"Folks with observable 2nd non-fatal MACE",
#"Folks with observable 2nd fatal MACE",
#"Folks with censored time to 2nd mace due to observed nonCV death",
#"Folks with 3rd non-fatal MACE",
#"Folks with 3rd fatal MACE",
#"Folks with censored time to 3rd mace due to observed nonCV death",
#"Folks with 4th non-fatal MACE",
#"Folks with 4th fatal MACE",
#"Folks with censored time to 4th mace due to observed nonCV death",
#"Folks with 5th fatal MACE",
#"Folks with censored time to 4th mace due to observed nonCV death" )
#
#trans.mat <- (is.na(tmat)==F)*1 + (is.na(tmat)==T)*0
#temp      <- as.vector(t(trans.mat))
#count     <- 0
#for(i in 1:length(temp)){
#  if(temp[i] == 1){
#    count <- count + 1
#    temp[i] <- trans.n[count]}}
#trans.mat <- t(matrix(temp,nrow=7,ncol=7,byrow=T))
#
#pp<-plotmat(trans.mat,pos=cbind(c(1,3,5,7,9,15,15)*.05,c(10,8,6,4,2,10,6)*.085),
#  name=c("Rand", "1st NFMACE", "2nd NFMACE", "3rd NFMACE", "4th NFMACE", "CVD", "Non-CVD"),lwd=1,box.lwd=2,
#  cex.txt=.75,box.size=0.04, shadow.size=0,dtext=1.1, box.cex = .8,
#  arr.lcol=rbind((trans.mat[1,]>0)*0,(trans.mat[2,]>0)*1,(trans.mat[3,]>0)*1,(trans.mat[4,]>0)*1,(trans.mat[5,]>0)*1,(trans.mat[6,]>0)*3,(trans.mat[7,]>0)*2),
#  arr.col=rbind((trans.mat[1,]>0)*0,(trans.mat[2,]>0)*1,(trans.mat[3,]>0)*1,(trans.mat[4,]>0)*1,(trans.mat[5,]>0)*1,(trans.mat[6,]>0)*3,(trans.mat[7,]>0)*2))
#}
