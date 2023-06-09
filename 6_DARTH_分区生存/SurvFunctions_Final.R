##################################################################################
####    Survival analysis functions for decision modeling    ####
##################################################################################

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (1)	
# M.G. Myriam Hunink, MD, PhD (2,3)
# Hawre J. Jalal, MD, PhD (4) 
# Eline M. Krijkamp, MSc (2)	
# Petros Pechlivanoglou, PhD (5) 
# David Rios, MA (5)
# Alan Yang, MSc (5)

# In collaboration of: 		
# 1 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 2 Erasmus MC, Rotterdam, The Netherlands
# 3 Harvard T.H. Chan School of Public Health, Boston, USA
# 4 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 5 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

#####################################################################################
# Please cite our publications when using this code
# Jalal H, et al. An Overview of R in Health Decision Sciences. Med. Decis. Making. 2017; 37(3): 735-746. 
# Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial. Med. Decis. Making. 2018; (in press). 

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR Sick CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property are exclusively owned by THE HOSPITAL FOR Sick CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, distributed or adapted in any way without written permission.



################# Fit different distributions to the data #####################

fit.fun <- function(time = "pfs", event = "event.pfs", trt = NULL, data = data ,add=FALSE)  
{
  data$time  = data[, time]  # was there an event of PFS?
  data$event = data[, event]  # was there an event of PFS?
  if (is.null(trt) == F ) data = subset(data, trt = trt)
  
  # Progression free survival  
  KM.fit     <-     survfit(Surv(time, event) ~ 1, data = data)                   # fit Kaplan-Meier curve 
  fit.llogis <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "llogis" ) # fit model with loglogistic distribution
  fit.weib   <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "weibull") # fit model with Weibull distribution
  fit.lnorm  <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "lnorm"  ) # fit model with lognormal distribution
  fit.gamma  <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "gamma"  ) # fit model with gamma distribution 
  fit.exp    <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "exp"    ) # fit model with exponential distribution 
 
  # extarapolate all models beyond the KM curve
  if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0,n.t), conf.int= F)}
  if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0,n.t), conf.int= F)}
  lines(fit.llogis,   t = times, col = 2, ci = F)
  lines(fit.weib,     t = times, col = 3, ci = F)
  lines(fit.lnorm,    t = times, col = 4, ci = F)
  lines(fit.gamma,    t = times, col = 5, ci = F)
  lines(fit.exp,      t = times, col = 6, ci = F)
  legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gamma", "Exponential"), col = 1:6, lty = rep(1, 6), bty="n")

  # compare AIC values for the PFS
  AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
               Weibull     = AIC(fit.weib), 
               Lognormal   = AIC(fit.lnorm), 
               Gamma       = AIC(fit.gamma),
               Exponentail = AIC(fit.exp))
  AIC= round(AIC,3)
  
  # compare BIC values for the PFS
  BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
               Weibull     = BIC(fit.weib), 
               Lognormal   = BIC(fit.lnorm), 
               Gamma       = BIC(fit.gamma),
               Exponential = BIC(fit.exp))
 
   BIC <- round(BIC,3)
  
   res <- list(Loglogistic = fit.llogis,
               Weibull     = fit.weib,
               Lognormal   = fit.lnorm, 
               Gamma       = fit.gamma,
               Exponential = fit.exp, 
               AIC         = AIC,
               BIC         = BIC)
  res
}

#### Fit partitioned survival model #### 

partsurv <- function(fit.pfs, fit.os, title = "trt", time = times){
  # Input
    # fit.pfs: flexsurv obj fitting pfs
    # fit.os: flexsurv obj fitting os
    # title:
    # time = numeric vector of time to estimate probabilities
  # output:
    #  res a list w/ one entry of a data frame w/ probabilities associated w/ stable ,prog and dead.
  
  pfs.surv <- summary(fit.pfs, t = time, ci = F)[[1]]$est
  os.surv <- summary(fit.os, t = time, ci = F)[[1]]$est
  prog                 <- os.surv - pfs.surv          # estimate the probability of remaining in the progressed state
  prog[prog < 0]       <- 0                           # in cases where the probability is negative replace with zero
  stable               <- pfs.surv                    # probability of remaining stable
  dead                 <- 1 - os.surv                 # probability of being dead
  trace <- data.frame(stable = stable, prog = prog, dead = dead)
  res   <- list(trace = trace)
  res   <- list(trace = trace)
  return(res)
}
