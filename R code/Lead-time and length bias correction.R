setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(dplyr)
rm(list=ls())

ayear <- 365.2425


## Markov model by days ####
library(msm)

# model
statetable.msm(state, id , data=screen)
Q <- rbind(c(0.9999, 0.00005, 0.00005),
           c(0, 0.999, 0.001),
           c(0, 0, 0))
rownames(Q) <- colnames(Q) <- c("disease-free", "preclinical", "clinical")
E <- rbind(c(0, 0.1, 0),
           c(0.1, 0, 0),
           c(0, 0, 0))
rownames(E) <- colnames(E) <- c("disease-free", "preclinical", "clinical")
msm.model <- msm(state ~ time, subject=id,
                 covariates = ~ sex + age,
                 ematrix=E, obstrue=obsture,
                 qmatrix=Q, deathexact=3, data=screen)

# Extract transition intensity matrix and sojourn time by sex and age
strata <- expand.grid(sex=sort(unique(screen$sex)), age=sort(unique(screen$age)))
sojourn <- NULL
for(i in 1:nrow(strata)) {
  sex.i <- strata$sex[i]
  age.i <- strata$age[i]
  sojourn1 <- qmatrix.msm(msm.model, ci="delta", sojourn=T, covariates=list(sex1=sex.i, ragersk=age.i))
  sojourn1 <- data.frame(sex=sex.i,
                         age=age.i,
                         trans=sojourn1$estimates["preclinical","clinical"],
                         lower=sojourn1$L["preclinical","clinical"],
                         upper=sojourn1$U["preclinical","clinical"],
                         sojourn=sojourn1$sojourn["preclinical"],
                         L=sojourn1$sojournL["preclinical"],
                         U=sojourn1$sojournU["preclinical"])
  sojourn <- rbind(sojourn, sojourn1)
}
sojourn <- as.data.frame(sojourn)
write.csv(sojourn, "sojourn time.csv", row.names=F)



## Lead-time bias and length bias correction ####
library(survival)

# Screening data
names(dat) <- c("id","sex","age","inc.date","inc.ICD","status",
                 "mor.date","mor.ICD", "time", "group")
dat$time.adj <- dat$time * ayear #convert year to day

# Lead-time bias correction
dat <- left_join(dat, sojourn, by=c("sex","age"))
dat$lead.time.death <- # for death cases
  (1 - exp(-dat$trans*dat$time.adj) - dat$trans*dat$time.adj*exp(-dat$trans*dat$time.adj)) /
  dat$trans * (1 - exp(-dat$trans*dat$time.adj))
dat$lead.time.live <- # for survival cases
  (1 - exp(-dat$trans*dat$time.adj)) / dat$trans
dat$lead.time <- ifelse(dat$status==1, dat$lead.time.live, dat$lead.time.death)
dat$time.adj <- (dat$time.adj - dat$lead.time) / ayear #convert day to year
dat <- as.data.frame(dat)
fit.COX <- summary(coxph(Surv(time.adj, status) ~ group + sex + age, data=dat, ties="breslow"))
HR.adj.leadtime <- round(fit.COX$conf.int[,c(1,3,4)], 2) #HR
print(round(HR.adj.leadtime, 2))

# Length bias correction
theta <- 0.9 #relative risk (θ) of death from slow-growing tumors vs faster-growing tumors (<1)
q <- 0.8 #the proportion of patients with faster-growing tumors (q)
p1 <- 1-surv[surv$group=="group=population", "surv"] #observed probability of cancer death for symptomatic cancer
p2 <- 1-surv[surv$group!="group=population", "surv"] #lead-time bias corrected probability of cancer death for screen-detected cancer
p3 <- sum(n.num[c("prevalent","incident")])/sum(n.num[c("prevalent","incident","interval")]) #observed probability of screen detection
phi <- p2*(q*theta+1-q) / (theta*p) #true relative risk (φ) of death from the cancer for screen-detected vs symptomatic tumors
pp <- (p1*(theta*q+1-q)*(1-p3))/((theta*q+1-q)*(theta+q*(1-theta))-p3*theta)
HR.adj.length <- log(1-phi*pp)/log(1-pp) #HR
print(round(HR.adj.length, 2))
