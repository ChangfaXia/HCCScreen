setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(dplyr)
library(haven)

rm(list=ls())
amonth <- 30.437
ayear <- 365.2425


## survival models ####
library(survival)

## screen-detected HCCs
dat1 <- inc.mor[,-c("name1","round")]
names(dat1) <- c("id","sex","age","inc.date","inc.ICD","status",
                 "mor.date","mor.ICD", "time", "group")
dat1$time.adj <- dat1$time * ayear


## lead-time bias
lead.para <- fread("sojourn time.csv")
dat1 <- left_join(dat1, lead.para, by=c("sex","age"))

dat1$lead.time.death <- # death
  (1 - exp(-dat1$trans*dat1$time.adj) - dat1$trans*dat1$time.adj*exp(-dat1$trans*dat1$time.adj)) /
  dat1$trans * (1 - exp(-dat1$trans*dat1$time.adj))
dat1$lead.time.live <-# alive
  (1 - exp(-dat1$trans*dat1$time.adj)) / dat1$trans
dat1$lead.time <- ifelse(dat1$status==1, dat1$lead.time.live, dat1$lead.time.death)
dat1$time.adj <- (dat1$time.adj - dat1$lead.time) / ayear

dat1 <- as.data.frame(dat1)
dat1 <- dat1[,-c(which(names(dat1)=="trans"):which(names(dat1)=="lead.time"))]


## data curation
dat1.all <- dat1
dat1.all$group <- "all"

dat2 <- interval[,-c("name1","round")]
dat2$group <- "interval"
names(dat2) <- c("id","sex","age","inc.date","inc.ICD","status",
                 "mor.date","mor.ICD", "time", "group")
dat2$time.adj <- dat2$time

dat3 <- control
dat3$group <- "population"
names(dat3) <- c("id","sex","age","inc.date","inc.ICD","status",
                 "mor.date","mor.ICD", "time", "group")
dat3$time.adj <- dat3$time

dat <- rbind(dat1.all, dat1, dat2, dat3)

# competing risk of death
dat$status1 <- ifelse(grepl("C22|c22",dat$mor.ICD)&dat$status==2, 2, 1)


## HRs by group
groups <- c("all","prevalent","incidence","interval")
i=1
for(i in 1:length(groups)){

  group1 <- groups[i]
  print(group1)
  dat.i <- subset(dat, group%in%c("population", group1))


  ## basecase survival model
  dat.i$group <- factor(dat.i$group, levels=c("population",group1))
  fit2 <- summary(coxph(Surv(time, status) ~ group, data=dat.i, ties="breslow"))
  aHR <- round(fit2$conf.int[,c(1,3,4)], 2)
  print(aHR)
  fit1 <- summary(survfit(Surv(time, status) ~ group, data=dat.i), times=1:3)
  n.surv <- data.frame(group=fit1$strata,
                       time=fit1$time,
                       surv=fit1$surv,
                       lower=fit1$lower,
                       upper=fit1$upper)
  n.surv <- n.surv[order(n.surv$group),]
  n.surv$rate <- paste0(sprintf("%.1f",n.surv$surv*100)," (",
                        sprintf("%.1f",n.surv$lower*100),"-",
                        sprintf("%.1f",n.surv$upper*100), ")")
  print(n.surv[,c("group","time","rate")])

  # adjusted for sex and age
  fit2 <- summary(coxph(Surv(time, status) ~ group + age + sex, data=dat.i, ties="breslow"))
  aHR <- round(fit2$conf.int[,c(1,3,4)], 2)
  print(aHR)

  # adjusted for sex and age (competing risk of death)
  fit2 <- summary(coxph(Surv(time, status1) ~ group + age + sex, data=dat.i, ties="breslow"))
  aHR <- round(fit2$conf.int[,c(1,3,4)], 2)
  print(aHR)


  ## survival model with lead-time correction
  dat.i$group <- factor(dat.i$group, levels=c("population",group1))

  fit2 <- summary(coxph(Surv(time.adj, status) ~ group, data=dat.i, ties="breslow"))
  aHR <- round(fit2$conf.int[,c(1,3,4)], 2)
  print(aHR)
  fit1 <- summary(survfit(Surv(time.adj, status) ~ group, data=dat.i), times=1:3)
  n.surv <- data.frame(group=fit1$strata,
                        time=fit1$time,
                        surv=fit1$surv,
                        lower=fit1$lower,
                        upper=fit1$upper)
  n.surv <- n.surv[order(n.surv$group),]
  n.surv$rate <- paste0(sprintf("%.1f",n.surv$surv*100)," (",
                        sprintf("%.1f",n.surv$lower*100),"-",
                        sprintf("%.1f",n.surv$upper*100), ")")
  print(n.surv[,c("group","time","rate")])

}

