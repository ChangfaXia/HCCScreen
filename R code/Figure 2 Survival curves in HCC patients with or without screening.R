setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(dplyr)
library(haven)

rm(list=ls())
amonth <- 30.437
ayear <- 365.2425



## ggplot Theme ####
library(Cairo)
library(ggplot2)
library(survminer)
library(survival)
library(ggsci)


cols1 <- pal_npg()(7)

font1 <- "Arial"

theme1 <-
  theme(panel.background=element_blank(),
        plot.background=element_blank(),
        plot.title=element_blank(),
        panel.grid.major=element_blank(),
        # panel.grid.major.x=element_blank(),
        # panel.grid.major.y=element_line(color="#d4d4d5",size=0.25),
        panel.grid.minor=element_blank(),
        panel.border= element_blank(),
        axis.line.x=element_line(color="black", size=0.3,colour="black"),
        axis.line.y=element_line(color="black", size=0.3,colour="black"),
        axis.ticks=element_line(size=0.3,colour="black"),
        axis.ticks.length=unit(0.2, "cm"),
        axis.text=element_text(colour="black",size=10,family=font1),
        axis.title.x=element_text(margin=unit(c(0.3,0,0,0),"cm"),size=11,family=font1),
        axis.title.y=element_text(margin=unit(c(0,-4.5,0,0),"cm"),size=11,family=font1),
        legend.position=c(0.62,0.15),
        legend.justification="right",
        legend.box="horizontal",
        legend.text=element_text(colour="black",size=10,family=font1,
                                 margin=margin(l=-0.1, unit="cm")),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm"),
        strip.text=element_text(face="bold",size=11,family=font1),
        strip.background=element_rect(fill=NA,color=NA,size=1))



## Survival plot ####

## data
dat1 <- inc.mor[,-c("name1","round")]
names(dat1) <- c("id","sex","age","inc.date","inc.ICD","status",
                 "mor.date","mor.ICD", "time", "group")
dat2 <- interval[,-c("name1","round")]
dat2$group <- "interval"
names(dat2) <- c("id","sex","age","inc.date","inc.ICD","status",
                 "mor.date","mor.ICD", "time", "group")
dat3 <- control
dat3$group <- "population"
names(dat3) <- c("id","sex","age","inc.date","inc.ICD","status",
                 "mor.date","mor.ICD", "time", "group")
dat <- rbind(dat1, dat2, dat3)
dat$time <- ifelse(dat$status==2,
                   as.Date(dat$mor.date)-as.Date(dat$inc.date),
                   as.Date("2021-10-01",origin="1970-01-01")-as.Date(dat$inc.date))
dat$time <- ifelse(dat$time<0, 0, dat$time/ayear)


## Survival models
dat$group <- factor(dat$group, levels=c("population","incidence","prevalent","interval"))

fit1 <- survfit(Surv(time, status) ~ group, data=dat)
conf <- paste0("Incident round vs. No-screened, corrected HR=0.74 (0.50-0.93)\n",
               "Prevalent round vs. No-screened, corrected HR=0.75 (0.60-0.86)\n",
               "Interval cancer vs. No-screened, crude HR=0.71 (0.59-0.85)")
# p.val <- fit2$sctest[3]
# p.val <- ifelse(p.val<0.001, "P<0.001", paste0("P=",sprintf("%.3f",p.val)))

# as tables
fit1.res <- data.frame(time=fit1$time,
                       surv=fit1$surv,
                       lower=fit1$lower,
                       upper=fit1$upper,
                       event=fit1$n.event,
                       risk=fit1$n.risk,
                       censor=fit1$n.censor,
                       group=rep(names(fit1$strata), fit1$strata))

fit1.res.min <-
  fit1.res %>%
  group_by(group) %>%
  summarise(across(c(time), ~min(.x, na.rm=T)), .groups="drop") %>%
  left_join(fit1.res, by=c("time","group")) %>%
  as.data.frame()
fit1.res.min$time <- 0
fit1.res.min$surv <- 1
fit1.res.min$lower <- 1
fit1.res.min$upper <- 1
fit1.res <- rbind(fit1.res.min, fit1.res)
fit1.res$group <- recode(fit1.res$group,
                         `group=incidence`="Incident rounds",
                         `group=interval`="Interval HCCs",
                         `group=population`="Non-screened",
                         `group=prevalent`="Prevalent round")
fit1.res <- fit1.res[order(fit1.res$group,decreasing=T),]

# censors
censor <- subset(fit1.res, censor!=0)[,1:3]

# Number at risk
num <- summary(fit1, times=seq(0,4,0.5))
num <- data.frame(group=num$strata,
                     time=num$time,
                     risk=num$n.risk)
num$group <- recode(num$group,
                    `group=incidence`="Incident rounds",
                    `group=interval`="Interval HCCs",
                    `group=population`="Non-screened",
                    `group=prevalent`="Prevalent round")


## plot
fit1.res <- fit1.res[order(fit1.res$time),]
fit1.res <- fit1.res[order(fit1.res$group),]
# K-M plot
aaa <- fit1.res
aaa[,c("surv","lower","upper")] <- aaa[c(1,1:(nrow(aaa)-1)), c("surv","lower","upper")]
aaa <- aaa[aaa$time!=0, ]
fit1.res <- rbind(fit1.res, aaa)
# reorder
fit1.res <- fit1.res[order(fit1.res$surv, decreasing=T),]
fit1.res <- fit1.res[order(fit1.res$time),]
fit1.res <- fit1.res[order(fit1.res$group),]

fit1.res$group <- factor(fit1.res$group,
                         levels=c("Incident rounds","Prevalent round","Interval HCCs","Non-screened"))

p1 <-
  ggplot(fit1.res, aes(x=time, y=100*surv, group=group, color=group)) +
  geom_ribbon(aes(ymin=100*upper,ymax=100*lower,fill=group),alpha=0.15,color=NA)+
  geom_line(size=0.5) +
  geom_point(data=censor, aes(x=time, y=100*surv, group=group, color=group),
             size=0.7,shape=3) +
  # annotate("text", label=paste0(conf),
  #          0.08, 0, label="Top-right", hjust=0, vjust=-0.25,size=3.4, family=font1)+
  xlab("Time Since Diagnosis (yr)") +
  ylab("Overall Survival (%)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,3.5), breaks=seq(0,3.5,0.5),
                     labels=c("0","0.5","1","1.5","2","2.5","3","3.5"))+
  scale_y_continuous(expand=c(0,0), limits=c(0,100.25))+
  scale_color_manual(values=cols1, labels=c("Screen-detected HCCs in the incident rounds",
                                            "Screen-detected HCCs in the prevalent round",
                                            "Interval HCCs",
                                            "Non-screened HCCs"))+
  scale_fill_manual(values=cols1, labels=c("Screen-detected HCCs in the incident rounds",
                                           "Screen-detected HCCs in the prevalent round",
                                           "Interval HCCs",
                                           "Non-screened HCCs"))+
  theme1

# Number at risk
num$group <- factor(num$group, levels=c("Incident rounds","Prevalent round","Interval HCCs","Non-screened"))
num$risk <- prettyNum(num$risk, big.mark=",", scientific=F)
p.num <-
  ggplot(num, aes(x=time, y=group, label=risk)) +
  geom_text(size=3.4, family=font1) +
  scale_x_continuous(expand=c(0,0), limits=c(0,3.5))+
  scale_y_discrete(name="",expand=c(0,0.5),limits=rev(levels(num$group)))+
  ggtitle("No. at Risk")+
  # theme1+
  theme_void()+
  theme(plot.margin=unit(c(0,0,0,0), "cm"),
        plot.title=element_text(colour="black",size=10,family=font1,
                                vjust=0, hjust=-0.2, margin=margin(0,0,5,0)),
        axis.text.y=element_text(colour="black",size=10,family=font1,
                                 hjust=1, margin=margin(0,18,0,0)))
p.num <- ggplot_gtable(ggplot_build(p.num))
p.num$layout$clip <- "off"

p <-
  ggarrange(p1, p.num,
            # labels=NA, vjust=0,hjust=0,
            font.label=list(size=9, family=font1),
            heights=c(0.8,0.2), align="v", ncol=1, nrow=2)+
  theme(plot.margin=margin(2.54,1.91,2.54,1.91, "cm"))

Cairo(file="Figure survival plot.pdf",type="pdf",width=21,height=18,units="cm")
print(p)
dev.off()
