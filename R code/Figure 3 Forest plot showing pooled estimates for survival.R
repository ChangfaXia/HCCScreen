setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(dplyr)
library(haven)

rm(list=ls())
amonth <- 30.437
ayear <- 365.2425



## ggplot theme ####
library(Cairo)
library(ggplot2)
library(ggsci)
library(survival)
library(ggpubr)

cols1 <- pal_npg()(7)
font1 <- "Arial"

theme3 <-
  theme(panel.background=element_blank(),
        plot.background=element_blank(),
        plot.title=element_blank(),
        panel.grid.major=element_blank(),
        # panel.grid.major.x=element_blank(),
        # panel.grid.major.y=element_line(color="#d4d4d5",size=0.25),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line.x=element_line(color="black", size=0.3,colour="black"),
        axis.line.y=element_blank(),
        axis.ticks.x=element_line(size=0.3,colour="black"),
        axis.ticks.y=element_blank(),
        axis.ticks.length=unit(0.2, "cm"),
        axis.text.x=element_text(colour="black",size=8,family=font1),
        axis.text.y=element_blank(),
        axis.title.x=element_text(margin=unit(c(0.2,0,0,0),"cm"),hjust=0.6,size=9,family=font1),
        axis.title.y=element_blank(),
        legend.position=c(0.34,0.9),
        legend.justification="right",
        legend.box="horizontal",
        legend.text=element_text(colour="black",size=8,family=font1,
                                 margin=margin(l=-0.1, unit="cm")),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        plot.margin=unit(c(0.8,0,0,0), "cm"),
        strip.text=element_text(face="bold",size=9,family=font1),
        strip.background=element_rect(fill=NA,color=NA,size=1))



## Mate forest plot ####
library(readxl)

dat1 <-
  read_xlsx("survival_HR_Fig3A.xlsx",sheet=1, col_names=T, guess_max=1e5) %>%
  as.data.frame()
dat2 <-
  read_xlsx("survival_HR_leadtime_Fig3B.xlsx",sheet=1, col_names=T, guess_max=1e5) %>%
  as.data.frame()
dat3 <-
  read_xlsx("survival_HR_leadtime_length_Fig3C.xlsx",sheet=1, col_names=T, guess_max=1e5) %>%
  as.data.frame()


## panel A
dat1$total <- ifelse(is.na(dat1$`weight (%)`),2,1)
dat1$weight <- sprintf("%.1f", dat1$`weight (%)`)
dat1$weight[dat1$weight=="NA"] <- NA
dat1$`weight (%)`[dat1$total==2] <- mean(dat1$`weight (%)`, na.rm=T)
dat1$HR1 <-
  paste0(sprintf("%.2f",dat1$HR)," (",
         sprintf("%.2f",dat1$cil),"-",
         sprintf("%.2f",dat1$ciu), ")")
dat1 <- rbind(dat1[c(1,1),], dat1)
dat1[c(1,2),] <- NA
dat1[1,"study"] <- "Study"
dat1[1,"HR1"] <- "HR"
dat1[1,"weight"] <- "Weight (%)"
dat1[nrow(dat1),"study"] <- NA
dat1$no <- nrow(dat1):1
dat1$total <- factor(dat1$total)

p1 <-
  ggplot(dat1, aes(no, HR, fill=total, shape=total))+
  geom_hline(aes(yintercept=1),linetype="22", size=0.3)+
  geom_errorbar(aes(ymin=cil, ymax=ciu), size=0.3, width=0)+
  geom_point(aes(size=`weight (%)`),stroke=0.1,color="black") +
  geom_text(aes(y=0.1,label=study), size=3, hjust=0) +
  annotate("text", label=paste0("bold(Overall)~(Heterogeneity~italic(I)^2==",deparse(sprintf("%.1f",100*max(dat1$异质性I方,na.rm=T))), "*\'%\')"), #空格在parse=T时要用~表示
           x=1,y=0.1, parse=T, size=3,color="black",family=font1,hjust=0)+
  geom_text(aes(y=2.7,label=HR1), size=3) +
  geom_text(aes(y=4.5,label=weight), size=3) +
  annotate("text",label="Favour screening",
           x=31,y=0.96,size=3.0,hjust=1,color="black",family=font1,fontface="bold")+
  annotate("text",label="Does not favour\nscreening",
           x=30.5,y=1.04,size=3.0,hjust=0,color="black",family=font1,fontface="bold")+
  coord_flip(clip="off") +
  scale_x_continuous(expand=c(0,0),limits=c(0.25,31.25))+
  scale_y_log10(expand=c(0,0), limits=c(0.1,5),
                breaks=c(0.2,0.5,1,1.5,2), labels=c("0.2","0.5","1.0","1.5","2"))+
  scale_size_continuous(range=c(2,3.5), guide="none") +
  scale_fill_manual(values=cols1, guide="none") +
  scale_shape_manual(values=c(22,23), guide="none") +
  annotation_logticks(side="b",size=0.3,outside=T,
                      short=unit(0.1, "cm"),
                      mid=unit(0.15, "cm"),
                      long=unit(0.2, "cm"),)+
  ylab("Hazard Ratio (95% CI)") +
  xlab(NULL) +
  theme3


## panel B
dat2$total <- ifelse(is.na(dat2$`weight (%)`),2,1)
dat2$weight <- sprintf("%.1f", dat2$`weight (%)`)
dat2$weight[dat2$weight=="NA"] <- NA
dat2$`weight (%)`[dat2$total==2] <- mean(dat2$`weight (%)`, na.rm=T)
dat2$HR1 <-
  paste0(sprintf("%.2f",dat2$HR)," (",
         sprintf("%.2f",dat2$cil),"-",
         sprintf("%.2f",dat2$ciu), ")")
dat2 <- rbind(dat2[c(1,1),], dat2)
dat2[c(1,2),] <- NA
dat2[1,"study"] <- "Study"
dat2[1,"HR1"] <- "HR"
dat2[1,"weight"] <- "Weight (%)"
dat2[nrow(dat2),"study"] <- NA
dat2$no <- nrow(dat2):1
dat2$total <- factor(dat2$total)

p2 <-
  ggplot(dat2, aes(no, HR, fill=total, shape=total))+
  geom_hline(aes(yintercept=1),linetype="22", size=0.3)+
  geom_errorbar(aes(ymin=cil, ymax=ciu), size=0.3, width=0)+
  geom_point(aes(size=`weight (%)`),stroke=0.1,color="black") +
  geom_text(aes(y=0.1,label=study), size=3, hjust=0) +
  annotate("text", label=paste0("bold(Overall)~(Heterogeneity~italic(I)^2==",deparse(sprintf("%.1f",100*max(dat2$异质性I方,na.rm=T))), "*\'%\')"), #空格在parse=T时要用~表示
           x=1,y=0.1, parse=T, size=3,color="black",family=font1,hjust=0)+
  geom_text(aes(y=2.7,label=HR1), size=3) +
  geom_text(aes(y=4.5,label=weight), size=3) +
  annotate("text",label="Favour screening",
           x=17,y=0.96,size=3.0,hjust=1,color="black",family=font1,fontface="bold")+
  annotate("text",label="Does not favour\nscreening",
           x=16.5,y=1.04,size=3.0,hjust=0,color="black",family=font1,fontface="bold")+
  coord_flip(clip="off") +
  scale_x_continuous(expand=c(0,0),limits=c(0.25,17.25))+
  scale_y_log10(expand=c(0,0), limits=c(0.1,5),
                breaks=c(0.2,0.5,1,1.5,2), labels=c("0.2","0.5","1.0","1.5","2"))+
  scale_size_continuous(range=c(2,3.5), guide="none") +
  scale_fill_manual(values=cols1, guide="none") +
  scale_shape_manual(values=c(22,23), guide="none") +
  annotation_logticks(side="b",size=0.3,outside=T,
                      short=unit(0.1, "cm"),
                      mid=unit(0.15, "cm"),
                      long=unit(0.2, "cm"),)+
  ylab("Hazard Ratio (95% CI)") +
  xlab(NULL) +
  theme3


## panel C
dat3$total <- ifelse(is.na(dat3$`weight (%)`),2,1)
dat3$weight <- sprintf("%.1f", dat3$`weight (%)`)
dat3$weight[dat3$weight=="NA"] <- NA
dat3$`weight (%)`[dat3$total==2] <- mean(dat3$`weight (%)`, na.rm=T)
dat3$HR1 <-
  paste0(sprintf("%.2f",dat3$HR)," (",
         sprintf("%.2f",dat3$cil),"-",
         sprintf("%.2f",dat3$ciu), ")")
dat3 <- rbind(dat3[c(1,1),], dat3)
dat3[c(1,2),] <- NA
dat3[1,"study"] <- "Study"
dat3[1,"HR1"] <- "HR"
dat3[1,"weight"] <- "Weight (%)"
dat3[nrow(dat3),"study"] <- NA
dat3$no <- nrow(dat3):1
dat3$total <- factor(dat3$total)

p3 <-
  ggplot(dat3, aes(no, HR, fill=total, shape=total))+
  geom_hline(aes(yintercept=1),linetype="22", size=0.3)+
  geom_errorbar(aes(ymin=cil, ymax=ciu), size=0.3, width=0)+
  geom_point(aes(size=`weight (%)`),stroke=0.1,color="black") +
  geom_text(aes(y=0.1,label=study), size=3, hjust=0) +
  annotate("text", label=paste0("bold(Overall)~(Heterogeneity~italic(I)^2==",deparse(sprintf("%.1f",100*max(dat3$异质性I方,na.rm=T))), "*\'%\')"), #空格在parse=T时要用~表示
           x=1,y=0.1, parse=T, size=3,color="black",family=font1,hjust=0)+
  geom_text(aes(y=2.7,label=HR1), size=3) +
  geom_text(aes(y=4.5,label=weight), size=3) +
  annotate("text",label="Favour screening",
           x=5,y=0.96,size=3.0,hjust=1,color="black",family=font1,fontface="bold")+
  annotate("text",label="Does not favour\nscreening",
           x=4.5,y=1.04,size=3.0,hjust=0,color="black",family=font1,fontface="bold")+
  coord_flip(clip="off") +
  scale_x_continuous(expand=c(0,0),limits=c(0.25,5.25))+
  scale_y_log10(expand=c(0,0), limits=c(0.1,5),
                breaks=c(0.2,0.5,1,1.5,2), labels=c("0.2","0.5","1.0","1.5","2"))+
  scale_size_continuous(range=c(3,3.5), guide="none") +
  scale_fill_manual(values=cols1, guide="none") +
  scale_shape_manual(values=c(22,23), guide="none") +
  annotation_logticks(side="b",size=0.3,outside=T,
                      short=unit(0.1, "cm"),
                      mid=unit(0.15, "cm"),
                      long=unit(0.2, "cm"),)+
  ylab("Hazard Ratio (95% CI)") +
  xlab(NULL) +
  theme3


## out-print ggplot
p <-
  ggarrange(p1,p2,#p3,
            labels=c("(A)",
                     "(B)"),
            vjust=2.2,hjust=-0.4,
            font.label=list(size=8.5, family=font1), #,face="plain"
            heights=c(0.62,0.38), align="v", ncol=1, nrow=2)+
  theme(plot.margin=margin(2.14,1.91,2.54,1.61, "cm"))

Cairo(file="Figure meta plot",type="pdf",width=21,height=29.7,units="cm")
print(p)
dev.off()
