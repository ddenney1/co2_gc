library(tidyverse)
library(ggplot2)
library(car)
library(lme4)
library(effects)
require(lmerTest)
require(MuMIn)
library(performance)
require(visreg)
library(emmeans)
library(coxme)
library(DHARMa)
library(performance)
library(gridExtra)
library(broom.mixed)

#setwd("~/Documents/personnel/Denney")
setwd("C:/Users/derek/OneDrive/Documents/Boechera/Growth_Chambers")

fitness <- read.csv("full_summary_fixed_3May23.csv", header = T, fileEncoding = "UTF-8-BOM")
sapply(fitness,class)

fitness$Treatment<-as.factor(fitness$Treatment)
fitness$Temperature<-as.factor(fitness$Temperature)
fitness$CO2<-as.factor(fitness$CO2)
fitness$Genotype <-as.factor(fitness$Genotype)
fitness$Population <-as.factor(fitness$Population)
fitness$Tray <-as.factor(fitness$Tray)
fitness$Block<-paste(fitness$Tray, fitness$Treatment, sep="_")
fitness$Block<-as.factor(fitness$Block)

fitness$S_elev<-scale(fitness$Elevation,center=TRUE, scale=TRUE)
fitness$day_transplant<-scale(fitness$OD_transplant,center=TRUE, scale=TRUE)

fitness$Elev <- fitness$Elevation/1000
sapply(fitness,class)

###Vernalization ended 4/10/22, ordinal day 100. All dates in the experiment startd with day of transplanting (OD_transplant), which was in 2011, so we need to include (365-OD_transplant) to account for the time alive in 2021
fitness$Days_to_Flower<- fitness$Date_PostVern_Flowering-100-(365-fitness$OD_transplant)
fitness$Flowering_Duration<-fitness$Date_silique-fitness$Date_PostVern_Flowering
fitness$Height_flowering <- rowSums(fitness[,c("Height1_flowering", "Height2_flowering","Height3_flowering", "Height4_flowering","Height5_flowering", "Height6_flowering","Height7_flowering", "Height8_flowering")], na.rm=TRUE)
fitness$Root_to_Shoot<-fitness$Root_DW/(fitness$Flower_DW+fitness$Bolt_Leaves_DW+fitness$Stem_DW+fitness$Rosette_DW)
fitness$instant_WUE<- fitness$Licor_A/fitness$Licor_E
fitness$intrinsic_WUE<- fitness$Licor_A/fitness$Licor_gsw
fitness$Root_to_Shoot[is.infinite(fitness$Root_to_Shoot)] <- NA

# Cohort information by filter (excel has exact cohorts)
fitness$Transplant_Date <- as.Date(fitness$Transplant_Date, format = "%m/%d/%Y")
fitness$Cohort <- ifelse(fitness$Transplant_Date < as.Date("2021-09-04"), 1, 2)
fitness$Cohort<-as.factor(fitness$Cohort)

hist(fitness$Days_to_Flower)
hist(fitness$Flowering_Duration)
hist(fitness$Height_flowering)
fitness["Height_flowering"][fitness["Height_flowering"] == 0] <- NA


##############################################
######Transpiration ##########
##############################################
fitness$Licor_E_mmol <- (fitness$Licor_E *1000)
# With cohort
hist(fitness$Licor_E)
transpiration <- lmer(Licor_E_mmol ~Cohort+ Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(transpiration,type="III")

emtrends(transpiration, specs = c("Temperature", "CO2"), var = "S_elev")
emtrends(transpiration, specs = "CO2", var = "S_elev")

plot(predictorEffects(transpiration, ~ S_elev), partial.residuals=TRUE)

y <- expression("Transpiration rate (mmol m"^-2*" s"^-1*")")
c <- c('High'="High CO2", 'Low' = "Low CO2")

fitness$Temperature <- factor(fitness$Temperature, levels = rev(levels(fitness$Temperature)))
fitness$CO2 <- factor(fitness$CO2, levels = rev(levels(fitness$CO2)))

gr <- ref_grid(transpiration, cov.keep= c('S_elev', 'CO2'))
emm <- emmeans(gr, spec= c('S_elev', 'CO2'), level= 0.95)
gg <- ggplot(data= fitness, aes(x= S_elev, y= Licor_E_mmol)) +
  geom_ribbon(data= data.frame(emm), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80') +
  geom_line(data= data.frame(emm), aes(y= emmean)) +
  geom_point() +
  theme_bw()+theme(text = element_text(size=20),
                   axis.line.x = element_line(colour = "black"), 
                   axis.line.y = element_line(colour = "black"), 
                   panel.border = element_blank(),
                   panel.grid.major =element_blank(), 
                   panel.grid.minor=element_blank(),legend.position = "bottom")+
  labs(x="Source elevation (standardized)", y=y) +
  facet_grid(. ~ CO2, labeller = as_labeller(c))

gg

# Rerun model without temperature to verify slopes of CO2 treatments
transpiration_CO2 <- lmer(Licor_E_mmol ~ CO2* S_elev + (1|Block) +(1| Genotype) , 
                          data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(transpiration_CO2,type="III")
plot(predictorEffects(transpiration_CO2, ~ S_elev), partial.residuals=TRUE)

## Extract confidence intervals
# CI for S_elev under low temperature
cc <- confint(transpiration,parm="beta_")  ## slow (~ 11 seconds)
ctab <- cbind(est=fixef(transpiration),cc)
ctab<-as.data.frame(ctab)
names(ctab)[2] <- "lower"
names(ctab)[3] <- "upper"


# CI for S_elev under higher temperature
fitness $co2_low<-factor(fitness $CO2, levels = c("Low","High"))
transpiration_low <- lmer(Licor_E_mmol ~ Temperature* co2_low* S_elev + (1|Block) +(1| Genotype) , 
                          data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
cc_low <- confint(transpiration_low,parm="beta_")  ## slow (~ 11 seconds)
ctab_low <- cbind(est=fixef(transpiration_low), cc_low)
ctab_low <-as.data.frame(ctab_low)
names(ctab_low)[2] <- "lower"
names(ctab_low)[3] <- "upper"

ctab
ctab_low


tidy(transpiration, conf.int=T, effects="fixed")

##############################################
######Photosynthesis ##########
##############################################

# With cohort in model
photo <- lmer(Licor_A ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(photo,type="III")
x <- expression("Levels of CO"[2])
y <- expression("Assimilation rate (??mol m"^-2*" s"^-1*")")
emmip(photo,~ CO2, type="response", CIs=TRUE)+theme_bw()+
  theme(text = element_text(size=20),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), 
        panel.border = element_blank(),
        panel.grid.major =element_blank(),
        panel.grid.minor=element_blank())+geom_point(size=3)+ 
  ylab(y)+ xlab(x)+
  scale_colour_grey()+ theme(legend.position = c(0.2, 0.5))


library(vioplot)
fitness$CO2 <- factor(fitness$CO2, levels = rev(levels(fitness$CO2)))
vioplot(Licor_A ~ CO2, data= fitness, plotCentre = "point",  pchMed = 23,  
        horizontal=FALSE,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00"),
        col=c("#56B4E9","#D55E00"), 
        ylab=y, xlab=x)+ 
  stripchart(Licor_A ~ CO2, data= fitness,  
             method = "jitter",col = "black", 
             vertical = T, pch = 1, add = TRUE)

Photo= ggplot(data = fitness, aes(x= Elevation, y= Licor_A, group=Temperature, shape=Temperature, color= CO2))+    geom_point(size = 4, position=position_dodge(width=0.8)) +    xlab("Elevation of origin (m)")+ ylab("Photosynthesis") 
Photo +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", se=FALSE,  size=1.6)+ facet_grid(~ Temperature:CO2)

##############################################
######Stomatal Conductance ##########
##############################################
fitness$licor_gsw_mmol <- fitness$Licor_gsw*1000
hist(na.omit(fitness$licor_gsw_mmol))
# With cohort
hist(fitness$Licor_gsw)
condu <- lmer(licor_gsw_mmol ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(condu,type="III")
summary(condu)
plot(predictorEffects(condu, ~ S_elev), partial.residuals=TRUE)


cond= ggplot(data = fitness, aes(x= Elevation, y= Licor_gsw, group=Temperature, shape=Temperature, color= CO2))+    geom_point(size = 4, position=position_dodge(width=0.8)) +    xlab("Elevation of origin (m)")+ ylab("Stomatal conductance") 
cond +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", se=FALSE,  size=1.6)+ facet_grid(~ Temperature:CO2)

y <- expression("Stomotal conductance (mmol m"^-2*" s"^-1*")")
c <- c('High' = "High CO2", 'Low' = "Low CO2")
gr <- ref_grid(condu, cov.keep= c('S_elev', 'CO2'))
emm <- emmeans(gr, spec= c('S_elev', 'CO2'), level= 0.95)
gg <- ggplot(data= fitness, aes(x= S_elev, y= licor_gsw_mmol)) +
  geom_ribbon(data= data.frame(emm), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80') +
  geom_line(data= data.frame(emm), aes(y= emmean)) +
  geom_point() +
  theme_bw()+theme(text = element_text(size=20),
                   axis.line.x = element_line(colour = "black"), 
                   axis.line.y = element_line(colour = "black"), 
                   panel.border = element_blank(),
                   panel.grid.major =element_blank(), 
                   panel.grid.minor=element_blank(),legend.position = "bottom")+
  labs(x="Source elevation (standardized)", y=y) +
  facet_grid(. ~ CO2, labeller = as_labeller(c))

gg


# Rerun model without temperature to verify slopes of CO2 treatments
condu_CO2 <- lmer(licor_gsw_mmol ~CO2* S_elev + (1|Block) +(1| Genotype) , 
                  data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(condu_CO2,type="III")
plot(predictorEffects(condu_CO2, ~ S_elev), partial.residuals=TRUE)

## Extract confidence intervals
# CI for S_elev under low temperature
cc <- confint(condu,parm="beta_")  ## slow (~ 11 seconds)
ctab <- cbind(est=fixef(condu),cc)
ctab<-as.data.frame(ctab)
names(ctab)[2] <- "lower"
names(ctab)[3] <- "upper"


# CI for S_elev under higher temperature
fitness $co2_low<-factor(fitness $CO2, levels = c("Low","High"))
condu_low <- lmer(licor_gsw_mmol ~Temperature* co2_low* S_elev + (1|Block) +(1| Genotype) , 
                  data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
cc_low <- confint(condu_low,parm="beta_")  ## slow (~ 11 seconds)
ctab_low <- cbind(est=fixef(condu_low), cc_low)
ctab_low <-as.data.frame(ctab_low)
names(ctab_low)[2] <- "lower"
names(ctab_low)[3] <- "upper"
ctab
ctab_low

### Evaluating CO2 and T on their own to visualize effects
x <- expression("Levels of CO"[2])
# CO2
emmip(condu,~ CO2, type="response", CIs=TRUE)+
  theme_bw()+theme(text = element_text(size=20),
                   axis.line.x = element_line(colour = "black"),
                   axis.line.y = element_line(colour = "black"),
                   panel.border = element_blank(),panel.grid.major =element_blank(),
                   panel.grid.minor=element_blank())+geom_point(size=3)+  
  xlab(x) +
  ylab("Stomatal conductance (gsw)")+
  scale_colour_grey()+ theme(legend.position = c(0.2, 0.5))

# temperature
x <- "Temperature"
emmip(condu,~ Temperature, type="response", CIs=TRUE)+theme_bw()+
  theme(text = element_text(size=20),
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black"), 
        panel.border = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor=element_blank())+geom_point(size=3)+  
  xlab(x) + ylab("Stomatal conductance (gsw)")+ 
  scale_colour_grey()+ theme(legend.position = c(0.2, 0.5))



fitness$quartile_group <- cut(fitness$S_elev, 
                              breaks = quantile(fitness$S_elev, 
                                                probs = c(0, 0.25, 0.5, 0.75, 1)), 
                              labels = FALSE, include.lowest = T)

# Assign values 1 through 4 to quartile groups
fitness$quartile_group <- match(fitness$quartile_group, unique(fitness$quartile_group))


gs <- fitness %>% group_by(CO2, quartile_group) %>% summarise(mean_gsw = mean(Licor_gsw, na.rm= T))
summary(gs)

plot(gs$CO2, gs$mean_gsw, col = gs$quartile_group)

##############################################
######Intrinsic WUE ##########
##############################################
hist(fitness$intrinsic_WUE)
# With cohort
intrinsicWUE <- lmer(intrinsic_WUE ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype), data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(intrinsicWUE,type="III")
summary(intrinsicWUE)
plot(predictorEffects(intrinsicWUE, ~ S_elev), partial.residuals=TRUE)

x.1 <- expression("Levels of CO"[2])
# CO2
c.iwue <- emmip(intrinsicWUE,~ CO2, type="response", CIs=TRUE)+
  theme_bw()+theme(text = element_text(size=20),
                   axis.line.x = element_line(colour = "black"),
                   axis.line.y = element_line(colour = "black"),
                   panel.border = element_blank(),panel.grid.major =element_blank(),
                   panel.grid.minor=element_blank())+geom_point(size=3)+  
  xlab(x.1) +
  ylab("Intrinsic Water-Use Efficiency (A/gsw)")+
  scale_colour_grey()+ theme(legend.position = c(0.2, 0.5))

# temperature
x <- "Temperature"
t.iwue <- emmip(intrinsicWUE,~ Temperature , type="response", CIs=TRUE)+theme_bw()+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = "black"), 
        #axis.line.y = element_line(colour = "black"), 
        panel.border = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor=element_blank())+geom_point(size=3)+ 
  xlab(x)+
  scale_colour_grey()+ theme(legend.position = c(0.2, 0.5))


grid.arrange(c.iwue, t.iwue, ncol = 2)

WUE= ggplot(data = fitness, aes(x= Elevation, y= intrinsic_WUE, group=Temperature, shape=Temperature, color= CO2))+    geom_point(size = 4, position=position_dodge(width=0.8)) +    xlab("Elevation of origin (m)")+ ylab("Instrinsic WUE ") 
WUE +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", se=FALSE,  size=1.6)+ facet_grid(~ Temperature:CO2)


gr <- ref_grid(intrinsicWUE, cov.keep= c('S_elev', 'CO2'))
emm <- emmeans(gr, spec= c('S_elev', 'CO2'), level= 0.95)
gg <- ggplot(data= fitness, aes(x= S_elev, y= intrinsic_WUE)) +
  geom_ribbon(data= data.frame(emm), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80') +
  geom_line(data= data.frame(emm), aes(y= emmean)) +
  geom_point() +
  theme_bw()+theme(text = element_text(size=20),
                   axis.line.x = element_line(colour = "black"), 
                   axis.line.y = element_line(colour = "black"), 
                   panel.border = element_blank(),
                   panel.grid.major =element_blank(), 
                   panel.grid.minor=element_blank(),legend.position = "bottom")+
  labs(x="Source elevation (standardized)", y=y) +
  facet_grid(. ~ CO2, labeller = as_labeller(c))

gg
library(vioplot)
# Violin plot
fitness$CO2 <- factor(fitness$CO2, levels = rev(levels(fitness$CO2)))
vioplot(intrinsic_WUE ~  CO2  , data= fitness, plotCentre = "point",  pchMed = 23,  
        vertical=T,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00","#E1BE6A","#40B0A6"),
        col=c("#56B4E9","#D55E00", "#E1BE6A","#40B0A6"), 
        ylab="iWUE", xlab="CO2")
fitness$Temperature <- factor(fitness$Temperature, levels = rev(levels(fitness$Temperature)))
vioplot(intrinsic_WUE ~  Temperature  , data= fitness, plotCentre = "point",  pchMed = 23,  
        vertical=T,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00","#E1BE6A","#40B0A6"),
        col=c("#E1BE6A","#40B0A6","#E1BE6A","#40B0A6"), 
        ylab="iWUE", xlab="Temperature")
stripchart(intrinsic_WUE ~ Temperature, data= fitness,  
           method = "jitter",col = "gray", 
           vertical = T, pch = 1, add = TRUE, side = "left")


# Violin with all iwue treatments together
vioplot(intrinsic_WUE ~  CO2 +Temperature  , data= fitness, plotCentre = "point",  pchMed = 23,  
        vertical=T,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00","#E1BE6A","#40B0A6"),
        col=c("#56B4E9","#D55E00", "#E1BE6A","#40B0A6"), 
        ylab="iWUE", xlab="CO2 and Temperature")
stripchart(intrinsicWUE ~ CO2 + Temperature, data= fitness,  
           method = "jitter",col = "gray", 
           vertical = T, pch = 1, add = TRUE)

model <- lm(intrinsic_WUE ~ CO2, data = fitness)

# Create a reference grid for CO2 only
gr_co2 <- ref_grid(simple_model, cov.keep = 'CO2')
emm_co2 <- emmeans(gr_co2, "CO2", level = 0.95, adjust = "tukey")
emm_co2

model <- lm(intrinsic_WUE ~ Temperature, data = fitness)
gr_temp <- ref_grid(model, cov.keep = 'Temperature')
emm_temp <- emmeans(gr_temp, "Temperature", level = 0.95, adjust = "tukey")
emm_temp

emms <- emmeans(intrinsicWUE, ~ CO2 + Temperature)
pairwise_comparisons <- pairs(emms)
pairwise_comparisons

emmeans_grid <- emmeans(intrinsicWUE, ~ CO2 + Temperature)
emmeans_grid


# Rerun model without temperature to verify slopes of CO2 treatments
intrinsicWUE_CO2 <- lmer(intrinsic_WUE ~CO2* S_elev + (1|Block) +(1| Genotype) , 
                         data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(intrinsicWUE_CO2,type="III")
plot(predictorEffects(intrinsicWUE_CO2, ~ S_elev), partial.residuals=TRUE)

# CI for S_elev under low temperature
cc <- confint(intrinsicWUE,parm="beta_")  ## slow (~ 11 seconds)
ctab <- cbind(est=fixef(intrinsicWUE),cc)
ctab<-as.data.frame(ctab)
names(ctab)[2] <- "lower"
names(ctab)[3] <- "upper"


# CI for S_elev under higher temperature
fitness $co2_low<-factor(fitness $CO2, levels = c("Low","High"))
intrinsicWUE_low <- lmer(intrinsic_WUE ~Temperature* co2_low* S_elev + (1|Block) +(1| Genotype) , 
                         data = fitness, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
cc_low <- confint(intrinsicWUE_low,parm="beta_")  ## slow (~ 11 seconds)
ctab_low <- cbind(est=fixef(intrinsicWUE_low), cc_low)
ctab_low <-as.data.frame(ctab_low)
names(ctab_low)[2] <- "lower"
names(ctab_low)[3] <- "upper"
ctab
ctab_low

##########################
### Field Licor Data #####
###########################

sco <- read.delim("From Jill/Sco_ecophys.txt", header = T)
sco $S_elev <- (sco $elevation..m. - mean(sco $elevation..m., na.rm = TRUE)) / sd(sco $elevation..m.,na.rm = TRUE)
head(sco)

# Photosynthesis
plot(sco$photosynthesis)

sco_photo <- lmer(photosynthesis ~ elevation..m. + (1|Quad) +(1| genotype), data = sco, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(sco_photo, type="III")
summary(sco_photo)
plot(predictorEffects(sco_photo, ~ elevation..m.), partial.residuals=T)

gr <- ref_grid(sco_photo, cov.keep= "elevation..m.")
emm <- emmeans(gr, spec= "elevation..m.", level= 0.95)
gg <- ggplot(data= sco, aes(x= elevation..m., y= photosynthesis)) +
  #geom_ribbon(data= data.frame(emm), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80') +
  geom_line(data= data.frame(emm), aes(y= emmean)) +
  geom_point() +
  theme_bw()+theme(text = element_text(size=20),
                   axis.line.x = element_line(colour = "black"), 
                   axis.line.y = element_line(colour = "black"), 
                   panel.border = element_blank(),
                   panel.grid.major =element_blank(), 
                   panel.grid.minor=element_blank(),legend.position = "bottom")+
  labs(x="Source elevation (m)", y="Assimilation rate")

gg

cc <- confint(sco_photo,parm="beta_")  ## slow (~ 11 seconds)
ctab <- cbind(est=fixef(sco_photo),cc)
ctab<-as.data.frame(ctab)
names(ctab)[2] <- "lower"
names(ctab)[3] <- "upper"
ctab