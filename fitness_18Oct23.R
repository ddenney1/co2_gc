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
library(glmmTMB)
library(ggeffects)
library(dplyr)


setwd("C:/Users/derek/OneDrive/Documents/Boechera/Growth_Chambers")
#setwd("C:/Users/dd66718/Downloads/Growth_Chambers")
#setwd("~/Documents/personnel/Denney")


fit <- read.csv("From Jill/fitness.csv", header = T, fileEncoding = "UTF-8-BOM")
sapply(fit,class)

fit$Treatment<-as.factor(fit$Treatment)
fit$Temperature<-as.factor(fit$Temperature)
fit$CO2<-as.factor(fit$CO2)
fit$Genotype <-as.factor(fit$Genotype)
fit$Population <-as.factor(fit$Population)
fit$Tray <-as.factor(fit$Tray)
fit$Block<-paste(fit$Tray, fit$Treatment, sep="_")
fit$Block<-as.factor(fit$Block)

fit $S_elev <- (fit $Elevation - mean(fit $Elevation, na.rm = TRUE)) / sd(fit $Elevation,na.rm = TRUE)
fit$day_transplant<-scale(fit$OD_transplant,center=TRUE, scale=TRUE)
fit$Elev <- fit$Elevation/1000
sapply(fit,class)


# Cohort information by filter (excel has exact cohorts)
fit$Transplant_Date <- as.Date(fit$Transplant_Date, format = "%m/%d/%Y")
fit$Cohort <- ifelse(fit$OD_transplant < 247, 1, 2)
fit$Cohort<-as.factor(fit$Cohort)

fit$Selev2<-I(fit$S_elev^2)



plot(fit$total_length_siliques ~fit$total_silique_number)

##Read in climatic data
climate <- read.csv("From Jill/source_populations_climates.csv", header = T)
climate$Population<-as.factor(climate$Population)
##Mean snowmelt timinigs and GDD
current_snow <-aggregate(climate $SnowMeltDOY, list(climate $Population), FUN=mean) 
current_snow <- current_snow %>% 
  rename("Population" = "Group.1",   "SnowMeltDOY"=  "x")

current_GDD <-aggregate(climate $GDD_nosnow, list(climate $Population), FUN=mean) 
current_GDD <- current_GDD %>% 
  rename("Population" = "Group.1",   "GDD_nosnow"=  "x")

currentclimi <-merge(current_snow, current_GDD,by=c("Population"))



##merge the datasets
fitness<-merge(currentclimi, fit,by=c("Population"))
#write.csv(fitness, "fitness_climate.csv")
fitness $snow_melt<-scale(fitness $SnowMeltDOY,center=TRUE, scale=TRUE)
fitness $sGDD <- (fitness $GDD_nosnow - mean(fitness $GDD_nosnow, na.rm = TRUE)) / sd(fitness $GDD_nosnow, na.rm = TRUE)

hist(fitness $snow_melt)
hist(fitness $sGDD)
plot(fitness$SnowMeltDOY~fitness$Elevation)
plot(fitness$GDD_nosnow ~fitness$Elevation)

fitness $sGDD2<-I(fitness $sGDD ^2)

mod1<-lm(Elevation~S_elev, data= fitness)
summary(mod1)

mod2<-lm(GDD_nosnow ~sGDD, data= fitness)
summary(mod2)


##############################################
######Survival: Logistic regression ##########
##############################################

survival_model <- glmer(Season_survival ~ Cohort+Temperature*CO2 *S_elev+I(S_elev ^2)*Temperature*CO2+ (1|Block) +(1| Genotype), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(survival_model, type="III")
plot(predictorEffects(survival_model, ~ S_elev), partial.residuals=TRUE)

visreg(survival_model,"S_elev",  by= "Temperature",cond=list(CO2 ="Low"), overlay= TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Probabilty of survival", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9"))) 
visreg(survival_model,"S_elev", by= "Temperature", cond=list(CO2 ="High"), overlay=TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Probabilty of survival", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9")))

# Obtain slopes for each treatment combination. These slopes need to be exponentiated to calcuate odds ratios
linear <- emtrends(survival_model, specs = c("Temperature", "CO2"), var = "S_elev")
# pairwise comparisons of these slopes
pairs(linear)

quadratic <- emtrends(survival_model, specs = c("Temperature", "CO2"), var = "I(S_elev ^2)")
# pairwise comparisons of these slopes
pairs(quadratic)

#For comparison 
library(broom.mixed)
##High temperature, High CO2
tidy(survival_model,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(survival_model,conf.int=TRUE,exponentiate=FALSE,effects="fixed")


surv <- glmer(Season_survival ~ Cohort+Temperature*CO2 *S_elev+Selev2*Temperature*CO2+ (1|Block) +(1| Genotype), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(surv, type="III")
emtrends(surv, specs = c("Temperature", "CO2"), var = "S_elev")
emtrends(surv, specs = c("Temperature", "CO2"), var = "Selev2")

##Test the random effects
survival_nogeno <- glmer(Season_survival ~ Cohort+Temperature*CO2 *S_elev+I(S_elev ^2)*Temperature*CO2+ (1|Block), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(survival_model, survival_nogeno)


survival_noblock <- glmer(Season_survival ~ Cohort+Temperature*CO2 *S_elev+I(S_elev ^2)*Temperature*CO2+ (1| Genotype), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(survival_model, survival_noblock)


mydf <- ggpredict(survival_model, c("S_elev[all]",  "Temperature[Low, High]","CO2[Low,High]"), type = "re")
plot(mydf, show_data=TRUE,facet=TRUE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of survival")
plot(mydf, show_data=TRUE,facet=TRUE,show_ci = FALSE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of survival")






###############################
##With growing degree days
###############################


survival_model_GDD<-glmer(Season_survival ~Cohort+ sGDD*Temperature*CO2+ I(sGDD^2)*Temperature*CO2+(1|Genotype)+(1|Block), data= fitness,family=binomial(link="logit"))
Anova(survival_model_GDD, type="III")
visreg(survival_model_GDD,"sGDD",  by= "Temperature",cond=list(CO2 ="Low"), overlay= TRUE,  scale = "response", xlab="Standardized Growing Degree Days", ylab="Probabilty of survival", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9"))) 
visreg(survival_model_GDD,"sGDD", by= "Temperature", cond=list(CO2 ="High"), overlay=TRUE,  scale = "response", xlab="Standardized Growing Degree Days", ylab="Probabilty of survival", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9")))


survival_model_GDDcoef<-glmer(Season_survival ~Cohort+ sGDD*Temperature*CO2+ sGDD2*Temperature*CO2+(1|Genotype)+(1|Block), data= fitness,family=binomial(link="logit"))
emtrends(survival_model_GDDcoef, specs = c("Temperature", "CO2"), var = "sGDD")
emtrends(survival_model_GDDcoef, specs = c("Temperature", "CO2"), var = "sGDD2")


##Test the random effects
survival_nogenoGDD <- glmer(Season_survival ~ Cohort+Temperature*CO2 * sGDD +I(sGDD ^2)*Temperature*CO2+ (1|Block), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(survival_model_GDD, survival_nogenoGDD)


survival_noblockGDD <- glmer(Season_survival ~ Cohort+Temperature*CO2 * sGDD +I(sGDD ^2)*Temperature*CO2+ (1| Genotype), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(survival_model_GDD, survival_noblockGDD)


survGDD <- ggpredict(survival_model_GDD, c("sGDD[all]",  "Temperature[Low, High]","CO2[Low,High]"), type = "re")
plot(survGDD, show_data=TRUE,facet=TRUE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of survival")
plot(survGDD, show_data=TRUE,facet=TRUE,show_ci = FALSE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of survival")


##############################################
######Survival: Days until mortality ##########
##############################################
Mod_days_mort<-coxme(Surv(Date_first_death_observed, Status )~  Cohort +S_elev*Temperature*CO2+I(S_elev^2)*Temperature*CO2+(1|Genotype)+(1|Block), data = fitness, na.action=na.exclude)
#Model coefficients
sum_mort<-summary(Mod_days_mort) 
#Anova table. This model takes a moment to run
anova_mort<-Anova(Mod_days_mort)  
anova_mort


##############################################
######Probability of reproduction: Logistic regression ##########
##############################################

repro_model <- glmer(Lifetime_Fruited ~ Cohort+Temperature*CO2 *S_elev+ (1|Block) +(1| Genotype), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(repro_model, type="III")
plot(predictorEffects(repro_model, ~ S_elev), partial.residuals=TRUE)

visreg(repro_model,"S_elev",  by= "Temperature",cond=list(CO2 ="Low"), overlay= TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Probabilty of repro", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9"))) 
visreg(repro_model,"S_elev", by= "Temperature", cond=list(CO2 ="High"), overlay=TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Probabilty of repro", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9")))

# Obtain slopes for each treatment combination. These slopes need to be exponentiated to calcuate odds ratios
linear <- emtrends(repro_model, specs = c("Temperature", "CO2"), var = "S_elev")
# pairwise comparisons of these slopes
pairs(linear)

#For comparison 
library(broom.mixed)
##High temperature, High CO2
tidy(repro_model,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(repro_model,conf.int=TRUE,exponentiate=FALSE,effects="fixed")


##Test the random effects
repro_nogeno <- glmer(Lifetime_Fruited ~ Cohort+Temperature*CO2 *S_elev+ (1|Block), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(repro_model, repro_nogeno)


repro_noblock <- glmer(Lifetime_Fruited ~ Cohort+Temperature*CO2 *S_elev+ (1| Genotype), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(repro_model, repro_noblock)


mydf <- ggpredict(repro_model, c("S_elev[all]",  "Temperature[Low, High]","CO2[Low,High]"), type = "re")
plot(mydf, show_data=TRUE,facet=TRUE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of repro")
plot(mydf, show_data=TRUE,facet=TRUE,show_ci = FALSE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of repro")



###############################
##With growing degree days
###############################


repro_model_GDD<-glmer(Lifetime_Fruited ~Cohort+ sGDD*Temperature*CO2+(1|Genotype)+(1|Block), data= fitness,family=binomial(link="logit"))
Anova(repro_model_GDD, type="III")
visreg(repro_model_GDD,"sGDD",  by= "Temperature",cond=list(CO2 ="Low"), overlay= TRUE,  scale = "response", xlab="Standardized Growing Degree Days", ylab="Probabilty of repro", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9"))) 
visreg(repro_model_GDD,"sGDD", by= "Temperature", cond=list(CO2 ="High"), overlay=TRUE,  scale = "response", xlab="Standardized Growing Degree Days", ylab="Probabilty of repro", ylim=c(0,1), partial=FALSE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9")))

emtrends(repro_model_GDD, specs = c("Temperature", "CO2"), var = "sGDD")


##Test the random effects
repro_nogenoGDD <- glmer(Lifetime_Fruited ~ Cohort+Temperature*CO2 * sGDD + (1|Block), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(repro_model_GDD, repro_nogenoGDD)


repro_noblockGDD <- glmer(Lifetime_Fruited ~ Cohort+Temperature*CO2 * sGDD + (1| Genotype), data = fitness, family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(repro_model_GDD, repro_noblockGDD)


reproGDD <- ggpredict(repro_model_GDD, c("sGDD[all]",  "Temperature[Low, High]","CO2[Low,High]"), type = "re")
plot(reproGDD, show_data=TRUE,facet=TRUE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of repro")
plot(reproGDD, show_data=TRUE,facet=TRUE,show_ci = FALSE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Probablity of repro")


##############################################
######Fecundity: Silique number ##########
##############################################

fruited<-subset(fitness, Lifetime_Fruited=="1")
plot(fruited$Mature_silique_number~fruited$Mature_length_siliques)


Number_siliques_elev <- glmer(Mature_silique_number ~ Cohort+Temperature*CO2 *S_elev+ (1|Block) +(1| Genotype), data = fitness, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(Number_siliques_elev, type="III")
visreg(Number_siliques_elev,"S_elev",  by= "Temperature",cond=list(CO2 ="Low"), overlay= TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Fecundity (number of fruits)", partial=TRUE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9"))) 
visreg(Number_siliques_elev,"S_elev", by= "Temperature", cond=list(CO2 ="High"), overlay=TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Fecundity (number of fruits)",  partial=TRUE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9")))

fecelev <- ggpredict(Number_siliques_elev, c("S_elev",  "Temperature[Low, High]","CO2[Low,High]"), type = "re")
plot(fecelev, show_data=TRUE,facet=TRUE, show_ci = FALSE,colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Fecundity (number of fruits)")
# plot(fecelev, show_data=TRUE,facet=TRUE,show_ci = FALSE, 
#      line=list(lty=1:2,col=c("#D55E00","#56B4E9")), 
#      points=list(col=c("#D55E00","#56B4E9")) +
#        labs(x="Source elevation (standardized)", y="Fecundity (number of fruits)"))
emtrends(Number_siliques_elev, specs = c("Temperature", "CO2"), var = "S_elev")


##Test the random effects
#Doesn't converge
Number_siliques_elev_nogeno <- glmer(Mature_silique_number ~ Cohort+Temperature*CO2 *S_elev+ (1|Block), data = fitness, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(Number_siliques_elev, Number_siliques_elev_nogeno)

Number_siliques_elev_noblock <- glmer(Mature_silique_number ~ Cohort+Temperature*CO2 *S_elev+(1| Genotype), data = fitness, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(Number_siliques_elev, Number_siliques_elev_noblock)
Number_siliques_elev_nogeno <- glm(Mature_silique_number ~ Cohort+Temperature*CO2 *S_elev, data = fitness, family=Gamma(link="log"))
anova(Number_siliques_elev_noblock, Number_siliques_elev_nogeno)




##With growing degree days
Number_siliques_GDD <- glmer(Mature_silique_number ~ Cohort+Temperature*CO2 *sGDD+ (1|Block) +(1| Genotype), data = fitness, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(Number_siliques_GDD, type="III")
visreg(Number_siliques_GDD,"S_elev",  by= "Temperature",cond=list(CO2 ="Low"), overlay= TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Fecundity (number of fruits)", partial=TRUE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9"))) 
visreg(Number_siliques_GDD,"S_elev", by= "Temperature", cond=list(CO2 ="High"), overlay=TRUE,  scale = "response", xlab="Standardized source elevation", ylab="Fecundity (number of fruits)",  partial=TRUE, type="conditional",line=list(lty=1:2,col=c("#D55E00","#56B4E9")), points=list(col=c("#D55E00","#56B4E9")))


fecGDD <- ggpredict(Number_siliques_GDD, c("sGDD",  "Temperature[Low, High]","CO2[Low,High]"), type = "re")
plot(fecGDD, show_data=TRUE,facet=TRUE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Fecundity (number of fruits)")
plot(fecGDD, show_data=TRUE,facet=TRUE,show_ci = FALSE, colors = c("#D55E00","#56B4E9")) + labs(x="Source elevation (standardized)", y="Fecundity (number of fruits)")

emtrends(Number_siliques_GDD, specs = c("Temperature", "CO2"), var = "sGDD")


##Test the random effects
#Doesn't converge
Number_siliques_GDDnogeno <- glmer(Mature_silique_number ~ Cohort+Temperature*CO2 * sGDD + (1|Block), data = fitness, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(Number_siliques_GDD, Number_siliques_GDDnogeno)

Number_siliques_GDDnoblock <- glmer(Mature_silique_number ~ Cohort+Temperature*CO2 * sGDD +(1| Genotype), data = fitness, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(Number_siliques_GDD, Number_siliques_GDDnoblock)
Number_siliques_GDD_nogeno <- glm(Mature_silique_number ~ Cohort+Temperature*CO2 * sGDD, data = fitness, family=Gamma(link="log"))
anova(Number_siliques_GDDnoblock, Number_siliques_GDD_nogeno)




