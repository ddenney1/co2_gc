library(tidyverse)
library(ggplot2)
library(car)
library(emmeans)
library(DHARMa)
library(glmmTMB)
library(ggeffects)


setwd("C:/Users/derek/OneDrive/Documents/Boechera/Growth_Chambers")
#setwd("C:/Users/dd66718/Downloads/Growth_Chambers")
#setwd("~/Documents/personnel/Denney")


fit <- read.csv("From Jill/fit.csv", header = T, stringsAsFactors=TRUE)
sapply(fit,class)

fit$Population <-as.factor(fit$Population)
fit$Tray <-as.factor(fit$Tray)
fit$Block<-paste(fit$Tray, fit$Treatment, sep="_")
fit$Block<-as.factor(fit$Block)
fit$Cohort<-as.factor(fit$Cohort)

fit $S_elev <- (fit $Elevation - mean(fit $Elevation, na.rm = TRUE)) / sd(fit $Elevation,na.rm = TRUE)
fit$day_transplant<-scale(fit$OD_transplant,center=TRUE, scale=TRUE)
fit$Elev <- fit$Elevation/1000
sapply(fit,class)

##Change baseline to improve plotting
fit$Temperature<-factor(fit$Temperature, levels = c("Low","High"))
fit$CO2<-factor(fit$CO2, levels = c("Low","High"))
fit$day_transplant<-scale(fit$OD_transplant,center=TRUE, scale=TRUE)

# Switch back
fit$Temperature<-factor(fit$Temperature, levels = c("High","Low"))
fit$CO2<-factor(fit$CO2, levels = c("High","Low"))

##############################################
######Survival: Logistic regression ##########
##############################################

survival_model <- glmmTMB(Season_survival ~Cohort+ Temperature*CO2 *S_elev+I(S_elev^2)*Temperature*CO2+ (1|Block) +(1| Genotype), data = fit, family=binomial(link="logit"))
Anova(survival_model, type="III")
summary(survival_model)
##Assess residuals
simulationOutput <- simulateResiduals(fittedModel= survival_model, plot = T, re.form = NULL)

cols=c("#56B4E9","#D55E00")
pred <- ggpredict(survival_model, terms = c("S_elev[all]", "Temperature", "CO2"), type = "re", interval="confidence")
plot(pred, show_data=TRUE, colors = cols, facet=TRUE,jitter=0.05)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Probability of survival")+ aes(linetype=group)



##Test the random effects
survival_nogeno <- glmmTMB(Season_survival ~ Cohort+Temperature*CO2 *S_elev+I(S_elev ^2)*Temperature*CO2+ (1|Block), data = fit, family=binomial(link="logit"))
anova(survival_model, survival_nogeno)

survival_noblock <- glmmTMB(Season_survival ~ Cohort+Temperature*CO2 *S_elev+I(S_elev ^2)*Temperature*CO2+ (1| Genotype), data = fit, family=binomial(link="logit"))
anova(survival_model, survival_noblock)




# Obtain slopes for each treatment combination. These slopes need to be exponentiated to calcuate odds ratios
linear <- emtrends(survival_model, specs = c("Temperature", "CO2"), var = "S_elev")
# pairwise comparisons of these slopes
pairs(linear)

quadratic <- emtrends(survival_model, specs = c("Temperature", "CO2"), var = "I(S_elev^2)")
# pairwise comparisons of these slopes
pairs(quadratic)

#For comparison 
library(broom.mixed)
##High temperature, High CO2
tidy(survival_model,conf.int=TRUE,exponentiate=TRUE,effects="fixed")
tidy(survival_model,conf.int=TRUE,exponentiate=FALSE,effects="fixed")


##############################################
######Probability of reproduction: Logistic regression ##########
##############################################
repro_model <- glmmTMB(Lifetime_Fruited ~ Cohort+Temperature*CO2 *S_elev+ (1|Block) +(1| Genotype), data = fit, family=binomial(link="logit"))
Anova(repro_model, type="III")
simulationOutput <- simulateResiduals(fittedModel= repro_model, plot = T, re.form = NULL)

cols=c("#56B4E9","#D55E00")
pred <- ggpredict(repro_model, terms = c("S_elev[all]", "Temperature", "CO2"), type = "re", interval="confidence")
plot(pred, show_data=TRUE, colors = cols, facet=TRUE,jitter=0.05)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Probability of reproduction")+ aes(linetype=group)



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
repro_nogeno <- glmmTMB(Lifetime_Fruited ~ Cohort +Temperature*CO2 *S_elev+ (1|Block), data = fit, family=binomial(link="logit"))
anova(repro_model, repro_nogeno)


repro_noblock <- glmmTMB(Lifetime_Fruited ~ Cohort +Temperature*CO2 *S_elev+ (1| Genotype), data = fit, family=binomial(link="logit"))
anova(repro_model, repro_noblock)


##############################################
######Fecundity: Silique number ##########
##############################################

fruited<-subset(fit, Lifetime_Fruited=="1")
fruited $S_elev <- (fruited $Elevation - mean(fruited $Elevation, na.rm = TRUE)) / sd(fruited $Elevation,na.rm = TRUE)

Number_siliques_elev <- glmmTMB(Mature_silique_number ~ Cohort+Temperature*CO2 *S_elev+ (1|Block) +(1| Genotype), data = fruited, family=Gamma(link="log"))
Anova(Number_siliques_elev, type="III")
simulationOutput <- simulateResiduals(fittedModel= Number_siliques_elev, plot = T, re.form = NULL)

cols=c("#56B4E9","#D55E00")
pred <- ggpredict(Number_siliques_elev, terms = c("S_elev[all]", "Temperature", "CO2"), type = "re", interval="confidence")
plot(pred, show_data=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Fecundity (number of fruits)")+ aes(linetype=group)


emtrends(Number_siliques_elev, specs = c("Temperature", "CO2"), var = "S_elev")


##Test the random effects
Number_siliques_elev_nogeno <- glmmTMB(Mature_silique_number ~ Cohort +Temperature*CO2 *S_elev+ (1|Block), data = fruited, family=Gamma(link="log"))
anova(Number_siliques_elev, Number_siliques_elev_nogeno)

Number_siliques_elev_noblock <- glmmTMB(Mature_silique_number ~ Cohort +Temperature*CO2 *S_elev+(1| Genotype), data = fruited, family=Gamma(link="log"))
anova(Number_siliques_elev, Number_siliques_elev_noblock)

## Test cohorts separately 
fruit1 <- subset(fruited, Cohort == "1")
fruit2 <- subset(fruited, Cohort == "2")
Number_siliques_elev1 <- glmer(Mature_silique_number ~ Temperature*CO2 *S_elev+ (1|Block) +(1| Genotype), data = fruit1, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(Number_siliques_elev1, type="III")
Number_siliques_elev2 <- glmer(Mature_silique_number ~ Temperature*CO2 *S_elev+ (1|Block) +(1| Genotype), data = fruit2, family=Gamma(link="log"),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(Number_siliques_elev2, type="III")

fecelev1 <- ggpredict(Number_siliques_elev1, c("S_elev[all",  "Temperature","CO2"), type = "re", interval = "confidence")
plot(fecelev1, show_data=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Fecundity (number of fruits)")+ aes(linetype=group)

fecelev2 <- ggpredict(Number_siliques_elev2, c("S_elev[all",  "Temperature","CO2"), type = "re", interval = "confidence")
plot(fecelev2, show_data=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Fecundity (number of fruits)")+ aes(linetype=group)
