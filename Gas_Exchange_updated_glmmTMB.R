library(tidyverse)
library(ggplot2)
library(car)
library(lme4)
library(effects)
require(glmmTMBTest)
require(MuMIn)
library(performance)
require(visreg)
library(emmeans)
library(coxme)
library(DHARMa)
library(gridExtra)
library(broom.mixed)
library(glmmTMB)

#setwd("~/Documents/personnel/Denney")
setwd("C:/Users/derek/OneDrive/Documents/Boechera/Growth_Chambers")

physiology <- read.csv("From Jill/phys.csv", header = T, stringsAsFactors=TRUE)
sapply(physiology,class)
physiology$Population <-as.factor(physiology$Population)
physiology$Tray <-as.factor(physiology$Tray)
physiology$Block<-paste(physiology$Tray, physiology$Treatment, sep="_")
physiology$Block<-as.factor(physiology$Block)
physiology$Cohort<-as.factor(physiology$Cohort)

physiology$S_elev<-scale(physiology$Elevation,center=TRUE, scale=TRUE)
physiology$Elev <- physiology$Elevation/1000
physiology$instant_WUE<- physiology$Licor_A/physiology$Licor_E
physiology$intrinsic_WUE<- physiology$Licor_A/physiology$Licor_gsw

sapply(physiology,class)

##Change baseline to improve plotting
physiology$Temperature<-factor(physiology$Temperature, levels = c("Low","High"))
physiology$CO2<-factor(physiology$CO2, levels = c("Low","High"))

physiology$Temperature<-factor(physiology$Temperature, levels = c("High","Low"))
physiology$CO2<-factor(physiology$CO2, levels = c("High","Low"))


##############################################
######Photosynthesis ##########
##############################################
photo <- glmmTMB(Licor_A ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = physiology)
Anova(photo,type="III")
#Use the DHARMa package to examine the residuals, which are good
simulationOutput <- simulateResiduals(fittedModel= photo, plot = T, re.form = NULL)


#Test the random effects
photo_no_gen <- glmmTMB(Licor_A ~ Cohort+Temperature*CO2* S_elev + (1|Block) , data = physiology)
photo_no_block <- glmmTMB(Licor_A ~ Cohort+Temperature*CO2* S_elev + (1| Genotype) , data = physiology)
anova(photo, photo_no_gen)
anova(photo, photo_no_block)

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
physiology$CO2 <- factor(physiology$CO2, levels = rev(levels(physiology$CO2)))
vioplot(Licor_A ~ CO2, data= physiology, plotCentre = "point",  pchMed = 23,  
        horizontal=FALSE,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00"),
        col=c("#56B4E9","#D55E00"), 
        ylab=y, xlab=x)+ 
  stripchart(Licor_A ~ CO2, data= physiology,  
             method = "jitter",col = "black", 
             vertical = T, pch = 1, add = TRUE)

##############################################
######Intrinsic WUE ##########
##############################################
hist(physiology$intrinsic_WUE)
#Square root transform intrinsic water-use efficiency to improve residuals
physiology $iWUE_sq <- sqrt(physiology $intrinsic_WUE)
hist(physiology $iWUE_sq)

intrinsicWUE <- glmmTMB(iWUE_sq ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype), data = physiology)
Anova(intrinsicWUE,type="III")
summary(intrinsicWUE)

#Test the random effects
WUE_no_gen <- glmmTMB(iWUE_sq ~ Cohort+Temperature*CO2* S_elev + (1|Block) , data = physiology)
WUE_no_block <- glmmTMB(iWUE_sq ~ Cohort+Temperature*CO2* S_elev + (1| Genotype) , data = physiology)
anova(intrinsicWUE, WUE_no_gen)
anova(intrinsicWUE, WUE_no_block)



plot(predictorEffects(intrinsicWUE, ~ S_elev), partial.residuals=TRUE)
#Use the DHARMa package to examine the residuals, which are ok
simulationOutput <- simulateResiduals(fittedModel= intrinsicWUE, plot = T, re.form = NULL)

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

WUE= ggplot(data = physiology, aes(x= Elevation, y= intrinsic_WUE, group=Temperature, shape=Temperature, color= CO2))+    geom_point(size = 4, position=position_dodge(width=0.8)) +    xlab("Elevation of origin (m)")+ ylab("Instrinsic WUE ") 
WUE +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", se=FALSE,  size=1.6)+ facet_grid(~ Temperature:CO2)


gr <- ref_grid(intrinsicWUE, cov.keep= c('S_elev', 'CO2'))
emm <- emmeans(gr, spec= c('S_elev', 'CO2'), level= 0.95)
gg <- ggplot(data= physiology, aes(x= S_elev, y= intrinsic_WUE)) +
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
physiology$CO2 <- factor(physiology$CO2, levels = rev(levels(physiology$CO2)))
vioplot(intrinsic_WUE ~  CO2  , data= physiology, plotCentre = "point",  pchMed = 23,  
        vertical=T,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00","#E1BE6A","#40B0A6"),
        col=c("#56B4E9","#D55E00", "#E1BE6A","#40B0A6"), 
        ylab="iWUE", xlab="CO2")
physiology$Temperature <- factor(physiology$Temperature, levels = rev(levels(physiology$Temperature)))
vioplot(intrinsic_WUE ~  Temperature  , data= physiology, plotCentre = "point",  pchMed = 23,  
        vertical=T,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00","#E1BE6A","#40B0A6"),
        col=c("#E1BE6A","#40B0A6","#E1BE6A","#40B0A6"), 
        ylab="iWUE", xlab="Temperature")
stripchart(intrinsic_WUE ~ Temperature, data= physiology,  
           method = "jitter",col = "gray", 
           vertical = T, pch = 1, add = TRUE, side = "left")


# Violin with all iwue treatments together
vioplot(intrinsic_WUE ~  CO2 +Temperature  , data= physiology, plotCentre = "point",  pchMed = 23,  
        vertical=T,
        colMed = "black",
        colMed2 = c("#56B4E9","#D55E00","#E1BE6A","#40B0A6"),
        col=c("#56B4E9","#D55E00", "#E1BE6A","#40B0A6"), 
        ylab="iWUE", xlab="CO2 and Temperature")
stripchart(intrinsic_WUE ~ CO2 + Temperature, data= physiology,  
           method = "jitter",col = "gray", 
           vertical = T, pch = 1, add = TRUE)

simple_model <- lm(intrinsic_WUE ~ CO2, data = physiology)

# Create a reference grid for CO2 only
gr_co2 <- ref_grid(simple_model, cov.keep = 'CO2')
emm_co2 <- emmeans(gr_co2, "CO2", level = 0.95, adjust = "tukey")
emm_co2

model <- lm(intrinsic_WUE ~ Temperature, data = physiology)
gr_temp <- ref_grid(model, cov.keep = 'Temperature')
emm_temp <- emmeans(gr_temp, "Temperature", level = 0.95, adjust = "tukey")
emm_temp

emms <- emmeans(intrinsicWUE, ~ CO2 + Temperature)
pairwise_comparisons <- pairs(emms)
pairwise_comparisons

emmeans_grid <- emmeans(intrinsicWUE, ~ CO2 + Temperature)
emmeans_grid


# Rerun model without temperature to verify slopes of CO2 treatments
intrinsicWUE_CO2 <- glmmTMB(intrinsic_WUE ~CO2* S_elev + (1|Block) +(1| Genotype) , data = physiology)
Anova(intrinsicWUE_CO2,type="III")
plot(predictorEffects(intrinsicWUE_CO2, ~ S_elev), partial.residuals=TRUE)

# CI for S_elev under low temperature
cc <- confint(intrinsicWUE,parm="beta_")  ## slow (~ 11 seconds)
ctab <- cbind(est=fixef(intrinsicWUE),cc)
ctab<-as.data.frame(ctab)
names(ctab)[2] <- "lower"
names(ctab)[3] <- "upper"


# CI for S_elev under higher temperature
physiology $co2_low<-factor(physiology $CO2, levels = c("Low","High"))
intrinsicWUE_low <- glmmTMB(intrinsic_WUE ~Temperature* co2_low* S_elev + (1|Block) +(1| Genotype) , data = physiology)
cc_low <- confint(intrinsicWUE_low,parm="beta_")  ## slow (~ 11 seconds)
ctab_low <- cbind(est=fixef(intrinsicWUE_low), cc_low)
ctab_low <-as.data.frame(ctab_low)
names(ctab_low)[2] <- "lower"
names(ctab_low)[3] <- "upper"
ctab
ctab_low

##############################################
######Transpiration ##########
##############################################
physiology$Licor_E_mmol <- (physiology$Licor_E *1000)

# Transform to improve residuals
physiology $Trans_sq <- sqrt(physiology $Licor_E_mmol)

hist(physiology$Licor_E_mmol)
hist(physiology$Trans_sq)

##model transpiration using the square root of transpiration
trans <- glmmTMB(Trans_sq ~Cohort+ Temperature*CO2* S_elev+ (1|Block) +(1|Genotype) , data = physiology)
Anova(trans,type="III")
emtrends(trans, specs = c("Temperature", "CO2"), var = "S_elev")
emtrends(trans, specs = "CO2", var = "S_elev")
plot(predictorEffects(trans, ~ S_elev), partial.residuals=TRUE)
simulationOutput <- simulateResiduals(fittedModel= trans, plot = T, re.form = NULL)
testOutliers(simulationOutput)
testDispersion(simulationOutput)

#Test the random effects
trans_no_gen <- glmmTMB(Trans_sq ~ Cohort+Temperature*CO2* S_elev + (1|Block) , data = physiology)
trans_no_block <- glmmTMB(Trans_sq ~ Cohort+Temperature*CO2* S_elev + (1| Genotype) , data = physiology)
anova(trans, trans_no_gen)
anova(trans, trans_no_block)

library(ggeffects)
cols=c("#D55E00","#56B4E9")
pred <- ggpredict(trans, terms = c("S_elev[all]", "CO2"), type = "re", interval="confidence")
plot(pred, show_residuals=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Transpiration rate (square root transformed)")

physiology $co2_low<-factor(physiology $CO2, levels = c("Low","High"))
transpiration_low <- glmmTMB(Trans_sq ~Cohort+ co2_low* S_elev+ (1|Block) +(1|Genotype) , data = physiology)
cc_low <- confint(transpiration_low,parm="beta_")  ## slow (~ 11 seconds)
ctab_low <- cbind(est=fixef(transpiration_low), cc_low)
ctab_low <-as.data.frame(ctab_low)
names(ctab_low)[2] <- "lower"
names(ctab_low)[3] <- "upper"

transpiration_high <- glmmTMB(Trans_sq ~Cohort+ CO2* S_elev+ (1|Block) +(1|Genotype) , data = physiology)
cc_high <- confint(transpiration_high,parm="beta_")  ## slow (~ 11 seconds)
ctab_high <- cbind(est=fixef(transpiration_high), cc_low)
ctab_high <-as.data.frame(ctab_high)
names(ctab_high)[2] <- "lower"
names(ctab_high)[3] <- "upper"


ctab_high
ctab_low

tidy(transpiration_low, conf.int=T, effects="fixed")


##To model with the raw data
transpiration <- glmmTMB(Licor_E_mmol ~Cohort+ Temperature*CO2* S_elev+ (1|Block) +(1|Genotype) , data = physiology)
Anova(transpiration,type="III")
emtrends(transpiration, specs = c("Temperature", "CO2"), var = "S_elev")
emtrends(transpiration, specs = "CO2", var = "S_elev")
plot(predictorEffects(transpiration, ~ S_elev), partial.residuals=TRUE)
#Use the DHARMa package to examine the residuals, which are ok
simulationOutput <- simulateResiduals(fittedModel= transpiration, plot = T, re.form = NULL)


y <- expression("Transpiration rate (mmol m"^-2*" s"^-1*")")
c <- c('High'="High CO2", 'Low' = "Low CO2")

physiology$Temperature <- factor(physiology$Temperature, levels = rev(levels(physiology$Temperature)))
physiology$CO2 <- factor(physiology$CO2, levels = rev(levels(physiology$CO2)))

gr <- ref_grid(transpiration, cov.keep= c('S_elev', 'CO2'))
emm <- emmeans(gr, spec= c('S_elev', 'CO2'), level= 0.95)
gg <- ggplot(data= physiology, aes(x= S_elev, y= Licor_E_mmol)) +
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
transpiration_CO2 <- glmmTMB(Licor_E_mmol ~ CO2* S_elev + (1|Block) +(1| Genotype) , 
                          data = physiology)
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
physiology $co2_low<-factor(physiology $CO2, levels = c("Low","High"))
transpiration_low <- glmmTMB(Licor_E_mmol ~ Temperature* co2_low* S_elev + (1|Block) +(1| Genotype) , 
                          data = physiology)
cc_low <- confint(transpiration_low,parm="beta_")  ## slow (~ 11 seconds)
ctab_low <- cbind(est=fixef(transpiration_low), cc_low)
ctab_low <-as.data.frame(ctab_low)
names(ctab_low)[2] <- "lower"
names(ctab_low)[3] <- "upper"

ctab
ctab_low


tidy(transpiration, conf.int=T, effects="fixed")


##############################################
######Stomatal Conductance ##########
##############################################
physiology$licor_gsw_mmol <- physiology$Licor_gsw*1000

#Square root transform stomatal conductance to acheive normality and homoscedasticity of resdiuals
physiology $gsw_sq <- sqrt(physiology $Licor_gsw)
hist(na.omit(physiology$licor_gsw_mmol))
hist(na.omit(physiology$gsw_sq))

hist(physiology$Licor_gsw)

condu <- glmmTMB(gsw_sq ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = physiology)
Anova(condu,type="III")
summary(condu)
plot(predictorEffects(condu, ~ S_elev), partial.residuals=TRUE)
#Use the DHARMa package to examine the residuals, which are ok
simulationOutput <- simulateResiduals(fittedModel= condu, plot = T, re.form = NULL)

#Test the random effects
cond_no_gen <- glmmTMB(gsw_sq ~ Cohort+Temperature*CO2* S_elev + (1|Block) , data = physiology)
cond_no_block <- glmmTMB(gsw_sq ~ Cohort+Temperature*CO2* S_elev + (1| Genotype) , data = physiology)
anova(condu, cond_no_block)
anova(condu, cond_no_gen)

cols=c("#D55E00","#56B4E9")
pred <- ggpredict(condu, terms = c("S_elev[all]", "CO2"), type = "re", interval="confidence")
plot(pred, show_residuals=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Stomatal conductance")




cond= ggplot(data = physiology, aes(x= Elevation, y= Licor_gsw, group=Temperature, shape=Temperature, color= CO2))+    geom_point(size = 4, position=position_dodge(width=0.8)) +    xlab("Elevation of origin (m)")+ ylab("Stomatal conductance") 
cond +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", se=FALSE,  size=1.6)+ facet_grid(~ Temperature:CO2)

y <- expression("Stomotal conductance (mmol m"^-2*" s"^-1*")")
c <- c('High' = "High CO2", 'Low' = "Low CO2")
gr <- ref_grid(condu, cov.keep= c('S_elev', 'CO2'))
emm <- emmeans(gr, spec= c('S_elev', 'CO2'), level= 0.95)
gg <- ggplot(data= physiology, aes(x= S_elev, y= licor_gsw_mmol)) +
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
condu_CO2 <- glmmTMB(licor_gsw_mmol ~CO2* S_elev + (1|Block) +(1| Genotype) , 
                  data = physiology, control=glmmTMBControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
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
physiology $co2_low<-factor(physiology $CO2, levels = c("Low","High"))
condu_low <- glmmTMB(licor_gsw_mmol ~Temperature* co2_low* S_elev + (1|Block) +(1| Genotype) , 
                  data = physiology, control=glmmTMBControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
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



physiology$quartile_group <- cut(physiology$S_elev, 
                              breaks = quantile(physiology$S_elev, 
                                                probs = c(0, 0.25, 0.5, 0.75, 1)), 
                              labels = FALSE, include.lowest = T)

# Assign values 1 through 4 to quartile groups
physiology$quartile_group <- match(physiology$quartile_group, unique(physiology$quartile_group))


gs <- physiology %>% group_by(CO2, quartile_group) %>% summarise(mean_gsw = mean(Licor_gsw, na.rm= T))
summary(gs)

plot(gs$CO2, gs$mean_gsw, col = gs$quartile_group)



##########################
### Field Licor Data #####
###########################
sco <- read.delim("From Jill/Sco_ecophys.txt", header = T)
sco $S_elev <- (sco $elevation..m. - mean(sco $elevation..m., na.rm = TRUE)) / sd(sco $elevation..m.,na.rm = TRUE)
head(sco)

# Photosynthesis
plot(sco$photosynthesis)

sco_photo <- glmmTMB(photosynthesis ~ elevation..m. + (1|Quad) +(1| genotype), data = sco)
Anova(sco_photo, type="III")
#Use the DHARMa package to examine the residuals, which are ok
simulationOutput <- simulateResiduals(fittedModel= sco_photo, plot = T, re.form = NULL)

#Check random effects
sco_photo_noblock <- glmmTMB(photosynthesis ~ elevation..m. + (1| genotype), data = sco)
sco_photo_nogeno <- glmmTMB(photosynthesis ~ elevation..m. + (1|Quad) , data = sco)
anova(sco_photo_noblock, sco_photo)
anova(sco_photo_nogeno, sco_photo)

summary(sco_photo)
plot(predictorEffects(sco_photo, ~ elevation..m.), partial.residuals=T)

pred_field <- ggpredict(sco_photo, terms = c("elevation..m."), type = "re", interval="confidence")
plot(pred_field, show_residuals=TRUE, )+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Assimilation rate")


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
