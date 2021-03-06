### Analysis of thesis data for manuscript 
### 
### Bird abundances versus habitat variables; 
### post-breeding habitat breadth and selectivity of birds


library(tidyverse)
library(lme4)
library(bbmle)
library(DHARMa)
# library(emmeans)
# library(effects)
library(glmmTMB)
library(vegan)



load("data/01_dataimport.RData")


##### Checking for correlation between habitat variables ############


# Strong correlation at high CCavg is a direct physical result because with higher
# density of trees, there are less gaps overall and hence the variability in
# cover is also lower.
cor.test(habvar$CCavg, habvar$TreeDens, method = "pearson")
cor.test(habvar$CCavgsd, habvar$CCavg, method = "pearson")
cor.test(habvar$CCavgsd, habvar$TreeDens, method = "pearson")

plot(habvar$CCavg, habvar$TreeDens)
plot(habvar$CCavgsd, habvar$CCavg)
plot(habvar$CCavgsd, habvar$TreeDens)
# CCavgsd does not correlate with TreeDens. 
# Hence, I can either include just CCavg in the models, which captures partly the 
# essences of CCavgsd and TreeDens, or I can include the latter two in the models 
# and exclude CCavg where the two variables would then capture a broader scope
# of variability than would have been possible with CCavg alone, but might miss
# some of the core essence of CCavg

cor.test(habvar$TreePropDeci, habvar$TreeDiv, method="spearman")
plot(habvar$TreePropDeci, habvar$TreeDiv)
# # devtools::install_github("ProcessMiner/nlcor")
# library(nlcor)
# nlcor(habvar$TreePropDeci, habvar$TreeDiv, plt=T)
# I think it was Salek 2016 who said proportion of deciduous trees was a better
# predictor than tree diversity, and since the (non-linear) correlation is very 
# strong here, I think it is better to include TreePropDeci and exclude TreeDiv.

# Using CCavgsd+TreeDens over using CCavg is not statistically significantly better
# but has 1.79 less residual deviance and only uses one extra df, so it might be
# the way to go.
# Using TreeDiv over TreePropDeci has, technically, lower resid. dev. and lower
# Cp, but the difference is negligible (0.25).


cor.test(habvar$TreePropDeci, habvar$TreeDens, method="pearson")
plot(habvar$TreePropDeci, habvar$TreeDens) # correlated
cor.test(habvar$TreePropDeci, habvar$CCavg, method="pearson")

cor.test(habvar$UndDens, habvar$TreeDens, method="pearson")
cor.test(habvar$UndDiv, habvar$TreeDens, method="pearson")
cor.test(habvar$UndRich, habvar$TreeDens, method="pearson")
cor.test(habvar$UndDens, habvar$CCavgsd, method="pearson")
cor.test(habvar$UndDiv, habvar$CCavgsd, method="pearson")
cor.test(habvar$UndRich, habvar$CCavgsd, method="pearson")



### ###

### Analysis 1: Building models for all birds ####

boxplot(m_all$All_Abun ~ m_all$Period)
plot(m_all$All_Abun ~ m_all$Days)


all.0 <- glmer(All_Abun ~ Observer + (1|Point) + (1|Days),
               data = mdata, family=poisson())
anova(all.0)


# adding Week and testing for polynomial 
all.1a <- update(all.0, .~. + Week)
all.1b <- update(all.0, .~. + poly(Week,2))
all.1c <- update(all.0, .~. + poly(Week,3))
anova(all.1a, all.1c) # 2nd and 3rd not better than linear
AICctab(all.1a, all.1b, all.1c)
rm(list = c("all.1a", "all.1b", "all.1c"))

all.1 <- update(all.0, .~. + Week)
summary(all.1)


# testing nuisance variables like time of count during day 
all.2a <- update(all.1, .~. + CoD)
all.2b <- update(all.1, .~. + poly(CoD,2))
anova(all.1, all.2a) # neither linear nor poly significant (AIC and chi-sq both)
all.2c <- update(all.1, .~. + Weather)
summary(all.2c)
anova(all.1, all.2c) # Weather not significant

# finally moving to predictors of interest
all.2d <- update(all.1, .~. + HabClass)
all.2e <- update(all.1, .~. + Moss)
AICctab(all.1, all.2d, all.2e)
anova(all.1, all.2d)
summary(all.2d)
# summary function shows that Interior is marginally significant (0.04) but other levels
# are not. AIC is not different at all. So will have to ignore.
# summary() p value cannot be relied on.
summary(all.2e)
anova(all.1, all.2e) # Moss ns too

all.2f <- update(all.1, .~. + CCavgscaled)
all.2g <- update(all.1, .~. + log(CCavgsd) + log(TreeDens))
AICctab(all.1, all.2f, all.2g) # as expected, sd+TreeDens better than simple CC
summary(all.2g) # dAIC 2.5 so add
all.2h <- update(all.2g, .~. + log(CCavgsd):log(TreeDens))
summary(all.2h) # drop interaction 
drop1(all.2g) # but keep main effects
all.2i <- update(all.1, .~. + poly(log(CCavgsd),2) + poly(log(TreeDens),2))
all.2j <- update(all.1, .~. + poly(log(CCavgsd),2) + log(TreeDens))
all.2k <- update(all.1, .~. + log(CCavgsd) + poly(log(TreeDens),2))
AICctab(all.2g, all.2i, all.2j, all.2k)
# two poly is worst. poly of TreeDens is technically better than poly of CCavgsd
# but both worse than simple linear

rm(list=c("all.2a","all.2b","all.2c","all.2d","all.2e","all.2f",
          "all.2g","all.2h","all.2i","all.2j","all.2k"))
all.2 <- update(all.1, .~. + log(CCavgsd) + log(TreeDens))
summary(all.2)


# interaction of predictors with Week (of pertinence to question)
all.3a <- update(all.2, .~. + Week:log(CCavgsd))
all.3b <- update(all.2, .~. + Week:log(TreeDens))
all.3c <- update(all.2, .~. + Week:(log(CCavgsd)+log(TreeDens)))
AICctab(all.2, all.3a, all.3b, all.3c) 
# TreeDens seems more important in general (3b over 3a) but both recommended by AIC

rm(list=c("all.3a", "all.3b", "all.3c"))
all.3 <- update(all.2, .~. + Week:(log(CCavgsd)+log(TreeDens)))
summary(all.3)


# testing TreePropDeci and perhaps dropping others for this
all.4a <- update(all.3, .~. + TPDscaled)
all.4b <- update(all.3, .~. + poly(TPDscaled,2))
anova(all.4a, all.4b) # poly not preferred
drop1(all.4a)
all.4c <- update(all.3, .~. + TPDscaled + Week:TPDscaled)
all.4d <- update(all.1, .~. + TPDscaled + Week:TPDscaled) # replacing other pred.
all.4e <- update(all.1, .~. + poly(TPDscaled,2))
AICctab(all.1, all.2, all.3, all.4a, all.4c, all.4d, all.4e)
anova(all.1, all.3)
# TreePropDeci not useful
# trying TreeDiv
all.4f <- update(all.3, .~. + TreeDiv) # failed to converge 0.0112 (tot 0.002)
all.4g <- update(all.3, .~. + poly(TreeDiv,2)) # failed to converge 0.0248 (tot 0.002)
AICctab(all.3, all.4f, all.4g)
# TreeDiv not useful

# trying UndDens
all.4h <- update(all.3, .~. + sqrt(UndDens))
all.4i <- update(all.3, .~. + poly(sqrt(UndDens),2))
all.4j <- update(all.1, .~. + poly(sqrt(UndDens),2))
all.4k <- update(all.1, .~. + Week*poly(sqrt(UndDens),2))
AICctab(all.3, all.4h, all.4i, all.4j, all.4k)
anova(all.3, all.4i)
# UndDens not useful
# trying DOM
all.4l <- update(all.3, .~. + DOM)
anova(all.3, all.4l)
all.4m <- update(all.3, .~. + DOM + Week:DOM)
AICctab(all.3, all.4l, all.4m)
# DOM main effect and interaction with Week


rm(list=c("all.4a","all.4b","all.4c","all.4d","all.4e","all.4f","all.4g",
          "all.4h","all.4i","all.4j","all.4k","all.4l","all.4m"))
all.4 <- update(all.3, .~. + DOM + Week:DOM)


summary(all.4) # failure to converge 0.0117 (tot 0.002)
hist(resid(all.4)) # very Gaussian(!), except for those two outliers
all.4sim <- simulateResiduals(all.4, plot=T, n=1000) # beautiful :)
plot(all.4sim)
testDispersion(all.4sim) # 1.2055 but p value sig because of sample size

hist(all.4sim) # uniform as should be
plotResiduals(all.4sim, form=mdata$Week)
plotResiduals(all.4sim, form=mdata$CCavgsd)
plotResiduals(all.4sim, form=mdata$TreeDens)
plotResiduals(all.4sim, form=mdata$DOM)
# all good

# exploring convergence issue 
all.4opt <- allFit(all.4)
summary(all.4opt)$which.OK # only 1/7 optimisers failed
summary(all.4opt)$llik # all pretty similar values so I think bobyqa good choice
summary(all.4opt)$fixef
# use bobyqa


all.5 <- update(all.4, control=glmerControl(optimizer = "bobyqa"))
summary(all.5)

all.5sim <- simulateResiduals(all.5, n=1000)
plot(all.5sim)
testDispersion(all.5sim)
hist(all.5sim)


# further selection

all.6a <- update(all.5, .~. - Week:(log(CCavgsd)+log(TreeDens)))
all.6b <- update(all.5, .~. - Week:(log(CCavgsd)+log(TreeDens)) - log(CCavgsd))
all.6c <- update(all.5, .~. - Week:(log(CCavgsd)+log(TreeDens)) - log(TreeDens))
all.6d <- update(all.5, .~. - Week*(log(CCavgsd)+log(TreeDens)))

AICctab(all.1, all.2, all.3, all.4, all.5, all.6a, all.6b, all.6c, all.6d,
        weights=T, base=T, logLik=T)
rm(list=c("all.6a","all.6b","all.6c","all.6d"))

all.6 <- update(all.5, .~. - Week:(log(CCavgsd)+log(TreeDens)))
AICctab(all.1, all.2, all.3, all.4, all.5, all.6, weights=T, base=T, logLik=T)
summary(all.6)

# interactions between predictors. CCavgsd:DOM sig but mostly for Moss (only 2 Points)
# all.6e <- update(all.6, .~. + DOM*log(CCavgsd))
# all.6f <- update(all.6, .~. + DOM*log(TreeDens))
# all.6g <- update(all.6, .~. + DOM*(log(CCavgsd)+log(TreeDens)))
# AICctab(all.6, all.6e, all.6f, all.6g)
# 
# all.6h <- update(all.6, .~. + DOM*log(CCavgsd) - log(TreeDens))
# AICctab(all.6, all.6e, all.6f, all.6g, all.6h)
# 
# summary(all.6e)
# plot(allEffects(all.6e))
# plot(allEffects(all.6h))
# rm(list=c("all.6e","all.6f","all.6g","all.6h"))



## Matching fruiting phenology to abundance ##


allfruit.0 <- glmmTMB(All_Abun ~ Observer + (1|Point) + (1|Days),
                      data=mdata_fr, family=poisson)
summary(allfruit.0)

allfruit.1 <- update(allfruit.0, .~. + Rubus + Vaccinium + Fran_alnu + Sorb_aucu)
allfruit.2 <- update(allfruit.0, .~. + poly(Rubus,2) + poly(Vaccinium,2) + 
                       poly(Fran_alnu,2) + Sorb_aucu)
allfruit.3 <- update(allfruit.0, .~. + poly(Rubus,3) + poly(Vaccinium,3) + 
                       poly(Fran_alnu,3) + Sorb_aucu)
allfruit.4 <- update(allfruit.0, .~. + poly(Rubus,4) + poly(Vaccinium,4) + 
                       poly(Fran_alnu,4) + Sorb_aucu)

AICctab(allfruit.0, allfruit.1, allfruit.2)

### ###

### Analysis 2: Building models for guilds ####

# invertebrate has 22.6% zeroes (94/416); omnivore 7.9%

plot(m_data_A1_inv$GuildAbun ~ m_data_A1_inv$Week)
plot(m_data_A1_omn$GuildAbun ~ m_data_A1_omn$Week)
# omni seem to show much clearer pattern (increase). migrants

### Analysis 2: Building models for guilds: invertebrate-feeding ####



# all the previous models I tried and got to dead ends with (in CodeBin):
# rm(list=c("inv.0min","inv.1min","inv.2min","inv.3min","inv.3minsim",
#           "inv.4min","inv.4min2","inv.4minNB","inv.4minnorm","inv.4minnormsim",
#           "inv.4minsim","inv.NB0","inv.NB1","invTMB.0","invTMB.1","invTMB.2",
#           "invTMB.3","invTMB.4","invTMB.4NB","invTMB.4NBsim","invTMB.4poi",
#           "invTMB.4poisim","invTMB.4sim"))


# finally decided on a mixed model with Poisson error using glmmTMB with no
# zero inflation. singular fit issues resolved by using glmmTMB over glmer, and estimated 
# parameters in both methods were pretty much the same.


inv.0 <- glmmTMB(GuildAbun ~ Observer + (1|Point),
                 data = m_data_A1_inv, family=poisson())
summary(inv.0)
inv.0sim <- simulateResiduals(inv.0)
plot(inv.0sim)
hist(inv.0sim)


# adding Week and testing poly
inv.1a <- update(inv.0, .~. + Week)
inv.1b <- update(inv.0, .~. + poly(Week,2))
AICctab(inv.0, inv.1a, inv.1b)
# main effect not significant
# nuisance variables like CoD and Weather
inv.1c <- update(inv.0, .~. + CoD)
inv.1d <- update(inv.0, .~. + poly(CoD,2))
inv.1e <- update(inv.0, .~. + Weather)
AICctab(inv.0, inv.1c, inv.1d, inv.1e)
# not significant

# predictors of interest
inv.1f <- update(inv.0, .~. + HabClass)
inv.1g <- update(inv.0, .~. + Moss)
AICctab(inv.0, inv.1f, inv.1g)
inv.1h <- update(inv.0, .~. + HabClass*Week)
AICctab(inv.0, inv.1f, inv.1g, inv.1h)
# habclass not useful at all. Moss is.
inv.1i <- update(inv.0, .~. + Moss*Week)
AICctab(inv.0, inv.1g, inv.1i)
summary(inv.1i)
# interaction has 1.4 less AIC on 2 df, so will not include, 
# even though it is pertinent to question. 
inv.1j <- update(inv.0, .~. + Moss*poly(Week,2))
AICctab(inv.0, inv.1g, inv.1i, inv.1j)

rm(list=c("inv.1a","inv.1b","inv.1c","inv.1d","inv.1e","inv.1f","inv.1g","inv.1h","inv.1i","inv.1j"))
inv.1 <- update(inv.0, .~. + Moss)
anova(inv.0, inv.1)


inv.2a <- update(inv.1, .~. + log(CCavgsd) + log(TreeDens))
inv.2b <- update(inv.1, .~. + poly(log(CCavgsd),2) + log(TreeDens))
inv.2c <- update(inv.1, .~. + log(CCavgsd) + poly(log(TreeDens),2))
inv.2d <- update(inv.1, .~. + poly(log(CCavgsd),2) + poly(log(TreeDens),2))
inv.2e <- update(inv.1, .~. + log(CCavgsd)*Week + log(TreeDens)*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e)
# although still ns, the interaction with Week has the closest AIC to inv.1
inv.2f <- update(inv.1, .~. + poly(log(CCavgsd),2)*Week + poly(log(TreeDens),2)*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f)

inv.2g <- update(inv.1, .~. + log(CCavgsd)*Week + log(TreeDens))
inv.2h <- update(inv.1, .~. + log(CCavgsd) + log(TreeDens)*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f, inv.2g, inv.2h)
inv.2i <- update(inv.1, .~. + log(CCavgsd)*Week + poly(log(TreeDens),2))
inv.2j <- update(inv.1, .~. + log(CCavgsd) + poly(log(TreeDens),2)*Week)
AICctab(inv.1, inv.2a, inv.2c, inv.2e, inv.2h, inv.2i, inv.2j)
inv.2k <- update(inv.1, .~. + log(TreeDens)*Week)
inv.2l <- update(inv.1, .~. + log(TreeDens))
inv.2m <- update(inv.1, .~. + log(CCavgsd)*Week)
inv.2n <- update(inv.1, .~. + log(CCavgsd))
AICctab(inv.1, inv.2l, inv.2m, inv.2n, inv.2k)
# not useful!

inv.2o <- update(inv.1, .~. + TPDscaled)
inv.2p <- update(inv.1, .~. + TPDscaled*Week)
inv.2q <- update(inv.1, .~. + poly(TPDscaled,2))
inv.2r <- update(inv.1, .~. + poly(TPDscaled,2)*Week)
AICctab(inv.1, inv.2p, inv.2q, inv.2r, inv.2o)
inv.2s <- update(inv.1, .~. + TreePropDeci)
inv.2t <- update(inv.1, .~. + poly(TreePropDeci,2))
AICctab(inv.1, inv.2o, inv.2p, inv.2q, inv.2t)
inv.2u <- update(inv.1, .~. + poly(TPDscaled,3))
inv.2v <- update(inv.1, .~. + poly(TPDscaled,4))
AICctab(inv.1, inv.2o, inv.2q, inv.2u, inv.2v)
# only linear main effect of TPD
rm(list=c("inv.2a","inv.2b","inv.2c","inv.2d","inv.2e","inv.2f","inv.2g","inv.2h","inv.2i",
          "inv.2j","inv.2k","inv.2l","inv.2m","inv.2n","inv.2o","inv.2p","inv.2q","inv.2r",
          "inv.2s","inv.2t","inv.2u","inv.2v"))
inv.2 <- update(inv.1, .~. + TPDscaled)


inv.3a <- update(inv.2, .~. - TPDscaled + TreeDiv)
inv.3b <- update(inv.2, .~. + TreeDiv)
inv.3c <- update(inv.2, .~. + poly(TreeDiv,2))
inv.3d <- update(inv.2, .~. + TreeDiv*Week)
AICctab(inv.2, inv.3a, inv.3b, inv.3c, inv.3d)
# TreeDiv not important
inv.3e <- update(inv.2, .~. + sqrt(UndDens))
inv.3f <- update(inv.2, .~. + poly(sqrt(UndDens),2))
inv.3g <- update(inv.2, .~. + sqrt(UndDens)*Week)
inv.3h <- update(inv.2, .~. + poly(sqrt(UndDens),2)*Week)
AICctab(inv.2, inv.3e, inv.3f, inv.3g, inv.3h)
# UndDens not important
inv.3i <- update(inv.2, .~. + TreeRich)
inv.3j <- update(inv.2, .~. + TreeRich*Week)
AICctab(inv.2, inv.3i, inv.3j)
inv.3k <- update(inv.2, .~. + poly(TreeRich,2))
AICctab(inv.2, inv.3i, inv.3k)
# interesting how abundance increases with TreeRich but decreases with TPD (and poly ns)

rm(list=c("inv.3a","inv.3b","inv.3c","inv.3d","inv.3e","inv.3f","inv.3g","inv.3h","inv.3i",
          "inv.3j","inv.3k"))
inv.3 <- update(inv.2, .~. + TreeRich)
summary(inv.3)


# finally, testing DOM and comparing with Moss
inv.4a <- update(inv.3, .~. + DOM) 
inv.4b <- update(inv.3, .~. + DOM*Week)
inv.4c <- update(inv.3, .~. + DOM - Moss) 
inv.4d <- update(inv.3, .~. + DOM*Week - Moss)
AICctab(inv.3, inv.4a, inv.4b, inv.4c, inv.4d) # 4c better only 0.1 dAICc on 1df
summary(inv.4c)
summary(inv.4a)
inv.4e <- update(inv.3, .~. + DOM - TreeRich)
AICctab(inv.3, inv.4a, inv.4c, inv.4e)
inv.4f <- update(inv.3, .~. + DOM - TreeRich - Moss)
AICctab(inv.3, inv.4a, inv.4c, inv.4e, inv.4f)
# effect of TreeRich prob. caused by STRONG effect of Carex. 
# All Carex Points have low TreeRich
rm(list=c("inv.4a","inv.4b","inv.4c","inv.4d","inv.4e","inv.4f"))
inv.4 <- update(inv.3, .~. + DOM - TreeRich - Moss)
summary(inv.4)


inv.4sim <- simulateResiduals(inv.4, n=500)
plot(inv.4sim)
hist(inv.4sim)
testDispersion(inv.4sim)
testZeroInflation(inv.4sim)



# glmer seems to be working fine with this final model too. The problem earlier might have
# been due to Moss and DOM together, because trying a glmer with inv.4 + Moss gives 
# singular fit again. Regardless, glmmTMB performs better overall.


### ###

### Analysis 2: Building models for guilds: omnivore-feeding ####


# all the previous models I tried and got to dead ends with (in CodeBin):
# rm(list=c("omnGLM.0","omnGLM.1","omnGLM.2","omnGLM.3","omnGLM.4","omnGLM.4sim"))



omn.0 <- glmmTMB(GuildAbun ~ Observer + (1|Point), data=m_data_A1_omn, family=poisson)
omn.0sim <- simulateResiduals(omn.0)
plot(omn.0sim)
# KS test, dispersion test, outlier test all sig
omn.0NB <- glmmTMB(GuildAbun ~ Observer + (1|Point), data=m_data_A1_omn, family = nbinom2)
omn.0NBsim <- simulateResiduals(omn.0NB)
plot(omn.0NBsim)


# adding Week and testing poly
omn.1a <- update(omn.0NB, .~. + Week)
omn.1b <- update(omn.0NB, .~. + poly(Week,2))
AICctab(omn.0, omn.0NB, omn.1a, omn.1b)
# linear effect of Week very good
rm(list=c("omn.1a","omn.1b"))
omn.1 <- update(omn.0NB, .~. + Week)


# nuisance variables like CoD and Weather
omn.2a <- update(omn.1, .~. + CoD)
omn.2b <- update(omn.1, .~. + poly(CoD,2))
omn.2c <- update(omn.1, .~. + Weather)
AICctab(omn.1, omn.2a, omn.2b, omn.2c)
# not significant
# predictors of interest
omn.2d <- update(omn.1, .~. + HabClass)
omn.2e <- update(omn.1, .~. + Moss)
AICctab(omn.1, omn.2d, omn.2e) # habclass sig
omn.2f <- update(omn.1, .~. + HabClass + Moss)
AICctab(omn.1, omn.2d, omn.2e, omn.2f) # Moss ns
omn.2g <- update(omn.1, .~. + HabClass*Week)
omn.2h <- update(omn.1, .~. + Moss*Week)
AICctab(omn.1, omn.2d, omn.2g, omn.2e, omn.2h)
# only main effect of HabClass
rm(list=c("omn.2a","omn.2b","omn.2c","omn.2d","omn.2e","omn.2f","omn.2g","omn.2h"))
omn.2 <- update(omn.1, .~. + HabClass)
anova(omn.1, omn.2)


omn.3a <- update(omn.2, .~. + log(CCavgsd) + log(TreeDens))
omn.3b <- update(omn.2, .~. + poly(log(CCavgsd),2) + log(TreeDens))
omn.3c <- update(omn.2, .~. + log(CCavgsd) + poly(log(TreeDens),2))
omn.3d <- update(omn.2, .~. + poly(log(CCavgsd),2) + poly(log(TreeDens),2))
omn.3e <- update(omn.2, .~. + log(CCavgsd)*Week + log(TreeDens)*Week)
AICctab(omn.2, omn.3a, omn.3b, omn.3c, omn.3d, omn.3e)
# no poly
omn.3f <- update(omn.2, .~. + log(CCavgsd)*Week + log(TreeDens))
omn.3g <- update(omn.2, .~. + log(CCavgsd) + log(TreeDens)*Week)
AICctab(omn.2, omn.3a, omn.3b, omn.3c, omn.3d, omn.3e, omn.3f, omn.3g)
# interaction CCavgsd:Week + TD is best
omn.3h <- update(omn.2, .~. + log(CCavgsd)*Week) 
omn.3i <- update(omn.2, .~. + log(CCavgsd))
AICctab(omn.2, omn.3a, omn.3f, omn.3h, omn.3i)

rm(list=c("omn.3a","omn.3b","omn.3c","omn.3d","omn.3e","omn.3f","omn.3g","omn.3h","omn.3i"))
omn.3 <- update(omn.2, .~. + log(CCavgsd)*Week + log(TreeDens))


omn.4a <- update(omn.3, .~. + TPDscaled - log(TreeDens))
omn.4b <- update(omn.3, .~. + TPDscaled) 
AICctab(omn.3, omn.4a, omn.4b) # keep
omn.4c <- update(omn.3, .~. + TPDscaled*Week)
omn.4d <- update(omn.3, .~. + poly(TPDscaled,2))
omn.4e <- update(omn.3, .~. + poly(TPDscaled,2)*Week)
AICctab(omn.3, omn.4a, omn.4b, omn.4c, omn.4d)
# not significant
omn.4f <- update(omn.3, .~. + TreeDiv)
omn.4g <- update(omn.3, .~. + poly(TreeDiv,2))
omn.4h <- update(omn.3, .~. + TreeDiv*Week)
omn.4i <- update(omn.3, .~. + poly(TreeDiv,2)*Week)
AICctab(omn.3, omn.4f, omn.4g, omn.4h, omn.4i)
# not significant
omn.4j <- update(omn.3, .~. + UndDensscaled)
omn.4k <- update(omn.3, .~. + poly(UndDensscaled,2))
omn.4l <- update(omn.3, .~. + UndDensscaled*Week)
omn.4m <- update(omn.3, .~. + poly(UndDensscaled,2)*Week)
AICctab(omn.3, omn.4j, omn.4k, omn.4l, omn.4m)
# not significant
omn.4n <- update(omn.3, .~. + TreeRich)
omn.4o <- update(omn.3, .~. + TreeRich*Week)
omn.4p <- update(omn.3, .~. + poly(TreeRich,2))
omn.4q <- update(omn.3, .~. + poly(TreeRich,2)*Week)
AICctab(omn.3, omn.4n, omn.4o, omn.4p, omn.4q)
# not significant
omn.4r <- update(omn.3, .~. + DOM) 
omn.4s <- update(omn.3, .~. + DOM*Week)
omn.4t <- update(omn.3, .~. + DOM - HabClass)
AICctab(omn.3, omn.4r, omn.4s, omn.4t)
# DOM has 3.4 dAICc but 4 df ##
rm(list=c("omn.4a","omn.4b","omn.4c","omn.4d","omn.4e","omn.4f","omn.4g","omn.4h","omn.4i",
          "omn.4j","omn.4k","omn.4l","omn.4m","omn.4n","omn.4o","omn.4p","omn.4q","omn.4r",
          "omn.4s","omn.4t"))
# no selection. stick with omn.3

summary(omn.3)
omn.3sim <- simulateResiduals(omn.3) 
plot(omn.3sim)
hist(omn.3sim) 
testDispersion(omn.3sim)
testZeroInflation(omn.3sim)




# glmer is not working for this final model, but glmmTMB took time for each model



### ###


### Analysis 2: Building models for migration classification ####

plot(m_data_A1_NM$Abun ~ m_data_A1_NM$Week)
plot(m_data_A1_PM$Abun ~ m_data_A1_PM$Week)
plot(m_data_A1_CM$Abun ~ m_data_A1_CM$Week)


### ###

### Analysis 2: Building models for migration classification: Non-migrants ####


NM.0 <- glmmTMB(Abun ~ Observer + (1|Point) + (1|Days),
                data = m_data_A1_NM, family=poisson())
summary(NM.0)
NM.0sim <- simulateResiduals(NM.0)
plot(NM.0sim)
hist(NM.0sim)


# adding Week and testing poly
NM.1a <- update(NM.0, .~. + Week)
NM.1b <- update(NM.0, .~. + poly(Week,2))
AICctab(NM.0, NM.1a, NM.1b)

rm(list=c("NM.1a","NM.1b"))
NM.1 <- update(NM.0, .~. + Week)
anova(NM.0, NM.1)


# nuisance variables like CoD and Weather
NM.2a <- update(NM.1, .~. + CoD)
NM.2b <- update(NM.1, .~. + poly(CoD,2))
NM.2c <- update(NM.1, .~. + Weather)
NM.2d <- update(NM.1, .~. + CoD + Weather)
AICctab(NM.1, NM.2a, NM.2b, NM.2c, NM.2d)
summary(NM.2d)

rm(list=c("NM.2a", "NM.2b","NM.2c","NM.2d"))
NM.2 <- update(NM.1, .~. + CoD + Weather)


# predictors of interest
NM.3a <- update(NM.2, .~. + HabClass)
NM.3b <- update(NM.2, .~. + Moss)
AICctab(NM.2, NM.3a, NM.3b)
NM.3c <- update(NM.2, .~. + HabClass*Week)
NM.3d <- update(NM.2, .~. + Moss*Week)
# HabClass not useful but Moss*Week is
AICctab(NM.2, NM.3a, NM.3b, NM.3c, NM.3d)
summary(NM.3d)
NM.3e <- update(NM.2, .~. + Moss*poly(Week,2)) # convergence issue
AICctab(NM.2, NM.3b, NM.3d, NM.3e)

rm(list=c("NM.3a","NM.3b","NM.3c","NM.3d","NM.3e"))
NM.3 <- update(NM.2, .~. + Moss*Week)


NM.4a <- update(NM.3, .~. + log(CCavgsd) + log(TreeDens))
NM.4b <- update(NM.3, .~. + poly(log(CCavgsd),2) + log(TreeDens))
NM.4c <- update(NM.3, .~. + log(CCavgsd) + poly(log(TreeDens),2))
NM.4d <- update(NM.3, .~. + poly(log(CCavgsd),2) + poly(log(TreeDens),2))
NM.4e <- update(NM.3, .~. + log(CCavgsd)*Week + log(TreeDens)*Week)
AICctab(NM.3, NM.4a, NM.4b, NM.4c, NM.4d, NM.4e)
# although still ns, the interaction with Week has the closest AIC to NM.3
NM.4f <- update(NM.3, .~. + poly(log(CCavgsd),2)*Week + poly(log(TreeDens),2)*Week)
AICctab(NM.3, NM.4a, NM.4b, NM.4c, NM.4d, NM.4e, NM.4f)

NM.4g <- update(NM.3, .~. + log(CCavgsd)*Week + log(TreeDens))
NM.4h <- update(NM.3, .~. + log(CCavgsd) + log(TreeDens)*Week)
AICctab(NM.3, NM.4a, NM.4b, NM.4c, NM.4d, NM.4e, NM.4f, NM.4g, NM.4h)
NM.4i <- update(NM.3, .~. + log(CCavgsd)*Week + poly(log(TreeDens),2))
NM.4j <- update(NM.3, .~. + log(CCavgsd) + poly(log(TreeDens),2)*Week)
AICctab(NM.3, NM.4a, NM.4c, NM.4e, NM.4h, NM.4i, NM.4j)
NM.4k <- update(NM.3, .~. + log(TreeDens)*Week)
NM.4l <- update(NM.3, .~. + log(TreeDens))
NM.4m <- update(NM.3, .~. + log(CCavgsd)*Week)
NM.4n <- update(NM.3, .~. + log(CCavgsd))
AICctab(NM.3, NM.4l, NM.4m, NM.4n, NM.4k)
# not useful!
NM.4o <- update(NM.3, .~. + TPDscaled)
NM.4p <- update(NM.3, .~. + TPDscaled*Week)
NM.4q <- update(NM.3, .~. + poly(TPDscaled,2))
NM.4r <- update(NM.3, .~. + poly(TPDscaled,2)*Week)
AICctab(NM.3, NM.4p, NM.4q, NM.4r, NM.4o)
NM.4s <- update(NM.3, .~. + TreePropDeci)
NM.4t <- update(NM.3, .~. + poly(TreePropDeci,2))
AICctab(NM.3, NM.4o, NM.4p, NM.4q, NM.4t)
NM.4u <- update(NM.3, .~. + poly(TPDscaled,3))
NM.4v <- update(NM.3, .~. + poly(TPDscaled,4))
AICctab(NM.3, NM.4o, NM.4q, NM.4u, NM.4v)
NM.4w <- update(NM.3, .~. + TreeDiv)
NM.4x <- update(NM.3, .~. + poly(TreeDiv,2))
NM.4y <- update(NM.3, .~. + TreeDiv*Week)
NM.4z <- update(NM.3, .~. + poly(TreeDiv,3))
AICctab(NM.3, NM.4w, NM.4x, NM.4y, NM.4z)
# not useful
rm(list=c("NM.4a","NM.4b","NM.4c","NM.4d","NM.4e","NM.4f","NM.4g","NM.4h","NM.4i",
          "NM.4j","NM.4k","NM.4l","NM.4m","NM.4n","NM.4o","NM.4p","NM.4q","NM.4r",
          "NM.4s","NM.4t","NM.4u","NM.4v","NM.4w","NM.4x","NM.4y","NM.4z"))
NM.4 <- NM.3


NM.5a <- update(NM.4, .~. + sqrt(UndDens))
NM.5b <- update(NM.4, .~. + poly(sqrt(UndDens),2))
NM.5c <- update(NM.4, .~. + sqrt(UndDens)*Week)
NM.5d <- update(NM.4, .~. + poly(sqrt(UndDens),2)*Week)
AICctab(NM.4, NM.5a, NM.5b, NM.5c, NM.5d)
# UndDens not important
NM.5e <- update(NM.4, .~. + TreeRich)
NM.5f <- update(NM.4, .~. + TreeRich*Week)
NM.5g <- update(NM.4, .~. + poly(TreeRich,2))
AICctab(NM.4, NM.5e, NM.5f, NM.5g)
# TreeRich not important
NM.5h <- update(NM.4, .~. + DOM) 
NM.5i <- update(NM.4, .~. + DOM*Week)
NM.5j <- update(NM.4, .~. + DOM - Moss - Moss:Week) 
NM.5k <- update(NM.4, .~. + DOM*Week - Moss - Moss:Week)
AICctab(NM.4, NM.5h, NM.5i, NM.5j, NM.5k)
# interaction DOM*Week too complex and doesn't give enough dAICc
NM.5l <- update(NM.4, .~. + DOM - Weather)
NM.5m <- update(NM.4, .~. + poly(TreeDiv,2) - Weather)
AICctab(NM.4, NM.5l, NM.5m)
rm(list=c("NM.5a","NM.5b","NM.5c","NM.5d","NM.5e","NM.5f","NM.5g","NM.5h","NM.5i",
          "NM.5j","NM.5k", "NM.5l","NM.5m"))
NM.5 <- update(NM.4, .~. + DOM - Weather)
summary(NM.5)


NM.5sim <- simulateResiduals(NM.5, n=500)
plot(NM.5sim)
hist(NM.5sim)
testDispersion(NM.5sim)
testZeroInflation(NM.5sim)


### ###

### Analysis 2: ## Building models for migration classification: Partial migrants ####


PM.0 <- glmmTMB(Abun ~ Observer + (1|Point) + (1|Days),
                data = m_data_A1_PM, family=poisson())
summary(PM.0)
PM.0sim <- simulateResiduals(PM.0)
plot(PM.0sim)
hist(PM.0sim)
testDispersion(PM.0sim)
testZeroInflation(PM.0sim)
testGeneric(PM.0sim, summary = countOnes, alternative = "greater") # 1-inflation

PM.0nb <- update(PM.0, family=nbinom2)
PM.0nbsim <- simulateResiduals(PM.0nb)
plot(PM.0nbsim)


# adding Week and testing poly
PM.1a <- update(PM.0, .~. + Week)
PM.1b <- update(PM.0, .~. + poly(Week,2))
AICctab(PM.0, PM.1a, PM.1b)

rm(list=c("PM.1a","PM.1b"))
PM.1 <- update(PM.0, .~. + Week)
anova(PM.0, PM.1)


# nuisance variables like CoD and Weather
PM.2a <- update(PM.1, .~. + CoD)
PM.2b <- update(PM.1, .~. + poly(CoD,2))
PM.2c <- update(PM.1, .~. + Weather)
PM.2d <- update(PM.1, .~. + CoD + Weather)
AICctab(PM.1, PM.2a, PM.2b, PM.2c, PM.2d)
summary(PM.2d)

rm(list=c("PM.2a", "PM.2b","PM.2c","PM.2d"))
PM.2 <- update(PM.1, .~. + CoD + Weather)


# predictors of interest
PM.3a <- update(PM.2, .~. + HabClass)
PM.3b <- update(PM.2, .~. + Moss)
AICctab(PM.2, PM.3a, PM.3b)
PM.3c <- update(PM.2, .~. + HabClass*Week)
PM.3d <- update(PM.2, .~. + Moss*Week)
# HabClass not useful but Moss*Week is
AICctab(PM.2, PM.3a, PM.3b, PM.3c, PM.3d)
summary(PM.3d)
PM.3e <- update(PM.2, .~. + Moss*poly(Week,2)) # convergence issue
AICctab(PM.2, PM.3b, PM.3d, PM.3e)

rm(list=c("PM.3a","PM.3b","PM.3c","PM.3d","PM.3e"))
PM.3 <- update(PM.2, .~. + Moss*Week)


PM.4a <- update(PM.3, .~. + log(CCavgsd) + log(TreeDens))
PM.4b <- update(PM.3, .~. + poly(log(CCavgsd),2) + log(TreeDens))
PM.4c <- update(PM.3, .~. + log(CCavgsd) + poly(log(TreeDens),2))
PM.4d <- update(PM.3, .~. + poly(log(CCavgsd),2) + poly(log(TreeDens),2))
PM.4e <- update(PM.3, .~. + log(CCavgsd)*Week + log(TreeDens)*Week)
AICctab(PM.3, PM.4a, PM.4b, PM.4c, PM.4d, PM.4e)
# although still ns, the interaction with Week has the closest AIC to PM.3
PM.4f <- update(PM.3, .~. + poly(log(CCavgsd),2)*Week + poly(log(TreeDens),2)*Week)
AICctab(PM.3, PM.4a, PM.4b, PM.4c, PM.4d, PM.4e, PM.4f)

PM.4g <- update(PM.3, .~. + log(CCavgsd)*Week + log(TreeDens))
PM.4h <- update(PM.3, .~. + log(CCavgsd) + log(TreeDens)*Week)
AICctab(PM.3, PM.4a, PM.4b, PM.4c, PM.4d, PM.4e, PM.4f, PM.4g, PM.4h)
PM.4i <- update(PM.3, .~. + log(CCavgsd)*Week + poly(log(TreeDens),2))
PM.4j <- update(PM.3, .~. + log(CCavgsd) + poly(log(TreeDens),2)*Week)
AICctab(PM.3, PM.4a, PM.4c, PM.4e, PM.4h, PM.4i, PM.4j)
PM.4k <- update(PM.3, .~. + log(TreeDens)*Week)
PM.4l <- update(PM.3, .~. + log(TreeDens))
PM.4m <- update(PM.3, .~. + log(CCavgsd)*Week)
PM.4n <- update(PM.3, .~. + log(CCavgsd))
AICctab(PM.3, PM.4l, PM.4m, PM.4n, PM.4k)
# not useful!
PM.4o <- update(PM.3, .~. + TPDscaled)
PM.4p <- update(PM.3, .~. + TPDscaled*Week)
PM.4q <- update(PM.3, .~. + poly(TPDscaled,2))
PM.4r <- update(PM.3, .~. + poly(TPDscaled,2)*Week)
AICctab(PM.3, PM.4p, PM.4q, PM.4r, PM.4o)
PM.4s <- update(PM.3, .~. + TreePropDeci)
PM.4t <- update(PM.3, .~. + poly(TreePropDeci,2))
AICctab(PM.3, PM.4o, PM.4p, PM.4q, PM.4t)
PM.4u <- update(PM.3, .~. + poly(TPDscaled,3))
PM.4v <- update(PM.3, .~. + poly(TPDscaled,4))
AICctab(PM.3, PM.4o, PM.4q, PM.4u, PM.4v)
PM.4w <- update(PM.3, .~. + TreeDiv)
PM.4x <- update(PM.3, .~. + poly(TreeDiv,2))
PM.4y <- update(PM.3, .~. + TreeDiv*Week)
PM.4z <- update(PM.3, .~. + poly(TreeDiv,3))
AICctab(PM.3, PM.4w, PM.4x, PM.4y, PM.4z)
# not useful
rm(list=c("PM.4a","PM.4b","PM.4c","PM.4d","PM.4e","PM.4f","PM.4g","PM.4h","PM.4i",
          "PM.4j","PM.4k","PM.4l","PM.4m","PM.4n","PM.4o","PM.4p","PM.4q","PM.4r",
          "PM.4s","PM.4t","PM.4u","PM.4v","PM.4w","PM.4x","PM.4y","PM.4z"))
PM.4 <- PM.3


PM.5a <- update(PM.4, .~. + sqrt(UndDens))
PM.5b <- update(PM.4, .~. + poly(sqrt(UndDens),2))
PM.5c <- update(PM.4, .~. + sqrt(UndDens)*Week)
PM.5d <- update(PM.4, .~. + poly(sqrt(UndDens),2)*Week)
AICctab(PM.4, PM.5a, PM.5b, PM.5c, PM.5d)
# UndDens not important
PM.5e <- update(PM.4, .~. + TreeRich)
PM.5f <- update(PM.4, .~. + TreeRich*Week)
PM.5g <- update(PM.4, .~. + poly(TreeRich,2))
AICctab(PM.4, PM.5e, PM.5f, PM.5g)
# TreeRich not important
PM.5h <- update(PM.4, .~. + DOM) 
PM.5i <- update(PM.4, .~. + DOM*Week)
PM.5j <- update(PM.4, .~. + DOM - Moss - Moss:Week) 
PM.5k <- update(PM.4, .~. + DOM*Week - Moss - Moss:Week)
AICctab(PM.4, PM.5h, PM.5i, PM.5j, PM.5k)
# interaction DOM*Week too complex and doesn't give enough dAICc
PM.5l <- update(PM.4, .~. + DOM - Weather)
PM.5m <- update(PM.4, .~. + poly(TreeDiv,2) - Weather)
AICctab(PM.4, PM.5l, PM.5m)
rm(list=c("PM.5a","PM.5b","PM.5c","PM.5d","PM.5e","PM.5f","PM.5g","PM.5h","PM.5i",
          "PM.5j","PM.5k"))
PM.5 <- update(PM.4, .~. + DOM - Weather)
summary(PM.5)


PM.5sim <- simulateResiduals(PM.5, n=500)
plot(PM.5sim)
hist(PM.5sim)
testDispersion(PM.5sim)
testZeroInflation(PM.5sim)

### ###

### Analysis 2: ## Building models for migration classification: Complete migrants ####


CM.0 <- glmmTMB(Abun ~ Observer + (1|Point) + (1|Days),
                data = m_data_A1_CM, family=poisson())
summary(CM.0)
CM.0sim <- simulateResiduals(CM.0)
plot(CM.0sim)
hist(CM.0sim)
testDispersion(CM.0sim)
testZeroInflation(CM.0sim)
testGeneric(CM.0sim, summary = countOnes, alternative = "greater") # 1-inflation


# adding Week and testing poly
CM.1a <- update(CM.0, .~. + Week)
CM.1b <- update(CM.0, .~. + poly(Week,2))
AICctab(CM.0, CM.1a, CM.1b)

rm(list=c("CM.1a","CM.1b"))
CM.1 <- update(CM.0, .~. + Week)
anova(CM.0, CM.1)


# nuisance variables like CoD and Weather
CM.2a <- update(CM.1, .~. + CoD)
CM.2b <- update(CM.1, .~. + poly(CoD,2))
CM.2c <- update(CM.1, .~. + Weather)
CM.2d <- update(CM.1, .~. + CoD + Weather)
AICctab(CM.1, CM.2a, CM.2b, CM.2c, CM.2d)
summary(CM.2d)

rm(list=c("CM.2a", "CM.2b","CM.2c","CM.2d"))
CM.2 <- update(CM.1, .~. + CoD + Weather)


# predictors of interest
CM.3a <- update(CM.2, .~. + HabClass)
CM.3b <- update(CM.2, .~. + Moss)
AICctab(CM.2, CM.3a, CM.3b)
CM.3c <- update(CM.2, .~. + HabClass*Week)
CM.3d <- update(CM.2, .~. + Moss*Week)
# HabClass not useful but Moss*Week is
AICctab(CM.2, CM.3a, CM.3b, CM.3c, CM.3d)
summary(CM.3d)
CM.3e <- update(CM.2, .~. + Moss*poly(Week,2)) # convergence issue
AICctab(CM.2, CM.3b, CM.3d, CM.3e)

rm(list=c("CM.3a","CM.3b","CM.3c","CM.3d","CM.3e"))
CM.3 <- update(CM.2, .~. + Moss*Week)


CM.4a <- update(CM.3, .~. + log(CCavgsd) + log(TreeDens))
CM.4b <- update(CM.3, .~. + poly(log(CCavgsd),2) + log(TreeDens))
CM.4c <- update(CM.3, .~. + log(CCavgsd) + poly(log(TreeDens),2))
CM.4d <- update(CM.3, .~. + poly(log(CCavgsd),2) + poly(log(TreeDens),2))
CM.4e <- update(CM.3, .~. + log(CCavgsd)*Week + log(TreeDens)*Week)
AICctab(CM.3, CM.4a, CM.4b, CM.4c, CM.4d, CM.4e)
# although still ns, the interaction with Week has the closest AIC to CM.3
CM.4f <- update(CM.3, .~. + poly(log(CCavgsd),2)*Week + poly(log(TreeDens),2)*Week)
AICctab(CM.3, CM.4a, CM.4b, CM.4c, CM.4d, CM.4e, CM.4f)

CM.4g <- update(CM.3, .~. + log(CCavgsd)*Week + log(TreeDens))
CM.4h <- update(CM.3, .~. + log(CCavgsd) + log(TreeDens)*Week)
AICctab(CM.3, CM.4a, CM.4b, CM.4c, CM.4d, CM.4e, CM.4f, CM.4g, CM.4h)
CM.4i <- update(CM.3, .~. + log(CCavgsd)*Week + poly(log(TreeDens),2))
CM.4j <- update(CM.3, .~. + log(CCavgsd) + poly(log(TreeDens),2)*Week)
AICctab(CM.3, CM.4a, CM.4c, CM.4e, CM.4h, CM.4i, CM.4j)
CM.4k <- update(CM.3, .~. + log(TreeDens)*Week)
CM.4l <- update(CM.3, .~. + log(TreeDens))
CM.4m <- update(CM.3, .~. + log(CCavgsd)*Week)
CM.4n <- update(CM.3, .~. + log(CCavgsd))
AICctab(CM.3, CM.4l, CM.4m, CM.4n, CM.4k)
# not useful!
CM.4o <- update(CM.3, .~. + TPDscaled)
CM.4p <- update(CM.3, .~. + TPDscaled*Week)
CM.4q <- update(CM.3, .~. + poly(TPDscaled,2))
CM.4r <- update(CM.3, .~. + poly(TPDscaled,2)*Week)
AICctab(CM.3, CM.4p, CM.4q, CM.4r, CM.4o)
CM.4s <- update(CM.3, .~. + TreePropDeci)
CM.4t <- update(CM.3, .~. + poly(TreePropDeci,2))
AICctab(CM.3, CM.4o, CM.4p, CM.4q, CM.4t)
CM.4u <- update(CM.3, .~. + poly(TPDscaled,3))
CM.4v <- update(CM.3, .~. + poly(TPDscaled,4))
AICctab(CM.3, CM.4o, CM.4q, CM.4u, CM.4v)
CM.4w <- update(CM.3, .~. + TreeDiv)
CM.4x <- update(CM.3, .~. + poly(TreeDiv,2))
CM.4y <- update(CM.3, .~. + TreeDiv*Week)
CM.4z <- update(CM.3, .~. + poly(TreeDiv,3))
AICctab(CM.3, CM.4w, CM.4x, CM.4y, CM.4z)
# not useful
rm(list=c("CM.4a","CM.4b","CM.4c","CM.4d","CM.4e","CM.4f","CM.4g","CM.4h","CM.4i",
          "CM.4j","CM.4k","CM.4l","CM.4m","CM.4n","CM.4o","CM.4p","CM.4q","CM.4r",
          "CM.4s","CM.4t","CM.4u","CM.4v","CM.4w","CM.4x","CM.4y","CM.4z"))
CM.4 <- CM.3


CM.5a <- update(CM.4, .~. + sqrt(UndDens))
CM.5b <- update(CM.4, .~. + poly(sqrt(UndDens),2))
CM.5c <- update(CM.4, .~. + sqrt(UndDens)*Week)
CM.5d <- update(CM.4, .~. + poly(sqrt(UndDens),2)*Week)
AICctab(CM.4, CM.5a, CM.5b, CM.5c, CM.5d)
# UndDens not important
CM.5e <- update(CM.4, .~. + TreeRich)
CM.5f <- update(CM.4, .~. + TreeRich*Week)
CM.5g <- update(CM.4, .~. + poly(TreeRich,2))
AICctab(CM.4, CM.5e, CM.5f, CM.5g)
# TreeRich not important
CM.5h <- update(CM.4, .~. + DOM) 
CM.5i <- update(CM.4, .~. + DOM*Week)
CM.5j <- update(CM.4, .~. + DOM - Moss - Moss:Week) 
CM.5k <- update(CM.4, .~. + DOM*Week - Moss - Moss:Week)
AICctab(CM.4, CM.5h, CM.5i, CM.5j, CM.5k)
# interaction DOM*Week too complex and doesn't give enough dAICc
CM.5l <- update(CM.4, .~. + DOM - Weather)
CM.5m <- update(CM.4, .~. + poly(TreeDiv,2) - Weather)
AICctab(CM.4, CM.5l, CM.5m)
rm(list=c("CM.5a","CM.5b","CM.5c","CM.5d","CM.5e","CM.5f","CM.5g","CM.5h","CM.5i",
          "CM.5j","CM.5k"))
CM.5 <- update(CM.4, .~. + DOM - Weather)
summary(CM.5)


CM.5sim <- simulateResiduals(CM.5, n=500)
plot(CM.5sim)
hist(CM.5sim)
testDispersion(CM.5sim)
testZeroInflation(CM.5sim)


### ###

### Analysis 3: ## Building models for species ####

plot(m_data_A1_EGT$Abun ~ m_data_A1_EGT$Week)
plot(m_data_A1_EBT$Abun ~ m_data_A1_EBT$Week)
plot(m_data_A1_Nut$Abun ~ m_data_A1_Nut$Week)
plot(m_data_A1_Phy$Abun ~ m_data_A1_Phy$Week)
plot(m_data_A1_Rob$Abun ~ m_data_A1_Rob$Week)
plot(m_data_A1_Cer$Abun ~ m_data_A1_Cer$Week)


### Analysis 3: ## Building models for species: European Great Tit ####

# tried poisson, nb and zi for EGT, EBT and Nut. none fit well, all seem to be 1-inflated

rm(list=c("EGT.0","EGT.0sim","EGT.x","EGT.xsim","EGT.xgp","EGT.xgpsim","EGT.xcmp","EGT.xcmpsim",
          "EGT.xnb","EGT.xnbsim","EGT.xnb1","EGT.xnb1sim","EGT.xqp","EBT.0","EBT.0sim",
          "Nut.0","Phy.0","Rob.0","Rob.x","Nut.0sim","Phy.0sim","Rob.0sim","Rob.xsim"))

EGT.0 <- glmmTMB(Abun ~ Observer + (1|Point) + (1|Days),
                 data = m_data_A1_EGT, family=poisson())
summary(EGT.0)
EGT.0sim <- simulateResiduals(EGT.0, n=1000)
plot(EGT.0sim)
hist(EGT.0sim)
testDispersion(EGT.0sim)
testZeroInflation(EGT.0sim)

countOnes <- function(x) sum(x == 1)  # testing for number of 1s
testGeneric(EGT.0sim, summary = countOnes, alternative = "greater") # 1-inflation
countTwos <- function(x) sum(x == 2)  # testing for number of 2s
testGeneric(EGT.0sim, summary = countTwos, alternative = "greater") # 2-inflation
countThrees <- function(x) sum(x == 3)  # testing for number of 3s
testGeneric(EGT.0sim, summary = countThrees, alternative = "greater") # 3-inflation


EGT.x <- update(EGT.0, .~. + HabClass + log(CCavgsd)*Week + log(TreeDens))
summary(EGT.x)
EGT.xsim <- simulateResiduals(EGT.x)
plot(EGT.xsim)
testGeneric(EGT.xsim, summary = countOnes, alternative = "greater") # 1-inflation

EGT.xgp <- update(EGT.0, family=genpois)
summary(EGT.xcmp)
EGT.xgpsim <- simulateResiduals(EGT.xgp)
plot(EGT.xgpsim)
testGeneric(EGT.xgpsim, summary = countOnes, alternative = "greater") # 1-inflation


EGT.xcmp <- update(EGT.0, family=compois)
summary(EGT.xcmp)
EGT.xcmpsim <- simulateResiduals(EGT.xcmp)
plot(EGT.xcmpsim)
testGeneric(EGT.xcmpsim, summary = countOnes, alternative = "greater") # 1-inflation

# still one-inflated



EGT.xnb <- update(EGT.0, family=nbinom2)
summary(EGT.xnb)
EGT.xnbsim <- simulateResiduals(EGT.xnb)
plot(EGT.xnbsim)
testGeneric(EGT.xnbsim, summary = countOnes, alternative = "greater") # 1-inflation

EGT.xnb1 <- update(EGT.0, family=nbinom1)
summary(EGT.xnb1)
EGT.xnb1sim <- simulateResiduals(EGT.xnb1)
plot(EGT.xnb1sim)
testGeneric(EGT.xnb1sim, summary = countOnes, alternative = "greater") # 1-inflation



EGT.xqp <- gamm(Abun ~ Observer + HabClass + log(CCavgsd)*Week + log(TreeDens),
                random = list(Point=~1, Days=~1),
                data = m_data_A1_EGT, family=quasipoisson())
summary(EGT.xqp$gam)
plot(EGT.xqp$gam)


EBT.0 <- glmmTMB(Abun ~ Observer + (1|Point) + (1|Days),
                 data = m_data_A1_EBT, family=poisson())
summary(EBT.0)
EBT.0sim <- simulateResiduals(EBT.0, n=1000)
plot(EBT.0sim)
hist(EBT.0sim)
testDispersion(EBT.0sim)
testZeroInflation(EBT.0sim)
testGeneric(EBT.0sim, summary = countOnes, alternative = "greater") # 1-inflation


Nut.0 <- glmmTMB(Abun ~ Observer + (1|Point),
                 data = m_data_A1_Nut, family=poisson())
summary(Nut.0)
Nut.0sim <- simulateResiduals(Nut.0, n=1000)
plot(Nut.0sim)
hist(Nut.0sim)
testDispersion(Nut.0sim)
testZeroInflation(Nut.0sim)
testGeneric(Nut.0sim, summary = countOnes, alternative = "greater") # 1-inflation


Phy.0 <- glmmTMB(Abun ~ Observer + (1|Point),
                 data = m_data_A1_Phy, family=poisson())
summary(Phy.0)
Phy.0sim <- simulateResiduals(Phy.0, n=1000)
plot(Phy.0sim)
hist(Phy.0sim)
testDispersion(Phy.0sim)
testZeroInflation(Phy.0sim)
testGeneric(Phy.0sim, summary = countOnes, alternative = "greater") # 1-inflation


Rob.0 <- glmmTMB(Abun ~ Observer + (1|Point),
                 data = m_data_A1_Rob, family=poisson())
summary(Rob.0)
Rob.0sim <- simulateResiduals(Rob.0, n=1000)
plot(Rob.0sim)
hist(Rob.0sim)
testDispersion(Rob.0sim)
testZeroInflation(Rob.0sim)
testGeneric(Rob.0sim, summary = countOnes, alternative = "greater") # 1-inflation

Rob.x <- update(Rob.0, .~. + HabClass + log(CCavgsd)*Week + log(TreeDens))
summary(Rob.x)
Rob.xsim <- simulateResiduals(Rob.x)
plot(Rob.xsim)
testGeneric(Rob.xsim, summary = countOnes, alternative = "greater") # 1-inflation



### ###

