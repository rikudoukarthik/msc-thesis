### Analysis of thesis data for manuscript ### 

### Bird abundances versus habitat variables; 
### post-breeding habitat breadth and selectivity of birds

library(tidyverse)
library(lme4)
library(bbmle)
library(DHARMa)
library(glmmTMB)
library(vegan)
# devtools::install_github("ProcessMiner/nlcor")
library(nlcor)

load("data/01_dataimport.RData")

### Checking for correlation between habitat variables ############

cor.test(habvar$CC, habvar$TDens, method = "pearson")
cor.test(habvar$CH, habvar$CC, method = "pearson")
cor.test(habvar$CH, habvar$TDens, method = "pearson")

plot(habvar$CC, habvar$TDens)
plot(habvar$CH, habvar$CC)
plot(habvar$CH, habvar$TDens)

# Strong correlation of CC and CH at high CC is a direct physical result 
# because with higher density of trees, there are less gaps overall and hence the 
# variability/heterogeneity in cover is also lower.

# CH does not correlate with TDens. 

# Hence, I can either include just CC in the models, which captures partly the 
# essences of CH and TDens, or I can include the latter two in the models 
# and exclude CC where the two variables would then capture a broader scope
# of variability than would have been possible with CC alone, but might miss
# some of the core essence of CC

cor.test(habvar$TPD, habvar$TDiv, method = "spearman")
plot(habvar$TPD, habvar$TDiv)
# should test non-linear correlation

nlcor(habvar$TPD, habvar$TDiv, plt = T)
# from the plot, it is clear that there is a hump-shaped correlation, but the function
# is unable to pick it up, probably due to the low sample size (it was working during thesis)

# I think it was Salek 2016 who said proportion of deciduous trees was a better
# predictor than tree diversity, and since the (non-linear) correlation is very
# strong here, I think it is better to include TPD and exclude TDiv.

# Using CH+TDens over using CC is not statistically significantly better
# but has 1.79 less residual deviance and only uses one extra df, so it might be
# the way to go.
# Using TDiv over TPD has, technically, lower resid. dev. and lower
# Cp, but the difference is negligible (0.25).


cor.test(habvar$TPD, habvar$TDens, method="pearson")
plot(habvar$TPD, habvar$TDens) # correlated

cor.test(habvar$UDiv, habvar$TDens, method="pearson") # close to sig.
plot(habvar$UDiv, habvar$TDens) 

# not significant
cor.test(habvar$TPD, habvar$CC, method="pearson")
cor.test(habvar$UDens, habvar$TDens, method="pearson")
cor.test(habvar$UR, habvar$TDens, method="pearson")
cor.test(habvar$UDens, habvar$CH, method="pearson")
cor.test(habvar$UDiv, habvar$CH, method="pearson")
cor.test(habvar$UR, habvar$CH, method="pearson")

### ###

### All birds ####

boxplot(m_all$Point_Abun ~ m_all$Period)
plot(m_all$Point_Abun ~ m_all$Days)


all.0 <- glmer(Point_Abun ~ Observer + (1|Point) + (1|Days),
               data = m_all, family = poisson())
summary(all.0)


# adding Week and testing for polynomial 
all.1a <- update(all.0, . ~ . + Week)
all.1b <- update(all.0, . ~ . + poly(Week, 2))
all.1c <- update(all.0, . ~ . + poly(Week, 3))
anova(all.1a, all.1b) 
AICctab(all.1a, all.1b, all.1c) # returns best model on top and delta calculated from that
# 2nd and 3rd not better than linear
rm(list = c("all.1a", "all.1b", "all.1c"))

all.1 <- update(all.0, .~. + Week)
summary(all.1)


# testing nuisance variables like time of count during day 
all.2a <- update(all.1, . ~ . + CoD)
all.2b <- update(all.1, . ~ . + poly(CoD, 2))
anova(all.1, all.2a) 
AICctab(all.1, all.2a, all.2b)
# both significant, poly more
all.2c <- update(all.2b, . ~ . + Weather)
AICctab(all.1, all.2b, all.2c) 
# Weather significant

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

all.2f <- update(all.1, .~. + scaleCC)
all.2g <- update(all.1, .~. + log(CH) + log(TDens))
AICctab(all.1, all.2f, all.2g) # as expected, sd+TDens better than simple CC
summary(all.2g) 
# dAIC 2.5 so add
all.2h <- update(all.2g, .~. + log(CH):log(TDens))
summary(all.2h) 
# drop interaction 
drop1(all.2g) # drop main effect of het.
all.2i <- update(all.1, .~. + poly(log(CH),2) + poly(log(TDens),2))
all.2j <- update(all.1, .~. + poly(log(CH),2) + log(TDens))
all.2k <- update(all.1, .~. + log(CH) + poly(log(TDens),2))
AICctab(all.2g, all.2i, all.2j, all.2k)
# two poly is worst. poly of TDens is technically better than poly of CH
# but both worse than simple linear

rm(list=c("all.2a","all.2b","all.2c","all.2d","all.2e","all.2f",
          "all.2g","all.2h","all.2i","all.2j","all.2k"))
all.2 <- update(all.1, .~. + log(CH) + log(TDens))
summary(all.2)


# interaction of predictors with Week (of pertinence to question)
all.3a <- update(all.2, .~. + Week:log(CH))
all.3b <- update(all.2, .~. + Week:log(TDens))
all.3c <- update(all.2, .~. + Week:(log(CH)+log(TDens)))
AICctab(all.2, all.3a, all.3b, all.3c) 
# TDens seems more important in general (3b over 3a) but both recommended by AIC

rm(list=c("all.3a", "all.3b", "all.3c"))
all.3 <- update(all.2, .~. + Week:(log(CH)+log(TDens)))
summary(all.3)


# testing TPD and perhaps dropping others for this
all.4a <- update(all.3, .~. + scaleTPD)
all.4b <- update(all.3, .~. + poly(scaleTPD,2))
anova(all.4a, all.4b) # poly not preferred
drop1(all.4a)
all.4c <- update(all.3, .~. + scaleTPD + Week:scaleTPD)
all.4d <- update(all.1, .~. + scaleTPD + Week:scaleTPD) # replacing other pred.
all.4e <- update(all.1, .~. + poly(scaleTPD,2))
AICctab(all.1, all.2, all.3, all.4a, all.4c, all.4d, all.4e)
anova(all.1, all.3)
# TPD not useful
# trying TDiv
all.4f <- update(all.3, .~. + TDiv) # failed to converge 0.0112 (tot 0.002)
all.4g <- update(all.3, .~. + poly(TDiv,2)) # failed to converge 0.0248 (tot 0.002)
AICctab(all.3, all.4f, all.4g)
# TDiv not useful

# trying UDens
all.4h <- update(all.3, .~. + sqrt(UDens))
all.4i <- update(all.3, .~. + poly(sqrt(UDens),2))
all.4j <- update(all.1, .~. + poly(sqrt(UDens),2))
all.4k <- update(all.1, .~. + Week*poly(sqrt(UDens),2))
AICctab(all.3, all.4h, all.4i, all.4j, all.4k)
anova(all.3, all.4i)
# UDens not useful
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
plotResiduals(all.4sim, form=mdata$CH)
plotResiduals(all.4sim, form=mdata$TDens)
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

all.6a <- update(all.5, .~. - Week:(log(CH)+log(TDens)))
all.6b <- update(all.5, .~. - Week:(log(CH)+log(TDens)) - log(CH))
all.6c <- update(all.5, .~. - Week:(log(CH)+log(TDens)) - log(TDens))
all.6d <- update(all.5, .~. - Week*(log(CH)+log(TDens)))

AICctab(all.1, all.2, all.3, all.4, all.5, all.6a, all.6b, all.6c, all.6d,
        weights=T, base=T, logLik=T)
rm(list=c("all.6a","all.6b","all.6c","all.6d"))

all.6 <- update(all.5, .~. - Week:(log(CH)+log(TDens)))
AICctab(all.1, all.2, all.3, all.4, all.5, all.6, weights=T, base=T, logLik=T)
summary(all.6)

# interactions between predictors. CH:DOM sig but mostly for Moss (only 2 Points)
# all.6e <- update(all.6, .~. + DOM*log(CH))
# all.6f <- update(all.6, .~. + DOM*log(TDens))
# all.6g <- update(all.6, .~. + DOM*(log(CH)+log(TDens)))
# AICctab(all.6, all.6e, all.6f, all.6g)
# 
# all.6h <- update(all.6, .~. + DOM*log(CH) - log(TDens))
# AICctab(all.6, all.6e, all.6f, all.6g, all.6h)
# 
# summary(all.6e)
# plot(allEffects(all.6e))
# plot(allEffects(all.6h))
# rm(list=c("all.6e","all.6f","all.6g","all.6h"))



## Matching fruiting phenology to abundance ##


allfruit.0 <- glmmTMB(Point_Abun ~ Observer + (1|Point) + (1|Days),
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

### Guilds: invertebrate-feeding ####

# invertebrate has 22.6% zeroes (94/416); omnivore 7.9%

plot(m_data_A1_inv$GuildAbun ~ m_data_A1_inv$Week)

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


inv.2a <- update(inv.1, .~. + log(CH) + log(TDens))
inv.2b <- update(inv.1, .~. + poly(log(CH),2) + log(TDens))
inv.2c <- update(inv.1, .~. + log(CH) + poly(log(TDens),2))
inv.2d <- update(inv.1, .~. + poly(log(CH),2) + poly(log(TDens),2))
inv.2e <- update(inv.1, .~. + log(CH)*Week + log(TDens)*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e)
# although still ns, the interaction with Week has the closest AIC to inv.1
inv.2f <- update(inv.1, .~. + poly(log(CH),2)*Week + poly(log(TDens),2)*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f)

inv.2g <- update(inv.1, .~. + log(CH)*Week + log(TDens))
inv.2h <- update(inv.1, .~. + log(CH) + log(TDens)*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f, inv.2g, inv.2h)
inv.2i <- update(inv.1, .~. + log(CH)*Week + poly(log(TDens),2))
inv.2j <- update(inv.1, .~. + log(CH) + poly(log(TDens),2)*Week)
AICctab(inv.1, inv.2a, inv.2c, inv.2e, inv.2h, inv.2i, inv.2j)
inv.2k <- update(inv.1, .~. + log(TDens)*Week)
inv.2l <- update(inv.1, .~. + log(TDens))
inv.2m <- update(inv.1, .~. + log(CH)*Week)
inv.2n <- update(inv.1, .~. + log(CH))
AICctab(inv.1, inv.2l, inv.2m, inv.2n, inv.2k)
# not useful!

inv.2o <- update(inv.1, .~. + scaleTPD)
inv.2p <- update(inv.1, .~. + scaleTPD*Week)
inv.2q <- update(inv.1, .~. + poly(scaleTPD,2))
inv.2r <- update(inv.1, .~. + poly(scaleTPD,2)*Week)
AICctab(inv.1, inv.2p, inv.2q, inv.2r, inv.2o)
inv.2s <- update(inv.1, .~. + TPD)
inv.2t <- update(inv.1, .~. + poly(TPD,2))
AICctab(inv.1, inv.2o, inv.2p, inv.2q, inv.2t)
inv.2u <- update(inv.1, .~. + poly(scaleTPD,3))
inv.2v <- update(inv.1, .~. + poly(scaleTPD,4))
AICctab(inv.1, inv.2o, inv.2q, inv.2u, inv.2v)
# only linear main effect of TPD
rm(list=c("inv.2a","inv.2b","inv.2c","inv.2d","inv.2e","inv.2f","inv.2g","inv.2h","inv.2i",
          "inv.2j","inv.2k","inv.2l","inv.2m","inv.2n","inv.2o","inv.2p","inv.2q","inv.2r",
          "inv.2s","inv.2t","inv.2u","inv.2v"))
inv.2 <- update(inv.1, .~. + scaleTPD)


inv.3a <- update(inv.2, .~. - scaleTPD + TDiv)
inv.3b <- update(inv.2, .~. + TDiv)
inv.3c <- update(inv.2, .~. + poly(TDiv,2))
inv.3d <- update(inv.2, .~. + TDiv*Week)
AICctab(inv.2, inv.3a, inv.3b, inv.3c, inv.3d)
# TDiv not important
inv.3e <- update(inv.2, .~. + sqrt(UDens))
inv.3f <- update(inv.2, .~. + poly(sqrt(UDens),2))
inv.3g <- update(inv.2, .~. + sqrt(UDens)*Week)
inv.3h <- update(inv.2, .~. + poly(sqrt(UDens),2)*Week)
AICctab(inv.2, inv.3e, inv.3f, inv.3g, inv.3h)
# UDens not important
inv.3i <- update(inv.2, .~. + TR)
inv.3j <- update(inv.2, .~. + TR*Week)
AICctab(inv.2, inv.3i, inv.3j)
inv.3k <- update(inv.2, .~. + poly(TR,2))
AICctab(inv.2, inv.3i, inv.3k)
# interesting how abundance increases with TR but decreases with TPD (and poly ns)

rm(list=c("inv.3a","inv.3b","inv.3c","inv.3d","inv.3e","inv.3f","inv.3g","inv.3h","inv.3i",
          "inv.3j","inv.3k"))
inv.3 <- update(inv.2, .~. + TR)
summary(inv.3)


# finally, testing DOM and comparing with Moss
inv.4a <- update(inv.3, .~. + DOM) 
inv.4b <- update(inv.3, .~. + DOM*Week)
inv.4c <- update(inv.3, .~. + DOM - Moss) 
inv.4d <- update(inv.3, .~. + DOM*Week - Moss)
AICctab(inv.3, inv.4a, inv.4b, inv.4c, inv.4d) # 4c better only 0.1 dAICc on 1df
summary(inv.4c)
summary(inv.4a)
inv.4e <- update(inv.3, .~. + DOM - TR)
AICctab(inv.3, inv.4a, inv.4c, inv.4e)
inv.4f <- update(inv.3, .~. + DOM - TR - Moss)
AICctab(inv.3, inv.4a, inv.4c, inv.4e, inv.4f)
# effect of TR prob. caused by STRONG effect of Carex. 
# All Carex Points have low TR
rm(list=c("inv.4a","inv.4b","inv.4c","inv.4d","inv.4e","inv.4f"))
inv.4 <- update(inv.3, .~. + DOM - TR - Moss)
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

### Guilds: omnivore-feeding ####

plot(m_data_A1_omn$GuildAbun ~ m_data_A1_omn$Week)
# omni seem to show much clearer pattern (increase). migrants

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


omn.3a <- update(omn.2, .~. + log(CH) + log(TDens))
omn.3b <- update(omn.2, .~. + poly(log(CH),2) + log(TDens))
omn.3c <- update(omn.2, .~. + log(CH) + poly(log(TDens),2))
omn.3d <- update(omn.2, .~. + poly(log(CH),2) + poly(log(TDens),2))
omn.3e <- update(omn.2, .~. + log(CH)*Week + log(TDens)*Week)
AICctab(omn.2, omn.3a, omn.3b, omn.3c, omn.3d, omn.3e)
# no poly
omn.3f <- update(omn.2, .~. + log(CH)*Week + log(TDens))
omn.3g <- update(omn.2, .~. + log(CH) + log(TDens)*Week)
AICctab(omn.2, omn.3a, omn.3b, omn.3c, omn.3d, omn.3e, omn.3f, omn.3g)
# interaction CH:Week + TD is best
omn.3h <- update(omn.2, .~. + log(CH)*Week) 
omn.3i <- update(omn.2, .~. + log(CH))
AICctab(omn.2, omn.3a, omn.3f, omn.3h, omn.3i)

rm(list=c("omn.3a","omn.3b","omn.3c","omn.3d","omn.3e","omn.3f","omn.3g","omn.3h","omn.3i"))
omn.3 <- update(omn.2, .~. + log(CH)*Week + log(TDens))


omn.4a <- update(omn.3, .~. + scaleTPD - log(TDens))
omn.4b <- update(omn.3, .~. + scaleTPD) 
AICctab(omn.3, omn.4a, omn.4b) # keep
omn.4c <- update(omn.3, .~. + scaleTPD*Week)
omn.4d <- update(omn.3, .~. + poly(scaleTPD,2))
omn.4e <- update(omn.3, .~. + poly(scaleTPD,2)*Week)
AICctab(omn.3, omn.4a, omn.4b, omn.4c, omn.4d)
# not significant
omn.4f <- update(omn.3, .~. + TDiv)
omn.4g <- update(omn.3, .~. + poly(TDiv,2))
omn.4h <- update(omn.3, .~. + TDiv*Week)
omn.4i <- update(omn.3, .~. + poly(TDiv,2)*Week)
AICctab(omn.3, omn.4f, omn.4g, omn.4h, omn.4i)
# not significant
omn.4j <- update(omn.3, .~. + scaleUDens)
omn.4k <- update(omn.3, .~. + poly(scaleUDens,2))
omn.4l <- update(omn.3, .~. + scaleUDens*Week)
omn.4m <- update(omn.3, .~. + poly(scaleUDens,2)*Week)
AICctab(omn.3, omn.4j, omn.4k, omn.4l, omn.4m)
# not significant
omn.4n <- update(omn.3, .~. + TR)
omn.4o <- update(omn.3, .~. + TR*Week)
omn.4p <- update(omn.3, .~. + poly(TR,2))
omn.4q <- update(omn.3, .~. + poly(TR,2)*Week)
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


