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


all0 <- glmer(Point_Abun ~ Observer + (1|Point) + (1|Days),
               data = m_all, family = poisson())
summary(all0)


# adding Week and testing for polynomial 
all1a <- update(all0, . ~ . + Week)
all1b <- update(all0, . ~ . + poly(Week, 2))
all1c <- update(all0, . ~ . + poly(Week, 3))
anova(all1a, all1b) 
AICctab(all1a, all1b, all1c) # returns best model on top and delta calculated from that
# 2nd and 3rd not better than linear
rm(list = c("all1a", "all1b", "all1c"))

all1 <- update(all0, .~. + Week)
summary(all1)


# testing nuisance variables like time of count during day 
all2a <- update(all1, . ~ . + CoD)
all2b <- update(all1, . ~ . + poly(CoD, 2))
anova(all1, all2a) 
AICctab(all1, all2a, all2b)
# neither linear nor poly significant
all2c <- update(all2b, . ~ . + Weather)
AICctab(all1, all2b, all2c) 
# Weather not significant

# finally moving to predictors of interest
all2d <- update(all1, .~. + HabClass)
all2e <- update(all1, .~. + Moss)
AICctab(all1, all2d, all2e)
anova(all1, all2d)
summary(all2d)
# summary function shows that Interior is marginally significant (0.04) but other levels
# are not. AIC is not different at all So will have to ignore.
# summary() p value cannot be relied on.
summary(all2e)
anova(all1, all2e) # Moss ns too

all2f <- update(all1, .~. + scaleCC)
all2g <- update(all1, .~. + log(CH) + log(TDens))
AICctab(all1, all2f, all2g) # as expected, sd+TDens better than simple CC
summary(all2g) 
# dAIC 2.5 so add
all2h <- update(all2g, .~. + log(CH):log(TDens))
summary(all2h) 
# drop interaction 
drop1(all2g) # drop main effect of het.
all2i <- update(all1, .~. + poly(log(CH),2) + poly(log(TDens),2))
all2j <- update(all1, .~. + poly(log(CH),2) + log(TDens))
all2k <- update(all1, .~. + log(CH) + poly(log(TDens),2))
AICctab(all2g, all2i, all2j, all2k)
# two poly is worst. poly of TDens is technically better than poly of CH
# but both worse than simple linear

rm(list=c("all2a","all2b","all2c","all2d","all2e","all2f",
          "all2g","all2h","all2i","all2j","all2k"))
all2 <- update(all1, .~. + log(CH) + log(TDens))
summary(all2)


# interaction of predictors with Week (of pertinence to question)
all3a <- update(all2, .~. + Week:log(CH))
all3b <- update(all2, .~. + Week:log(TDens))
all3c <- update(all2, .~. + Week:(log(CH)+log(TDens)))
AICctab(all2, all3a, all3b, all3c) 
# TDens seems more important in general (3b over 3a) but both recommended by AIC

rm(list=c("all3a", "all3b", "all3c"))
all3 <- update(all2, .~. + Week:(log(CH)+log(TDens)))
summary(all3)


# testing TPD and perhaps dropping others for this
all4a <- update(all3, .~. + scaleTPD)
all4b <- update(all3, .~. + poly(scaleTPD,2))
anova(all4a, all4b) # poly not preferred
drop1(all4a)
all4c <- update(all3, .~. + scaleTPD + Week:scaleTPD)
all4d <- update(all1, .~. + scaleTPD + Week:scaleTPD) # replacing other pred.
all4e <- update(all1, .~. + poly(scaleTPD,2))
AICctab(all1, all2, all3, all4a, all4c, all4d, all4e)
anova(all1, all3)
# TPD not useful
# trying TDiv
all4f <- update(all3, .~. + TDiv) # failed to converge 0.0112 (tot 0.002)
all4g <- update(all3, .~. + poly(TDiv,2)) # failed to converge 0.0248 (tot 0.002)
AICctab(all3, all4f, all4g)
# TDiv not useful

# trying UDens
all4h <- update(all3, .~. + sqrt(UDens))
all4i <- update(all3, .~. + poly(sqrt(UDens),2))
all4j <- update(all1, .~. + poly(sqrt(UDens),2))
all4k <- update(all1, .~. + Week*poly(sqrt(UDens),2))
AICctab(all3, all4h, all4i, all4j, all4k)
anova(all3, all4i)
# UDens not useful
# trying DOM
all4l <- update(all3, .~. + DOM)
anova(all3, all4l)
all4m <- update(all3, .~. + DOM + Week:DOM)
AICctab(all3, all4l, all4m)
# DOM main effect and interaction with Week


rm(list=c("all4a","all4b","all4c","all4d","all4e","all4f","all4g",
          "all4h","all4i","all4j","all4k","all4l","all4m"))
all4 <- update(all3, .~. + DOM + Week:DOM)


summary(all4) # failure to converge 0.0117 (tot 0.002)
hist(resid(all4)) # very Gaussian(!), except for those two outliers
all4sim <- simulateResiduals(all4, plot=T, n=1000) # beautiful :)
plot(all4sim)
testDispersion(all4sim) # 1.2055 but p value sig because of sample size

hist(all4sim) # uniform as should be
plotResiduals(all4sim, form=mdata$Week)
plotResiduals(all4sim, form=mdata$CH)
plotResiduals(all4sim, form=mdata$TDens)
plotResiduals(all4sim, form=mdata$DOM)
# all good

# exploring convergence issue 
all4opt <- allFit(all4)
summary(all4opt)$which.OK # only 1/7 optimisers failed
summary(all4opt)$llik # all pretty similar values so I think bobyqa good choice
summary(all4opt)$fixef
# use bobyqa


all5 <- update(all4, control=glmerControl(optimizer = "bobyqa"))
summary(all5)

all5sim <- simulateResiduals(all5, n=1000)
plot(all5sim)
testDispersion(all5sim)
hist(all5sim)


# further selection

all6a <- update(all5, .~. - Week:(log(CH)+log(TDens)))
all6b <- update(all5, .~. - Week:(log(CH)+log(TDens)) - log(CH))
all6c <- update(all5, .~. - Week:(log(CH)+log(TDens)) - log(TDens))
all6d <- update(all5, .~. - Week*(log(CH)+log(TDens)))

AICctab(all1, all2, all3, all4, all5, all6a, all6b, all6c, all6d,
        weights=T, base=T, logLik=T)
rm(list=c("all6a","all6b","all6c","all6d"))

all6 <- update(all5, .~. - Week:(log(CH)+log(TDens)))
AICctab(all1, all2, all3, all4, all5, all6, weights=T, base=T, logLik=T)
summary(all6)

# interactions between predictors. CH:DOM sig but mostly for Moss (only 2 Points)
# all6e <- update(all6, .~. + DOM*log(CH))
# all6f <- update(all6, .~. + DOM*log(TDens))
# all6g <- update(all6, .~. + DOM*(log(CH)+log(TDens)))
# AICctab(all6, all6e, all6f, all6g)
# 
# all6h <- update(all6, .~. + DOM*log(CH) - log(TDens))
# AICctab(all6, all6e, all6f, all6g, all6h)
# 
# summary(all6e)
# plot(allEffects(all6e))
# plot(allEffects(all6h))
# rm(list=c("all6e","all6f","all6g","all6h"))



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
# rm(list=c("inv0min","inv1min","inv2min","inv3min","inv3minsim",
#           "inv4min","inv4min2","inv4minNB","inv4minnorm","inv4minnormsim",
#           "inv4minsim","invNB0","invNB1","invTMB.0","invTMB.1","invTMB.2",
#           "invTMB.3","invTMB.4","invTMB.4NB","invTMB.4NBsim","invTMB.4poi",
#           "invTMB.4poisim","invTMB.4sim"))


# finally decided on a mixed model with Poisson error using glmmTMB with no
# zero inflation. singular fit issues resolved by using glmmTMB over glmer, and estimated 
# parameters in both methods were pretty much the same.


inv0 <- glmmTMB(GuildAbun ~ Observer + (1|Point),
                 data = m_data_A1_inv, family=poisson())
summary(inv0)
inv0sim <- simulateResiduals(inv0)
plot(inv0sim)
hist(inv0sim)


# adding Week and testing poly
inv1a <- update(inv0, .~. + Week)
inv1b <- update(inv0, .~. + poly(Week,2))
AICctab(inv0, inv1a, inv1b)
# main effect not significant
# nuisance variables like CoD and Weather
inv1c <- update(inv0, .~. + CoD)
inv1d <- update(inv0, .~. + poly(CoD,2))
inv1e <- update(inv0, .~. + Weather)
AICctab(inv0, inv1c, inv1d, inv1e)
# not significant

# predictors of interest
inv1f <- update(inv0, .~. + HabClass)
inv1g <- update(inv0, .~. + Moss)
AICctab(inv0, inv1f, inv1g)
inv1h <- update(inv0, .~. + HabClass*Week)
AICctab(inv0, inv1f, inv1g, inv1h)
# habclass not useful at all Moss is.
inv1i <- update(inv0, .~. + Moss*Week)
AICctab(inv0, inv1g, inv1i)
summary(inv1i)
# interaction has 1.4 less AIC on 2 df, so will not include, 
# even though it is pertinent to question. 
inv1j <- update(inv0, .~. + Moss*poly(Week,2))
AICctab(inv0, inv1g, inv1i, inv1j)

rm(list=c("inv1a","inv1b","inv1c","inv1d","inv1e","inv1f","inv1g","inv1h","inv1i","inv1j"))
inv1 <- update(inv0, .~. + Moss)
anova(inv0, inv1)


inv2a <- update(inv1, .~. + log(CH) + log(TDens))
inv2b <- update(inv1, .~. + poly(log(CH),2) + log(TDens))
inv2c <- update(inv1, .~. + log(CH) + poly(log(TDens),2))
inv2d <- update(inv1, .~. + poly(log(CH),2) + poly(log(TDens),2))
inv2e <- update(inv1, .~. + log(CH)*Week + log(TDens)*Week)
AICctab(inv1, inv2a, inv2b, inv2c, inv2d, inv2e)
# although still ns, the interaction with Week has the closest AIC to inv1
inv2f <- update(inv1, .~. + poly(log(CH),2)*Week + poly(log(TDens),2)*Week)
AICctab(inv1, inv2a, inv2b, inv2c, inv2d, inv2e, inv2f)

inv2g <- update(inv1, .~. + log(CH)*Week + log(TDens))
inv2h <- update(inv1, .~. + log(CH) + log(TDens)*Week)
AICctab(inv1, inv2a, inv2b, inv2c, inv2d, inv2e, inv2f, inv2g, inv2h)
inv2i <- update(inv1, .~. + log(CH)*Week + poly(log(TDens),2))
inv2j <- update(inv1, .~. + log(CH) + poly(log(TDens),2)*Week)
AICctab(inv1, inv2a, inv2c, inv2e, inv2h, inv2i, inv2j)
inv2k <- update(inv1, .~. + log(TDens)*Week)
inv2l <- update(inv1, .~. + log(TDens))
inv2m <- update(inv1, .~. + log(CH)*Week)
inv2n <- update(inv1, .~. + log(CH))
AICctab(inv1, inv2l, inv2m, inv2n, inv2k)
# not useful!

inv2o <- update(inv1, .~. + scaleTPD)
inv2p <- update(inv1, .~. + scaleTPD*Week)
inv2q <- update(inv1, .~. + poly(scaleTPD,2))
inv2r <- update(inv1, .~. + poly(scaleTPD,2)*Week)
AICctab(inv1, inv2p, inv2q, inv2r, inv2o)
inv2s <- update(inv1, .~. + TPD)
inv2t <- update(inv1, .~. + poly(TPD,2))
AICctab(inv1, inv2o, inv2p, inv2q, inv2t)
inv2u <- update(inv1, .~. + poly(scaleTPD,3))
inv2v <- update(inv1, .~. + poly(scaleTPD,4))
AICctab(inv1, inv2o, inv2q, inv2u, inv2v)
# only linear main effect of TPD
rm(list=c("inv2a","inv2b","inv2c","inv2d","inv2e","inv2f","inv2g","inv2h","inv2i",
          "inv2j","inv2k","inv2l","inv2m","inv2n","inv2o","inv2p","inv2q","inv2r",
          "inv2s","inv2t","inv2u","inv2v"))
inv2 <- update(inv1, .~. + scaleTPD)


inv3a <- update(inv2, .~. - scaleTPD + TDiv)
inv3b <- update(inv2, .~. + TDiv)
inv3c <- update(inv2, .~. + poly(TDiv,2))
inv3d <- update(inv2, .~. + TDiv*Week)
AICctab(inv2, inv3a, inv3b, inv3c, inv3d)
# TDiv not important
inv3e <- update(inv2, .~. + sqrt(UDens))
inv3f <- update(inv2, .~. + poly(sqrt(UDens),2))
inv3g <- update(inv2, .~. + sqrt(UDens)*Week)
inv3h <- update(inv2, .~. + poly(sqrt(UDens),2)*Week)
AICctab(inv2, inv3e, inv3f, inv3g, inv3h)
# UDens not important
inv3i <- update(inv2, .~. + TR)
inv3j <- update(inv2, .~. + TR*Week)
AICctab(inv2, inv3i, inv3j)
inv3k <- update(inv2, .~. + poly(TR,2))
AICctab(inv2, inv3i, inv3k)
# interesting how abundance increases with TR but decreases with TPD (and poly ns)

rm(list=c("inv3a","inv3b","inv3c","inv3d","inv3e","inv3f","inv3g","inv3h","inv3i",
          "inv3j","inv3k"))
inv3 <- update(inv2, .~. + TR)
summary(inv3)


# finally, testing DOM and comparing with Moss
inv4a <- update(inv3, .~. + DOM) 
inv4b <- update(inv3, .~. + DOM*Week)
inv4c <- update(inv3, .~. + DOM - Moss) 
inv4d <- update(inv3, .~. + DOM*Week - Moss)
AICctab(inv3, inv4a, inv4b, inv4c, inv4d) # 4c better only 0.1 dAICc on 1df
summary(inv4c)
summary(inv4a)
inv4e <- update(inv3, .~. + DOM - TR)
AICctab(inv3, inv4a, inv4c, inv4e)
inv4f <- update(inv3, .~. + DOM - TR - Moss)
AICctab(inv3, inv4a, inv4c, inv4e, inv4f)
# effect of TR prob. caused by STRONG effect of Carex. 
# All Carex Points have low TR
rm(list=c("inv4a","inv4b","inv4c","inv4d","inv4e","inv4f"))
inv4 <- update(inv3, .~. + DOM - TR - Moss)
summary(inv4)


inv4sim <- simulateResiduals(inv4, n=500)
plot(inv4sim)
hist(inv4sim)
testDispersion(inv4sim)
testZeroInflation(inv4sim)



# glmer seems to be working fine with this final model too. The problem earlier might have
# been due to Moss and DOM together, because trying a glmer with inv4 + Moss gives 
# singular fit again. Regardless, glmmTMB performs better overall


### ###

### Guilds: omnivore-feeding ####

plot(m_data_A1_omn$GuildAbun ~ m_data_A1_omn$Week)
# omni seem to show much clearer pattern (increase). migrants

# all the previous models I tried and got to dead ends with (in CodeBin):
# rm(list=c("omnGLM.0","omnGLM.1","omnGLM.2","omnGLM.3","omnGLM.4","omnGLM.4sim"))



omn0 <- glmmTMB(GuildAbun ~ Observer + (1|Point), data=m_data_A1_omn, family=poisson)
omn0sim <- simulateResiduals(omn0)
plot(omn0sim)
# KS test, dispersion test, outlier test all sig
omn0NB <- glmmTMB(GuildAbun ~ Observer + (1|Point), data=m_data_A1_omn, family = nbinom2)
omn0NBsim <- simulateResiduals(omn0NB)
plot(omn0NBsim)


# adding Week and testing poly
omn1a <- update(omn0NB, .~. + Week)
omn1b <- update(omn0NB, .~. + poly(Week,2))
AICctab(omn0, omn0NB, omn1a, omn1b)
# linear effect of Week very good
rm(list=c("omn1a","omn1b"))
omn1 <- update(omn0NB, .~. + Week)


# nuisance variables like CoD and Weather
omn2a <- update(omn1, .~. + CoD)
omn2b <- update(omn1, .~. + poly(CoD,2))
omn2c <- update(omn1, .~. + Weather)
AICctab(omn1, omn2a, omn2b, omn2c)
# not significant
# predictors of interest
omn2d <- update(omn1, .~. + HabClass)
omn2e <- update(omn1, .~. + Moss)
AICctab(omn1, omn2d, omn2e) # habclass sig
omn2f <- update(omn1, .~. + HabClass + Moss)
AICctab(omn1, omn2d, omn2e, omn2f) # Moss ns
omn2g <- update(omn1, .~. + HabClass*Week)
omn2h <- update(omn1, .~. + Moss*Week)
AICctab(omn1, omn2d, omn2g, omn2e, omn2h)
# only main effect of HabClass
rm(list=c("omn2a","omn2b","omn2c","omn2d","omn2e","omn2f","omn2g","omn2h"))
omn2 <- update(omn1, .~. + HabClass)
anova(omn1, omn2)


omn3a <- update(omn2, .~. + log(CH) + log(TDens))
omn3b <- update(omn2, .~. + poly(log(CH),2) + log(TDens))
omn3c <- update(omn2, .~. + log(CH) + poly(log(TDens),2))
omn3d <- update(omn2, .~. + poly(log(CH),2) + poly(log(TDens),2))
omn3e <- update(omn2, .~. + log(CH)*Week + log(TDens)*Week)
AICctab(omn2, omn3a, omn3b, omn3c, omn3d, omn3e)
# no poly
omn3f <- update(omn2, .~. + log(CH)*Week + log(TDens))
omn3g <- update(omn2, .~. + log(CH) + log(TDens)*Week)
AICctab(omn2, omn3a, omn3b, omn3c, omn3d, omn3e, omn3f, omn3g)
# interaction CH:Week + TD is best
omn3h <- update(omn2, .~. + log(CH)*Week) 
omn3i <- update(omn2, .~. + log(CH))
AICctab(omn2, omn3a, omn3f, omn3h, omn3i)

rm(list=c("omn3a","omn3b","omn3c","omn3d","omn3e","omn3f","omn3g","omn3h","omn3i"))
omn3 <- update(omn2, .~. + log(CH)*Week + log(TDens))


omn4a <- update(omn3, .~. + scaleTPD - log(TDens))
omn4b <- update(omn3, .~. + scaleTPD) 
AICctab(omn3, omn4a, omn4b) # keep
omn4c <- update(omn3, .~. + scaleTPD*Week)
omn4d <- update(omn3, .~. + poly(scaleTPD,2))
omn4e <- update(omn3, .~. + poly(scaleTPD,2)*Week)
AICctab(omn3, omn4a, omn4b, omn4c, omn4d)
# not significant
omn4f <- update(omn3, .~. + TDiv)
omn4g <- update(omn3, .~. + poly(TDiv,2))
omn4h <- update(omn3, .~. + TDiv*Week)
omn4i <- update(omn3, .~. + poly(TDiv,2)*Week)
AICctab(omn3, omn4f, omn4g, omn4h, omn4i)
# not significant
omn4j <- update(omn3, .~. + scaleUDens)
omn4k <- update(omn3, .~. + poly(scaleUDens,2))
omn4l <- update(omn3, .~. + scaleUDens*Week)
omn4m <- update(omn3, .~. + poly(scaleUDens,2)*Week)
AICctab(omn3, omn4j, omn4k, omn4l, omn4m)
# not significant
omn4n <- update(omn3, .~. + TR)
omn4o <- update(omn3, .~. + TR*Week)
omn4p <- update(omn3, .~. + poly(TR,2))
omn4q <- update(omn3, .~. + poly(TR,2)*Week)
AICctab(omn3, omn4n, omn4o, omn4p, omn4q)
# not significant
omn4r <- update(omn3, .~. + DOM) 
omn4s <- update(omn3, .~. + DOM*Week)
omn4t <- update(omn3, .~. + DOM - HabClass)
AICctab(omn3, omn4r, omn4s, omn4t)
# DOM has 3.4 dAICc but 4 df ##
rm(list=c("omn4a","omn4b","omn4c","omn4d","omn4e","omn4f","omn4g","omn4h","omn4i",
          "omn4j","omn4k","omn4l","omn4m","omn4n","omn4o","omn4p","omn4q","omn4r",
          "omn4s","omn4t"))
# no selection. stick with omn3

summary(omn3)
omn3sim <- simulateResiduals(omn3) 
plot(omn3sim)
hist(omn3sim) 
testDispersion(omn3sim)
testZeroInflation(omn3sim)




# glmer is not working for this final model, but glmmTMB took time for each model



### ###


