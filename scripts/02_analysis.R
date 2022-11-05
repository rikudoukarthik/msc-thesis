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
a <- update(all0, . ~ . + Week)
b <- update(all0, . ~ . + poly(Week, 2))
c <- update(all0, . ~ . + poly(Week, 3))
anova(a, b) 
AICctab(a, b, c) # returns best model on top and delta calculated from that
# 2nd and 3rd not better than linear
rm(list = c("a", "b", "c"))

all1 <- update(all0, .~. + Week)
summary(all1)


# testing nuisance variables like time of count during day 
a <- update(all1, . ~ . + CoD)
b <- update(all1, . ~ . + poly(CoD, 2))
anova(all1, a) 
AICctab(all1, a, b)
# neither linear nor poly significant
c <- update(b, . ~ . + Weather)
AICctab(all1, b, c) 
# Weather not significant

# finally moving to predictors of interest
d <- update(all1, .~. + HabClass)
e <- update(all1, .~. + Moss)
AICctab(all1, d, e)
anova(all1, d)
summary(d)
# summary function shows that Interior is marginally significant (0.04) but other levels
# are not. AIC is not different at all So will have to ignore.
# summary() p value cannot be relied on.
summary(e)
anova(all1, e) # Moss ns too

f <- update(all1, .~. + scaleCC)
g <- update(all1, .~. + log(CH) + log(TDens))
AICctab(all1, f, g) # as expected, sd+TDens better than simple CC
summary(g) 
# dAIC 2.5 so add
h <- update(g, .~. + log(CH):log(TDens))
summary(h) 
# drop interaction 
drop1(g) # drop main effect of het.
i <- update(all1, .~. + poly(log(CH),2) + poly(log(TDens),2))
j <- update(all1, .~. + poly(log(CH),2) + log(TDens))
k <- update(all1, .~. + log(CH) + poly(log(TDens),2))
AICctab(g, i, j, k)
# two poly is worst. poly of TDens is technically better than poly of CH
# but both worse than simple linear

rm(list=c("a","b","c","d","e","f",
          "g","h","i","j","k"))
all2 <- update(all1, .~. + log(CH) + log(TDens))
summary(all2)


# interaction of predictors with Week (of pertinence to question)
a <- update(all2, .~. + Week:log(CH))
b <- update(all2, .~. + Week:log(TDens))
c <- update(all2, .~. + Week:(log(CH)+log(TDens)))
AICctab(all2, a, b, c) 
# TDens seems more important in general (3b over 3a) but both recommended by AIC

rm(list=c("a", "b", "c"))
all3 <- update(all2, .~. + Week:(log(CH)+log(TDens)))
summary(all3)


# testing TPD and perhaps dropping others for this
a <- update(all3, .~. + scaleTPD)
b <- update(all3, .~. + poly(scaleTPD,2))
anova(a, b) # poly not preferred
drop1(a)
c <- update(all3, .~. + scaleTPD + Week:scaleTPD)
d <- update(all1, .~. + scaleTPD + Week:scaleTPD) # replacing other pred.
e <- update(all1, .~. + poly(scaleTPD,2))
AICctab(all1, all2, all3, a, c, d, e)
anova(all1, all3)
# TPD not useful
# trying TDiv
f <- update(all3, .~. + TDiv) # failed to converge 0.0112 (tot 0.002)
g <- update(all3, .~. + poly(TDiv,2)) # failed to converge 0.0248 (tot 0.002)
AICctab(all3, f, g)
# TDiv not useful

# trying UDens
h <- update(all3, .~. + sqrt(UDens))
i <- update(all3, .~. + poly(sqrt(UDens),2))
j <- update(all1, .~. + poly(sqrt(UDens),2))
k <- update(all1, .~. + Week*poly(sqrt(UDens),2))
AICctab(all3, h, i, j, k)
anova(all3, i)
# UDens not useful
# trying DOM
l <- update(all3, .~. + DOM)
anova(all3, l)
m <- update(all3, .~. + DOM + Week:DOM)
AICctab(all3, l, m)
# DOM main effect and interaction with Week


rm(list=c("a","b","c","d","e","f","g",
          "h","i","j","k","l","m"))
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

a <- update(all5, .~. - Week:(log(CH)+log(TDens)))
b <- update(all5, .~. - Week:(log(CH)+log(TDens)) - log(CH))
c <- update(all5, .~. - Week:(log(CH)+log(TDens)) - log(TDens))
d <- update(all5, .~. - Week*(log(CH)+log(TDens)))

AICctab(all1, all2, all3, all4, all5, a, b, c, d,
        weights=T, base=T, logLik=T)
rm(list=c("a","b","c","d"))

all6 <- update(all5, .~. - Week:(log(CH)+log(TDens)))
AICctab(all1, all2, all3, all4, all5, all6, weights=T, base=T, logLik=T)
summary(all6)

# interactions between predictors. CH:DOM sig but mostly for Moss (only 2 Points)
# e <- update(all6, .~. + DOM*log(CH))
# f <- update(all6, .~. + DOM*log(TDens))
# g <- update(all6, .~. + DOM*(log(CH)+log(TDens)))
# AICctab(all6, e, f, g)
# 
# h <- update(all6, .~. + DOM*log(CH) - log(TDens))
# AICctab(all6, e, f, g, h)
# 
# summary(e)
# plot(allEffects(e))
# plot(allEffects(h))
# rm(list=c("e","f","g","h"))



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
a <- update(inv0, .~. + Week)
b <- update(inv0, .~. + poly(Week,2))
AICctab(inv0, a, b)
# main effect not significant
# nuisance variables like CoD and Weather
c <- update(inv0, .~. + CoD)
d <- update(inv0, .~. + poly(CoD,2))
e <- update(inv0, .~. + Weather)
AICctab(inv0, c, d, e)
# not significant

# predictors of interest
f <- update(inv0, .~. + HabClass)
g <- update(inv0, .~. + Moss)
AICctab(inv0, f, g)
h <- update(inv0, .~. + HabClass*Week)
AICctab(inv0, f, g, h)
# habclass not useful at all Moss is.
i <- update(inv0, .~. + Moss*Week)
AICctab(inv0, g, i)
summary(i)
# interaction has 1.4 less AIC on 2 df, so will not include, 
# even though it is pertinent to question. 
j <- update(inv0, .~. + Moss*poly(Week,2))
AICctab(inv0, g, i, j)

rm(list=c("a","b","c","d","e","f","g","h","i","j"))
inv1 <- update(inv0, .~. + Moss)
anova(inv0, inv1)


a <- update(inv1, .~. + log(CH) + log(TDens))
b <- update(inv1, .~. + poly(log(CH),2) + log(TDens))
c <- update(inv1, .~. + log(CH) + poly(log(TDens),2))
d <- update(inv1, .~. + poly(log(CH),2) + poly(log(TDens),2))
e <- update(inv1, .~. + log(CH)*Week + log(TDens)*Week)
AICctab(inv1, a, b, c, d, e)
# although still ns, the interaction with Week has the closest AIC to inv1
f <- update(inv1, .~. + poly(log(CH),2)*Week + poly(log(TDens),2)*Week)
AICctab(inv1, a, b, c, d, e, f)

g <- update(inv1, .~. + log(CH)*Week + log(TDens))
h <- update(inv1, .~. + log(CH) + log(TDens)*Week)
AICctab(inv1, a, b, c, d, e, f, g, h)
i <- update(inv1, .~. + log(CH)*Week + poly(log(TDens),2))
j <- update(inv1, .~. + log(CH) + poly(log(TDens),2)*Week)
AICctab(inv1, a, c, e, h, i, j)
k <- update(inv1, .~. + log(TDens)*Week)
l <- update(inv1, .~. + log(TDens))
m <- update(inv1, .~. + log(CH)*Week)
n <- update(inv1, .~. + log(CH))
AICctab(inv1, l, m, n, k)
# not useful!

o <- update(inv1, .~. + scaleTPD)
p <- update(inv1, .~. + scaleTPD*Week)
q <- update(inv1, .~. + poly(scaleTPD,2))
r <- update(inv1, .~. + poly(scaleTPD,2)*Week)
AICctab(inv1, p, q, r, o)
s <- update(inv1, .~. + TPD)
t <- update(inv1, .~. + poly(TPD,2))
AICctab(inv1, o, p, q, t)
u <- update(inv1, .~. + poly(scaleTPD,3))
v <- update(inv1, .~. + poly(scaleTPD,4))
AICctab(inv1, o, q, u, v)
# only linear main effect of TPD
rm(list=c("a","b","c","d","e","f","g","h","i",
          "j","k","l","m","n","o","p","q","r",
          "s","t","u","v"))
inv2 <- update(inv1, .~. + scaleTPD)


a <- update(inv2, .~. - scaleTPD + TDiv)
b <- update(inv2, .~. + TDiv)
c <- update(inv2, .~. + poly(TDiv,2))
d <- update(inv2, .~. + TDiv*Week)
AICctab(inv2, a, b, c, d)
# TDiv not important
e <- update(inv2, .~. + sqrt(UDens))
f <- update(inv2, .~. + poly(sqrt(UDens),2))
g <- update(inv2, .~. + sqrt(UDens)*Week)
h <- update(inv2, .~. + poly(sqrt(UDens),2)*Week)
AICctab(inv2, e, f, g, h)
# UDens not important
i <- update(inv2, .~. + TR)
j <- update(inv2, .~. + TR*Week)
AICctab(inv2, i, j)
k <- update(inv2, .~. + poly(TR,2))
AICctab(inv2, i, k)
# interesting how abundance increases with TR but decreases with TPD (and poly ns)

rm(list=c("a","b","c","d","e","f","g","h","i",
          "j","k"))
inv3 <- update(inv2, .~. + TR)
summary(inv3)


# finally, testing DOM and comparing with Moss
a <- update(inv3, .~. + DOM) 
b <- update(inv3, .~. + DOM*Week)
c <- update(inv3, .~. + DOM - Moss) 
d <- update(inv3, .~. + DOM*Week - Moss)
AICctab(inv3, a, b, c, d) # 4c better only 0.1 dAICc on 1df
summary(c)
summary(a)
e <- update(inv3, .~. + DOM - TR)
AICctab(inv3, a, c, e)
f <- update(inv3, .~. + DOM - TR - Moss)
AICctab(inv3, a, c, e, f)
# effect of TR prob. caused by STRONG effect of Carex. 
# All Carex Points have low TR
rm(list=c("a","b","c","d","e","f"))
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
a <- update(omn0NB, .~. + Week)
b <- update(omn0NB, .~. + poly(Week,2))
AICctab(omn0, omn0NB, a, b)
# linear effect of Week very good
rm(list=c("a","b"))
omn1 <- update(omn0NB, .~. + Week)


# nuisance variables like CoD and Weather
a <- update(omn1, .~. + CoD)
b <- update(omn1, .~. + poly(CoD,2))
c <- update(omn1, .~. + Weather)
AICctab(omn1, a, b, c)
# not significant
# predictors of interest
d <- update(omn1, .~. + HabClass)
e <- update(omn1, .~. + Moss)
AICctab(omn1, d, e) # habclass sig
f <- update(omn1, .~. + HabClass + Moss)
AICctab(omn1, d, e, f) # Moss ns
g <- update(omn1, .~. + HabClass*Week)
h <- update(omn1, .~. + Moss*Week)
AICctab(omn1, d, g, e, h)
# only main effect of HabClass
rm(list=c("a","b","c","d","e","f","g","h"))
omn2 <- update(omn1, .~. + HabClass)
anova(omn1, omn2)


a <- update(omn2, .~. + log(CH) + log(TDens))
b <- update(omn2, .~. + poly(log(CH),2) + log(TDens))
c <- update(omn2, .~. + log(CH) + poly(log(TDens),2))
d <- update(omn2, .~. + poly(log(CH),2) + poly(log(TDens),2))
e <- update(omn2, .~. + log(CH)*Week + log(TDens)*Week)
AICctab(omn2, a, b, c, d, e)
# no poly
f <- update(omn2, .~. + log(CH)*Week + log(TDens))
g <- update(omn2, .~. + log(CH) + log(TDens)*Week)
AICctab(omn2, a, b, c, d, e, f, g)
# interaction CH:Week + TD is best
h <- update(omn2, .~. + log(CH)*Week) 
i <- update(omn2, .~. + log(CH))
AICctab(omn2, a, f, h, i)

rm(list=c("a","b","c","d","e","f","g","h","i"))
omn3 <- update(omn2, .~. + log(CH)*Week + log(TDens))


a <- update(omn3, .~. + scaleTPD - log(TDens))
b <- update(omn3, .~. + scaleTPD) 
AICctab(omn3, a, b) # keep
c <- update(omn3, .~. + scaleTPD*Week)
d <- update(omn3, .~. + poly(scaleTPD,2))
e <- update(omn3, .~. + poly(scaleTPD,2)*Week)
AICctab(omn3, a, b, c, d)
# not significant
f <- update(omn3, .~. + TDiv)
g <- update(omn3, .~. + poly(TDiv,2))
h <- update(omn3, .~. + TDiv*Week)
i <- update(omn3, .~. + poly(TDiv,2)*Week)
AICctab(omn3, f, g, h, i)
# not significant
j <- update(omn3, .~. + scaleUDens)
k <- update(omn3, .~. + poly(scaleUDens,2))
l <- update(omn3, .~. + scaleUDens*Week)
m <- update(omn3, .~. + poly(scaleUDens,2)*Week)
AICctab(omn3, j, k, l, m)
# not significant
n <- update(omn3, .~. + TR)
o <- update(omn3, .~. + TR*Week)
p <- update(omn3, .~. + poly(TR,2))
q <- update(omn3, .~. + poly(TR,2)*Week)
AICctab(omn3, n, o, p, q)
# not significant
r <- update(omn3, .~. + DOM) 
s <- update(omn3, .~. + DOM*Week)
t <- update(omn3, .~. + DOM - HabClass)
AICctab(omn3, r, s, t)
# DOM has 3.4 dAICc but 4 df ##
rm(list=c("a","b","c","d","e","f","g","h","i",
          "j","k","l","m","n","o","p","q","r",
          "s","t"))
# no selection. stick with omn3

summary(omn3)
omn3sim <- simulateResiduals(omn3) 
plot(omn3sim)
hist(omn3sim) 
testDispersion(omn3sim)
testZeroInflation(omn3sim)




# glmer is not working for this final model, but glmmTMB took time for each model



### ###


