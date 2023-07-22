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
              data = m_all, family = poisson(), 
              control = glmerControl(optimizer = "bobyqa"))
summary(all0)
all0sim <- simulateResiduals(all0)
plot(all0sim)

all0NB <- glmmTMB(Point_Abun ~ Observer + (1|Point) + (1|Days),
                  data = m_all, family = nbinom2)
all0NBsim <- simulateResiduals(all0NB)
plot(all0NBsim)

AICctab(all0, all0NB) # returns best model on top and delta calculated from that

# # choose optimiser to avoid convergence issue:
# 
# # install.packages("optimx")
# # install.packages("dfoptim")
# opt <- allFit(model_not_converging)
# summary(opt)$which.OK # only 1/7 optimisers failed
# summary(opt)$llik # all pretty similar values so I think bobyqa good choice
# summary(opt)$fixef
# # use bobyqa


# adding Week and testing for polynomial 
a <- update(all0NB, . ~ . + Week)
b <- update(all0NB, . ~ . + poly(Week, 2))
c <- update(all0NB, . ~ . + poly(Week, 3))
anova(a, b) 
AICctab(all0NB, a, b, c) # returns best model on top and delta calculated from that
# 2nd and 3rd not better than linear

all1 <- a
summary(all1)


# testing nuisance variables like time of count during day 
a <- update(all1, . ~ . + CoD)
b <- update(all1, . ~ . + poly(CoD, 2))
# anova(all1, a) 
AICctab(all1, a, b)
# neither linear nor poly significant
c <- update(all1, . ~ . + Weather)
AICctab(all1, a, b, c) 
# Weather not significant


# finally moving to predictors of interest
a <- update(all1, .~. + HabClass)
b <- update(all1, .~. + Moss)
AICctab(all1, a, b) # both ns
anova(all1, a)
summary(a)
c <- update(all1, .~. + scaleCC)
d <- update(all1, .~. + logCH) # log better logically and statistically than scale
e <- update(all1, .~. + logTDens)
f <- update(all1, .~. + scaleTPD)
AICctab(all1, a, b, c, d, e, f) 

all2 <- e
summary(all2)

a <- update(all2, .~. + TPD) # failed to converge, need to stick with scaled var
b <- update(all2, .~. + TDiv)
c <- update(all2, .~. + scaleUDens)
d <- update(all2, .~. + UDens)
e <- update(all2, .~. + DOM)
AICctab(all2, a, b, c, d, e) 

all3 <- e
summary(all3)

a <- update(all3, .~. + logCH)
b <- update(all3, .~. + logCH - logTDens)
c <- update(all3, .~. + logCH - DOM)
d <- update(all3, .~. + scaleTPD)
e <- update(all3, .~. + scaleTPD - logTDens)
f <- update(all3, .~. + scaleCC)
g <- update(all3, .~. + scaleCC - logTDens)
AICctab(all3, a, b, c, d, e, f, g) 
# most variables not significant, but testing poly before discarding
h <- update(all3, .~. + poly(scaleCC, 2))
i <- update(all3, .~. + poly(logCH, 2))
j <- update(all3, .~. + poly(scaleTPD, 2))
k <- update(all3, .~. + poly(TPD, 2))
l <- update(all3, .~. + poly(TDiv, 2))
m <- update(all3, .~. + poly(scaleUDens, 2))
n <- update(all3, .~. + poly(UDens, 2))
AICctab(all3, h, i, j, k, l, m, n) 
# poly variables also don't add info

a <- update(all2, .~. + logCH + logCH:logTDens)
b <- update(all3, .~. + logCH + logCH:logTDens)
AICctab(all2, all3, a, b) 
# b is the best
# but trying poly before making decision
c <- update(all2, .~. + poly(logCH, 2))
d <- update(all1, .~. + poly(logCH, 2))
e <- update(all3, .~. + poly(logCH, 2))
f <- update(all3, .~. + poly(logCH, 2) + poly(logCH, 2):logTDens)
g <- update(all1, .~. + poly(logTDens, 2))
AICctab(all3, a, b, c, d, e, f, g)
# poly not better than linear

a <- b
# interaction of predictors with Week (of pertinence to question)
b <- update(a, .~. + Week:logCH)
c <- update(a, .~. + Week:logTDens)
d <- update(a, .~. + Week:(logCH + logTDens))
e <- update(a, .~. + Week:logCH - logCH:logTDens)
f <- update(a, .~. + Week:logTDens - logCH:logTDens)
g <- update(a, .~. + Week:logCH + Week:logTDens - logCH:logTDens)
AICctab(all3, a, b, c, d, e, f, g)
# interaction is poor
h <- update(a, .~. + Week:DOM)
i <- update(a, .~. + Week:DOM - logCH:logTDens)
j <- update(a, .~. + Week:DOM - logCH*logTDens)
AICctab(all3, a, h, i, j)


all4 <- a
# final model is all4
summary(all4)

hist(resid(all4)) 
simres <- simulateResiduals(all4, plot = T, n = 1000) # beautiful :)
testDispersion(simres) 
testZeroInflation(simres) 

hist(simres) # uniform as should be
plotResiduals(simres, form = m_all$Week)
plotResiduals(simres, form = m_all$logTDens)
plotResiduals(simres, form = m_all$logCH)
plotResiduals(simres, form = m_all$DOM)
# all good, only CH tests sig but the plot shows no major deviance


# further selection
drop1(all4)

# can't drop any existing variables

### Guilds: invertebrate-feeding ####

m_guild_inv <- m_guild %>% filter(GuildFeed == "Invertebrate")

# invertebrate has 22.6% zeroes (94/416); omnivore 7.9%


# finally decided on a mixed model with Poisson error using glmmTMB with no
# zero inflation. singular fit issues resolved by using glmmTMB over glmer, and estimated 
# parameters in both methods were pretty much the same.


inv0 <- glmmTMB(Point_Abun ~ Observer + (1|Point),
                data = m_guild_inv, family = poisson())
summary(inv0)
inv0sim <- simulateResiduals(inv0)
plot(inv0sim)
hist(inv0sim)


# adding Week and testing poly
a <- update(inv0, .~. + Week)
b <- update(inv0, .~. + poly(Week, 2))
AICctab(inv0, a, b)
# main effect not significant
# nuisance variables like CoD and Weather
c <- update(inv0, .~. + CoD)
d <- update(inv0, .~. + poly(CoD, 2))
e <- update(inv0, .~. + Weather)
AICctab(inv0, c, d, e) # CoD important

inv1 <- c
summary(inv1)


# predictors of interest
a <- update(inv1, .~. + HabClass)
b <- update(inv1, .~. + Moss)
c <- update(inv1, .~. + HabClass*Week)
AICctab(inv1, a, b, c)
# habclass not useful at all, Moss is.
d <- update(inv1, .~. + Moss*Week)
e <- update(inv1, .~. + Moss:Week)
# CoD is nuisance variable, so will retain (won't test keeping moss and removing CoD)
AICctab(inv1, b, d, e)
f <- update(inv1, .~. + Moss*poly(Week, 2))
g <- update(inv1, .~. + Moss:poly(Week, 2))
AICctab(inv1, d, e, f, g)
# Moss:Week could be captured in DOM
h <- update(inv1, .~. + DOM*Week)
i <- update(inv1, .~. + DOM:Week)
AICctab(inv1, d, e, h, i) # it's not

AICctab(inv1, d, e)
# e only has interaction, without main effects, and d only has 1.6 AICc higher for 1 df,
# meaning it is not significantly worse than e
inv2 <- d
anova(inv1, inv2)


a <- update(inv2, .~. + scaleCC)
b <- update(inv2, .~. + logCH) # log better logically and statistically than scale
c <- update(inv2, .~. + logTDens)
d <- update(inv2, .~. + scaleTPD)
e <- update(inv2, .~. + TPD)
f <- update(inv2, .~. + TDiv)
g <- update(inv2, .~. + scaleUDens)
h <- update(inv2, .~. + UDens)
i <- update(inv2, .~. + DOM)
AICctab(inv2, a, b, c, d, e, f, g, h, i) 
# most variables not significant, but testing poly before discarding
j <- update(inv2, .~. + poly(scaleCC, 2))
k <- update(inv2, .~. + poly(logCH, 2)) # log better logically and statistically than scale
l <- update(inv2, .~. + poly(logTDens, 2))
m <- update(inv2, .~. + poly(scaleTPD, 2))
n <- update(inv2, .~. + poly(TPD, 2))
o <- update(inv2, .~. + poly(TDiv, 2))
p <- update(inv2, .~. + poly(scaleUDens, 2))
q <- update(inv2, .~. + poly(UDens, 2))
AICctab(inv2, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q) 
# poly is even worse for those variables
# testing interaction with week for CH and TDens
r <- update(inv2, .~. + logCH:Week) # log better logically and statistically than scale
s <- update(inv2, .~. + logTDens:Week)
AICctab(inv2, a, b, c, d, e, f, g, h, i, r, s) 
t <- update(inv2, .~. + logCH + logTDens)
AICctab(inv2, a, b, c, d, e, f, g, h, i, t) 
# not useful!
u <- update(inv2, .~. + scaleTPD*Week)
v <- update(inv2, .~. + poly(scaleTPD, 3))
w <- update(inv2, .~. + poly(scaleTPD, 4))
AICctab(inv2, a, b, c, d, e, f, g, h, i, u, v, w)
y <- update(inv2, .~. + DOM*Week)
AICctab(inv2, i, y)

# only DOM
inv3 <- i


a <- update(inv3, .~. - DOM + TDiv)
b <- update(inv3, .~. + TDiv)
c <- update(inv3, .~. + poly(TDiv, 2))
d <- update(inv3, .~. + TDiv*Week)
AICctab(inv3, a, b, c, d)
# TDiv not important
e <- update(inv3, .~. + sqrt(UDens))
f <- update(inv3, .~. + poly(sqrt(UDens), 2))
g <- update(inv3, .~. + sqrt(UDens)*Week)
h <- update(inv3, .~. + poly(sqrt(UDens), 2)*Week)
AICctab(inv3, e, f, g, h)
# UDens not important
i <- update(inv3, .~. + TR)
j <- update(inv3, .~. + TR*Week)
AICctab(inv3, i, j)
k <- update(inv3, .~. + poly(TR, 2))
AICctab(inv3, i, k)

# finally, testing DOM and comparing with Moss
l <- update(inv3, .~. + DOM:Week)
m <- update(inv3, .~. + DOM:Week - Moss)
n <- update(inv3, .~. - Moss:Week)
AICctab(inv3, l, m, n)

# inv3 is the final model
summary(inv3)


hist(resid(inv3)) # very Gaussian(!), except for those two outliers
simres <- simulateResiduals(inv3, plot = T, n = 1000) 
testDispersion(simres) 
testZeroInflation(simres)

hist(simres) # uniform as should be
plotResiduals(simres, form = m_guild_inv$Week)
plotResiduals(simres, form = m_guild_inv$Moss)
plotResiduals(simres, form = m_guild_inv$CoD)
plotResiduals(simres, form = m_guild_inv$DOM)
# all good
# only Moss in DOM has slightly higher mean so sig. Levene test, but huge variance


# further selection
drop1(inv3)

# can't drop any existing variables

### Guilds: omnivore-feeding ####

m_guild_omn <- m_guild %>% filter(GuildFeed == "Omnivore")
# 7.9% zeroes (33/416)

plot(m_guild_omn$Point_Abun ~ m_guild_omn$Week)
# omni seem to show much clearer pattern (increase). migrants


omn0 <- glmmTMB(Point_Abun ~ Observer + (1|Point), 
                data = m_guild_omn, family = poisson)
omn0sim <- simulateResiduals(omn0)
plot(omn0sim)
# KS test, dispersion test, outlier test all sig
omn0NB <- glmmTMB(Point_Abun ~ Observer + (1|Point), 
                  data = m_guild_omn, family = nbinom2)
omn0NBsim <- simulateResiduals(omn0NB)
plot(omn0NBsim)
# corrected for overdispersion using negbin


# adding Week and testing poly
a <- update(omn0NB, .~. + Week)
b <- update(omn0NB, .~. + poly(Week, 2))
AICctab(omn0, omn0NB, a, b)
# linear effect of Week very good

omn1 <- a

# nuisance variables like CoD and Weather
a <- update(omn1, .~. + CoD)
b <- update(omn1, .~. + poly(CoD, 2))
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

omn2 <- d
anova(omn1, omn2)


a <- update(omn2, .~. + scaleCC)
b <- update(omn2, .~. + logCH) # log better logically and statistically than scale
c <- update(omn2, .~. + logTDens)
d <- update(omn2, .~. + scaleTPD)
e <- update(omn2, .~. + TPD)
f <- update(omn2, .~. + TDiv)
g <- update(omn2, .~. + scaleUDens)
h <- update(omn2, .~. + UDens)
i <- update(omn2, .~. + DOM)
AICctab(omn2, a, b, c, d, e, f, g, h, i) 
# logCH best, logTDens close
j <- update(omn2, .~. + logCH + logTDens)
k <- update(omn2, .~. + logCH*logTDens)
AICctab(omn2, b, j, k) # j is best
# testing poly before discarding
l <- update(omn2, .~. + poly(scaleCC, 2))
m <- update(omn2, .~. + poly(logCH, 2)) # log better logically and statistically than scale
n <- update(omn2, .~. + poly(logTDens, 2))
o <- update(omn2, .~. + poly(scaleTPD, 2))
p <- update(omn2, .~. + poly(TPD, 2))
q <- update(omn2, .~. + poly(TDiv, 2))
r <- update(omn2, .~. + poly(scaleUDens, 2))
s <- update(omn2, .~. + poly(UDens, 2))
AICctab(omn2, j, l, m, n, o, p, q, r, s) 
# poly is even worse for those variables

a <- j
# testing interaction with week for CH and TDens
b <- update(a, .~. + logCH:Week) 
c <- update(a, .~. + logTDens:Week)
d <- update(a, .~. + logCH:Week - logCH) 
e <- update(a, .~. + logTDens:Week - logTDens)
f <- update(a, .~. + logCH:Week - logCH - logTDens) 
g <- update(a, .~. + logTDens:Week - logTDens - logCH)
AICctab(omn2, a, b, c, d, e, f, g) 
# TDens interaction with week is not important, but main effect is important
AICctab(a, b, d) 
# interaction of CH with week without its main effect is clearly the best; even though 
# it's only interaction, dAICc and df can't be argued with


a <- d
b <- update(omn2, .~. + scaleTPD*Week)
c <- update(omn2, .~. + poly(scaleTPD, 3))
d <- update(omn2, .~. + poly(scaleTPD, 4))
e <- update(omn2, .~. + logCH + logTDens*Week)
AICctab(omn2, a, b, c, d, e) 

omn3 <- a


a <- update(omn3, .~. + scaleTPD - logCH:Week)
b <- update(omn3, .~. + scaleTPD) 
AICctab(omn3, a, b) 
c <- update(omn3, .~. + scaleTPD*Week)
d <- update(omn3, .~. + poly(scaleTPD,2))
e <- update(omn3, .~. + poly(scaleTPD,2)*Week)
AICctab(omn3, a, b, c, d, e)
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
# too df-costly


# final model is omn3

summary(omn3)

hist(resid(omn3)) 
simres <- simulateResiduals(omn3, plot = T, n = 1000) # beautiful :)
testDispersion(simres)
testZeroInflation(simres)

hist(simres) # uniform as should be
plotResiduals(simres, form = m_guild_omn$Week)
plotResiduals(simres, form = m_guild_omn$HabClass)
plotResiduals(simres, form = m_guild_omn$logTDens)
plotResiduals(simres, form = m_guild_omn$logCH)
# all good!
# only CH tests sig but the plot shows no major deviance (just overfitted)


### Saving data ####

save(all0, all1, all2, all3, all4,
     inv0, inv1, inv2, inv3,
     omn0, omn1, omn2, omn3,
     file = "data/02_analysis.RData", pos = ".GlobalEnv")
