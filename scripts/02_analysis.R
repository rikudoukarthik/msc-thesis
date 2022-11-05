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
a <- update(all0, . ~ . + Week)
b <- update(all0, . ~ . + poly(Week, 2))
c <- update(all0, . ~ . + poly(Week, 3))
anova(a, b) 
AICctab(all0, a, b, c) # returns best model on top and delta calculated from that
# 2nd and 3rd not better than linear

all1 <- a
summary(all1)


# testing nuisance variables like time of count during day 
a <- update(all1, . ~ . + CoD)
b <- update(all1, . ~ . + poly(CoD, 2))
anova(all1, a) 
AICctab(all1, a, b)
# neither linear nor poly significant
c <- update(all1, . ~ . + Weather)
AICctab(all1, a, b, c) 
# Weather not significant


# finally moving to predictors of interest
a <- update(all1, .~. + HabClass)
b <- update(all1, .~. + Moss)
AICctab(all1, a, b) # Moss ns
anova(all1, a)
summary(a)
# summary function shows that Interior is marginally significant (0.04) but other levels
# are not. AIC is not different at all So will have to ignore.
# summary() p value cannot be relied on.

c <- update(all1, .~. + scaleCC)
d <- update(all1, .~. + logCH) # log better logically and statistically than scale
e <- update(all1, .~. + logTDens)
f <- update(all1, .~. + scaleTPD)
g <- update(all1, .~. + TPD)
h <- update(all1, .~. + TDiv)
i <- update(all1, .~. + scaleUDens)
j <- update(all1, .~. + UDens)
k <- update(all1, .~. + DOM)
AICctab(all1, c, d, e, f, g, h, i, j, k) 
# most variables not signficant, but testing poly before discarding
f <- update(all1, .~. + poly(scaleTPD, 2))
g <- update(all1, .~. + poly(TPD, 2))
h <- update(all1, .~. + poly(TDiv, 2))
i <- update(all1, .~. + poly(scaleUDens, 2))
j <- update(all1, .~. + poly(UDens, 2))
AICctab(all1, c, d, e, f, g, h, i, j, k) 
# poly is ever worse for those variables

# now, moving to variables which have an effect:
all2 <- k
# DOM clearly important, so trying other variables in addition to DOM
a <- update(all1, .~. + scaleCC)
b <- update(all1, .~. + logCH) # log better logically and statistically than scale
c <- update(all1, .~. + logTDens)

d <- update(all2, .~. + scaleCC)
e <- update(all2, .~. + logCH)
f <- update(all2, .~. + logTDens)
AICctab(all1, all2, a, b, c, d, e, f) 
# as expected, CH & TDens better than simple CC so dropping CC (models a and d)
# model with only TDens (c) does not pass dAICc threshold so drop
# contest between models b, e and f---both e and f with DOM are better so drop b

a <- update(all2, .~. + logCH)
b <- update(all2, .~. + logTDens)
c <- update(all2, .~. + logCH + logTDens)
d <- update(all2, .~. + logCH + logTDens + logCH:logTDens)
AICctab(all2, a, b, c, d) 
# both together (or with interaction) does not improve model (delta = 0.3 for 1 df) 
# among a and b, TDens better than CH, so retain only CH

# but trying poly before making decision
e <- update(a, .~. + poly(logCH, 2))
f <- update(a, .~. + poly(logTDens, 2))
AICctab(all2, a, b, e, f)
# poly very poor, both linear better than both poly

all3 <- b
summary(all3)


# interaction of predictors with Week (of pertinence to question)
a <- update(all3, .~. + Week:DOM)
b <- update(all3, .~. + Week:logTDens)
c <- update(all3, .~. + Week:(DOM + logTDens))
AICctab(all3, a, b, c) 
# interaction with both is poor, so drop c
drop1(c)
# dropping interaction with TDens is better than with DOM

all4 <- a


summary(all4)

hist(resid(all4)) # very Gaussian(!), except for those two outliers
simres <- simulateResiduals(all4, plot = T, n = 1000) # beautiful :)
plot(simres)
testDispersion(simres) # 1.4445 but p value sig because of sample size

hist(simres) # uniform as should be
plotResiduals(simres, form = m_all$Week)
plotResiduals(simres, form = m_all$TDens)
plotResiduals(simres, form = m_all$DOM)
# all good


# further selection

drop1(all4)
# can't drop any existing variables

# interactions between predictors. CH:DOM sig but mostly for Moss (only 2 Points)
a <- update(all4, .~. + DOM:logTDens)
AICctab(all4, a)
# close to 2dAICc but better to not include for parsimony
b <- update(all4, .~. - logTDens + DOM:logTDens)
AICctab(all4, a, b)
summary(a)
summary(b)
# same again, although a and b are rather similar wrt AIC and df!

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


a <- update(inv1, .~. + logCH + logTDens)
b <- update(inv1, .~. + poly(logCH,2) + logTDens)
c <- update(inv1, .~. + logCH + poly(logTDens,2))
d <- update(inv1, .~. + poly(logCH,2) + poly(logTDens,2))
e <- update(inv1, .~. + logCH*Week + logTDens*Week)
AICctab(inv1, a, b, c, d, e)
# although still ns, the interaction with Week has the closest AIC to inv1
f <- update(inv1, .~. + poly(logCH,2)*Week + poly(logTDens,2)*Week)
AICctab(inv1, a, b, c, d, e, f)

g <- update(inv1, .~. + logCH*Week + logTDens)
h <- update(inv1, .~. + logCH + logTDens*Week)
AICctab(inv1, a, b, c, d, e, f, g, h)
i <- update(inv1, .~. + logCH*Week + poly(logTDens,2))
j <- update(inv1, .~. + logCH + poly(logTDens,2)*Week)
AICctab(inv1, a, c, e, h, i, j)
k <- update(inv1, .~. + logTDens*Week)
l <- update(inv1, .~. + logTDens)
m <- update(inv1, .~. + logCH*Week)
n <- update(inv1, .~. + logCH)
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


a <- update(omn2, .~. + logCH + logTDens)
b <- update(omn2, .~. + poly(logCH,2) + logTDens)
c <- update(omn2, .~. + logCH + poly(logTDens,2))
d <- update(omn2, .~. + poly(logCH,2) + poly(logTDens,2))
e <- update(omn2, .~. + logCH*Week + logTDens*Week)
AICctab(omn2, a, b, c, d, e)
# no poly
f <- update(omn2, .~. + logCH*Week + logTDens)
g <- update(omn2, .~. + logCH + logTDens*Week)
AICctab(omn2, a, b, c, d, e, f, g)
# interaction CH:Week + TD is best
h <- update(omn2, .~. + logCH*Week) 
i <- update(omn2, .~. + logCH)
AICctab(omn2, a, f, h, i)

rm(list=c("a","b","c","d","e","f","g","h","i"))
omn3 <- update(omn2, .~. + logCH*Week + logTDens)


a <- update(omn3, .~. + scaleTPD - logTDens)
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


