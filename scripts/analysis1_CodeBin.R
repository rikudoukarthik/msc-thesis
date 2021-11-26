# Import datasets ####
# first definition of bird data used for analysis
b_data_A1 <- b_data %>%  
  group_by(Point, Week, Date, Observer, Weather, Wind, Visibility, Spec_code) %>% 
  summarise(Abundance=sum(Number)) %>% 
  ungroup() %>% 
  mutate(Days = as.integer((as.Date(Date, format="%d-%m-%y") -
                              as.Date("07-07-2020", format="%d-%m-%y"))), .keep = "unused") %>% 
  mutate(Period = as.factor(case_when(Week %in% 1:4 ~ 1,
                                      Week %in% 5:8 ~ 2,
                                      Week %in% 9:13 ~ 3))) %>% 
  group_by(Week, Spec_code) %>% 
  mutate(RelAbun=round(Abundance/sum(Abundance),4)) %>% 
  inner_join(b_speccode, by="Spec_code")


## for ordination ##

# combining information from Weeks 1 and 2 for first ordination
b_data_ordW1 <- b_data_A1 %>% ungroup() %>% 
  mutate(Spec_code = fct_collapse(Spec_code, 
                                  Cer_sp = c("Cer_bra","Cer_fam","Cer_sp"),
                                  Phy_sp = c("Phy_col","Phy_tro","Phy_sp"),
                                  Reg_sp = c("Reg_reg","Reg_ign","Reg_sp"),
                                  Den_sp = c("Den_maj","Den_med","Den_sp"))) %>% 
  filter(Spec_code=="Par_maj"|
           Spec_code=="Cya_cae"|
           Spec_code=="Sit_eur"|
           Spec_code=="Eri_rub"|
           Spec_code=="Cer_sp"|
           Spec_code=="Phy_sp"|
           Spec_code=="Reg_sp"|
           Spec_code=="Den_sp") %>% 
  filter(Week==1|Week==2) %>%
  group_by(Point, Spec_code) %>% summarise(Abundance = sum(Abundance)) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")


# combining information from Weeks 12 and 13 for second ordination
b_data_ordW13 <- b_data_A1 %>% ungroup() %>% 
  mutate(Spec_code = fct_collapse(Spec_code, 
                                  Cer_sp = c("Cer_bra","Cer_fam","Cer_sp"),
                                  Phy_sp = c("Phy_col","Phy_tro","Phy_sp"),
                                  Reg_sp = c("Reg_reg","Reg_ign","Reg_sp"),
                                  Den_sp = c("Den_maj","Den_med","Den_sp"))) %>% 
  filter(Spec_code=="Par_maj"|
           Spec_code=="Cya_cae"|
           Spec_code=="Sit_eur"|
           Spec_code=="Eri_rub"|
           Spec_code=="Cer_sp"|
           Spec_code=="Phy_sp"|
           Spec_code=="Reg_sp"|
           Spec_code=="Den_sp") %>% 
  filter(Week==12|Week==13) %>%
  group_by(Point, Spec_code) %>% summarise(Abundance = sum(Abundance)) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")



# using mean rather than sum #
# combining information from Weeks 1:3 for first ordination
b_data_ordPBmean <- b_data_A1 %>% ungroup() %>% 
  mutate(Spec_code = fct_collapse(Spec_code, 
                                  Cer_sp = c("Cer_bra","Cer_fam","Cer_sp"),
                                  Phy_sp = c("Phy_col","Phy_tro","Phy_sp"),
                                  Reg_sp = c("Reg_reg","Reg_ign","Reg_sp"),
                                  Den_sp = c("Den_maj","Den_med","Den_sp"))) %>% 
  filter(Spec_code=="Par_maj"|
           Spec_code=="Cya_cae"|
           Spec_code=="Sit_eur"|
           Spec_code=="Eri_rub"|
           Spec_code=="Cer_sp"|
           Spec_code=="Phy_sp"|
           Spec_code=="Reg_sp"|
           Spec_code=="Den_sp") %>% 
  filter(Week==1|Week==2|Week==3) %>%
  group_by(Point, Spec_code) %>% summarise(Abundance = ceiling(mean(Abundance))) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")

# combining information from Weeks 11:13 for second ordination
b_data_ordPMmean <- b_data_A1 %>% ungroup() %>% 
  mutate(Spec_code = fct_collapse(Spec_code, 
                                  Cer_sp = c("Cer_bra","Cer_fam","Cer_sp"),
                                  Phy_sp = c("Phy_col","Phy_tro","Phy_sp"),
                                  Reg_sp = c("Reg_reg","Reg_ign","Reg_sp"),
                                  Den_sp = c("Den_maj","Den_med","Den_sp"))) %>% 
  filter(Spec_code=="Par_maj"|
           Spec_code=="Cya_cae"|
           Spec_code=="Sit_eur"|
           Spec_code=="Eri_rub"|
           Spec_code=="Cer_sp"|
           Spec_code=="Phy_sp"|
           Spec_code=="Reg_sp"|
           Spec_code=="Den_sp") %>% 
  filter(Week==11|Week==12|Week==13) %>%
  group_by(Point, Spec_code) %>% summarise(Abundance = ceiling(mean(Abundance))) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")




# Code to clean environment ####

# the models that are precursors to the final chosen models
rm(list=c("egtA.glm1","egtA.gam1","egtA.gam2","egtA.gam3","egtA.gam4","egtA.gam5",
          "egtRA.glm1","egtRA.glm2","egtRA.gam1","egtRA.gam2","egtRA.gam3","egtRA.gam4",
          "egtRA.gam5"))
rm(list=c("allA.glm0","allA.glm1","allA.glm2","allA.gam1","allA.gam2","allA.gam3"))
rm(list=c("allRA.glm0","allRA.gam0","allRA.gam1","allRA.gam2","allRA.gam4"))
rm(list=c("gldA.glm0","gldA.glm1","gldA.glm2"))

# the models that need not be created to get chosen models, but can be removed if
# they are created anyway
rm(list=c("egtA.glm3","egtA.glm4","egtA.gam7","allA.gam5","allRA.glm2","allRA.gam3",
          "gldA.glm4","gldA.gam1"))
### Analysis 1: Building models for all birds #1 ####

### GLM using Abundance
allA.glm1 <- step(allA.glm0, scope = ~Period+log(CCavgsd)+log(TreeDens)+TreePropDeci+
                    DOM+Moss+HabClass+sqrt(UndDens))
# checking if CCavg is better than CCavgsd+TreeDens
allA.glm2 <- step(allA.glm0, scope = ~Period+log(CCavg)+TreePropDeci+
                    DOM+Moss+HabClass+sqrt(UndDens))
anova(allA.glm1,allA.glm2,test="Chisq")
# no. using the two variables is better than just CCavg.


# interactions
allA.glm3 <- step(allA.glm1, scope = ~ . + Period:(DOM + HabClass + log(CCavgsd) + 
                                                     log(TreeDens) + TreePropDeci))
summary(allA.glm3)
# why is it not possible to have only interaction without main effect?

# Go with allA.glm3 ###

allA.glm4 <- update(allA.glm3, .~. + Observer + Weather)
summary(allA.glm4)


# visualising predicted and observed
#
# All.pred <- predict(allA.glm3qp, type="response", se.fit = T)
# All.predgrid <- m_data_A1_all
# All.predgrid$BirdAbun.fit <- All.pred$fit
# All.predgrid$BirdAbun.fitse <- All.pred$se.fit
# 
# ggplot(aes(x=log(TreeDens), y=BirdAbun, col=Period), data=All.predgrid) + 
#   theme_bw() + geom_point(size=3) + 
#   geom_line(aes(y=BirdAbun.fit), size=3) 


# checking quasipoisson. allA.glm3 is very overdispersed (1198 null dev. on 400df)

allA.glm3qp <- update(allA.glm3, family=quasipoisson())
anova(allA.glm3qp, test="F")
anova(allA.glm3qp, test="Cp")
summary(allA.glm3qp)

allA.glm3nb <- glm.nb(BirdAbun ~ Period + DOM + log(CCavgsd) + HabClass + 
                        TreePropDeci + log(TreeDens) + Period:DOM + Period:log(TreeDens) + 
                        Period:TreePropDeci, data = m_data_A1_all)
summary(allA.glm3nb)
drop1(allA.glm3nb)
# nb has better Deviance Goodness of Fit
anova(allA.glm3nb, test="Chisq")
anova(allA.glm3, allA.glm3qp, allA.glm3nb, test="Cp")

allA.glm3nb2 <- update(allA.glm3nb, .~ Period + DOM)
summary(allA.glm3nb2)
anova(allA.glm3, allA.glm3qp, allA.glm3nb, allA.glm3nb2, test="Cp")




### GLM using RelAbun

allRA.glm0 <- glm(BirdRelAbun ~ Observer + Weather + Period, 
                  family=Gamma(), data=m_data_A1_all)
summary(allRA.glm0)

allRA.glm1 <- step(allRA.glm0, scope = ~.+log(CCavgsd)+log(TreeDens)+
                     TreePropDeci+DOM+Moss+HabClass+sqrt(UndDens))
summary(allRA.glm1) # 38% explained deviance, Period removed

allRA.glm2 <- step(allRA.glm1, scope= ~ . + Period:(log(CCavgsd)+log(TreeDens)+
                                                      TreePropDeci+DOM+Moss+HabClass+sqrt(UndDens)))
# none of the interaction terms significant, and even overall model is poor
# Go with allRA.glm1 ###


### GAM using RelAbun

allRA.gam0 <- gam(BirdRelAbun ~ Observer + Period, 
                  family=Gamma(), data=m_data_A1_all)

allRA.gam1 <- update(allRA.gam0, .~. + s(log(CCavgsd)) + s(log(TreeDens)) +
                       s(TreePropDeci) + DOM + Moss + HabClass + s(sqrt(UndDens)),
                     select=T)
summary(allRA.gam1)
# remove s(TreePropDeci), Moss; Period still ns but need it for interaction

allRA.gam2 <- update(allRA.gam1, .~. - s(TreePropDeci) - Moss)
summary(allRA.gam2)

# no interaction is significant. takes time to run, so don't.
# allRA.gam3 <- update(allRA.gam2, .~. + Period:(DOM + Moss + HabClass) +
#                       s(TreePropDeci, by=Period) + s(log(CCavgsd), by=Period) + 
#                       s(log(TreeDens), by=Period) + s(sqrt(UndDens), by=Period),
#                      select=T)
# summary(allRA.gam3)


plot(allRA.gam2, shade=T, seWithMean = T, pages=1)


allRA.gam4 <- update(allRA.gam2, .~. + s(log(CCavgsd), by = Period) +
                       s(log(TreeDens), by = Period) + s(sqrt(UndDens), by = Period),
                     select = T)
summary(allRA.gam4)

allRA.gam5 <- update(allRA.gam4, .~. - Period - s(log(CCavgsd)) - s(log(TreeDens)) -
                       s(sqrt(UndDens), by = Period))
summary(allRA.gam5) # Observer 32%

AIC(allRA.gam2, allRA.gam4, allRA.gam5)
plot(allRA.gam5, shade=T, seWithMean = T, pages=1)

# Go with allRA.gam5 ###

### Analysis 1: Building models for all birds #2 ####
#2 
### Using Abundance ###

allA.0 <- gam(BirdAbun ~ Observer + CoD + Weather + Period, data=m_data_A1_all,
              family=nb())
summary(allA.0)
# 4.5% explained with just CoD and Period. 39.2% on adding Observer.
# dispersion is high. 
# tried s(CoD) but still not significant and edf=1.001 so linear term better
# CoD:Period also not significant

allA.0qp <- gam(BirdAbun ~ Observer + CoD + Weather + Period, data=m_data_A1_all,
                family=quasipoisson())
summary(allA.0qp)
# qp explains 0.1% more than nb

allA.0mix <- glmer.nb(BirdAbun ~ Observer + CoD + Weather + Period + (1|Point), 
                      data=m_data_A1_all)
summary(allA.0mix)


allA.1 <- update(allA.0, .~. - CoD + s(log(CCavgsd)) + s(log(TreeDens)) + 
                   s(TreePropDeci) + DOM + Moss + HabClass + s(sqrt(UndDens)),
                 select=T)
summary(allA.1)
# Moss and s(TreePropdeci) not significant. Weather became ns in this model.
# edf of s(log(TreeDens)) = 0.8 so linear term better



allA.2 <- update(allA.1, .~. - Moss - s(TreePropDeci) - s(log(TreeDens)) +
                   log(TreeDens), select=T)
summary(allA.2)
# weather still borderline significant so will keep in model

allA.2qp <- update(allA.2, family=quasipoisson)
summary(allA.2qp)
# qp explains 0.4% more. Weather significant


allA.3 <- update(allA.2, .~. + Period:(DOM + HabClass + log(TreeDens)) + 
                   s(log(CCavgsd), by=Period) + s(sqrt(UndDens), by=Period),
                 select=T)
summary(allA.3)
# interactions Period:HabClass, Period:log(TreeDens), Period:s(log(CCavgsd)) and 
# Period:s(sqrt(UndDens)) not significant. Weather and DOM became ns.

allA.4 <- update(allA.2, .~. + Period:DOM - DOM - Weather)
summary(allA.4)
# a lot of the significant variables become ns on removing Observer, inc. smooth terms
# model explains 12.5% without Observer and 50.1% with.

plot(allA.4, shade=T, seWithMean = T, pages=1)


allA.4qp <- update(allA.4, family=quasipoisson())
summary(allA.4qp)


allA.4mix <- gamm(BirdAbun ~ Observer + Period + s(log(CCavgsd)) + HabClass + 
                    s(sqrt(UndDens)) + log(TreeDens) + Period:DOM,
                  random=list(Point=~1),
                  family=quasipoisson(), data=m_data_A1_all)
summary(allA.4mix$lme)
summary(allA.4mix$gam)


allA.4glmm <- glmer.nb(BirdAbun ~ Observer + Period + poly(log(CCavgsd),2) + HabClass + 
                         log(TreeDens) + Period:DOM + (1|Point),
                       data=m_data_A1_all)
AIC((allA.4glmm))



### Using RelAbun ###

allRA.0 <- gam(BirdRelAbun ~ Observer + Weather + Period, 
               family=Gamma(), data=m_data_A1_all)

allRA.1 <- update(allRA.0, .~. + s(log(CCavgsd)) + s(log(TreeDens)) +
                    s(TreePropDeci) + DOM + Moss + HabClass + s(sqrt(UndDens)),
                  select=T)
summary(allRA.1)
# Weather, Period, Moss, s(TreePropDeci) not sig

allRA.2 <- update(allRA.1, .~. - Weather - Period - Moss - s(TreePropDeci))
summary(allRA.2)

allRA.3 <- update(allRA.2, .~. + s(log(CCavgsd), by = Period) +
                    s(log(TreeDens), by = Period) + s(sqrt(UndDens), by = Period),
                  select = T)
summary(allRA.3)
# interactions don't seem to be significant

allRA.4 <- update(allRA.3, .~. - s(log(CCavgsd), by=Period) - s(log(TreeDens)) -
                    s(sqrt(UndDens), by = Period))
summary(allRA.4)


allRA.4mix <- glmer(BirdRelAbun ~ Observer + DOM + HabClass  
                    + poly(log(CCavgsd),2)
                    + poly(sqrt(UndDens),2)
                    + poly(log(TreeDens),2) + (1|Point),
                    family=Gamma(), data=m_data_A1_all)
summary(allRA.4mix)


# with habsel2

allHS.4 <- gam(BirdHabSel2 ~ Observer + s(log(CCavgsd)) + DOM + HabClass + 
                 s(sqrt(UndDens)) + s(log(TreeDens), by = Period),
               family = Gamma(link = log), data=m_data_A1_all)
summary(allHS.4)

plot(allHS.4, seWithMean = T, shade=T, pages=1)

allHS.4gamm <- gamm(BirdHabSel2 ~ Observer + s(log(CCavgsd)) + DOM + HabClass + 
                      s(sqrt(UndDens)) + s(log(TreeDens), by = Period),
                    random=list(Point=~1),
                    family = Gamma(link = log), data=m_data_A1_all)
summary(allHS.4gamm$gam)
plot(allHS.4gamm$gam, seWithMean = T, shade=T, pages=1)
# log(TreeDens) linear


allHS.4glmm <- glmer(BirdHabSel2 ~ Observer + DOM + HabClass + poly(log(CCavgsd),2) +  
                       poly(sqrt(UndDens),2) + log(TreeDens):Period +
                       (1|Point),
                     family = Gamma(link = log), data=m_data_A1_all)
summary(allHS.4glmm)
plot(allHS.4glmm$fitted)
AIC(allHS.4glmm, allHS.4gamm$lme)
confint(allHS.4glmm)

p <- predict(allHS.4glmm, type="response")
plot(log(m_data_A1_all$CCavgsd), p)

### ###

#3

all.0i <- lme4::lmer(BirdAbun ~ Observer + (1|Point) + (1|Days),
                     data = m_data_A1_all)
all.0ii <- glmer(BirdAbun ~ Observer + (1|Point) + (1|Days),
                 data = m_data_A1_all, family=poisson())
summary(all.0i)
summary(all.0ii)

plot(all.0i, type=c("p","smooth"), col.line=1) # fit vs resid w/ smooth line
plot(all.0i, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres
lattice::qqmath(all.0i)

plot(all.0ii, type=c("p","smooth"), col.line=1) # fit vs resid w/ smooth line
plot(all.0ii, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres
lattice::qqmath(all.0ii)


all.1a <- lme4::lmer(BirdAbun ~ Observer + (1|Point) + (1|Days) + Week,
                     data = m_data_A1_all)
all.1b <- glmer(BirdAbun ~ Observer + (1|Point) + (1|Days) + Week,
                data = m_data_A1_all, family=poisson())
summary(all.1a)
summary(all.1b)
plot(all.1a, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres
lattice::qqmath(all.1a)
plot(all.1b, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres
lattice::qqmath(all.1b)


all.1ax <- update(all.1a, .~. + poly(log(CCavgsd),2) + poly(sqrt(UndDens),2) +
                    DOM + HabClass, REML=FALSE) # singular fit
summary(all.1ax)
plot(all.1ax, type=c("p","smooth"), col.line=1) # fit vs resid w/ smooth line
plot(all.1ax, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres
lattice::qqmath(all.1ax)

hist(resid(all.1ax))
all.1axSim <- simulateResiduals(fittedModel = all.1ax, plot=T, n=1000)
plot(all.1axSim)
testDispersion(all.1axSim)


all.1bx <- update(all.1b, .~. + poly(log(CCavgsd),2) + poly(sqrt(UndDens),2) +
                    DOM + HabClass)
summary(all.1bx)
plot(all.1bx, type=c("p","smooth"), col.line=1) # fit vs resid w/ smooth line
plot(all.1bx, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres
lattice::qqmath(all.1bx)
hist(resid(all.1bx))
all.1bxSim <- simulateResiduals(fittedModel = all.1bx, plot=T, n=1000)
plot(all.1bxSim)
testDispersion(all.1bxSim)
# lower n for testDispersion
all.1bxSim2 <- simulateResiduals(all.1bx, plot=T, n=100)
testDispersion(all.1bxSim2)
# testTemporalAutocorrelation(all.1bxSim2, m_data_A1_all$Week) cannot be done, as 
# there are repeated obs. per time point so would have to manually test for each Point
# separately. 
testZeroInflation(all.1bxSim)


# normal and poisson won't work

all.1t <- lme4::lmer(log(BirdAbun) ~ Observer + (1|Point) + (1|Days) + Week,
                     data = m_data_A1_all)
summary(all.1t)
plot(all.1t, type=c("p","smooth"), col.line=1) # fit vs resid w/ smooth line
plot(all.1t, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres
lattice::qqmath(all.1t)
hist(resid(all.1t), breaks=15)



all.0 <- glmer.nb(BirdAbun ~ Observer + (1|Point) + (1|Days) + Week, 
                  data = m_data_A1_all)
summary(all.0)
lattice::qqmath(all.0)

all.0x <- update(all.0, .~. + poly(log(CCavgsd),2) + poly(sqrt(UndDens),2) +
                   DOM + HabClass)
summary(all.0x)
lattice::qqmath(all.0x)


m_data_A1_all <- m_data_A1_all %>% ungroup() %>%  mutate(ObsID = as.factor(1:401))

all.0p <- glmer(BirdAbun ~ Observer + (1|Point) + (1|Days) + Week + (1|ObsID),
                data = m_data_A1_all, family=poisson())
summary(all.0p)
lattice::qqmath(all.0p)
plot(all.0p, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres


all.0px <- update(all.0p, .~. + poly(log(CCavgsd),2) + poly(sqrt(UndDens),2) +
                    DOM + HabClass)
summary(all.0px)
lattice::qqmath(all.0px)
plot(all.0px, sqrt(abs(resid(.)))~fitted(.), type=c("p","smooth"), col.line=1) # fit v sqabsres



## ##

# p values 

?afex::mixed
afex::mixed(BirdAbun ~ Observer + (1 | Point) + (1 | Days) + Week + log(CCavgsd) +  
              log(TreeDens) + DOM + Week:log(CCavgsd) + Week:log(TreeDens) +  
              Week:DOM,
            data=m_data_A1_all, family=poisson(), method="LRT")
# something seems to be wrong, different from earlier AIC and anova results
car::Anova(all.4)
# slightly closer to earlier results, but still interaction of Week with habvar ns



### ###

### Analysis 2: Building models for guilds #1 ####


gldA.glm0 <- glm.nb(GuildAbun ~ Observer + Weather + Period,
                    data=m_data_A1_gld)
summary(gldA.glm0)

gldA.glm1 <- update(gldA.glm0, .~. + GuildFeed)
summary(gldA.glm1)

gldA.glm2 <- update(gldA.glm1, .~. + Period:GuildFeed)
summary(gldA.glm2)

gldA.glm3 <- step(gldA.glm1, scope= ~ . + DOM + Moss + HabClass + TreePropDeci +
                    log(CCavgsd) + log(TreeDens) + sqrt(UndDens))
anova(gldA.glm1, gldA.glm3)
summary(gldA.glm3) # 45.27 explained deviance

# interactions
gldA.glm4 <- update(gldA.glm3, .~. + GuildFeed:log(TreeDens))
# GuildFeed:DOM and :HabClass ns
summary(gldA.glm4)

gldA.gam1 <- gam(GuildAbun ~ Observer + Period + GuildFeed + 
                   DOM + HabClass + s(log(TreeDens)) + s(log(TreeDens), by=GuildFeed), 
                 data = m_data_A1_gld, family=nb)
summary(gldA.gam1)                    
AIC(gldA.glm3, gldA.gam1)
# interaction only sig for carnivore, and increase in explanation is not major

# Go with gldA.glm3 ###


# Days 

gldA.glm5 <- update(gldA.glm3, .~. - Period + Days + GuildFeed:Days)
summary(gldA.glm5)

gldA.gam2 <- gam(GuildAbun ~ Observer + s(Days) + GuildFeed + 
                   DOM + HabClass + s(log(TreeDens)), 
                 data = m_data_A1_gld, family=nb)
summary(gldA.gam2)     


### Analysis 2: Building models for guilds #2 inv ####

inv.0 <- glmer(GuildAbun ~ Observer + (1|Point) + (1|Days),
               data = m_data_A1_inv, family=poisson())
summary(inv.0)

# adding Week and testing poly
inv.1a <- update(inv.0, .~. + Week)
inv.1b <- update(inv.0, .~. + poly(Week,2))
AICctab(inv.0, inv.1a, inv.1b)
# main effect not significant, but interaction with other predictors might be

# nuisance variables like CoD and Weather
inv.1c <- update(inv.0, .~. + CoD)
inv.1d <- update(inv.0, .~. + poly(CoD,2))
AICctab(inv.0, inv.1c, inv.1d)
anova(inv.0, inv.1c) 
inv.1e <- update(inv.0, .~. + Weather)
summary(inv.1e)
# not significant

# predictors of interest
inv.1f <- update(inv.0, .~. + HabClass)
inv.1g <- update(inv.0, .~. + Moss)
AICctab(inv.0, inv.1f, inv.1g)
anova(inv.0, inv.1f)
anova(inv.0, inv.1g)
inv.1h <- update(inv.0, .~. + HabClass*Week)
summary(inv.1h)
# habclass not useful at all. Moss is.
inv.1i <- update(inv.0, .~. + Moss*Week)
AICctab(inv.0, inv.1g, inv.1i)
summary(inv.1i)
# interaction has 1.4 less AIC, and although this is not much, will include because
# interaction is pertinent to question. 
# however, will have to compare with DOMmoss
# Week 0: 1.670883 (Moss0), 1.536729 (Moss1)
# Week 13: 1.494841 (Moss0), 2.37479 (Moss1)
# in the model without Week effect:
# 1.579613 (Moss0), 1.976406 (Moss1)
inv.1j <- update(inv.0, .~. + Moss*poly(Week,2))
summary(inv.1j)

rm(list=c("inv.1a","inv.1b","inv.1c","inv.1d","inv.1e","inv.1f","inv.1g","inv.1h","inv.1i","inv.1j"))
inv.1 <- update(inv.0, .~. + Moss*Week)
anova(inv.0, inv.1)


inv.2a <- update(inv.1, .~. + log(CCavgsd) + log(TreeDens))
inv.2b <- update(inv.1, .~. + poly(log(CCavgsd),2) + log(TreeDens))
inv.2c <- update(inv.1, .~. + log(CCavgsd) + poly(log(TreeDens),2))
inv.2d <- update(inv.1, .~. + poly(log(CCavgsd),2) + poly(log(TreeDens),2))
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d)
summary(inv.2a)
inv.2e <- update(inv.1, .~. + log(CCavgsd):Week + log(TreeDens):Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e)
# although still ns, the interaction with Week has the closest AIC to inv.1
inv.2f <- update(inv.1, .~. + poly(log(CCavgsd),2):Week + poly(log(TreeDens),2):Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f)
inv.2g <- update(inv.1, .~. + log(CCavgsd)*Week + log(TreeDens)*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f, inv.2g)
summary(inv.2e)
inv.2h <- update(inv.1, .~. + log(CCavgsd):Week + log(TreeDens))
inv.2i <- update(inv.1, .~. + log(CCavgsd) + log(TreeDens):Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f, inv.2g, inv.2h, inv.2i)
inv.2j <- update(inv.1, .~. + log(CCavgsd):Week + poly(log(TreeDens),2))
inv.2k <- update(inv.1, .~. + log(CCavgsd) + poly(log(TreeDens),2):Week)
AICctab(inv.1, inv.2a, inv.2c, inv.2e, inv.2h, inv.2i, inv.2j, inv.2k)
# not useful!

inv.2l <- update(inv.1, .~. + TPDscaled)
inv.2m <- update(inv.1, .~. + TPDscaled*Week)
inv.2n <- update(inv.1, .~. + poly(TPDscaled,2))
inv.2o <- update(inv.1, .~. + poly(TPDscaled,2)*Week)
AICctab(inv.1, inv.2l, inv.2m, inv.2n, inv.2o)
inv.2p <- update(inv.1, .~. + poly(TPDscaled,3))
AICctab(inv.1, inv.2l, inv.2n, inv.2p)
summary(inv.2l)
inv.2q <- update(inv.1, .~. + TreePropDeci)
inv.2r <- update(inv.1, .~. + poly(TreePropDeci,2))
AICctab(inv.1, inv.2l, inv.2q, inv.2n, inv.2r)
summary(inv.2q)
# only linear main effect of TPD

rm(list=c("inv.2a","inv.2b","inv.2c","inv.2d","inv.2e","inv.2f","inv.2g","inv.2h","inv.2i",
          "inv.2j","inv.2k","inv.2l","inv.2m","inv.2n","inv.2o","inv.2p","inv.2q","inv.2r"))
inv.2 <- update(inv.1, .~. + TPDscaled)


inv.3a <- update(inv.2, .~. - TPDscaled + TreeDiv)
inv.3b <- update(inv.2, .~. + TreeDiv)
AICctab(inv.2, inv.3a, inv.3b)
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
summary(inv.3i)
summary(inv.3j)
inv.3k <- update(inv.2, .~. + poly(TreeRich,2))
AICctab(inv.2, inv.3i, inv.3k)
# interesting how abundance increases with TreeRich but decreases with TPD (and poly ns)

rm(list=c("inv.3a","inv.3b","inv.3c","inv.3d","inv.3e","inv.3f","inv.3g","inv.3h","inv.3i",
          "inv.3j","inv.3k"))
inv.3 <- update(inv.2, .~. + TreeRich)


# finally, testing DOM and comparing with Moss
inv.4a <- update(inv.3, .~. + DOM)
inv.4b <- update(inv.3, .~. + DOM*Week)
inv.4c <- update(inv.3, .~. + DOM - Moss)
inv.4d <- update(inv.3, .~. + DOM*Week - Moss)
AICctab(inv.3, inv.4a, inv.4b, inv.4c, inv.4d)
summary(inv.4c)
inv.4e <- update(inv.3, .~. + DOM - Moss*Week)
AICctab(inv.3, inv.4a, inv.4c, inv.4e)
# moss:week important
inv.4f <- update(inv.3, .~. + DOM - Moss - Week)
AICctab(inv.3, inv.4a, inv.4c, inv.4f)
summary(inv.4f)
# don't know whether to keep only interaction (p-value) or to keep all (visualisation wiil
# give a good idea of the effect)
summary(inv.4a)
inv.4g <- update(inv.3, .~. + DOM - TreeRich)
AICctab(inv.3, inv.4a, inv.4g)
summary(inv.4g)
# effect of TreeRich prob. caused by STRONG effect of Carex. 
# All Carex Points have low TreeRich
inv.4h <- update(inv.3, .~. + DOM - Moss - TreeRich)
inv.4i <- update(inv.3, .~. + DOM - Moss - Week - TreeRich)
AICctab(inv.4g, inv.4h, inv.4i)
# for now, will keep main effects of moss and week but remove TreeRich. i.e., inv.4g


rm(list=c("inv.4a","inv.4b","inv.4c","inv.4d","inv.4e","inv.4f","inv.4g","inv.4h","inv.4i"))
inv.4 <- update(inv.3, .~. + DOM - TreeRich)

summary(inv.4)
hist(resid(inv.4)) # surprisingly Gaussian
inv.4sim <- simulateResiduals(inv.4, plot=T, n=1000) 
plot(inv.4sim) # here it is, the zero-inflation
testDispersion(inv.4sim) # underdispersed? but p<0.05 only cos of sample size
hist(inv.4sim) # not uniform! left-skewed
plotResiduals(inv.4sim, form=m_data_A1_inv$Week)
plotResiduals(inv.4sim, form=m_data_A1_inv$Moss)
plotResiduals(inv.4sim, form=m_data_A1_inv$TPDscaled)
plotResiduals(inv.4sim, form=m_data_A1_inv$DOM)
# Week has an issue
testZeroInflation(inv.4sim) # not seen, maybe cos no explicit zero in dataset?



# exploring convergence issue 
inv.4opt <- allFit(inv.4)
summary(inv.4opt)$which.OK # only 1/7 optimisers failed
summary(inv.4opt)$llik 
summary(inv.4opt)$fixef



library(glmmTMB)
inv.4ZI <- glmmTMB(GuildAbun ~ Observer + (1|Point) + (1|Days) + Moss + Week +  
                     TPDscaled + DOM + Moss:Week, data=m_data_A1_inv,
                   family=poisson, ziformula = ~1)
inv.4ZI2 <- glmmTMB(GuildAbun ~ Observer + (1|Point) + (1|Days) + Moss + Week +  
                      TPDscaled + DOM + Moss:Week, data=m_data_A1_inv,
                    family=poisson, ziformula = ~1+(1|Point)+(1|Days))
inv.4NB <- glmer.nb(GuildAbun ~ Observer + (1|Point) + (1|Days) + Moss + Week +  
                      TPDscaled + DOM + Moss:Week, data=m_data_A1_inv)
AICtab(inv.4, inv.4ZI, inv.4ZI2)

inv.4min <- update(inv.4, .~. - (1|Point))
summary(inv.4min)
inv.4min2 <- update(inv.4, .~. - (1|Days))
inv.4min3 <- glm(GuildAbun ~ Observer + Moss + Week + TPDscaled +  
                   DOM + Moss:Week, data=m_data_A1_inv, family=poisson)
# When Point as fixed effect, only the Points with Carex are significant. Hence, remove
summary(inv.4min2)
anova(inv.4, inv.4min, inv.4min2)
AICctab(inv.4, inv.4min, inv.4min2)


### Analysis 2: Building models for guilds #3 inv ####

# mixed model poisson error no zero-inflation (GLMER)
# (singular fit on adding DOM)
inv.0 <- glmer(GuildAbun ~ Observer + (1|Point),
               data = m_data_A1_inv, family=poisson())
summary(inv.0)

# adding Week and testing poly
inv.1a <- update(inv.0, .~. + Week)
inv.1b <- update(inv.0, .~. + poly(Week,2))
AICctab(inv.0, inv.1a, inv.1b)
# main effect not significant, but interaction with other predictors might be

# nuisance variables like CoD and Weather
inv.1c <- update(inv.0, .~. + CoD)
inv.1d <- update(inv.0, .~. + poly(CoD,2))
anova(inv.0, inv.1c) 
inv.1e <- update(inv.0, .~. + Weather)
AICctab(inv.0, inv.1c, inv.1d, inv.1e)
# not significant

# predictors of interest
inv.1f <- update(inv.0, .~. + HabClass)
inv.1g <- update(inv.0, .~. + Moss)
AICctab(inv.0, inv.1f, inv.1g)
inv.1h <- update(inv.0, .~. + HabClass*Week)
summary(inv.1h)
AICctab(inv.0, inv.1f, inv.1g, inv.1h)
# habclass not useful at all. Moss is.
inv.1i <- update(inv.0, .~. + Moss*Week)
AICctab(inv.0, inv.1g, inv.1i)
summary(inv.1i)
# interaction has 1.4 less AIC, and although this is not much, will include because
# interaction is pertinent to question. 
# however, will have to compare with DOMmoss
# Week 0: 1.670883 (Moss0), 1.536729 (Moss1)
# Week 13: 1.494841 (Moss0), 2.37479 (Moss1)
# in the model without Week effect:
# 1.579613 (Moss0), 1.976406 (Moss1)
inv.1j <- update(inv.0, .~. + Moss*poly(Week,2))
AICctab(inv.0, inv.1g, inv.1i, inv.1j)
summary(inv.1j)

rm(list=c("inv.1a","inv.1b","inv.1c","inv.1d","inv.1e","inv.1f","inv.1g","inv.1h","inv.1i","inv.1j"))
inv.1 <- update(inv.0, .~. + Moss*Week)
anova(inv.0, inv.1)


inv.2a <- update(inv.1, .~. + CCavgsdscaled + TreeDensscaled)
inv.2b <- update(inv.1, .~. + poly(CCavgsdscaled,2) + TreeDensscaled)
inv.2c <- update(inv.1, .~. + CCavgsdscaled + poly(TreeDensscaled,2))
inv.2d <- update(inv.1, .~. + poly(CCavgsdscaled,2) + poly(TreeDensscaled,2))
inv.2e <- update(inv.1, .~. + CCavgsdscaled:Week + TreeDensscaled:Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e)
# although still ns, the interaction with Week has the closest AIC to inv.1
inv.2f <- update(inv.1, .~. + poly(CCavgsdscaled,2):Week + poly(TreeDensscaled,2):Week)
inv.2g <- update(inv.1, .~. + CCavgsdscaled*Week + TreeDensscaled*Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f, inv.2g)
summary(inv.2e)
inv.2h <- update(inv.1, .~. + CCavgsdscaled:Week + TreeDensscaled)
inv.2i <- update(inv.1, .~. + CCavgsdscaled + TreeDensscaled:Week)
AICctab(inv.1, inv.2a, inv.2b, inv.2c, inv.2d, inv.2e, inv.2f, inv.2g, inv.2h, inv.2i)
inv.2j <- update(inv.1, .~. + CCavgsdscaled:Week + poly(TreeDensscaled,2))
inv.2k <- update(inv.1, .~. + CCavgsdscaled + poly(TreeDensscaled,2):Week)
AICctab(inv.1, inv.2a, inv.2c, inv.2e, inv.2h, inv.2i, inv.2j, inv.2k)
inv.2l <- update(inv.1, .~. + TreeDensscaled:Week)
inv.2m <- update(inv.1, .~. + TreeDensscaled)
inv.2n <- update(inv.1, .~. + CCavgsdscaled:Week)
inv.2o <- update(inv.1, .~. + CCavgsdscaled)
AICctab(inv.1, inv.2l, inv.2m, inv.2n, inv.2o)
# not useful!
inv.2p <- update(inv.1, .~. + TPDscaled)
inv.2q <- update(inv.1, .~. + TPDscaled*Week)
inv.2r <- update(inv.1, .~. + poly(TPDscaled,2))
inv.2s <- update(inv.1, .~. + poly(TPDscaled,2)*Week)
AICctab(inv.1, inv.2p, inv.2q, inv.2r, inv.2s)
inv.2t <- update(inv.1, .~. + TreePropDeci)
inv.2u <- update(inv.1, .~. + poly(TreePropDeci,2))
AICctab(inv.1, inv.2p, inv.2t, inv.2r, inv.2u)
summary(inv.2q)
# only linear main effect of TPD

rm(list=c("inv.2a","inv.2b","inv.2c","inv.2d","inv.2e","inv.2f","inv.2g","inv.2h","inv.2i",
          "inv.2j","inv.2k","inv.2l","inv.2m","inv.2n","inv.2o","inv.2p","inv.2q","inv.2r",
          "inv.2s","inv.2t","inv.2u"))
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
summary(inv.3i)
inv.3k <- update(inv.2, .~. + poly(TreeRich,2))
AICctab(inv.2, inv.3i, inv.3k)
# interesting how abundance increases with TreeRich but decreases with TPD (and poly ns)

rm(list=c("inv.3a","inv.3b","inv.3c","inv.3d","inv.3e","inv.3f","inv.3g","inv.3h","inv.3i",
          "inv.3j","inv.3k"))
inv.3 <- update(inv.2, .~. + TreeRich)
summary(inv.3)


# finally, testing DOM and comparing with Moss
inv.4a <- update(inv.3, .~. + DOM) # singular fit
inv.4b <- update(inv.3, .~. + DOM*Week)
inv.4c <- update(inv.3, .~. + DOM - Moss) # don't keep interaction only without main eff
inv.4d <- update(inv.3, .~. + DOM*Week - Moss)
AICctab(inv.3, inv.4a, inv.4b, inv.4c, inv.4d)
summary(inv.4c)
inv.4e <- update(inv.3, .~. + DOM - Moss*Week)
AICctab(inv.3, inv.4a, inv.4c, inv.4e)
# moss:week important
inv.4f <- update(inv.3, .~. + DOM - Moss - Week)
AICctab(inv.3, inv.4a, inv.4c, inv.4f)
summary(inv.4f)
# don't know whether to keep only interaction (p-value) or to keep all (visualisation wiil
# give a good idea of the effect)
summary(inv.4a)
inv.4g <- update(inv.3, .~. + DOM - TreeRich)
AICctab(inv.3, inv.4a, inv.4g)
summary(inv.4g)
# effect of TreeRich prob. caused by STRONG effect of Carex. 
# All Carex Points have low TreeRich
inv.4h <- update(inv.3, .~. + DOM - Moss - TreeRich)
inv.4i <- update(inv.3, .~. + DOM - Moss - Week - TreeRich)
AICctab(inv.4g, inv.4h, inv.4i)
# for now, will keep main effects of moss and week (inv.4a)
# but remove TreeRich (inv.4g)


rm(list=c("inv.4a","inv.4b","inv.4c","inv.4d","inv.4e","inv.4f","inv.4g","inv.4h","inv.4i"))
inv.4 <- update(inv.3, .~. + DOM - TreeRich)
summary(inv.4)


# mixed model poisson error no zero-inflation (TMB, not GLMER)
# (no singular fit like inv.4 with glmer, no convergence issue)
# (same issues as zero-inflated model, invTMB.4)
invTMB.4poi <- update(invTMB.4, ziformula=~0)
summary(invTMB.4poi)
hist(resid(invTMB.4poi))
invTMB.4poisim <- simulateResiduals(invTMB.4poi, n=322)
plot(invTMB.4poisim)
hist(invTMB.4poisim)
testDispersion(invTMB.4poisim, alternative = "less") # underdispersion not great
testZeroInflation(invTMB.4poisim, alternative = "less") 


# mixed model poisson error with zero inflation
# (no singular fit, no convergence issue)
# (histogram of residuals not uniform, but slightly left skewed; QQ slight underdisp but 
# cos of sample size; residuals vs. predicted slightly increasing to the right;
# zero inflation is 0)
invTMB.0 <- glmmTMB(GuildAbun ~ Observer + (1|Point),
                    data = m_data_A1_inv, family=poisson(), ziformula = ~1)
summary(invTMB.0)

# adding Week and testing poly
invTMB.1a <- update(invTMB.0, .~. + Week)
invTMB.1b <- update(invTMB.0, .~. + poly(Week,2))
AICctab(invTMB.0, invTMB.1a, invTMB.1b)
# main effect not significant, but interaction with other predictors might be

# nuisance variables like CoD and Weather
invTMB.1c <- update(invTMB.0, .~. + CoD)
invTMB.1d <- update(invTMB.0, .~. + poly(CoD,2))
invTMB.1e <- update(invTMB.0, .~. + Weather)
AICctab(invTMB.0, invTMB.1c, invTMB.1d, invTMB.1e)
# not significant

# predictors of interest
invTMB.1f <- update(invTMB.0, .~. + HabClass)
invTMB.1g <- update(invTMB.0, .~. + Moss)
AICctab(invTMB.0, invTMB.1f, invTMB.1g)
invTMB.1h <- update(invTMB.0, .~. + HabClass*Week)
AICctab(invTMB.0, invTMB.1f, invTMB.1g, invTMB.1h)
# habclass not useful at all. Moss is.
invTMB.1i <- update(invTMB.0, .~. + Moss*Week)
AICctab(invTMB.0, invTMB.1g, invTMB.1i)
# interaction has 1.4 less AIC, and although this is not much, will include because
# interaction is pertinent to question. 
# however, will have to compare with DOMmoss
# Week 0: 1.670883 (Moss0), 1.536729 (Moss1)
# Week 13: 1.494841 (Moss0), 2.37479 (Moss1)
# in the model without Week effect:
# 1.579613 (Moss0), 1.976406 (Moss1)
invTMB.1j <- update(invTMB.0, .~. + Moss*poly(Week,2))
AICctab(invTMB.0, invTMB.1g, invTMB.1i, invTMB.1j)


rm(list=c("invTMB.1a","invTMB.1b","invTMB.1c","invTMB.1d","invTMB.1e","invTMB.1f","invTMB.1g","invTMB.1h","invTMB.1i","invTMB.1j"))
invTMB.1 <- update(invTMB.0, .~. + Moss*Week)
anova(invTMB.0, invTMB.1)


invTMB.2a <- update(invTMB.1, .~. + CCavgsdscaled + TreeDensscaled)
invTMB.2b <- update(invTMB.1, .~. + poly(CCavgsdscaled,2) + TreeDensscaled)
invTMB.2c <- update(invTMB.1, .~. + CCavgsdscaled + poly(TreeDensscaled,2))
invTMB.2d <- update(invTMB.1, .~. + poly(CCavgsdscaled,2) + poly(TreeDensscaled,2))
invTMB.2e <- update(invTMB.1, .~. + CCavgsdscaled:Week + TreeDensscaled:Week)
AICctab(invTMB.1, invTMB.2a, invTMB.2b, invTMB.2c, invTMB.2d, invTMB.2e)
# although still ns, the interaction with Week has the closest AIC to invTMB.1
invTMB.2f <- update(invTMB.1, .~. + poly(CCavgsdscaled,2):Week + poly(TreeDensscaled,2):Week)
invTMB.2g <- update(invTMB.1, .~. + CCavgsdscaled*Week + TreeDensscaled*Week)
AICctab(invTMB.1, invTMB.2a, invTMB.2b, invTMB.2c, invTMB.2d, invTMB.2e, invTMB.2f, invTMB.2g)
summary(invTMB.2e)

invTMB.2h <- update(invTMB.1, .~. + CCavgsdscaled:Week + TreeDensscaled)
invTMB.2i <- update(invTMB.1, .~. + CCavgsdscaled + TreeDensscaled:Week)
AICctab(invTMB.1, invTMB.2a, invTMB.2b, invTMB.2c, invTMB.2d, invTMB.2e, invTMB.2f, invTMB.2g, invTMB.2h, invTMB.2i)

invTMB.2j <- update(invTMB.1, .~. + CCavgsdscaled:Week + poly(TreeDensscaled,2))
invTMB.2k <- update(invTMB.1, .~. + CCavgsdscaled + poly(TreeDensscaled,2):Week)
AICctab(invTMB.1, invTMB.2a, invTMB.2c, invTMB.2e, invTMB.2h, invTMB.2i, invTMB.2j, invTMB.2k)

invTMB.2l <- update(invTMB.1, .~. + TreeDensscaled:Week)
invTMB.2m <- update(invTMB.1, .~. + TreeDensscaled)
invTMB.2n <- update(invTMB.1, .~. + CCavgsdscaled:Week)
invTMB.2o <- update(invTMB.1, .~. + CCavgsdscaled)
AICctab(invTMB.1, invTMB.2l, invTMB.2m, invTMB.2n, invTMB.2o)
# not useful!
invTMB.2p <- update(invTMB.1, .~. + TPDscaled)
invTMB.2q <- update(invTMB.1, .~. + TPDscaled*Week)
invTMB.2r <- update(invTMB.1, .~. + poly(TPDscaled,2))
invTMB.2s <- update(invTMB.1, .~. + poly(TPDscaled,2)*Week)
AICctab(invTMB.1, invTMB.2p, invTMB.2q, invTMB.2r, invTMB.2s)

invTMB.2t <- update(invTMB.1, .~. + TreePropDeci)
invTMB.2u <- update(invTMB.1, .~. + poly(TreePropDeci,2))
AICctab(invTMB.1, invTMB.2p, invTMB.2t, invTMB.2r, invTMB.2u)
# only linear main effect of TPD

rm(list=c("invTMB.2a","invTMB.2b","invTMB.2c","invTMB.2d","invTMB.2e","invTMB.2f","invTMB.2g","invTMB.2h","invTMB.2i",
          "invTMB.2j","invTMB.2k","invTMB.2l","invTMB.2m","invTMB.2n","invTMB.2o","invTMB.2p","invTMB.2q","invTMB.2r",
          "invTMB.2s","invTMB.2t","invTMB.2u"))
invTMB.2 <- update(invTMB.1, .~. + TPDscaled)


invTMB.3a <- update(invTMB.2, .~. - TPDscaled + TreeDiv)
invTMB.3b <- update(invTMB.2, .~. + TreeDiv)
AICctab(invTMB.2, invTMB.3a, invTMB.3b)
invTMB.3c <- update(invTMB.2, .~. + poly(TreeDiv,2))
invTMB.3d <- update(invTMB.2, .~. + TreeDiv*Week)
AICctab(invTMB.2, invTMB.3a, invTMB.3b, invTMB.3c, invTMB.3d)
# TreeDiv not important

invTMB.3e <- update(invTMB.2, .~. + sqrt(UndDens))
invTMB.3f <- update(invTMB.2, .~. + poly(sqrt(UndDens),2))
invTMB.3g <- update(invTMB.2, .~. + sqrt(UndDens)*Week)
invTMB.3h <- update(invTMB.2, .~. + poly(sqrt(UndDens),2)*Week)
AICctab(invTMB.2, invTMB.3e, invTMB.3f, invTMB.3g, invTMB.3h)
# UndDens not important
invTMB.3i <- update(invTMB.2, .~. + TreeRich)
invTMB.3j <- update(invTMB.2, .~. + TreeRich*Week)
AICctab(invTMB.2, invTMB.3i, invTMB.3j)
summary(invTMB.3i)
invTMB.3k <- update(invTMB.2, .~. + poly(TreeRich,2))
AICctab(invTMB.2, invTMB.3i, invTMB.3k)
# interesting how abundance increases with TreeRich but decreases with TPD (and poly ns)

rm(list=c("invTMB.3a","invTMB.3b","invTMB.3c","invTMB.3d","invTMB.3e","invTMB.3f",
          "invTMB.3g","invTMB.3h","invTMB.3i", "invTMB.3j","invTMB.3k"))
invTMB.3 <- update(invTMB.2, .~. + TreeRich)
summary(invTMB.3)


# finally, testing DOM and comparing with Moss
invTMB.4a <- update(invTMB.3, .~. + DOM) # singular fit
invTMB.4b <- update(invTMB.3, .~. + DOM*Week)
invTMB.4c <- update(invTMB.3, .~. + DOM - Moss) # don't keep interaction only without main eff
invTMB.4d <- update(invTMB.3, .~. + DOM*Week - Moss)
AICctab(invTMB.3, invTMB.4a, invTMB.4b, invTMB.4c, invTMB.4d)
summary(invTMB.4c)
invTMB.4e <- update(invTMB.3, .~. + DOM - Moss*Week)
AICctab(invTMB.3, invTMB.4a, invTMB.4c, invTMB.4e)
# moss:week important
invTMB.4f <- update(invTMB.3, .~. + DOM - Moss - Week)
AICctab(invTMB.3, invTMB.4a, invTMB.4c, invTMB.4f)
summary(invTMB.4f)
# don't know whether to keep only interaction (p-value) or to keep all (visualisation will
# give a good idea of the effect)
summary(invTMB.4a)
invTMB.4g <- update(invTMB.3, .~. + DOM - TreeRich)
AICctab(invTMB.3, invTMB.4a, invTMB.4g)
summary(invTMB.4g)
# effect of TreeRich prob. caused by STRONG effect of Carex. 
# All Carex Points have low TreeRich
invTMB.4h <- update(invTMB.3, .~. + DOM - Moss - TreeRich)
invTMB.4i <- update(invTMB.3, .~. + DOM - Moss - Week - TreeRich)
AICctab(invTMB.4g, invTMB.4h, invTMB.4i)
# for now, will keep main effects of moss and week (invTMB.4a)
# but remove TreeRich (invTMB.4g)

rm(list=c("invTMB.4a","invTMB.4b","invTMB.4c","invTMB.4d","invTMB.4e","invTMB.4f",
          "invTMB.4g","invTMB.4h","invTMB.4i"))
invTMB.4 <- update(invTMB.3, .~. + DOM - TreeRich)
summary(invTMB.4)

hist(resid(invTMB.4))
invTMB.4sim <- simulateResiduals(invTMB.4, n=322)
plot(invTMB.4sim)
hist(invTMB.4sim)
testDispersion(invTMB.4sim, alternative = "less") # underdispersion not great
plotResiduals(invTMB.4sim, form=m_data_A1_inv$Week)
plotResiduals(invTMB.4sim, form=m_data_A1_inv$Moss)
plotResiduals(invTMB.4sim, form=m_data_A1_inv$TPDscaled)
plotResiduals(invTMB.4sim, form=m_data_A1_inv$DOM)
testZeroInflation(invTMB.4sim, alternative = "less") # not seen, maybe cos no explicit zero in dataset?



# GLM with poisson error, without random effect (Point as fixed)
# (produces singularities and NAs in predictor estimates)
inv.4min2 <- glm(GuildAbun ~ Observer + Point + Moss + Week + TPDscaled +  
                   DOM + Moss:Week, data=m_data_A1_inv, family=poisson)
AICctab(inv.3, inv.4, inv.4min, inv.4min2)
summary(inv.4min2)


# GLM with poisson error, without random effect 
# (no fitting issue)
# (slight underdispersion, not uniform, slight left skew in histogram, zero inflation 0)
inv.0min <- glm(GuildAbun ~ Observer, data = m_data_A1_inv, family=poisson())
inv.1min <- update(inv.0min, .~. + Moss*Week)
inv.2min <- update(inv.1min, .~. + TPDscaled)
inv.3min <- update(inv.2min, .~. + TreeRich)
inv.3minsim <- simulateResiduals(inv.3min)
plot(inv.3minsim) # same issue
inv.4min <- update(inv.3min, .~. + DOM - TreeRich)
summary(inv.4min)
AICctab(inv.0min, inv.1min, inv.2min, inv.3min, inv.4min)

hist(resid(inv.4min)) 
inv.4minsim <- simulateResiduals(inv.4min, n=1000) 
plot(inv.4minsim) 
testDispersion(inv.4minsim) 
hist(inv.4minsim) 
testZeroInflation(inv.4sim) 


# inv.4minnorm <- update(inv.4min, family=gaussian())
# inv.4minnormsim <- simulateResiduals(inv.4minnorm)
# plot(inv.4minnormsim)
# so it is really is not normal. must be underdispersed


# GLM with negbin error, without random effect
# (theta too high, too many iterations)
inv.4minNB <- glm.nb(GuildAbun ~ Observer + Moss + Week + TPDscaled +  
                       DOM + Moss:Week, data=m_data_A1_inv)
summary(inv.4minNB)
plot(inv.4minNB)


# mixed model negbin error no zero-inflation
# (convergence issue: non-positive definite Hessian matrix & false convergence)
invTMB.4NB <- update(invTMB.4poi, family=nbinom2)
summary(invTMB.4NB) # large overdispersion
hist(resid(invTMB.4NB)) 
invTMB.4NBsim <- simulateResiduals(invTMB.4NB, n=322) 
plot(invTMB.4NBsim) 
testDispersion(invTMB.4NBsim, alternative = "less") # underdispersion not great
hist(invTMB.4NBsim) # not uniform! left-skewed
plotResiduals(invTMB.4NBsim, form=m_data_A1_inv$Week)
plotResiduals(invTMB.4NBsim, form=m_data_A1_inv$Moss)
plotResiduals(invTMB.4NBsim, form=m_data_A1_inv$TPDscaled)
plotResiduals(invTMB.4NBsim, form=m_data_A1_inv$DOM)
# Week has an issue
testZeroInflation(invTMB.4NBsim, alternative = "less") 


# mixed model negbin error 
inv.NB0 <- glmer.nb(GuildAbun ~ Observer + (1|Point), data = m_data_A1_inv)
inv.NB0 <- glmmTMB(GuildAbun ~ Observer + (1|Point), data = m_data_A1_inv, family=nbinom2)
inv.NB1 <- update(inv.NB0, .~. + Moss*Week) # false convergence
### ### ### ###





### ###


### Analysis 2: Building models for guilds #4 omn ####

omn.0 <- glmer(GuildAbun ~ Observer + (1|Point), data=m_data_A1_omn, family=poisson)
omn.0sim <- simulateResiduals(omn.0, plot=T)

omn.0NB <- glmer.nb(GuildAbun ~ Observer + (1|Point), data=m_data_A1_omn)
omn.0NBsim <- simulateResiduals(omn.0NB, plot=T)

# going with nb for now

# adding Week and testing poly
omn.1a <- update(omn.0NB, .~. + Week)
omn.1b <- update(omn.0NB, .~. + poly(Week,2))
AICctab(omn.0, omn.1a, omn.1b)
summary(omn.1a)
# linear effect of Week HIGHLY significant, 123 dAICc
rm(list=c("omn.1a","omn.1b"))
omn.1 <- update(omn.0NB, .~. + Week)

# nuisance variables like CoD and Weather
omn.2a <- update(omn.1, .~. + CoD)
omn.2b <- update(omn.1, .~. + poly(CoD,2))
AICctab(omn.1, omn.2a, omn.2b) 
omn.2c <- update(omn.1, .~. + Weather)
AICctab(omn.1, omn.2a, omn.2b, omn.2c)
# not significant
# predictors of interest
omn.2d <- update(omn.1, .~. + HabClass)
omn.2e <- update(omn.1, .~. + Moss)
AICctab(omn.1, omn.2d, omn.2e) # habclass
omn.2f <- update(omn.1, .~. + HabClass + Moss)
AICctab(omn.1, omn.2d, omn.2e, omn.2f) # Moss ns
omn.2g <- update(omn.1, .~. + HabClass*Week)
omn.2h <- update(omn.1, .~. + Moss*Week)
AICctab(omn.1, omn.2d, omn.2g, omn.2e, omn.2h)
# only main effect of HabClass

rm(list=c("omn.2a","omn.2b","omn.2c","omn.2d","omn.2e","omn.2f","omn.2g","omn.2h"))
omn.2 <- update(omn.1, .~. + HabClass)
anova(omn.1, omn.2)

omn.3a <- update(omn.2, .~. + CCavgsdscaled + TreeDensscaled)
omn.3b <- update(omn.2, .~. + poly(CCavgsdscaled,2) + TreeDensscaled)
omn.3c <- update(omn.2, .~. + CCavgsdscaled + poly(TreeDensscaled,2))
omn.3d <- update(omn.2, .~. + poly(CCavgsdscaled,2) + poly(TreeDensscaled,2))
omn.3e <- update(omn.2, .~. + CCavgsdscaled*Week + TreeDensscaled*Week)
AICctab(omn.2, omn.3a, omn.3b, omn.3c, omn.3d, omn.3e)
summary(omn.3a)
# singular fit, but CC and TD significant
omn.3f <- update(omn.2, .~. + CCavgsdscaled + TreeDensscaled - HabClass) # no singular fit
AICctab(omn.2, omn.3a, omn.3f)
summary(omn.3f)
omn.3g <- update(omn.2, .~. + poly(CCavgsdscaled,2):Week + poly(TreeDensscaled,2):Week)
AICctab(omn.2, omn.3a, omn.3b, omn.3c, omn.3d, omn.3e, omn.3f, omn.3g)
# adding CC and TD after HabClass is very good, but singular fit. 
# next best is HabClass without CC and TD
omn.3h <- update(omn.2, .~. + TreeDensscaled:Week)
omn.3i <- update(omn.2, .~. + TreeDensscaled)
omn.3j <- update(omn.2, .~. + CCavgsdscaled:Week) # better than CC + TD
omn.3k <- update(omn.2, .~. + CCavgsdscaled)
AICctab(omn.2, omn.3a, omn.3f, omn.3h, omn.3i, omn.3j, omn.3k)
omn.3l <- update(omn.2, .~. + CCavgsdscaled:Week - HabClass)
AICctab(omn.2, omn.3a, omn.3f, omn.3j, omn.3l)
# still prefer HabClass because of singular fit

omn.3m <- update(omn.2, .~. + TPDscaled) 
omn.3n <- update(omn.2, .~. + TPDscaled*Week)
omn.3o <- update(omn.2, .~. + poly(TPDscaled,2))
omn.3p <- update(omn.2, .~. + poly(TPDscaled,2)*Week)
AICctab(omn.2, omn.3m, omn.3n, omn.3o, omn.3p)
# no singular fit but not significant

omn.3q <- update(omn.2, .~. + TreeDiv)
omn.3r <- update(omn.2, .~. + poly(TreeDiv,2))
omn.3s <- update(omn.2, .~. + TreeDiv*Week)
omn.3t <- update(omn.2, .~. + poly(TreeDiv,2)*Week)
AICctab(omn.2, omn.3q, omn.3r, omn.3s, omn.3t)
# no singular fit but not significant

omn.3u <- update(omn.2, .~. + UndDensscaled)
omn.3v <- update(omn.2, .~. + poly(UndDensscaled,2))
omn.3w <- update(omn.2, .~. + UndDensscaled*Week)
omn.3x <- update(omn.2, .~. + poly(UndDensscaled,2)*Week)
AICctab(omn.2, omn.3u, omn.3v, omn.3w, omn.3x)
# no singular fit but not significant

omn.3y <- update(omn.2, .~. + TreeRich)
omn.3z <- update(omn.2, .~. + TreeRich*Week)
AICctab(omn.2, omn.3y, omn.3z)
omn.3z <- update(omn.2, .~. + poly(TreeRich,2))
# no singular fit but not significant

rm(list=c("omn.3a","omn.3b","omn.3c","omn.3d","omn.3e","omn.3f","omn.3g","omn.3h","omn.3i",
          "omn.3j","omn.3k","omn.3l","omn.3m","omn.3n","omn.3o","omn.3p","omn.3q","omn.3r",
          "omn.3s","omn.3t","omn.3u","omn.3v","omn.3w","omn.3x","omn.3y","omn.3z"))

omn.3 <- update(omn.2, .~. + CCavgsdscaled + TreeDensscaled)
summary(omn.3)


# finally, testing DOM and comparing with Moss
omn.4a <- update(omn.3, .~. + DOM) # singular fit
omn.4b <- update(omn.3, .~. + DOM*Week)
omn.4c <- update(omn.3, .~. + DOM - HabClass) # no singular fit
omn.4d <- update(omn.3, .~. + DOM*Week - HabClass)
AICctab(omn.3, omn.4a, omn.4b, omn.4c, omn.4d)
# main effect of DOM could have been added if not for singular fit


rm(list=c("omn.4a","omn.4b","omn.4c","omn.4d"))
omn.4 <- update(omn.3, .~. + DOM - CCavgsdscaled) # without CC better
summary(omn.4)
omn.4sim <- simulateResiduals(omn.4, plot=T) # QQ and residuals fine
hist(omn.4sim) # but hist of residuals very slight left skew, but can see uniform



## trying with GLM to remove singular fit issue ##

omnGLM.0 <-  glm.nb(GuildAbun ~ Observer, data=m_data_A1_omn)
omnGLM.1 <- update(omnGLM.0, .~. + Week)
summary(omnGLM.1)
omnGLM.2 <- update(omnGLM.1, .~. + HabClass)
summary(omnGLM.2)


omnGLM.3a <- update(omnGLM.2, .~. + CCavgsdscaled)
omnGLM.3b <- update(omnGLM.2, .~. + TreeDensscaled)
omnGLM.3c <- update(omnGLM.2, .~. + CCavgsdscaled + TreeDensscaled)
AICctab(omnGLM.2, omnGLM.3a, omnGLM.3b, omnGLM.3c) # both
omnGLM.3d <- update(omnGLM.2, .~. + CCavgsdscaled + poly(TreeDensscaled,2))
omnGLM.3e <- update(omnGLM.2, .~. + TreeDensscaled + poly(CCavgsdscaled,2))
omnGLM.3f <- update(omnGLM.2, .~. + poly(CCavgsdscaled,2) + poly(TreeDensscaled,2))
AICctab(omnGLM.2, omnGLM.3c, omnGLM.3d, omnGLM.3e, omnGLM.3f) # linear better
summary(omnGLM.3c)
omnGLM.3g <- update(omnGLM.2, .~. + Week*(CCavgsdscaled + TreeDensscaled))
omnGLM.3h <- update(omnGLM.2, .~. + Week*CCavgsdscaled + TreeDensscaled)
omnGLM.3i <- update(omnGLM.2, .~. + CCavgsdscaled + Week*TreeDensscaled)
AICctab(omnGLM.2, omnGLM.3c, omnGLM.3g, omnGLM.3h, omnGLM.3i) 
# interaction is not much better than main effect
rm(list=c("omnGLM.3a","omnGLM.3b","omnGLM.3c","omnGLM.3d","omnGLM.3e","omnGLM.3f","omnGLM.3g"
          ,"omnGLM.3h", "omnGLM.3i"))
omnGLM.3 <- update(omnGLM.2, .~. + CCavgsdscaled + TreeDensscaled)
summary(omnGLM.3)


omnGLM.4a <- update(omnGLM.3, .~. + TPDscaled)
omnGLM.4b <- update(omnGLM.3, .~. + poly(TPDscaled,2))
omnGLM.4c <- update(omnGLM.3, .~. + TPDscaled*Week)
omnGLM.4d <- update(omnGLM.3, .~. + poly(TPDscaled,2)*Week)
AICctab(omnGLM.3, omnGLM.4a, omnGLM.4b, omnGLM.4c, omnGLM.4d) # not important
omnGLM.4e <- update(omnGLM.3, .~. + TreeDiv)
omnGLM.4f <- update(omnGLM.3, .~. + poly(TreeDiv,2))
omnGLM.4g <- update(omnGLM.3, .~. + TreeDiv*Week)
omnGLM.4h <- update(omnGLM.3, .~. + poly(TreeDiv,2)*Week)
AICctab(omnGLM.3, omnGLM.4e, omnGLM.4f, omnGLM.4g, omnGLM.4h) # not important
omnGLM.4i <- update(omnGLM.3, .~. + UndDensscaled)
omnGLM.4j <- update(omnGLM.3, .~. + poly(UndDensscaled,2))
omnGLM.4k <- update(omnGLM.3, .~. + UndDensscaled*Week)
omnGLM.4l <- update(omnGLM.3, .~. + poly(UndDensscaled,2)*Week)
AICctab(omnGLM.3, omnGLM.4i, omnGLM.4j, omnGLM.4k, omnGLM.4l) # not important
omnGLM.4m <- update(omnGLM.3, .~. + TreeRich)
omnGLM.4n <- update(omnGLM.3, .~. + poly(TreeRich,2))
omnGLM.4o <- update(omnGLM.3, .~. + TreeRich*Week)
omnGLM.4p <- update(omnGLM.3, .~. + poly(TreeRich,2)*Week)
AICctab(omnGLM.3, omnGLM.4m, omnGLM.4n, omnGLM.4o, omnGLM.4p) # not important
omnGLM.4q <- update(omnGLM.3, .~. + DOM)
omnGLM.4r <- update(omnGLM.3, .~. + DOM*Week)
AICctab(omnGLM.3, omnGLM.4q, omnGLM.4r) # dAICc 3.5 with main effect, but df 4
summary(omnGLM.4q) # CC now ns
omnGLM.4s <- update(omnGLM.3, .~. + DOM - CCavgsdscaled)
AICctab(omnGLM.3, omnGLM.4q, omnGLM.4s) # better without CC
summary(omnGLM.4s)

rm(list=c("omnGLM.4a","omnGLM.4b","omnGLM.4c","omnGLM.4d","omnGLM.4e","omnGLM.4f","omnGLM.4g",
          "omnGLM.4h","omnGLM.4i","omnGLM.4j","omnGLM.4k","omnGLM.4l","omnGLM.4m","omnGLM.4n",
          "omnGLM.4o","omnGLM.4p","omnGLM.4q","omnGLM.4r","omnGLM.4s"))
omnGLM.4 <- update(omnGLM.3, .~. + DOM - CCavgsdscaled)
AICctab(omnGLM.0, omnGLM.1, omnGLM.2, omnGLM.3, omnGLM.4)
omnGLM.4sim <- simulateResiduals(omnGLM.4, plot=T) # one quantile a bit different but overall ns
hist(omnGLM.4sim) # much better than mixed model

summary(omnGLM.4)
AICctab(omn.2, omnGLM.2)
anova(omn.2, omnGLM.2)
# to see if estimates of fixed effects are same without Point
summary(omn.2)
summary(omnGLM.2)
AICctab(omn.2, omnGLM.2, omnGLM.3, omnGLM.4)


# probably best to go with this model, without Point
plot(allEffects(omnGLM.4))

### Analysis 3: Building models for EGT ####

### GLM using Abundance

egtA.glm1 <- glm(Abundance ~ Period, family=poisson, data=m_data_A1_EGT)
summary(egtA.glm1)

egtA.glm2 <- step(egtA.glm1, scope = ~Period+log(CCavgsd)+log(TreeDens)+TreePropDeci+
                    DOM+Moss+HabClass+sqrt(UndDens))
summary(egtA.glm2) # 18% explained variation

egtA.glm3 <- update(egtA.glm2, . ~ . + log(TreeDens):Period)
summary(egtA.glm3)
anova(egtA.glm2, egtA.glm3, test="Chisq")
drop1(egtA.glm3) # interaction not significant

egtA.glm4 <- glm(Abundance ~ log(TreeDens) + log(TreeDens):Period, 
                 family=poisson, data=m_data_A1_EGT) # without main effect of Period
summary(egtA.glm4)
anova(egtA.glm4, egtA.glm2, test = "Chisq")

# Go with egtA.glm2 ###

# visualising the model and comparing observe vs. predicted values
# 
# EGT.predlist <- list(Period = c("1", "2", "3"),
#                  TreeDens = seq(1, 15, by=0.01))
# EGT.predgrid <- expand.grid(EGT.predlist)
# EGT.pred <- predict(egtA.glm2, EGT.predgrid, type="response", se.fit = T)
# 
# EGT.predgrid$Abundance <- EGT.pred$fit
# EGT.predgrid$Abundance.se <- EGT.pred$se.fit
# plot(EGT.predgrid$Abundance ~ EGT.predgrid$TreeDens)
# 
# ggplot(aes(x=TreeDens, y=Abundance, col=Period), data=m_data_A1_EGT) + 
#   theme_bw() + geom_point(size=3) + 
#   geom_line(data=EGT.predgrid, size=3) 


### GAM using Abundance


egtA.gam1 <- gam(Abundance ~ Period + s(log(TreeDens)) + s(log(CCavgsd)) + 
                   s(TreePropDeci) + DOM + Moss + HabClass + s(sqrt(UndDens)), 
                 family=poisson, data=m_data_A1_EGT, select=T)
summary(egtA.gam1)

egtA.gam2 <- update(egtA.gam1, .~ Period + s(log(TreeDens)))
summary(egtA.gam2) # 16% explained deviance, but smooth term is essentially linear,
# and the GLM explains more than GAM.

layout(matrix(1:2, ncol = 2))
plot(egtA.gam2, shade=T, seWithMean = T)
acf(egtA.gam2$residuals, main="ACF")
pacf(egtA.gam2$residuals, main="pACF")
layout(1)
# this suggests that there is no significant autocorrelation
# tried using corAR1 of lag 1, but didn't change anything so not useful

egtA.gam3 <- update(egtA.gam2, .~. + s(log(TreeDens), by=Period, k=3))
summary(egtA.gam3) # 18.1% explained deviance

egtA.gam4 <- update(egtA.gam3, .~ Period + s(log(TreeDens), by=Period, k=3))
summary(egtA.gam4)
AIC(egtA.gam2, egtA.gam3, egtA.gam4)
# main effect of log(TreeDens) plus its interaction with Period seems better!
# Go with egtA.gam3 ###

layout(matrix(1:4, ncol = 2))
gam.check(egtA.gam3, rep=500)
plot.gam(egtA.gam3, shade=T, seWithMean = T, pages=1)


AIC(egtA.gam3, egtA.glm2) # GAM has lower AIC and df is not too compromised

# visualisation 
# EGT.gampred <- data.frame(predict.gam(egtA.gam5, se.fit = T))
# EGT.gampred <- transform(EGT.gampred, upper=fit+2*se.fit, lower=fit-2*se.fit)
# EGT.gampredt <- cbind(m_data_A1_EGT, EGT.gampred)
# 
# ggplot(EGT.gampredt, aes(x=log(TreeDens), y=fit, col=Period)) + theme_bw() +
#   scale_colour_viridis_d() +
#   geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey", alpha=0.2) +
#   geom_line() + geom_point(aes(y=log(Abundance), col=Period))


# testing the non-focal predictors like Observer
egtA.gam5 <- update(egtA.gam3, .~. + Observer + Weather + Wind + Visibility, 
                    select=T)

egtA.gam6 <- update(egtA.gam3, .~. + Observer + Weather)
summary(egtA.gam6)

egtA.gam7 <- update(egtA.gam6, .~. + s(TreePropDeci), select=T)
# s(TreePropDeci) not significant in the approximate significance test in the summary
# function, but deviance explained increases.

AIC(egtA.gam3, egtA.gam6, egtA.gam7)
AIC(egtA.glm2, egtA.gam6)
# Go with egtA.gam6 ###

# testing for quasipoisson
egtA.gam6qp <- update(egtA.gam6, family=quasipoisson)
summary(egtA.gam6qp) # when Observer removed, only 6% less deviance explained


### GLM using RelAbun

egtRA.glm1 <- glm(RelAbun ~ Period, family=Gamma(link="inverse"), 
                  data=m_data_A1_EGT)

egtRA.glm2 <- step(egtRA.glm1, scope = ~Period+log(CCavgsd)+log(TreeDens)+TreePropDeci+
                     DOM+Moss+HabClass+sqrt(UndDens))
# AIC after selection at each step is lower than the AIC predicted in previous step
summary(egtRA.glm2) # 20.27% explained deviance

egtRA.glm3 <- step(egtRA.glm2, scope = ~ . + Period:(sqrt(UndDens)+log(CCavgsd)))
summary(egtRA.glm3)
anova(egtRA.glm3, egtRA.glm2, test="Chisq")
AIC(egtRA.glm3, egtRA.glm2) # 23.5% explained deviance
# Go for egtRA.glm3 ###


### GAM using RelAbun

egtRA.gam1 <- gam(RelAbun ~ Period, family = Gamma(link="inverse"),
                  data=m_data_A1_EGT)

egtRA.gam2 <- update(egtRA.gam1, .~. + s(log(TreeDens)) + s(log(CCavgsd)) + 
                       s(TreePropDeci) + DOM + Moss + HabClass + s(sqrt(UndDens)), 
                     select=T)
summary(egtRA.gam2)

egtRA.gam3 <- update(egtRA.gam2, .~ Period + DOM + s(log(CCavgsd)) + s(TreePropDeci) +
                       s(sqrt(UndDens)))
summary(egtRA.gam3)


egtRA.gam4 <- gam(RelAbun ~ Period + s(log(CCavgsd)) + s(sqrt(UndDens)), 
                  family = Gamma(link="inverse"), data=m_data_A1_EGT)
# only main effects
summary(egtRA.gam4)
plot.gam(egtRA.gam4, shade=T, seWithMean = T, pages=1) # explains more

egtRA.gam5 <- update(egtRA.gam4, .~Period + s(log(CCavgsd), by=Period) +
                       s(sqrt(UndDens), by=Period), select=T)
# with interactions
summary(egtRA.gam5)


AIC(egtRA.gam4, egtRA.gam5)
AIC(egtRA.glm3, egtRA.gam5)
# Go with egtRA.gam5 ###

plot.gam(egtRA.gam5, shade=T, seWithMean = T, pages=1)

# 2.2% variation explained by Observer
egtRA.gam6 <- update(egtRA.gam5, .~. + Observer) # warnings about NAs
summary(egtRA.gam6)
AIC(egtRA.gam5, egtRA.gam6)
# Go with egtRA.gam6 (with observer effect) ###
plot.gam(egtRA.gam6, shade=T, seWithMean = T, pages=1)


### Analysis 1: Building models with Days rather than Period ####

DallA.gam1 <- gam(BirdAbun ~ Observer + Weather + s(Week), 
                  family=quasipoisson, data=m_data_A1_all)
summary(DallA.gam1)
plot(DallA.gam1, seWithMean = T, shade=T, pages=1, residuals=T)


DallA.gam2 <- gam(BirdAbun ~ Observer + s(log(CCavgsd)) + DOM + HabClass + 
                    s(sqrt(UndDens)) +  s(Week),
                  family=quasipoisson, data=m_data_A1_all)
summary(DallA.gam2)

anova(allA.gam4qp2, DallA.gam2, test="F")
# using continuous time may not be very useful

DallRA.gam1 <- gam(BirdRelAbun ~ Observer + DOM + HabClass + s(sqrt(UndDens)) + 
                     s(Week) + ti(log(CCavgsd), Week) + ti(log(TreeDens), Week),
                   family=quasipoisson, data=m_data_A1_all)
summary(DallRA.gam1)
anova(allRA.gam5, DallRA.gam1, test="F")
# don't know how resid. dev. dropped from 113 to 3



#  Analysis 2 ####
# using PCA loading scores as variables in models ###

# can try with all habvar and also with reduced hab data 

h_data_A2.1 <- h_data_A1 %>% select(-c(12,15)) %>% column_to_rownames("Point")
h_data_A2.2 <- h_data_A2.1 %>% select(-c(1:3))

h_spec <- read.delim("clipboard", row.names = 1) # A:N of "PointTreeSpec" sheet

biplot(prcomp(h_data_A2.1, scale. = T))
biplot(prcomp(h_data_A2.2, scale. = T))
biplot(prcomp(h_spec))

# decided not to do it for thesis 
### Analysis 3: Ordinations of species ####

## 2 weeks filtered species mean ##
ordW1 <- rda(b_data_ordW1, scale=T)
summary(ordW1)
summary(eigenvals(ordW1))
plot(ordW1, scaling="species")

ordW13 <- rda(b_data_ordW13, scale=T)
summary(ordW13)
summary(eigenvals(ordW13))
plot(ordW13, scaling="species")


set.seed(420)
ordgld <- b_speccode %>% select(c(1,4)) # for grouping in plots

ordW1.si <- inner_join(rownames_to_column(data.frame(ordW1$CA$u), "Point"),
                       h_data_A1, by="Point") # site scores
ordW1.sp <- inner_join(rownames_to_column(data.frame(ordW1$CA$v), "Spec_code"), 
                       ordgld, by="Spec_code") # species scores
ordW1.ef <- rownames_to_column(data.frame(scores(
  envfit(ordW1, h_data_ord, choices=1:2, scaling="species", 
         permutations = 1000), display="vectors")), "HabVar")

ordW13.si <- inner_join(rownames_to_column(data.frame(ordW13$CA$u), "Point"),
                        h_data_A1, by="Point") # site scores
ordW13.sp <- inner_join(rownames_to_column(data.frame(ordW13$CA$v), "Spec_code"), 
                        ordgld, by="Spec_code") # species scores
ordW13.ef <- rownames_to_column(data.frame(scores(
  envfit(ordW13, h_data_ord, choices=1:2, scaling="species", 
         permutations = 1000), display="vectors")), "HabVar")


ordplotW1 <- ggplot(ordW1.sp) + theme_bw() +
  coord_fixed() +
  geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
             position = position_dodge2(width=0.2)) +
  geom_text(aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
            position = position_dodge2(width=0.2)) + 
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = c("#e41a1c","#4daf4a")) +
  scale_shape_manual(values = c(22,23)) +
  geom_point(data=ordW1.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
  geom_segment(data=ordW1.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
               arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
  geom_text(data=ordW1.ef, aes(x=PC1-0.05, y=PC2-0.05, label=HabVar), 
            size=4, colour="black") +
  theme(axis.title.x = element_text(size=14), # enlarge x-axis labels
        axis.title.y = element_text(size=14), # enlarge y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  ggtitle("PCA: Bird abundances in sampling points (Post-breeding)")


ordplotW13 <- ggplot(ordW13.sp) + theme_bw() +
  coord_fixed() +
  geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
             position = position_dodge2(width=0.2)) +
  geom_text(aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
            position = position_dodge2(width=0.2)) + 
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = c("#e41a1c","#4daf4a")) +
  scale_shape_manual(values = c(22,23)) +
  geom_point(data=ordW13.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
  geom_segment(data=ordW13.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
               arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
  geom_text(data=ordW13.ef, aes(x=PC1-0.05, y=PC2-0.05, label=HabVar), 
            size=4, colour="black") +
  theme(axis.title.x = element_text(size=14), # enlarge x-axis labels
        axis.title.y = element_text(size=14), # enlarge y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  ggtitle("PCA: Bird abundances in sampling points (Pre-migratory)")



## 3 weeks mean filtered species ##

ordPBmean <- rda(b_data_ordPBmean, scale=T)
summary(ordPBmean)
summary(eigenvals(ordPBmean))
plot(ordPBmean, scaling="species")

ordPMmean <- rda(b_data_ordPMmean, scale=T)
summary(ordPMmean)
summary(eigenvals(ordPMmean))
plot(ordPMmean, scaling="species")


ordPBmean.si <- inner_join(rownames_to_column(data.frame(ordPBmean$CA$u), "Point"),
                           h_data_A1, by="Point") # site scores
ordPBmean.sp <- inner_join(rownames_to_column(data.frame(ordPBmean$CA$v), "Spec_code"), 
                           ordgld, by="Spec_code") # species scores
ordPBmean.ef <- rownames_to_column(data.frame(scores(
  envfit(ordPBmean, h_data_ord, choices=1:2, scaling="species", 
         permutations = 1000), display="vectors")), "HabVar")

ordPMmean.si <- inner_join(rownames_to_column(data.frame(ordPMmean$CA$u), "Point"),
                           h_data_A1, by="Point") # site scores
ordPMmean.sp <- inner_join(rownames_to_column(data.frame(ordPMmean$CA$v), "Spec_code"), 
                           ordgld, by="Spec_code") # species scores
ordPMmean.ef <- rownames_to_column(data.frame(scores(
  envfit(ordPMmean, h_data_ord, choices=1:2, scaling="species", 
         permutations = 1000), display="vectors")), "HabVar")


ordplotPBmean <- ggplot(ordPBmean.sp) + theme_bw() +
  coord_fixed() +
  geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
             position = position_dodge2(width=0.2)) +
  geom_text(aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
            position = position_dodge2(width=0.2)) + 
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = c("#e41a1c","#4daf4a")) +
  scale_shape_manual(values = c(22,23)) +
  geom_point(data=ordPBmean.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
  geom_segment(data=ordPBmean.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
               arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
  geom_text(data=ordPBmean.ef, aes(x=PC1-0.05, y=PC2-0.05, label=HabVar), 
            size=4, colour="black") +
  theme(axis.title.x = element_text(size=14), # enlarge x-axis labels
        axis.title.y = element_text(size=14), # enlarge y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  ggtitle("PCA: Bird abundances in sampling points (Post-breeding)")


ordplotPMmean <- ggplot(ordPMmean.sp) + theme_bw() +
  coord_fixed() +
  geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
             position = position_dodge2(width=0.2)) +
  geom_text(aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
            position = position_dodge2(width=0.2)) + 
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = c("#e41a1c","#4daf4a")) +
  scale_shape_manual(values = c(22,23)) +
  geom_point(data=ordPMmean.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
  geom_segment(data=ordPMmean.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
               arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
  geom_text(data=ordPMmean.ef, aes(x=PC1-0.05, y=PC2-0.05, label=HabVar), 
            size=4, colour="black") +
  theme(axis.title.x = element_text(size=14), # enlarge x-axis labels
        axis.title.y = element_text(size=14), # enlarge y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  ggtitle("PCA: Bird abundances in sampling points (Pre-migratory)")
### Analysis 4: Building models for caterpillar predation ####


cat.0 <- glm(cbind(Bird,OK) ~ Point + Week, data=cat_data, family=binomial)
summary(cat.0)
cat.1 <- update(cat.0, .~ Point + poly(Week,2))
summary(cat.1)
# further polynomials not significant


AICctab(cat.0, cat.1)


cat.m0 <- glmer(cbind(Bird,OK) ~ (1|Point) + poly(Week,2), data=cat_data, family=binomial)
summary(cat.m0)
# testing only main effects first until I get good predictor, because question is whether
# predation rates differ based on different habitat variables, and no prior reason to expect
# interaction with Week. Will check after selecting main effect.
# update: data not enough for interactions. tried with glmmTMB too. 
# "rank deficient" and convergence issue
cat.m1 <- update(cat.m0, .~. + CCavgscaled)
cat.m2 <- update(cat.m0, .~. + CCavgsd + log(TreeDens))
cat.m3 <- update(cat.m0, .~. + poly(CCavg,2))
cat.m4 <- update(cat.m0, .~. + poly(CCavgsd,2) + poly(log(TreeDens),2))
cat.m5 <- update(cat.m0, .~. + CCavgsd)
cat.m6 <- update(cat.m0, .~. + log(TreeDens))
AICctab(cat.m0, cat.m1, cat.m2, cat.m3, cat.m4, cat.m5, cat.m6)
# m6 and m1 closest to m0 so maybe interaction with week will help
cat.m7 <- update(cat.m0, .~. + log(TreeDens)*Week)
cat.m8 <- update(cat.m0, .~. + CCavgscaled*Week)
AICctab(cat.m0, cat.m1, cat.m6, cat.m7, cat.m8)

cat.m9 <- update(cat.m0, .~. + sqrt(UndDens))
cat.m10 <- update(cat.m0, .~. + TreeRich)
cat.m11 <- update(cat.m0, .~. + UndRich*Week)
cat.m12 <- update(cat.m0, .~. + TreeDiv)
cat.m13 <- update(cat.m0, .~. + UndDiv*Week)
AICctab(cat.m0, cat.m9, cat.m10, cat.m11, cat.m12, cat.m13)
cat.m14 <- update(cat.m0, .~. + TPDscaled)
cat.m15 <- update(cat.m0, .~. + HabClass)
cat.m16 <- update(cat.m0, .~. + DOM)
cat.m17 <- update(cat.m0, .~. + Moss)
AICctab(cat.m0, cat.m14, cat.m15, cat.m16, cat.m17)
summary(cat.m15)
cat.m18 <- update(cat.m15, .~. + HabClass:Week)
### Visualising Analysis 1: all birds ####

## p-values ##

# only need to check pairwise comparisons of the categorical predictor, DOM
plot(emmeans(all.6, pairwise~DOM))
# Carex lowest abundance, and significantly different from Vacc and 0 which are highest

plot(allEffects(all.6))

emmip(all.6, DOM ~ Week, cov.reduce=range)
# Moss has the highest increase in abundance from Week 1 to 13 (lowest to tied highest)
# Carex has smallest increase
