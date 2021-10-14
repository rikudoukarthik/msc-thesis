# ### Code Metadata ####
### 
### Analysis of thesis data | Part 1 (code by Karthik Thrikkadeeri) 
### 
### Bird abundances versus habitat variables; 
### post-breeding habitat breadth and selectivity of birds
### 
### Start date: 22/02/2021   End date: ~10/04/2021


# Packages ####
library(tidyverse)
library(lme4)
library(bbmle)
library(DHARMa)
# library(emmeans)
# library(effects)
library(glmmTMB)
library(vegan)
library(stargazer)
library(Hmisc)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

# Import and modify datasets ####

# Bird data from "data_Birds_2021Feb26.xlsx"; 
# Sheet "Birds"; columns A:V (all columns); rows 1:7368 (all rows)
b_rawdata <- read.delim("clipboard", as.is=F)

b_rawdata_all <- b_rawdata %>% group_by(Point, Week, Observer, Weather) %>% 
  summarise(BirdAbun = sum(Number))


# Bird species codes from "data_Birds_2021Feb26.xlsx"; 
# Sheet "Species"; columns A:E (all columns); rows 1:83 (all rows)
b_speccode <- read.delim("clipboard", as.is=F)
b_speccode <- b_speccode %>%  mutate(Migration75 = factor(Migration75))

# Modifying bird raw data:
# -Joining bird species info to Genus and Species names
# -Excluding UNID (implicit result of the innerjoin for species codes) (7323/7367 left)
# -Excluding Sky flyovers (6752/7323 rows left)
# -Excluding birds that did not enter 30m radius (3208/7323; 2840/6752 left)
# This altogether results in exclusion of 61.5% of observations, 
# i.e. it is 38.5% of b_rawdata
b_data <- inner_join(b_rawdata, b_speccode, by="Spec_code") %>% 
  filter(Flyover!="2", Within_30m==1)

# Habitat data from "data_HabVar_2021Feb26.xlsx";
# Sheet "Point Descriptions"; columns A:R (all columns); rows 1:33 (all rows)
h_rawdata <- read.delim("clipboard", as.is = F)
h_data <- h_rawdata %>% mutate(CCavgscaled=scale(CCavg), 
                               CCavgsdscaled=scale(CCavgsd),
                               TreeDensscaled=scale(TreeDens),
                               TPDscaled=scale(TreePropDeci),
                               UndDensscaled=scale(UndDens))

# Fruit phenology data from "data_HabVar_2021Feb26.xlsx";
# Sheet "Phenology"; columns A:K (all columns); rows 1:2614 (all rows)
h_phen <- read.delim("clipboard", as.is=F) %>% 
  filter(Species=="Rubus"|Species=="Vaccinium"|Species=="Fran_alnu"|Species=="Sorb_aucu") %>%
  select(-c(2,5:8,10,11)) %>% 
  mutate(Species = as.factor(Species)) %>% group_by(Point, Week, Species) %>% 
  summarise(Fruits = min(Fruits)) %>% 
  pivot_wider(names_from = "Species", values_from = "Fruits", values_fill=0)



### Modifications for analysis 1 ###

# Modifying bird data for analysis 1:
# -Removing unnecessary columns to leave only Week, Observer, Point, Weather, Wind,
#  Visibility, and Number
# -Adding column CoD "Count of Day" to test for effect of time of day and to help 
#  distinguish from true observer effect
# -Adding columns Period and Days to represent the time parameter in my study
# -Summarising Number column to Abundance per species (grouping by Species)
# -Calculating relative abundance of each species in each Point every week
#  in relation to total abundance of that species in all Points that week
b_data_A1 <- b_data %>% mutate(Date=as.Date(Date, format="%d-%m-%y")) %>% 
  group_by(Week, Date, Observer) %>% 
  arrange(StartTime, by_group=T) %>%
  mutate(CoD = as.integer(factor(StartTime))) %>% 
  arrange(Week, Date, Observer) %>% 
  group_by(Point, Week, Date, Observer, CoD, Weather, Wind, Visibility, Spec_code) %>% 
  summarise(Abundance=sum(Number)) %>% 
  ungroup() %>% 
  mutate(Days = as.factor(as.integer((Date) - as.Date("07-07-2020", format="%d-%m-%y"))), 
         .keep = "unused") %>% 
  mutate(Period = as.factor(case_when(Week %in% 1:4 ~ 1,
                                      Week %in% 5:8 ~ 2,
                                      Week %in% 9:13 ~ 3))) %>% 
  group_by(Week, Spec_code) %>% 
  mutate(RelAbun=round(Abundance/sum(Abundance),4),
         HabSel=round(Abundance/mean(Abundance),4), # mean across points present in
         HabSel2=round(Abundance/(sum(Abundance)/32),4)) %>% # mean across all points
  inner_join(b_speccode, by="Spec_code")
# HabSel2 makes more sense because it considers Points with 0 abundance, but it works
# only for the levels of All Birds and Guilds. For individual species, it might work
# for ones super abundant and widespread, but otherwise it degenerates into not such a
# meaningful measure. In such cases, using HabSel might be justified, as we already know
# the species is absent from many Points, then we can simply compare preferences among
# the Points it IS present in.



# Modifying habitat data for analysis 1:
# -Removing unnecessary columns to leave only Latitude, Longitude, CCavg, TreeDens,
#  UndDens, TreeRich, UndRich, TreeDiv, UndDiv, DOM, Moss
h_data_A1 <- h_data %>% select(-c(5, 7:8)) %>% group_by(HabClass)



# Merging bird and habitat variables into one tibble for analysis 1:
m_data_A1 <- inner_join(b_data_A1, h_data_A1, by="Point")



# Modifying for different levels of response:

# all species
b_data_A1_all <- b_data_A1 %>% 
  ungroup() %>% group_by(Point, Period, Week, Days, Observer, CoD, Weather) %>% 
  summarise(BirdAbun=sum(Abundance)) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(BirdRelAbun=round(BirdAbun/sum(BirdAbun),4),
         BirdHabSel=round(BirdAbun/mean(BirdAbun),4), # mean across points present in
         BirdHabSel2=round(BirdAbun/(sum(BirdAbun)/32),4)) # mean across all points
m_data_A1_all <- inner_join(b_data_A1_all, h_data_A1, by="Point")
m_data_A1_all_fr <- inner_join(b_data_A1_all, h_phen, by=c("Point","Week"))


# individual guilds
b_data_A1_gld <- b_data_A1 %>% 
  ungroup() %>% group_by(Point, Period, Week, Days, Observer, CoD, Weather, GuildFeed) %>% 
  summarise(GuildAbun=sum(Abundance)) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(GuildRelAbun=round(GuildAbun/sum(GuildAbun),4),
         GuildHabSel=round(GuildAbun/mean(GuildAbun),4), # mean across points present in
         GuildHabSel2=round(GuildAbun/(sum(GuildAbun)/32),4)) # mean across all points
m_data_A1_gld <- inner_join(b_data_A1_gld, h_data_A1, by="Point")

m_data_A1_omn <- m_data_A1_gld %>% filter(GuildFeed=="Omnivore")
m_data_A1_inv <- m_data_A1_gld %>% filter(GuildFeed=="Invertebrate")


# migrants and non-migrants
b_data_A1_migration <- b_data_A1 %>% 
  ungroup() %>% group_by(Point, Period, Week, Days, Observer, CoD, Weather, Migration75) %>% 
  summarise(Abun=sum(Abundance)) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(RelAbun=round(Abun/sum(Abun),4),
         HabSel=round(Abun/mean(Abun),4), # mean across points present in
         HabSel2=round(Abun/(sum(Abun)/32),4)) # mean across all points
m_data_A1_migration <- inner_join(b_data_A1_migration, h_data_A1, by="Point")

m_data_A1_NM <- m_data_A1_migration %>% filter(Migration75=="0")
m_data_A1_PM <- m_data_A1_migration %>% filter(Migration75=="1")
m_data_A1_CM <- m_data_A1_migration %>% filter(Migration75=="2")



# individual species
b_specdet <- b_data_A1 %>% ungroup() %>% group_by(Spec_code, Week) %>% 
  summarise(totpoints=n_distinct(Point), Abundance=Abundance) %>% 
  ungroup() %>% group_by(Spec_code, Week) %>% 
  summarise(totdec=sum(Abundance), totpoints=unique(totpoints))

b_data_A1_EGT <- b_data_A1 %>% filter(Spec_code=="Par_maj") %>% 
  ungroup() %>% group_by(Point, Period, Week, Days, Observer, CoD, Weather) %>% 
  summarise(Abun=Abundance) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(RelAbun=round(Abun/sum(Abun),4),
         HabSel=round(Abun/mean(Abun),4), 
         HabSel2=round(Abun/(sum(Abun)/32),4)) 
m_data_A1_EGT <- inner_join(b_data_A1_EGT, h_data_A1, by="Point")

b_data_A1_EBT <- b_data_A1 %>% filter(Spec_code=="Cya_cae") %>% 
  ungroup() %>% group_by(Point, Period, Week, Days, Observer, CoD, Weather) %>% 
  summarise(Abun=Abundance) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(RelAbun=round(Abun/sum(Abun),4),
         HabSel=round(Abun/mean(Abun),4), 
         HabSel2=round(Abun/(sum(Abun)/32),4)) 
m_data_A1_EBT <- inner_join(b_data_A1_EBT, h_data_A1, by="Point")

b_data_A1_Nut <- b_data_A1 %>% filter(Spec_code=="Sit_eur") %>% 
  ungroup() %>% group_by(Point, Period, Week, Days, Observer, CoD, Weather) %>% 
  summarise(Abun=Abundance) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(RelAbun=round(Abun/sum(Abun),4),
         HabSel=round(Abun/mean(Abun),4), 
         HabSel2=round(Abun/(sum(Abun)/32),4)) 
m_data_A1_Nut <- inner_join(b_data_A1_Nut, h_data_A1, by="Point")

b_data_A1_Phy <- b_data_A1 %>% filter(Genus=="Phylloscopus") %>% 
  ungroup() %>% group_by(Period, Week, Point, Days, Observer, CoD, Weather) %>% 
  summarise(Abun=sum(Abundance)) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(RelAbun=round(Abun/sum(Abun),4),
         HabSel=round(Abun/mean(Abun),4), 
         HabSel2=round(Abun/(sum(Abun)/32),4)) 
m_data_A1_Phy <- inner_join(b_data_A1_Phy, h_data_A1, by="Point")

b_data_A1_Rob <- b_data_A1 %>% filter(Spec_code=="Eri_rub") %>% 
  ungroup() %>% group_by(Period, Week, Point, Days, Observer, CoD, Weather) %>% 
  summarise(Abun=sum(Abundance)) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(RelAbun=round(Abun/sum(Abun),4),
         HabSel=round(Abun/mean(Abun),4), 
         HabSel2=round(Abun/(sum(Abun)/32),4)) 
m_data_A1_Rob <- inner_join(b_data_A1_Rob, h_data_A1, by="Point")

b_data_A1_Cer <- b_data_A1 %>% filter(Genus=="Certhia") %>% 
  ungroup() %>% group_by(Period, Week, Point, Days, Observer, CoD, Weather) %>% 
  summarise(Abun=sum(Abundance)) %>% 
  ungroup() %>% group_by(Week) %>% 
  mutate(RelAbun=round(Abun/sum(Abun),4),
         HabSel=round(Abun/mean(Abun),4), 
         HabSel2=round(Abun/(sum(Abun)/32),4)) 
m_data_A1_Cer <- inner_join(b_data_A1_Cer, h_data_A1, by="Point")




### for ordination ###
# removing previous:
# rm(list=c("b_data_ordW1","b_data_ordW13","b_data_ordPBmean","b_data_ordPMmean"))

# scaled predictors because using envfit and not constrained ordination
h_data_ord <- h_data_A1 %>% select(c(1, 5:11, 14)) %>% column_to_rownames("Point")


# combining information from Weeks 1:3 for first ordination
b_data_ordPB <- b_data_A1 %>% ungroup() %>% 
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
  group_by(Point, Spec_code) %>% summarise(Abundance = sum(Abundance)) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")

# combining information from Weeks 11:13 for second ordination
b_data_ordPM <- b_data_A1 %>% ungroup() %>% 
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
  group_by(Point, Spec_code) %>% summarise(Abundance = sum(Abundance)) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")



# all species: combining information from Weeks 1:3 for first ordination
b_data_ordPBall <- b_data_A1 %>% ungroup() %>% 
  mutate(Spec_code = fct_collapse(Spec_code, 
                                  Cer_sp = c("Cer_bra","Cer_fam","Cer_sp"),
                                  Phy_sp = c("Phy_col","Phy_tro","Phy_sp"),
                                  Reg_sp = c("Reg_reg","Reg_ign","Reg_sp"),
                                  Den_sp = c("Den_maj","Den_med","Den_sp"),
                                  Ant_sp = c("Ant_sp","Ant_tri"),
                                  Poe_sp = c("Poe_mon","Poe_pal","Poe_sp"))) %>% 
  filter(Week==1|Week==2|Week==3) %>%
  group_by(Point, Spec_code) %>% summarise(Abundance = sum(Abundance)) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")

# all species: combining information from Weeks 11:13 for second ordination
b_data_ordPMall <- b_data_A1 %>% ungroup() %>% 
  mutate(Spec_code = fct_collapse(Spec_code, 
                                  Cer_sp = c("Cer_bra","Cer_fam","Cer_sp"),
                                  Phy_sp = c("Phy_col","Phy_tro","Phy_sp"),
                                  Reg_sp = c("Reg_reg","Reg_ign","Reg_sp"),
                                  Den_sp = c("Den_maj","Den_med","Den_sp"),
                                  Ant_sp = c("Ant_sp","Ant_tri"),
                                  Poe_sp = c("Poe_mon","Poe_pal","Poe_sp"))) %>% 
  filter(Week==11|Week==12|Week==13) %>%
  group_by(Point, Spec_code) %>% summarise(Abundance = sum(Abundance)) %>% 
  pivot_wider(names_from = "Spec_code", values_from = "Abundance", values_fill = 0) %>% 
  column_to_rownames("Point")



### Caterpillar predation ###

cat_data <- read.delim("clipboard", as.is=F) 
cat_data <- cat_data %>% inner_join(h_data_A1, by="Point")





### ###

# Analysis ####
# Strong correlation at high CCavg is a direct physical result because with higher
# density of trees, there are less gaps overall and hence the variability in
# cover is also lower.
cor.test(h_data_A1$CCavg, h_data_A1$TreeDens, method="pearson")
plot(h_data_A1$CCavg, h_data_A1$TreeDens)
plot(h_data_A1$CCavgsd, h_data_A1$CCavg)
cor.test(h_data_A1$CCavgsd, h_data_A1$CCavg, method="pearson")
cor.test(h_data_A1$CCavgsd, h_data_A1$TreeDens, method="pearson")
plot(h_data_A1$CCavgsd, h_data_A1$TreeDens)
# CCavgsd does not correlate with TreeDens. 
# Hence, I can either include just CCavg in the models, which captures partly the 
# essences of CCavgsd and TreeDens, or I can include the latter two in the models 
# and exclude CCavg where the two variables would then capture a broader scope
# of variability than would have been possible with CCavg alone, but might miss
# some of the core essence of CCavg

cor.test(h_data_A1$TreePropDeci, h_data_A1$TreeDiv, method="spearman")
plot(h_data_A1$TreePropDeci, h_data_A1$TreeDiv)
# library(devtools)
# install_github("ProcessMiner/nlcor")
library(nlcor)
nlcor(h_data_A1$TreePropDeci, h_data_A1$TreeDiv, plt=T)
# I think it was Salek 2016 who said proportion of deciduous trees was a better
# predictor than tree diversity, and since the (non-linear) correlation is very 
# strong here, I think it is better to include TreePropDeci and exclude TreeDiv.

# Using CCavgsd+TreeDens over using CCavg is not statistically significantly better
# but has 1.79 less residual deviance and only uses one extra df, so it might be
# the way to go.
# Using TreeDiv over TreePropDeci has, technically, lower resid. dev. and lower
# Cp, but the difference is negligible (0.25).


cor.test(h_data_A1$TreePropDeci, h_data_A1$TreeDens, method="pearson")
plot(h_data_A1$TreePropDeci, h_data_A1$TreeDens)
# correlated
cor.test(h_data_A1$TreePropDeci, h_data_A1$CCavg, method="pearson")

cor.test(h_data_A1$UndDens, h_data_A1$TreeDens, method="pearson")
cor.test(h_data_A1$UndDiv, h_data_A1$TreeDens, method="pearson")
cor.test(h_data_A1$UndRich, h_data_A1$TreeDens, method="pearson")
cor.test(h_data_A1$UndDens, h_data_A1$CCavgsd, method="pearson")
cor.test(h_data_A1$UndDiv, h_data_A1$CCavgsd, method="pearson")
cor.test(h_data_A1$UndRich, h_data_A1$CCavgsd, method="pearson")




tr_data <- read.delim("clipboard") # tree species and counts

plot(h_data_A1$DOM, h_data_A1$TreeDens)
plot(h_data_A1$DOM, h_data_A1$CCavgsd)
plot(h_data_A1$DOM, h_data_A1$CCavg)
plot(h_data_A1$DOM, sqrt(h_data_A1$UndDens))
plot(h_data_A1$DOM, h_data_A1$TreeRich)
plot(h_data_A1$DOM, h_data_A1$TreeDiv)
plot(h_data_A1$DOM, h_data_A1$UndDiv)
plot(h_data_A1$DOM, h_data_A1$TreePropDeci)

plot(h_data_A1$DOM, tr_data$Pice_abie)
plot(h_data_A1$DOM, tr_data$Pinu_sylv)
plot(h_data_A1$DOM, tr_data$Quer_rope)


plot(h_data_A1$HabClass, sqrt(h_data_A1$UndDens))
plot(h_data_A1$HabClass, h_data_A1$UndDiv)
plot(h_data_A1$HabClass, h_data_A1$UndRich)
plot(h_data_A1$HabClass, h_data_A1$TreeDens)
plot(h_data_A1$HabClass, h_data_A1$CCavgsd)
plot(h_data_A1$HabClass, h_data_A1$CCavg)
plot(h_data_A1$HabClass, h_data_A1$TreeRich)
plot(h_data_A1$HabClass, h_data_A1$TreeDiv)
plot(h_data_A1$HabClass, h_data_A1$TreePropDeci)



### ###

### Analysis 1: Building models for all birds ####

boxplot(m_data_A1_all$BirdAbun ~ m_data_A1_all$Period)
boxplot(m_data_A1_all$BirdRelAbun ~ m_data_A1_all$Period)
plot(m_data_A1_all$BirdAbun ~ m_data_A1_all$Days)
plot(m_data_A1_all$BirdRelAbun ~ m_data_A1_all$Days)


### Models ###


all.0 <- glmer(BirdAbun ~ Observer + (1|Point) + (1|Days),
               data = m_data_A1_all, family=poisson())
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
plotResiduals(all.4sim, form=m_data_A1_all$Week)
plotResiduals(all.4sim, form=m_data_A1_all$CCavgsd)
plotResiduals(all.4sim, form=m_data_A1_all$TreeDens)
plotResiduals(all.4sim, form=m_data_A1_all$DOM)
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


allfruit.0 <- glmmTMB(BirdAbun ~ Observer + (1|Point) + (1|Days),
                      data=m_data_A1_all_fr, family=poisson)
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

### Analysis 3: Ordinations of species ####

rm(list=c("ordW1","ordW13","ordW1.si","ordW1.sp","ordW1.ef","ordW13.si","ordW13.sp",
          "ordW13.ef","ordplotW1","ordplotW13",
          "ordPBmean","ordPMmean","ordPBmean.si","ordPBmean.sp","ordPBmean.ef","ordPMmean.si",
          "ordPMmean.sp","ordPMmean.ef","ordplotPBmean","ordplotPMmean"))

set.seed(420)


# 3 weeks sum filtered species 

ordPB <- rda(b_data_ordPB, scale=T)
summary(ordPB)
summary(eigenvals(ordPB))
plot(ordPB, scaling="species")

ordPM <- rda(b_data_ordPM, scale=T)
summary(ordPM)
summary(eigenvals(ordPM))
plot(ordPM, scaling="species")


ordPB.si <- inner_join(rownames_to_column(data.frame(ordPB$CA$u), "Point"),
                       h_data_A1, by="Point") # site scores
ordPB.sp <- inner_join(rownames_to_column(data.frame(ordPB$CA$v), "Spec_code"), 
                       ordgld, by="Spec_code") # species scores
ordPB.ef <- rownames_to_column(data.frame(scores(
  envfit(ordPB, h_data_ord, choices=1:2, scaling="species", 
         permutations = 1000), display="vectors")), "HabVar")

ordPM.si <- inner_join(rownames_to_column(data.frame(ordPM$CA$u), "Point"),
                       h_data_A1, by="Point") # site scores
ordPM.sp <- inner_join(rownames_to_column(data.frame(ordPM$CA$v), "Spec_code"), 
                       ordgld, by="Spec_code") # species scores
ordPM.ef <- rownames_to_column(data.frame(scores(
  envfit(ordPM, h_data_ord, choices=1:2, scaling="species", 
         permutations = 1000), display="vectors")), "HabVar")





## 3 weeks sum all species (don't use)
# 
# ordPBall <- rda(b_data_ordPBall, scale=T)
# summary(ordPBall)
# summary(eigenvals(ordPBall))
# plot(ordPBall, scaling="species")
# 
# ordPMall <- rda(b_data_ordPMall, scale=T)
# summary(ordPMall)
# summary(eigenvals(ordPMall))
# plot(ordPMall, scaling="species")
# 
# 
# ordPBall.si <- inner_join(rownames_to_column(data.frame(ordPBall$CA$u), "Point"),
#                            h_data_A1, by="Point") # site scores
# ordPBall.sp <- inner_join(rownames_to_column(data.frame(ordPBall$CA$v), "Spec_code"), 
#                            ordgld, by="Spec_code") # species scores
# ordPBall.ef <- rownames_to_column(data.frame(scores(
#   envfit(ordPBall, h_data_ord, choices=1:2, scaling="species", 
#          permutations = 1000), display="vectors")), "HabVar")
# 
# ordPMall.si <- inner_join(rownames_to_column(data.frame(ordPMall$CA$u), "Point"),
#                            h_data_A1, by="Point") # site scores
# ordPMall.sp <- inner_join(rownames_to_column(data.frame(ordPMall$CA$v), "Spec_code"), 
#                            ordgld, by="Spec_code") # species scores
# ordPMall.ef <- rownames_to_column(data.frame(scores(
#   envfit(ordPMall, h_data_ord, choices=1:2, scaling="species", 
#          permutations = 1000), display="vectors")), "HabVar")
# 
# 
# ordplotPBall <- ggplot(ordPBall.sp) + theme_bw() +
#   coord_fixed() +
#   geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
#              position = position_dodge2(width=0.2)) +
#   geom_text(data = filter(ordPBall.sp, 
#                           Spec_code=="Par_maj"|
#                             Spec_code=="Cya_cae"|
#                             Spec_code=="Sit_eur"|
#                             Spec_code=="Eri_rub"|
#                             Spec_code=="Cer_sp"|
#                             Spec_code=="Phy_sp"|
#                             Spec_code=="Reg_sp"|
#                             Spec_code=="Den_sp"), 
#             aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
#             position = position_dodge2(width=0.2)) + 
#   scale_colour_manual(values = cbbPalette) +
#   scale_fill_manual(values = c("#377eb8","#e41a1c","#4daf4a")) +
#   scale_shape_manual(values = c(23,22,24)) +
#   geom_point(data=ordPBall.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
#   geom_segment(data=ordPBall.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
#                arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
#   geom_text(data=ordPBall.ef, aes(x=PC1-0.05, y=PC2-0.05, label=HabVar), 
#             size=4, colour="black") +
#   theme(axis.title.x = element_text(size=14), # enlarge x-axis labels
#         axis.title.y = element_text(size=14), # enlarge y-axis labels
#         panel.background = element_blank(), 
#         panel.grid.major = element_blank(),  #remove major-grid labels
#         panel.grid.minor = element_blank(),  #remove minor-grid labels
#         plot.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype="dotted") +
#   geom_vline(xintercept = 0, linetype="dotted") +
#   ggtitle("PCA: Bird abundances in sampling points (Post-breeding)")
# 
# 
# ordplotPMall <- ggplot(ordPMall.sp) + theme_bw() +
#   coord_fixed() +
#   geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
#              position = position_dodge2(width=0.2)) +
#   geom_text(data = filter(ordPBall.sp, 
#                           Spec_code=="Par_maj"|
#                             Spec_code=="Cya_cae"|
#                             Spec_code=="Sit_eur"|
#                             Spec_code=="Eri_rub"|
#                             Spec_code=="Cer_sp"|
#                             Spec_code=="Phy_sp"|
#                             Spec_code=="Reg_sp"|
#                             Spec_code=="Den_sp"), 
#             aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
#             position = position_dodge2(width=0.2)) + 
#   scale_colour_manual(values = cbbPalette) +
#   scale_fill_manual(values = c("#377eb8","#e41a1c","#4daf4a")) +
#   scale_shape_manual(values = c(23,22,24)) +
#   geom_point(data=ordPMall.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
#   geom_segment(data=ordPMall.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
#                arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
#   geom_text(data=ordPMall.ef, aes(x=PC1-0.05, y=PC2-0.05, label=HabVar), 
#             size=4, colour="black") +
#   theme(axis.title.x = element_text(size=14), # enlarge x-axis labels
#         axis.title.y = element_text(size=14), # enlarge y-axis labels
#         panel.background = element_blank(), 
#         panel.grid.major = element_blank(),  #remove major-grid labels
#         panel.grid.minor = element_blank(),  #remove minor-grid labels
#         plot.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype="dotted") +
#   geom_vline(xintercept = 0, linetype="dotted") +
#   ggtitle("PCA: Bird abundances in sampling points (Pre-migratory)")



### ###

### Analysis 4: Building models for caterpillar predation ####

# tested only main effects first until I get good predictor, because question is whether
# predation rates differ based on different habitat variables, and no prior reason to expect
# interaction with Week (apart from main eff of week). Will check after selecting main effect.
# update: data not enough for interactions. tried with glmmTMB too. 
# "rank deficient" but this not necessarily a problem. 
# rm(list=c("cat.m0","cat.m1","cat.m2","cat.m3","cat.m4","cat.m5","cat.m6","cat.m7","cat.m8",
# "cat.m9","cat.m10","cat.m11","cat.m12","cat.m13","cat.m14","cat.m15","cat.m16",
# "cat.m17","cat.m18"))
# cat.m1 <- update(cat.m0, .~. + CCavgscaled)
# cat.m2 <- update(cat.m0, .~. + CCavgsd + log(TreeDens))
# cat.m3 <- update(cat.m0, .~. + poly(CCavg,2))
# cat.m4 <- update(cat.m0, .~. + poly(CCavgsd,2) + poly(log(TreeDens),2))
# cat.m5 <- update(cat.m0, .~. + CCavgsd)
# cat.m6 <- update(cat.m0, .~. + log(TreeDens))
# AICctab(cat.m0, cat.m1, cat.m2, cat.m3, cat.m4, cat.m5, cat.m6)
# # m6 and m1 closest to m0 so maybe interaction with week will help
# cat.m7 <- update(cat.m0, .~. + log(TreeDens)*Week)
# cat.m8 <- update(cat.m0, .~. + CCavgscaled*Week)
# AICctab(cat.m0, cat.m1, cat.m6, cat.m7, cat.m8)
# 
# cat.m9 <- update(cat.m0, .~. + sqrt(UndDens))
# cat.m10 <- update(cat.m0, .~. + TreeRich)
# cat.m11 <- update(cat.m0, .~. + UndRich*Week)
# cat.m12 <- update(cat.m0, .~. + TreeDiv)
# cat.m13 <- update(cat.m0, .~. + UndDiv*Week)
# AICctab(cat.m0, cat.m9, cat.m10, cat.m11, cat.m12, cat.m13)
# cat.m14 <- update(cat.m0, .~. + TPDscaled)
# cat.m15 <- update(cat.m0, .~. + HabClass)
# cat.m16 <- update(cat.m0, .~. + DOM)
# cat.m17 <- update(cat.m0, .~. + Moss)
# AICctab(cat.m0, cat.m14, cat.m15, cat.m16, cat.m17)
# summary(cat.m15)
# cat.m18 <- update(cat.m15, .~. + HabClass:Week)


# most variables not significant, but UndDiv marginal
cat.0 <- glmer(cbind(Bird,OK) ~ (1|Point) + Week, data=cat_data, family=binomial)
cat.1 <- glmer(cbind(Bird,OK) ~ (1|Point) + poly(Week,2), data=cat_data, family=binomial)
cat.2a <- update(cat.1, .~. + UndRich)
cat.2b <- update(cat.1, .~. + UndDiv)
cat.2c <- update(cat.1, .~. + sqrt(UndDens))
cat.2d <- update(cat.1, .~. + UndRich*Week)
cat.2e <- update(cat.1, .~. + UndDiv*Week)
cat.2f <- update(cat.1, .~. + sqrt(UndDens)*Week)
AICctab(cat.0, cat.1, cat.2a, cat.2b, cat.2c, cat.2d, cat.2e, cat.2f)
cat.2g <- update(cat.1, .~. + UndDiv*poly(Week,2))
AICctab(cat.1, cat.2e, cat.2g)
cat.2h <- update(cat.1, .~. + UndDiv*Week + HabClass)
cat.2i <- update(cat.1, .~. + HabClass)
AICctab(cat.1, cat.2e, cat.2h, cat.2i)
rm(list = c("cat.2a","cat.2b","cat.2c","cat.2d","cat.2e","cat.2f","cat.2g","cat.2h","cat.2i"))
# rank deficient problem:
model.matrix(cbind(Bird,OK) ~ Week + UndDiv*Week, 
             cat_data)

cat.2 <- update(cat.1, .~. + HabClass)
summary(cat.2)
cat.3 <- update(cat.1, .~. + UndDiv*Week + HabClass)
AICctab(cat.0, cat.1, cat.2, cat.3)

simulateResiduals(cat.3, plot=T)



# ## considering Lost as also predated ##
# 
# catTOT.0 <- glmer(cbind(TotalPred,OK) ~ (1|Point) + Week, data=cat_data, family=binomial)
# catTOT.1 <- glmer(cbind(TotalPred,OK) ~ (1|Point) + poly(Week,2), data=cat_data, family=binomial)
# AICctab(catTOT.0, catTOT.1)
# 
# catTOT.2a <- update(catTOT.1, .~. + UndRich)
# catTOT.2b <- update(catTOT.1, .~. + UndDiv)
# catTOT.2c <- update(catTOT.1, .~. + sqrt(UndDens))
# catTOT.2d <- update(catTOT.1, .~. + UndRich*Week)
# catTOT.2e <- update(catTOT.1, .~. + UndDiv*Week)
# catTOT.2f <- update(catTOT.1, .~. + sqrt(UndDens)*Week)
# AICctab(catTOT.0, catTOT.1, catTOT.2a, catTOT.2b, catTOT.2c, catTOT.2d, catTOT.2e, catTOT.2f)
# catTOT.2g <- update(catTOT.1, .~. + UndDiv*poly(Week,2))
# AICctab(catTOT.1, catTOT.2c, catTOT.2g)
# catTOT.2h <- update(catTOT.1, .~. + UndDiv*Week + HabClass)
# catTOT.2i <- update(catTOT.1, .~. + HabClass)
# AICctab(catTOT.1, catTOT.2c, catTOT.2h, catTOT.2i)
# rm(list = c("catTOT.2a","catTOT.2b","catTOT.2c","catTOT.2d","catTOT.2e","catTOT.2f","catTOT.2g","catTOT.2h","catTOT.2i"))
# 
# catTOT.2a <- update(catTOT.1, .~. + CCavgsd + log(TreeDens))
# catTOT.2b <- update(catTOT.1, .~. + CCavgsd)
# catTOT.2c <- update(catTOT.1, .~. + log(TreeDens))
# catTOT.2d <- update(catTOT.1, .~. + CCavgsd*Week)
# catTOT.2e <- update(catTOT.1, .~. + log(TreeDens)*Week)
# catTOT.2f <- update(catTOT.1, .~. + poly(CCavgsd,2) + poly(log(TreeDens),2))
# AICctab(catTOT.0, catTOT.1, catTOT.2a, catTOT.2b, catTOT.2c, catTOT.2d, catTOT.2e, catTOT.2f)
# rm(list = c("catTOT.2a","catTOT.2b","catTOT.2c","catTOT.2d","catTOT.2e","catTOT.2f"))
# 
# catTOT.2a <- update(catTOT.1, .~. + TPDscaled)
# catTOT.2b <- update(catTOT.1, .~. + DOM)
# catTOT.2c <- update(catTOT.1, .~. + Moss)
# catTOT.2d <- update(catTOT.1, .~. + TreeRich)
# catTOT.2e <- update(catTOT.1, .~. + TreeDiv)
# AICctab(catTOT.0, catTOT.1, catTOT.2a, catTOT.2b, catTOT.2c, catTOT.2d, catTOT.2e)
# rm(list = c("catTOT.2a","catTOT.2b","catTOT.2c","catTOT.2d","catTOT.2e"))
# 
# summary(catTOT.1)
# 
# 
# # correlation with bird abundance
# 
# cat_corrTOT <- inner_join(select(cat_data, c(1,2,25:28)), 
#                        select(b_data_A1_all, c(1,3,8)), 
#                        by=c("Week","Point")) %>% 
#   mutate(Pred = TotalPred/(TotalPred+OK))
# plot(cat_corrTOT$BirdAbun, cat_corrTOT$Pred)
# cor.test(cat_corrTOT$BirdAbun, cat_corrTOT$Pred, method="pearson")
# 
# 
# cat_corrTOTplot <- ggplot(data=cat_corrTOT,
#                        mapping=aes(x=BirdAbun, y=Pred)) +
#   coord_cartesian(ylim=c(0,0.5)) +
#   scale_x_continuous(breaks=seq(0,30,5)) +
#   geom_point(position = "jitter", alpha=0.5, size=2, colour="goldenrod3")
# 
# ggsave("cat_corrTOTplot.png", cat_corrTOTplot, 
#        width = 4, height = 5, units = "in", dpi=300)
# 
# 
# # with only insectivores
# cat_corrTOT2 <- inner_join(select(cat_data, c(1,2,25:28)), 
#                         select(m_data_A1_inv, c(1,3,9)), 
#                         by=c("Week","Point")) %>% 
#   mutate(Pred = TotalPred/(TotalPred+OK))
# plot(cat_corrTOT2$GuildAbun, cat_corrTOT2$Pred)
# cor.test(cat_corrTOT2$GuildAbun, cat_corrTOT2$Pred, method="pearson")
# 
# cat_corrTOT2plot <- ggplot(data=cat_corrTOT2,
#                         mapping=aes(x=GuildAbun, y=Pred)) +
#   coord_cartesian(ylim=c(0,0.5)) +
#   scale_x_continuous(breaks=seq(0,30,5)) +
#   geom_point(position = "jitter", alpha=0.5, size=2, colour="goldenrod3")
# 
# ggsave("cat_corrTOT2plot.png", cat_corrTOT2plot, 
#        width = 4, height = 5, units = "in", dpi=300)



### ###


# Visualisations ####

cmult <- 1.96

# theme_simple <- theme_bw() + theme(panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(),
# panel.border = element_rect(size=1.2))
theme_simple <- theme_classic()
theme_set(theme_simple)

# purples: #A204B4 #D91EFA
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072A2", "#D55E00", "#CC79A7")
Set1 <- brewer.pal(9, name="Set1")
Set3 <- brewer.pal(12, name="Set3")

### ###

### Visualising results: summaries ####

# summary statistics etc. of some results


### detections of birds ###

# number of unique species observed in total = 69
n_distinct( (b_rawdata %>% 
               inner_join(b_speccode, by="Spec_code") %>% 
               filter(Species != is.na(Species)))$Spec_code)
# number of unique species observed in relevant data = 45
n_distinct( (b_data %>% 
               filter(Species != is.na(Species)))$Spec_code)
# number of observations of foraging
dim(filter(b_rawdata, Foraging==1)) # 570 individual birds foraging
dim(filter(b_data, Foraging==1)) # 494 individual birds foraging within 30 m
dim(filter(b_data, Fruits==1)) # 2 observations of birds eating fruits
# visual/auditory
dim(filter(b_rawdata, Heard==1 & Seen==0)) # 6398 = 86.9%
dim(filter(b_rawdata, Heard==0 & Seen==1)) # 361 = 4.9%
dim(filter(b_rawdata, Heard==1 & Seen==1)) # 608 = 8.25%
dim(filter(b_data, Heard==1 & Seen==0)) # 2209 = 77.78%
dim(filter(b_data, Heard==0 & Seen==1)) # 201 = 7.08%
dim(filter(b_data, Heard==1 & Seen==1)) # 430 = 15.14%
# observer effect
dim(filter(b_rawdata, Observer=="IF")) # 3371   
dim(filter(b_rawdata, Observer=="KT")) # 3996   
dim(filter(b_data, Observer=="IF")) # 834  
dim(filter(b_data, Observer=="KT")) # 2006 

### ###

### Visualising results: tables  ####


### Table: model for all birds ###
allAIC <- AICctab(all.1, all.2, all.3, all.4, all.5, all.6, weights=T, base=T, logLik=T)
summary(all.6)
# printing AICc table for import into Word in Methods section
stargazer(as.data.frame(allAIC), 
          summary = F,
          type="html", 
          out = "D:/COLLEGE/2019 - MSc Ecology/Thesis/4 - birds/2 - Data Analysis/Analysis1_Birds-HabVar/all6summary.html")

### Table: models for guilds ###
summary(inv.4)
summary(omn.3)
invAIC <- AICctab(inv.1, inv.2, inv.3, inv.4, weights=T, base=T, logLik=T)
omnAIC <- AICctab(omn.1, omn.2, omn.3, weights=T, base=T, logLik=T)
stargazer(rbind(as.data.frame(invAIC), 
                as.data.frame(omnAIC)), 
          summary = F,
          type="html", 
          out = "D:/COLLEGE/2019 - MSc Ecology/Thesis/4 - birds/2 - Data Analysis/
          Analysis1_Birds-HabVar/gldsummary.html")

### Table: models for caterpillars ###
catAIC <- AICctab(cat.0, cat.1, cat.2, cat.3, weights=T, base=T, logLik=T)
stargazer(as.data.frame(catAIC), 
          summary = F,
          type="html", 
          out = "D:/COLLEGE/2019 - MSc Ecology/Thesis/4 - birds/2 - Data Analysis/Analysis1_Birds-HabVar/catsummary.html")



### Table of species detected  ###

b_speclist30m <- inner_join(unique((b_data %>%
                                      filter(Species != is.na(Species)) %>%
                                      select("Spec_code"))),
                            b_speccode, by = "Spec_code") %>% 
  arrange(Spec_code) %>% rownames_to_column("S.No.")

stargazer(b_speclist30m, type="html", summary=F,
          out = "D:/COLLEGE/2019 - MSc Ecology/Thesis/4 - birds/2 - Data Analysis/Analysis1_Birds-HabVar/specieslist.html")



### Table of points and habitat variables  ###

stargazer(select(h_data_A1, -c(16:20)), type="html", summary=F,
          out = "D:/COLLEGE/2019 - MSc Ecology/Thesis/4 - birds/2 - Data Analysis/Analysis1_Birds-HabVar/pointhabvar.html")


### Table of species with guild info and overall abundances ###

b_speclistOverall <- b_data %>% ungroup() %>% 
  mutate(Spec_code = fct_collapse(Spec_code, 
                                  Acc_sp = c("Acc_nis","Acc_sp"),
                                  Ant_sp = c("Ant_pra","Ant_tri","Ant_sp"),
                                  Cer_sp = c("Cer_bra","Cer_fam","Cer_sp"),
                                  Cor_sp = c("Cor_corax","Cor_corone","Cor_monedula","Cor_sp"),
                                  Phy_sp = c("Phy_col","Phy_tro","Phy_sp"),
                                  Reg_sp = c("Reg_reg","Reg_ign","Reg_sp"),
                                  Den_sp = c("Den_maj","Den_med","Den_sp"),
                                  Fic_sp = c("Fic_sp"),
                                  Fal_sp = c("Fal_tin","Fal_sp"),
                                  Mot_sp = c("Mot_alb","Mot_sp"),
                                  Pas_sp = c("Pas_mon","Pas_sp"),
                                  Poe_sp = c("Poe_mon","Poe_pal","Poe_sp"),
                                  Tur_sp = c("Tur_mer","Tur_phi","Tur_pil","Tur_vis","Tur_sp")
  )) %>% 
  group_by(GuildFeed, Spec_code) %>% 
  summarise(TotDet = n()) %>% 
  inner_join(b_speccode[,1:3], by = "Spec_code") %>% 
  arrange(GuildFeed, desc(TotDet)) %>% rownames_to_column("S.No.")

write.csv(b_speclistOverall, 
          file="specieslistOverall.csv", 
          row.names = F)






### ###

### Visualising results: (Fig4) Analysis 1: all birds ####


## Week-DOM interaction ##

all.pred1 <- data.frame(Week = seq(1, 13, 0.01),
                        DOM = sample(unique(m_data_A1_all$DOM), 1201, replace = T),
                        Observer = factor("KT", levels = levels(m_data_A1_all$Observer)),
                        CCavgsd = mean(m_data_A1_all$CCavgsd),
                        TreeDens = mean(m_data_A1_all$TreeDens))
all.pred1$BirdAbun <- predict(all.6, all.pred1, type="response", re.form=NA)
# no random effects because interested in population level patterns

# Ben Bolker's vignette for plotting confidence intervals:
#(https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions)
all.pred1mm <- model.matrix(terms(all.6), all.pred1)
all.pred1pvar <- diag(all.pred1mm %*% tcrossprod(vcov(all.6), all.pred1mm)) # fixef uncertainty only
# all.pred1tvar <- all.pred1pvar + VarCorr(all.6$Point[1]) + VarCorr(all.6$Days[1]) random uncertainty
all.pred1 <- data.frame(all.pred1,
                        lo = all.pred1$BirdAbun - exp(cmult*sqrt(all.pred1pvar)),
                        hi = all.pred1$BirdAbun + exp(cmult*sqrt(all.pred1pvar)))
# plotting
all.pred1plot <- ggplot(aes(x=Week, y=BirdAbun, col=DOM, fill=DOM), data=all.pred1) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
  scale_colour_manual(values = cbbPalette,
                      labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
  geom_line(size=2) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.3, colour=NA) +
  scale_y_continuous(breaks = seq(0,40,4)) +
  scale_x_continuous(breaks = seq(1,13,1)) +
  coord_cartesian(ylim = c(2,18)) +
  labs(y = "Bird detections per point count") +
  guides(col=guide_legend(title = "Ground veg."),
         fill=guide_legend(title = "Ground veg."))


## Difference between DOM layers ##
# Week 1
all.pred2a <- data.frame(Week = 1,
                         DOM = sample(unique(m_data_A1_all$DOM), 1201, replace = T),
                         Observer = factor("KT", levels = levels(m_data_A1_all$Observer)),
                         CCavgsd = mean(m_data_A1_all$CCavgsd),
                         TreeDens = mean(m_data_A1_all$TreeDens))
all.pred2a$BirdAbun <- predict(all.6, all.pred2a, type="response", re.form=NA)
all.pred2amm <- model.matrix(terms(all.6), all.pred2a)
all.pred2apvar <- diag(all.pred2amm %*% tcrossprod(vcov(all.6), all.pred2amm)) # fixef uncertainty only
all.pred2a <- data.frame(all.pred2a,
                         lo = all.pred2a$BirdAbun - exp(cmult*sqrt(all.pred2apvar)),
                         hi = all.pred2a$BirdAbun + exp(cmult*sqrt(all.pred2apvar)))
# Week 13
all.pred2b <- data.frame(Week = 13,
                         DOM = sample(unique(m_data_A1_all$DOM), 1201, replace = T),
                         Observer = factor("KT", levels = levels(m_data_A1_all$Observer)),
                         CCavgsd = mean(m_data_A1_all$CCavgsd),
                         TreeDens = mean(m_data_A1_all$TreeDens))
all.pred2b$BirdAbun <- predict(all.6, all.pred2b, type="response", re.form=NA)
all.pred2bmm <- model.matrix(terms(all.6), all.pred2b)
all.pred2bpvar <- diag(all.pred2bmm %*% tcrossprod(vcov(all.6), all.pred2bmm)) # fixef uncertainty only
all.pred2b <- data.frame(all.pred2b,
                         lo = all.pred2b$BirdAbun - exp(cmult*sqrt(all.pred2bpvar)),
                         hi = all.pred2b$BirdAbun + exp(cmult*sqrt(all.pred2bpvar)))
# joining
all.pred2 <- rbind(all.pred2a, all.pred2b)
# plotting
all.pred2plot <- ggplot(data=all.pred2, mapping=aes(x=DOM, y=BirdAbun, col=factor(Week))) +
  guides(col = guide_legend(title = "Week", reverse = T)) +
  scale_colour_manual(values = cbbPalette[c(7, 6)]) +
  stat_summary(fun = "mean", geom = "point", size = 3.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, size = 1.5) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.3, colour=NA) +
  scale_y_continuous(breaks = seq(0,40,4)) +
  scale_x_discrete(labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
  coord_cartesian(ylim = c(2,18)) +
  labs(y = "Bird detections per point count",
       x = "Ground vegetation") +
  theme(axis.text.x = element_text(size=7))


## Abun - CCavgsd (for 0 DOM ground layer) ##
# Week 1
all.pred3a <- data.frame(Week = 1,
                         DOM = factor("0", levels = levels(m_data_A1_all$DOM)),
                         Observer = factor("KT", levels = levels(m_data_A1_all$Observer)),
                         CCavgsd = seq(1, 20, 0.001),
                         TreeDens = mean(m_data_A1_all$TreeDens))

all.pred3a$BirdAbun <- predict(all.6, all.pred3a, type="response", re.form=NA)
all.pred3amm <- model.matrix(terms(all.6), all.pred3a)
all.pred3apvar <- diag(all.pred3amm %*% tcrossprod(vcov(all.6), all.pred3amm))
all.pred3a <- data.frame(all.pred3a,
                         lo = all.pred3a$BirdAbun - exp(cmult*sqrt(all.pred3apvar)),
                         hi = all.pred3a$BirdAbun + exp(cmult*sqrt(all.pred3apvar)))
# Week 13
all.pred3b <- data.frame(Week = 13,
                         DOM = factor("0", levels = levels(m_data_A1_all$DOM)),
                         Observer = factor("KT", levels = levels(m_data_A1_all$Observer)),
                         CCavgsd = seq(1, 20, 0.001),
                         TreeDens = mean(m_data_A1_all$TreeDens))

all.pred3b$BirdAbun <- predict(all.6, all.pred3b, type="response", re.form=NA)
all.pred3bmm <- model.matrix(terms(all.6), all.pred3b)
all.pred3bpvar <- diag(all.pred3bmm %*% tcrossprod(vcov(all.6), all.pred3bmm))
all.pred3b <- data.frame(all.pred3b,
                         lo = all.pred3b$BirdAbun - exp(cmult*sqrt(all.pred3bpvar)),
                         hi = all.pred3b$BirdAbun + exp(cmult*sqrt(all.pred3bpvar)))
# joining
all.pred3 <- rbind(all.pred3a, all.pred3b)
# plotting
all.pred3plot <- ggplot(data=all.pred3, mapping=aes(x=CCavgsd, y=BirdAbun, 
                                                    col=factor(Week), fill=factor(Week))) +
  guides(col=guide_legend(title="Week", reverse = T),
         fill=guide_legend(title="Week", reverse = T)) +
  scale_colour_manual(values=cbbPalette[c(7,6)]) +
  scale_fill_manual(values=cbbPalette[c(7,6)]) +
  geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lo, ymax=hi), col=NA, alpha=0.3) +
  scale_y_continuous(breaks = seq(0,40,4)) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  coord_cartesian(ylim = c(2,20)) +
  labs(y = "Bird detections per point count",
       x = "Canopy heterogeneity") +
  theme(legend.justification = c(1,0.5))


## Abun - TreeDens (for 0 DOM ground layer) ##
# Week 1
all.pred4a <- data.frame(Week = 1,
                         DOM = factor("0", levels = levels(m_data_A1_all$DOM)),
                         Observer = factor("KT", levels = levels(m_data_A1_all$Observer)),
                         CCavgsd = mean(m_data_A1_all$CCavgsd),
                         TreeDens = seq(1,15,0.001))

all.pred4a$BirdAbun <- predict(all.6, all.pred4a, type="response", re.form=NA)
all.pred4amm <- model.matrix(terms(all.6), all.pred4a)
all.pred4apvar <- diag(all.pred4amm %*% tcrossprod(vcov(all.6), all.pred4amm))
all.pred4a <- data.frame(all.pred4a,
                         lo = all.pred4a$BirdAbun - exp(cmult*sqrt(all.pred4apvar)),
                         hi = all.pred4a$BirdAbun + exp(cmult*sqrt(all.pred4apvar)))
# Week 13
all.pred4b <- data.frame(Week = 13,
                         DOM = factor("0", levels = levels(m_data_A1_all$DOM)),
                         Observer = factor("KT", levels = levels(m_data_A1_all$Observer)),
                         CCavgsd = mean(m_data_A1_all$CCavgsd),
                         TreeDens = seq(1,15,0.001))

all.pred4b$BirdAbun <- predict(all.6, all.pred4b, type="response", re.form=NA)
all.pred4bmm <- model.matrix(terms(all.6), all.pred4b)
all.pred4bpvar <- diag(all.pred4bmm %*% tcrossprod(vcov(all.6), all.pred4bmm))
all.pred4b <- data.frame(all.pred4b,
                         lo = all.pred4b$BirdAbun - exp(cmult*sqrt(all.pred4bpvar)),
                         hi = all.pred4b$BirdAbun + exp(cmult*sqrt(all.pred4bpvar)))
# joining
all.pred4 <- rbind(all.pred4a, all.pred4b)
# plotting
all.pred4plot <- ggplot(data=all.pred4, mapping=aes(x=TreeDens, y=BirdAbun, 
                                                    col=factor(Week), fill=factor(Week))) +
  guides(col=guide_legend(title="Week", reverse = T),
         fill=guide_legend(title="Week", reverse = T)) +
  scale_colour_manual(values=cbbPalette[c(7,6)]) +
  scale_fill_manual(values=cbbPalette[c(7,6)]) +
  geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lo, ymax=hi), col=NA, alpha=0.3) +
  scale_y_continuous(breaks = seq(0,40,4)) +
  scale_x_continuous(breaks = seq(0,20,2)) +
  coord_cartesian(ylim = c(2,20)) +
  labs(y = "Bird detections per point count",
       x = expression(Tree~density~per~100~m^2)) +
  theme(legend.justification = c(1,0.5))



## Fig.4: analysis 1 ##

( ( all.pred1plot | (all.pred2plot + theme(axis.title.y = element_blank())) ) ) /
  ( ( all.pred3plot | (all.pred4plot + theme(axis.title.y = element_blank())) ) ) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) -> fig4

ggsave("Fig4.png", fig4, 
       width = 18, height = 24, units = "cm", dpi=300)



### ###


### Visualising results: (Fig5) Analysis 2: inv-feeders ####


## Difference between DOM layers ##
inv.pred1 <- data.frame(DOM = sample(unique(m_data_A1_inv$DOM), 1201, replace = T),
                        Observer = factor("KT", levels = levels(m_data_A1_inv$Observer)),
                        TPDscaled = mean(m_data_A1_inv$TPDscaled),
                        Point = sample(unique(m_data_A1_inv$Point), 1201, replace = T))
inv.pred1$GuildAbun <- predict(inv.4, inv.pred1, type="response", re.form=NA)
# don't know why it requires to specify Point even with re.form=NA
inv.pred1mm <- model.matrix(terms(inv.4), inv.pred1)
inv.pred1var <- diag(inv.pred1mm %*% vcov(inv.4)[["cond"]] %*% t(inv.pred1mm))
inv.pred1 <- data.frame(inv.pred1,
                        lo = inv.pred1$GuildAbun - exp(cmult*sqrt(inv.pred1var)),
                        hi = inv.pred1$GuildAbun + exp(cmult*sqrt(inv.pred1var)))
# plotting
inv.pred1plot <- ggplot(data=inv.pred1, mapping=aes(x=DOM, y=GuildAbun)) +
  geom_point(size=3.5, col="#655B1B") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, size = 1.5, col="#655B1B") +
  scale_y_continuous(breaks = seq(0,40,2)) +
  scale_x_discrete(labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
  coord_cartesian(ylim = c(0,10)) +
  labs(y = "Invertebrate-feeder detections per point count",
       x = "Ground vegetation") +
  theme(axis.text.x = element_text(size=7))


## Abun - TPDscaled (for 0 DOM ground layer) ##
inv.pred2 <- data.frame(DOM = factor("0", levels = levels(m_data_A1_inv$DOM)),
                        Observer = factor("KT", levels = levels(m_data_A1_inv$Observer)),
                        TPDscaled = seq(-1.5,1.5,0.001),
                        Point = sample(unique(m_data_A1_inv$Point), 3001, replace = T))
inv.pred2$GuildAbun <- predict(inv.4, inv.pred2, type="response", re.form=NA)
inv.pred2mm <- model.matrix(terms(inv.4), inv.pred2)
inv.pred2var <- diag(inv.pred2mm %*% vcov(inv.4)[["cond"]] %*% t(inv.pred2mm))
inv.pred2 <- data.frame(inv.pred2,
                        lo = inv.pred2$GuildAbun - exp(cmult*sqrt(inv.pred2var)),
                        hi = inv.pred2$GuildAbun + exp(cmult*sqrt(inv.pred2var)))
# plotting
inv.pred2plot <- ggplot(data=inv.pred2, mapping=aes(x=TPDscaled, y=GuildAbun)) +
  geom_line(size=1.5, col="#655B1B") +
  geom_ribbon(aes(ymin=lo, ymax=hi), col=NA, alpha=0.3, fill="#655B1B") +
  scale_y_continuous(breaks = seq(0,40,2)) +
  coord_cartesian(ylim=c(0,10)) +
  labs(y = "Invertebrate-feeder detections per point count",
       x = "Proportion of deciduous trees (scaled)") 


## fig5 ##

fig5 <- inv.pred1plot + (inv.pred2plot + theme(axis.title.y = element_blank())) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16),
                                            plot.tag.position = c(0,1)) 

ggsave("Fig5.png", fig5, 
       width = 18, height = 16, units = "cm", dpi=300)


### ###


### Visualising results: (Fig6) Analysis 2: omnivores ####


## Difference between HabClass layers ##
# Week 1
omn.pred1a <- data.frame(Week = 1,
                         HabClass = sample(unique(m_data_A1_omn$HabClass), 1201, replace = T),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = mean(m_data_A1_omn$TreeDens),
                         CCavgsd = mean(m_data_A1_omn$CCavgsd),
                         Point = sample(unique(m_data_A1_omn$Point), 1201, replace = T))
omn.pred1a$GuildAbun <- predict(omn.4, omn.pred1a, type="response", re.form=NA)
# don't know why it requires to specify Point even with re.form=NA
omn.pred1amm <- model.matrix(terms(omn.4), omn.pred1a)
omn.pred1avar <- diag(omn.pred1amm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred1amm))
omn.pred1a <- data.frame(omn.pred1a,
                         lo = omn.pred1a$GuildAbun - exp(cmult*sqrt(omn.pred1avar)),
                         hi = omn.pred1a$GuildAbun + exp(cmult*sqrt(omn.pred1avar)))
# Week 13
omn.pred1b <- data.frame(Week = 13,
                         HabClass = sample(unique(m_data_A1_omn$HabClass), 1201, replace = T),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = mean(m_data_A1_omn$TreeDens),
                         CCavgsd = mean(m_data_A1_omn$CCavgsd),
                         Point = sample(unique(m_data_A1_omn$Point), 1201, replace = T))
omn.pred1b$GuildAbun <- predict(omn.4, omn.pred1b, type="response", re.form=NA)
omn.pred1bmm <- model.matrix(terms(omn.4), omn.pred1b)
omn.pred1bvar <- diag(omn.pred1bmm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred1bmm))
omn.pred1b <- data.frame(omn.pred1b,
                         lo = omn.pred1b$GuildAbun - exp(cmult*sqrt(omn.pred1bvar)),
                         hi = omn.pred1b$GuildAbun + exp(cmult*sqrt(omn.pred1bvar)))
# joining
omn.pred1 <- rbind(omn.pred1a, omn.pred1b)
# plotting
omn.pred1plot <- ggplot(data=omn.pred1, 
                        mapping=aes(x=HabClass, y=GuildAbun, col=factor(Week))) +
  guides(col = guide_legend(title = "Week", reverse = T)) +
  scale_colour_manual(values = cbbPalette[c(7, 6)]) +
  stat_summary(fun = "mean", geom = "point", size = 3.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, size = 1.5) +
  scale_y_continuous(breaks = seq(0,40,2)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Omnivore detections per point count",
       x = "Habitat class") +
  theme(axis.text.x = element_text(size=8))


## Abun - TreeDens (for Edge) ##
# Week 1
omn.pred2a <- data.frame(Week = 1,
                         HabClass = factor("Edge", levels=levels(m_data_A1_omn$HabClass)),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = seq(1,15,0.01),
                         CCavgsd = mean(m_data_A1_omn$CCavgsd),
                         Point = sample(unique(m_data_A1_omn$Point), 1401, replace = T))
omn.pred2a$GuildAbun <- predict(omn.4, omn.pred2a, type="response", re.form=NA)
omn.pred2amm <- model.matrix(terms(omn.4), omn.pred2a)
omn.pred2avar <- diag(omn.pred2amm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred2amm))
omn.pred2a <- data.frame(omn.pred2a,
                         lo = omn.pred2a$GuildAbun - exp(cmult*sqrt(omn.pred2avar)),
                         hi = omn.pred2a$GuildAbun + exp(cmult*sqrt(omn.pred2avar)))
# Week 7
omn.pred2b <- data.frame(Week = 7,
                         HabClass = factor("Edge", levels=levels(m_data_A1_omn$HabClass)),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = seq(1,15,0.01),
                         CCavgsd = mean(m_data_A1_omn$CCavgsd),
                         Point = sample(unique(m_data_A1_omn$Point), 1401, replace = T))
omn.pred2b$GuildAbun <- predict(omn.4, omn.pred2b, type="response", re.form=NA)
omn.pred2bmm <- model.matrix(terms(omn.4), omn.pred2b)
omn.pred2bvar <- diag(omn.pred2bmm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred2bmm))
omn.pred2b <- data.frame(omn.pred2b,
                         lo = omn.pred2b$GuildAbun - exp(cmult*sqrt(omn.pred2bvar)),
                         hi = omn.pred2b$GuildAbun + exp(cmult*sqrt(omn.pred2bvar)))
# Week 13
omn.pred2c <- data.frame(Week = 13,
                         HabClass = factor("Edge", levels=levels(m_data_A1_omn$HabClass)),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = seq(1,15,0.01),
                         CCavgsd = mean(m_data_A1_omn$CCavgsd),
                         Point = sample(unique(m_data_A1_omn$Point), 1401, replace = T))
omn.pred2c$GuildAbun <- predict(omn.4, omn.pred2c, type="response", re.form=NA)
omn.pred2cmm <- model.matrix(terms(omn.4), omn.pred2c)
omn.pred2cvar <- diag(omn.pred2cmm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred2cmm))
omn.pred2c <- data.frame(omn.pred2c,
                         lo = omn.pred2c$GuildAbun - exp(cmult*sqrt(omn.pred2cvar)),
                         hi = omn.pred2c$GuildAbun + exp(cmult*sqrt(omn.pred2cvar)))
# joining
omn.pred2 <- rbind(omn.pred2a,omn.pred2b,omn.pred2c)
# plotting
omn.pred2plot <- ggplot(data=omn.pred2, mapping=aes(x=TreeDens, y=GuildAbun,
                                                    col=factor(Week), fill=factor(Week))) +
  guides(col=guide_legend(title = "Week", reverse = T),
         fill=guide_legend(title = "Week", reverse = T)) +
  scale_colour_manual(values = cbbPalette[c(7, 8, 6)]) +
  scale_fill_manual(values = cbbPalette[c(7, 8, 6)]) +
  geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lo, ymax=hi), col=NA, alpha=0.3) +
  scale_y_continuous(breaks = seq(0,40,2)) +
  scale_x_continuous(breaks = seq(0,20,2)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Omnivore detections per point count",
       x = expression(Tree~density~per~100~m^2)) 




## Abun - CCavgsd (for Edge) ##

# Week 1
omn.pred3a <- data.frame(Week = 1,
                         HabClass = factor("Edge", levels=levels(m_data_A1_omn$HabClass)),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = mean(m_data_A1_omn$TreeDens),
                         CCavgsd = seq(1,20,0.01),
                         Point = sample(unique(m_data_A1_omn$Point), 1901, replace = T))
omn.pred3a$GuildAbun <- predict(omn.4, omn.pred3a, type="response", re.form=NA)
omn.pred3amm <- model.matrix(terms(omn.4), omn.pred3a)
omn.pred3avar <- diag(omn.pred3amm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred3amm))
omn.pred3a <- data.frame(omn.pred3a,
                         lo = omn.pred3a$GuildAbun - exp(cmult*sqrt(omn.pred3avar)),
                         hi = omn.pred3a$GuildAbun + exp(cmult*sqrt(omn.pred3avar)))
# Week 7
omn.pred3b <- data.frame(Week = 7,
                         HabClass = factor("Edge", levels=levels(m_data_A1_omn$HabClass)),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = mean(m_data_A1_omn$TreeDens),
                         CCavgsd = seq(1,20,0.01),
                         Point = sample(unique(m_data_A1_omn$Point), 1901, replace = T))
omn.pred3b$GuildAbun <- predict(omn.4, omn.pred3b, type="response", re.form=NA)
omn.pred3bmm <- model.matrix(terms(omn.4), omn.pred3b)
omn.pred3bvar <- diag(omn.pred3bmm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred3bmm))
omn.pred3b <- data.frame(omn.pred3b,
                         lo = omn.pred3b$GuildAbun - exp(cmult*sqrt(omn.pred3bvar)),
                         hi = omn.pred3b$GuildAbun + exp(cmult*sqrt(omn.pred3bvar)))
# Week 13
omn.pred3c <- data.frame(Week = 13,
                         HabClass = factor("Edge", levels=levels(m_data_A1_omn$HabClass)),
                         Observer = factor("KT", levels = levels(m_data_A1_omn$Observer)),
                         TreeDens = mean(m_data_A1_omn$TreeDens),
                         CCavgsd = seq(1,20,0.01),
                         Point = sample(unique(m_data_A1_omn$Point), 1901, replace = T))
omn.pred3c$GuildAbun <- predict(omn.4, omn.pred3c, type="response", re.form=NA)
omn.pred3cmm <- model.matrix(terms(omn.4), omn.pred3c)
omn.pred3cvar <- diag(omn.pred3cmm %*% vcov(omn.4)[["cond"]] %*% t(omn.pred3cmm))
omn.pred3c <- data.frame(omn.pred3c,
                         lo = omn.pred3c$GuildAbun - exp(cmult*sqrt(omn.pred3cvar)),
                         hi = omn.pred3c$GuildAbun + exp(cmult*sqrt(omn.pred3cvar)))
# joining
omn.pred3 <- rbind(omn.pred3a,omn.pred3b,omn.pred3c)
# plotting
omn.pred3plot <- ggplot(data=omn.pred3, mapping=aes(x=CCavgsd, y=GuildAbun,
                                                    col=factor(Week), fill=factor(Week))) +
  guides(col=guide_legend(title = "Week", reverse = T),
         fill=guide_legend(title = "Week", reverse = T)) +
  scale_colour_manual(values = cbbPalette[c(7,8,6)]) +
  scale_fill_manual(values = cbbPalette[c(7,8,6)]) +
  geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lo, ymax=hi), col=NA, alpha=0.3) +
  scale_y_continuous(breaks = seq(0,40,2)) +
  scale_x_continuous(breaks = seq(0,20,4)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Omnivore detections per point count",
       x = "Canopy heterogeneity") 


## Fig6 ##


fig6 <- omn.pred1plot + 
  (omn.pred2plot + theme(axis.title.y = element_blank())) +
  (omn.pred3plot + theme(axis.title.y = element_blank())) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides="collect") & 
  theme(plot.tag = element_text(size = 16)) 

ggsave("Fig6.png", fig6, 
       width = 18, height = 16, units = "cm", dpi=300)



### ###


### Visualising results: (Fig7) Analysis 3: ordination ####

ordplotPB <- ggplot(ordPB.sp) + coord_fixed() +
  geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
             position = position_dodge2(width=0.2)) +
  geom_point(data=ordPB.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
  geom_segment(data=ordPB.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
               arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
  geom_text(data=ordPB.ef, aes(x=PC1-0.01, y=PC2-0.01, label=HabVar), 
            size=3, colour="black") +
  geom_text(aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
            size=3.5, position = position_dodge2(width=0.2)) + 
  scale_colour_manual(values = Set1) +
  scale_fill_manual(values = Set3[c(10,12)],
                    labels = c("Invert.-feeding","Omnivore")) +
  scale_shape_manual(values = c(22,23),
                     labels = c("Invert.-feeding","Omnivore")) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_x_continuous(breaks = seq(-1,1,0.25)) +
  guides(col=guide_legend(title = "Habitat class"),
         fill=guide_legend(title = "Feeding guild"),
         shape=guide_legend(title = "Feeding guild")) 

ordplotPM <- ggplot(ordPM.sp) + coord_fixed() +
  geom_point(aes(x=PC1, y=PC2, shape=GuildFeed, fill=GuildFeed), size=4,
             position = position_dodge2(width=0.2)) +
  geom_point(data=ordPM.si, aes(x=PC1, y=PC2, col=HabClass), size=2) +
  geom_segment(data=ordPM.ef, aes(x=0,xend=PC1,y=0,yend=PC2), size=1,
               arrow = arrow(length = unit(0.25,"cm")), colour="black") + 
  geom_text(data=ordPM.ef, aes(x=PC1+0.01, y=PC2+0.02, label=HabVar), 
            size=3, colour="black") +
  geom_text(aes(x=PC1, y=PC2, label=Spec_code), fontface = "bold", vjust=1.5, 
            size=3.5, position = position_dodge2(width=0.2)) + 
  scale_colour_manual(values = Set1) +
  scale_fill_manual(values = Set3[c(10,12)],
                    labels = c("Invert.-feeding","Omnivore")) +
  scale_shape_manual(values = c(22,23),
                     labels = c("Invert.-feeding","Omnivore")) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_x_continuous(breaks = seq(-1,1,0.25)) +
  guides(col=guide_legend(title = "Habitat class"),
         fill=guide_legend(title = "Feeding guild"),
         shape=guide_legend(title = "Feeding guild")) 


## Fig7 ##

fig7 <- ordplotPB / ordplotPM + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16)) 

ggsave("Fig7.png", fig7,
       width = 18, height = 24, units = "cm", dpi=300)


### ###

### Visualising results: (Fig8) Analysis 4: predation ####


cat.2pred1a <- data.frame(Week = seq(5,12,0.01),
                          HabClass = factor("Interior", levels = levels(cat_data$HabClass)))
cat.2pred1a$Pred <- predict(cat.2, cat.2pred1a, type="response", re.form=NA)

cat.2pred1b <- data.frame(Week = seq(5,12,0.01),
                          HabClass = factor("Road", levels = levels(cat_data$HabClass)))
cat.2pred1b$Pred <- predict(cat.2, cat.2pred1b, type="response", re.form=NA)

cat.2pred1c <- data.frame(Week = seq(5,12,0.01),
                          HabClass = factor("Edge", levels = levels(cat_data$HabClass)))
cat.2pred1c$Pred <- predict(cat.2, cat.2pred1c, type="response", re.form=NA)

cat.2pred1 <- rbind(cat.2pred1a,cat.2pred1b,cat.2pred1c)


cat.2pred1plot <- ggplot(data=cat.2pred1,
                         mapping=aes(x=Week, y=Pred,
                                     col=factor(HabClass), fill=factor(HabClass),
                                     group = HabClass)) +
  guides(col=guide_legend(title="Habitat class"),
         fill=guide_legend(title="Habitat class")) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_line(size=1.5) +
  geom_point(data=cat_data, aes(y = Bird/(Bird+OK)), position = "jitter", 
             alpha=0.25, size=3) +
  coord_cartesian(xlim=c(1,13), ylim=c(0,0.5)) +
  scale_x_continuous(breaks=1:13) +
  labs(y = "Predation rate")



cat.2pred2a <- data.frame(Week = 5,
                          HabClass = sample(unique(cat_data$HabClass), 1201, replace = T))
cat.2pred2a$Pred <- predict(cat.2, cat.2pred2a, type="response", re.form=NA)

cat.2pred2b <- data.frame(Week = 9,
                          HabClass = sample(unique(cat_data$HabClass), 1201, replace = T))
cat.2pred2b$Pred <- predict(cat.2, cat.2pred2b, type="response", re.form=NA)

cat.2pred2c <- data.frame(Week = 12,
                          HabClass = sample(unique(cat_data$HabClass), 1201, replace = T))
cat.2pred2c$Pred <- predict(cat.2, cat.2pred2c, type="response", re.form=NA)

cat.2pred2 <- rbind(cat.2pred2a,cat.2pred2b,cat.2pred2c)


cat.2pred2plot <- ggplot(data=cat.2pred2,
                         mapping=aes(x=HabClass, y=Pred,
                                     col=factor(Week), fill=factor(Week))) +
  guides(col=guide_legend(title="Week"),
         fill=guide_legend(title="Week")) +
  scale_colour_manual(values = cbbPalette[c(7,8,6)]) +
  scale_fill_manual(values = cbbPalette[c(7,8,6)]) +
  geom_point(fun = "mean", stat = "summary", 
             size = 4, stroke=1.5, shape=21, col="black", 
             position = position_dodge(width=0.4)) +
  geom_point(data=filter(cat_data, Week==5|Week==9|Week==12), 
             aes(y = Bird/(Bird+OK)), position = "jitter", 
             alpha=0.25, size=2) +
  coord_cartesian(ylim=c(0,0.5)) +
  labs(x = "Habitat class", y = "Predation rate")


cat.2pred1plot + 
  (cat.2pred2plot + theme(axis.title.y = element_blank())) + 
  plot_layout(guides = "collect", widths = c(3,2)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16)) -> fig8a


# this model is not much better than previous, but interaction UndDiv*Week is marginally
# significant so it can help us visualise the relationship that we can potentially unravel
# with a better dataset
#
# cat.3pred1a <- data.frame(Week = 5,
#                           HabClass = factor("Interior", levels = levels(cat_data$HabClass)),
#                           UndDiv = seq(0,1.5,0.0001))
# cat.3pred1a$Pred <- predict(cat.3, cat.3pred1a, type="response", re.form=NA)
# 
# cat.3pred1b <- data.frame(Week = 12,
#                           HabClass = factor("Interior", levels = levels(cat_data$HabClass)),
#                           UndDiv = seq(0,1.5,0.0001))
# cat.3pred1b$Pred <- predict(cat.3, cat.3pred1b, type="response", re.form=NA)
# 
# cat.3pred1 <- rbind(cat.3pred1a,cat.3pred1b)
# 
# cat.3pred1plot <- ggplot(data=cat.3pred1,
#                          mapping=aes(x=UndDiv, y=Pred,
#                                      col=factor(Week), fill=factor(Week))) +
#   guides(col=guide_legend(title="Week"),
#          fill=guide_legend(title="Week")) +
#   scale_colour_manual(values=cbbPalette[c(2,6)]) +
#   scale_fill_manual(values=cbbPalette[c(2,6)]) +
#   geom_line(size=1.5) +
#   coord_cartesian(ylim=c(0,0.5)) +
#   geom_point(data=filter(cat_data, Week==5 | Week==12), aes(y = Bird/(Bird+OK)),
#              position = "jitter", alpha=0.25, size=2) 
# 
# ggsave("cat.3pred.png", cat.3pred1plot, 
#        width = 4, height = 5, units = "in", dpi=300)


# correlation with bird abundance

cat_corr <- left_join(select(cat_data, c(1,2,25:28)), 
                      select(b_data_A1_all, c(1,3,8)), 
                      by=c("Week","Point")) %>% 
  mutate(Pred = Bird/(Bird+OK),
         BirdAbun = ifelse(is.na(BirdAbun), 0, BirdAbun))

plot(cat_corr$BirdAbun, cat_corr$Pred)
cor.test(cat_corr$BirdAbun, cat_corr$Pred, method="pearson")
# perhaps not linear relation?
cat_corrplot <- ggplot(data=cat_corr,
                       mapping=aes(x=BirdAbun, y=Pred)) +
  coord_cartesian(ylim=c(0,0.5)) +
  scale_x_continuous(breaks=seq(0,40,4)) +
  geom_point(alpha=0.5, size=2, colour="#655B1B", position = "jitter") +
  labs(x = "Bird abundance", y = "Predation rate")


# with only insectivores

cat_corr2 <- left_join(select(cat_data, c(1,2,25:28)), 
                       select(m_data_A1_inv, c(1,3,9)), 
                       by=c("Week","Point")) %>% 
  mutate(Pred = Bird/(Bird+OK),
         GuildAbun = ifelse(is.na(GuildAbun), 0, GuildAbun))

plot(cat_corr2$GuildAbun, cat_corr2$Pred)
cor.test(cat_corr2$GuildAbun, cat_corr2$Pred, method="pearson")

cat_corr2plot <- ggplot(data=cat_corr2,
                        mapping=aes(x=GuildAbun, y=Pred)) +
  coord_cartesian(ylim=c(0,0.5)) +
  scale_x_continuous(breaks=seq(0,40,4)) +
  geom_point(alpha=0.5, size=2, colour="#655B1B", position = "jitter") +
  labs(x = "Invert.-feeder abundance", y = "Predation rate")


# with only omnivores
cat_corr3 <- left_join(select(cat_data, c(1,2,25:28)), 
                       select(m_data_A1_omn, c(1,3,9)), 
                       by=c("Week","Point")) %>% 
  mutate(Pred = Bird/(Bird+OK),
         GuildAbun = ifelse(is.na(GuildAbun), 0, GuildAbun))

plot(cat_corr3$GuildAbun, cat_corr3$Pred)
cor.test(cat_corr3$GuildAbun, cat_corr3$Pred, method="pearson")

cat_corr3plot <- ggplot(data=cat_corr3,
                        mapping=aes(x=GuildAbun, y=Pred)) +
  coord_cartesian(ylim=c(0,0.5)) +
  scale_x_continuous(breaks=seq(0,40,4)) +
  geom_point(alpha=0.5, size=2, colour="#655B1B", position = "jitter") +
  labs(x = "Omnivore abundance", y = "Predation rate")



cat_corrplot +
  (cat_corr2plot + theme(axis.title.y = element_blank())) + 
  (cat_corr3plot + theme(axis.title.y = element_blank())) -> fig8b


## Fig8 ##

fig8 <- (fig8a) / (fig8b & theme(axis.title.x = element_text(size=9))) + 
  plot_layout(guides = "collect", heights = c(2,1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) 

ggsave("Fig8.png", fig8, 
       width = 18, height = 24, units = "cm", dpi=300)


### ###



### Visualising results: secondary graphs ####

### Fig.3: average bird detections per point count over the weeks and points ###

fig3 <- 
  (ggplot(b_data_A1_all, aes(x=Week, y=BirdAbun, group=Week)) +
     geom_point(col="#213B73", size=3.5, stat = "summary",
                fun.data = mean_se, fun.args = list(mult=2)) +
     stat_summary(col="#213B73", geom = "errorbar", size=1.5, width=0.2, 
                  fun.data = mean_se, fun.args = list(mult=2)) +
     scale_y_continuous(breaks = seq(0,40,4)) +
     scale_x_continuous(breaks = seq(1,13,1)) +
     coord_cartesian(ylim = c(0,16)) +
     labs(y = "Bird detections per point count")) /
  (ggplot(b_data_A1_all, aes(x=Point, y=BirdAbun)) + 
     geom_boxplot(fill="#2E799E", alpha=0.6, size=0.9, outlier.alpha = 1) +
     scale_y_continuous(breaks = seq(0,40,4)) +
     labs(y = "Bird detections per point count") +
     theme(axis.text.x = element_text(angle = 45, hjust = 0.75, size = 8))) + # plots median
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) 

ggsave("Fig3.png", fig3, 
       width = 18, height = 24, units = "cm", dpi=300)



### Fig.9: correlations between canopy measures ###

# ( ggplot(h_data_A1) + geom_point(aes(CCavg, TreeDens), size=2, colour="#655B1B") +
#   scale_y_continuous(breaks = seq(0,20,2)) + scale_x_continuous(breaks = seq(0,100,10)) +
#   labs(x = "Canopy cover", y = expression(Tree~density~per~100~m^2)) ) +
# ( ggplot(h_data_A1) + geom_point(aes(CCavg, CCavgsd), size=2, colour="#655B1B") +
#   scale_y_continuous(breaks = seq(0,24,4)) + scale_x_continuous(breaks = seq(0,100,10)) +
#   labs(x = "Canopy cover", y = "Canopy heterogeneity") ) +
# ( ggplot(h_data_A1) + geom_point(aes(CCavgsd, TreeDens), size=2, colour="#655B1B") +
#   scale_x_continuous(breaks = seq(0,24,4)) + scale_y_continuous(breaks = seq(0,20,2)) +
#   labs(x = "Canopy heterogeneity", y = expression(Tree~density~per~100~m^2)) ) +
# ( ggplot(h_data_A1) + geom_point(aes(TreePropDeci, TreeDens), size=2, colour="#655B1B") +
#   scale_x_continuous(breaks = seq(0,100,10)) + scale_y_continuous(breaks = seq(0,20,2)) +
#   labs(x = "Proportion of deciduous trees", y = expression(Tree~density~per~100~m^2)) ) +
# plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) -> fig9

# with correlation lines
( ggscatter(h_data_A1, x="CCavg", y="TreeDens", 
            cor.method="pearson", cor.coef=T,
            conf.int=T, add="reg.line", size=2, color="#655B1B") +
    scale_y_continuous(breaks = seq(0,20,2)) + scale_x_continuous(breaks = seq(0,100,10)) +
    labs(x = "Canopy cover", y = expression(Tree~density~per~100~m^2)) ) +
  ( ggscatter(h_data_A1, x="CCavg", y="CCavgsd", 
              cor.method="pearson", cor.coef=T,
              conf.int=T, add="reg.line", size=2, color="#655B1B") +
      scale_y_continuous(breaks = seq(0,24,4)) + scale_x_continuous(breaks = seq(0,100,10)) +
      labs(x = "Canopy cover", y = "Canopy heterogeneity") ) +
  ( ggscatter(h_data_A1, x="CCavgsd", y="TreeDens", 
              cor.method="pearson", cor.coef=T,
              conf.int=F, size=2, color="#655B1B") +
      scale_x_continuous(breaks = seq(0,24,4)) + scale_y_continuous(breaks = seq(0,20,2)) +
      labs(x = "Canopy heterogeneity", y = expression(Tree~density~per~100~m^2)) ) +
  ( ggscatter(h_data_A1, x="TreePropDeci", y="TreeDens", 
              cor.method="pearson", cor.coef=T,
              conf.int=T, add="reg.line", size=2, color="#655B1B") +
      scale_x_continuous(breaks = seq(0,100,10)) + scale_y_continuous(breaks = seq(0,20,2)) +
      labs(x = "Proportion of deciduous trees", y = expression(Tree~density~per~100~m^2)) ) +
  plot_annotation(tag_levels = "A") & theme_simple & 
  theme(plot.tag = element_text(size = 16)) -> fig9

ggsave("Fig9.png", fig9, 
       width = 18, height = 18, units = "cm", dpi=300)



### Fig.10: trends of habvar with DOM ###

h_spruce <- read.delim("clipboard") # "PointTreeSpec" sheet

h_data_DOMtrends <- inner_join(select(h_data_A1, c(1,5,6,12,14,15)),
                               select(h_spruce, c(1,8)),
                               by = "Point")


( ggplot(h_data_DOMtrends) + 
    geom_boxplot(aes(DOM, TreeDens), 
                 fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
    scale_y_continuous(breaks = seq(0,20,2)) +
    scale_x_discrete(labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
    labs(y = expression(Tree~density~per~100~m^2),
         x = "Ground vegetation") ) + 
  ( ggplot(h_data_DOMtrends) + 
      geom_boxplot(aes(DOM, CCavgsd), 
                   fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
      scale_y_continuous(breaks = seq(0,20,2)) +
      scale_x_discrete(labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
      labs(y = "Canopy heterogeneity",
           x = "Ground vegetation") ) +
  ( ggplot(h_data_DOMtrends) + 
      geom_boxplot(aes(DOM, TreePropDeci), 
                   fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
      scale_y_continuous(breaks = seq(0,100,10)) +
      scale_x_discrete(labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
      coord_cartesian(ylim = c(0,100)) +
      labs(y = "Proportion of deciduous trees",
           x = "Ground vegetation") ) +
  ( ggplot(h_data_DOMtrends) + 
      geom_boxplot(aes(DOM, Pice_abie), 
                   fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
      scale_y_continuous(breaks = seq(0,100,10)) +
      scale_x_discrete(labels = c("Bare","Graminoid","Moss","Rubus","Vaccinium")) +
      coord_cartesian(ylim = c(0,100)) +
      labs(y = "Number of spruce trees in canopy",
           x = "Ground vegetation") ) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16),
        axis.text.x = element_text(size=6)) -> fig10

ggsave("Fig10.png", fig10, 
       width = 16, height = 21, units = "cm", dpi=300)





### Fig.11&12: Detections by observer ###

( ggplot(m_data_A1_all, aes(Observer, BirdAbun)) +
    geom_boxplot(fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
    labs(y = "Bird detections per point count") +
    scale_y_continuous(breaks = seq(0,40,4)) +
    coord_cartesian(ylim=c(0,28)) ) +
  ( ggplot(b_rawdata_all, aes(Observer, BirdAbun)) +
      geom_boxplot(fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
      labs(y = "Bird detections per point count") +
      scale_y_continuous(breaks = seq(0,60,8)) +
      coord_cartesian(ylim=c(0,56)) +
      theme(axis.title.y = element_blank()) ) +
  ( ggplot(m_data_A1_all, aes(Observer, BirdAbun, group=Observer)) +
      geom_bar(stat="summary", fun="sum", fill="#655B1B", alpha=0.6, size=0.7, col="black") +
      labs(y = "Overall bird detections") +
      coord_cartesian(ylim=c(0,4000)) +
      scale_y_continuous(breaks = seq(0,4000,500)) +
      geom_hline(yintercept = sum(m_data_A1_all$BirdAbun), linetype="dashed", size=1.2) ) +
  ( ggplot(b_rawdata_all, aes(Observer, BirdAbun, group=Observer)) +
      geom_bar(stat="summary", fun="sum", fill="#655B1B", alpha=0.6, size=0.7, col="black") +
      labs(y = "Overall bird detections") +
      coord_cartesian(ylim=c(0,8000)) +
      scale_y_continuous(breaks = seq(0,9000,1000)) +
      geom_hline(yintercept = sum(b_rawdata_all$BirdAbun), linetype="dashed", size=1.2) +
      theme(axis.title.y = element_blank()) ) +
  plot_layout(ncol=4) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) -> fig11

ggsave("Fig11.png", fig11, 
       width = 20, height = 10, units = "cm", dpi=300)



( ggplot(m_data_A1_all) +
    geom_col(aes(Week, BirdAbun, fill=Observer)) +
    scale_fill_viridis_d(end=0.5, direction=-1) +
    scale_x_continuous(breaks = 1:13) +
    scale_y_continuous(breaks = seq(0,500,50)) +
    labs(y = "Number of detections") ) +
  ( ggplot(m_data_A1_all) +
      geom_col(aes(Point, BirdAbun, fill=Observer)) +
      scale_fill_viridis_d(end=0.5, direction=-1) +
      scale_y_continuous(breaks = seq(0,200,20)) +
      coord_cartesian(ylim = c(0,140)) +
      labs(y = "Number of detections") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) ) +
  ( ggplot(b_rawdata_all) +
      geom_col(aes(Week, BirdAbun, fill=Observer)) +
      scale_fill_viridis_d(end=0.5, direction=-1) +
      scale_x_continuous(breaks = 1:13) +
      scale_y_continuous(breaks = seq(0,1000,100)) +
      labs(y = "Number of detections") ) +
  ( ggplot(b_rawdata_all) +
      geom_col(aes(Point, BirdAbun, fill=Observer)) +
      scale_fill_viridis_d(end=0.5, direction=-1) +
      scale_y_continuous(breaks = seq(0,500,50)) +
      coord_cartesian(ylim = c(0,300)) +
      labs(y = "Number of detections") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) ) +
  plot_layout(nrow=2, widths = c(2,3), guides="collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16)) -> fig12

ggsave("Fig12.png", fig12, 
       width = 21, height = 19, units = "cm", dpi=300)



### Fig.13: Detections by visual/auditory ###

b_data_visaud <- b_data %>% 
  mutate(Detection = factor(ifelse(Seen==1 & Heard==0, "Visual", ifelse(Seen==1 & Heard==1,
                                                                        "Both", "Auditory")),
                            levels = c("Visual","Both","Auditory"))) %>% 
  select(c(1,5,11,12,16,17,27))

( ggplot(b_data_visaud) +
    geom_col(aes(Week, Number, fill=Detection)) +
    scale_fill_viridis_d(direction=-1) +
    scale_x_continuous(breaks = 1:13) +
    scale_y_continuous(breaks = seq(0,500,50)) +
    labs(y = "Number of detections") ) +
  ( ggplot(b_data_visaud) +
      geom_col(aes(Point, Number, fill=Detection)) +
      scale_fill_viridis_d(direction=-1) +
      scale_y_continuous(breaks = seq(0,200,20)) +
      coord_cartesian(ylim = c(0,140)) +
      labs(y = "Number of detections") +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.75, size = 8)) ) +
  plot_layout(nrow=2, guides="collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) -> fig13

ggsave("Fig13.png", fig13, 
       width = 18, height = 16, units = "cm", dpi=300)


### changes in abundances of species over weeks ###

(ggplot(filter(b_data_A1, 
               Genus=="Sitta"|
                 Genus=="Parus"), 
        aes(x=Week, y=Abundance, col=Genus), group=Week) +
    geom_point(size=3.5, stat = "summary",
               fun.data = mean_se, fun.args = list(mult=2)) +
    stat_summary(geom = "errorbar", size=1.5, width=0.3, 
                 fun.data = mean_se, fun.args = list(mult=2)) +
    scale_color_manual(values = cbbPalette[c(8,4)],
                       labels = c("P. major","S. europaea")) +
    guides(col = guide_legend(title = "Species")) +
    scale_y_continuous(breaks = seq(0,20,2)) +
    scale_x_continuous(breaks = seq(1,13,1)) +
    coord_cartesian(ylim = c(0,8), xlim = c(1,13)) +
    labs(y = "Detections per point count")) +
  (ggplot(filter(b_data_A1, 
                 Genus=="Fringilla"|
                   Genus=="Cyanistes"), 
          aes(x=Week, y=Abundance, col=Genus), group=Week) +
     geom_point(size=3.5, stat = "summary",
                fun.data = mean_se, fun.args = list(mult=2)) +
     stat_summary(geom = "errorbar", size=1.5, width=0.3, 
                  fun.data = mean_se, fun.args = list(mult=2)) +
     scale_color_manual(values = cbbPalette,
                        labels = c("C. caeruleus","F. coelebs")) +
     guides(col = guide_legend(title = "Species")) +
     scale_y_continuous(breaks = seq(0,20,2)) +
     scale_x_continuous(breaks = seq(1,13,1)) +
     coord_cartesian(ylim = c(0,8), xlim = c(1,13)) +
     labs(y = "Detections per point count") +
     theme(axis.title.y = element_blank())) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16),
        legend.text = element_text(face = "italic")) -> figSpAbWe


ggsave("FigSpAbWe.png", figSpAbWe, 
       width = 28, height = 18, units = "cm", dpi=300)

### ###


