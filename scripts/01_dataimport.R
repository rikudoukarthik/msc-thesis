### Data import and modification for analysis
### 
### Bird abundances versus habitat variables; 
### post-breeding habitat breadth and selectivity of birds
### 


library(tidyverse)



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