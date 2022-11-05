### Data import and modification for analysis ###

library(tidyverse)
library(lubridate)

# Birds -----------------------------------------------------------------------------

# Bird species codes from "data_Birds_2021Feb26.xlsx"; Sheet "Species"; 
# columns A:E (all columns); rows 1:83 (all rows)
birds_codes <- readxl::read_xlsx("data/data_Birds_2021Feb26.xlsx", "Species") %>% 
  mutate(Migration75 = factor(Migration75),
         # making true NAs
         Species = if_else(Species == "NA", NA_character_, Species))


# Bird data from "data_Birds_2021Feb26.xlsx"; Sheet "Birds"; columns A:V (all columns); 
# rows 1:7368 (all rows)

# Modifying bird raw data:
# - Joining Genus and Species names to bird data
# - Excluding UNID (implicit result of the innerjoin for species codes) (7323/7367 left)
# - Excluding Sky flyovers (6752/7323 rows left)
# - Excluding birds that did not enter 30m radius (3208/7323; 2840/6752 left)

# This altogether results in exclusion of 61.5% of observations, i.e., it is 38.5% of birds0

birds0 <- readxl::read_xlsx("data/data_Birds_2021Feb26.xlsx", "Birds")

birds1 <- birds0 %>% 
  inner_join(birds_codes, by = "Spec_code") %>% 
  filter(Flyover != "2", Within_30m == 1)

100 * length(birds1$Spec_code) / length(birds0$Spec_code)


# - Removing columns to retain Week, Observer, Point, Weather, Wind, Visibility, Number
# - Adding column CoD "Count of Day" to test for effect of time of day and to help 
#   distinguish from true observer effect
# - Summarising Number column to Abundance per species (grouping by Species)
# - Adding columns Period and Days to represent the time parameter in my study
# - Calculating relative abundance of each species in each Point every week
#   in relation to total abundance of that species in all Points that week
birds2 <- birds1 %>% 
  mutate(Date = as_date(Date),
         StartTime = (str_split(StartTime, " ", simplify = T)[, 2]),
         EndTime = (str_split(EndTime, " ", simplify = T)[, 2])) %>%
  # Count of Day
  group_by(Week, Date, Observer) %>% 
  arrange(StartTime, by_group = T) %>%
  mutate(CoD = as.integer(factor(StartTime))) %>% 
  # summarising per species to get total count instead of multiple presence rows
  arrange(Week, Date, Observer) %>% 
  group_by(Point, Week, Date, Observer, CoD, Weather, Wind, Visibility, Spec_code) %>% 
  summarise(GuildFeed = GuildFeed,
            Migration75 = Migration75,
            Spec_Abun = sum(Number)) %>%
  ungroup() %>% 
  # Days and Period
  mutate(Days = (Date - as_date("2020-07-07")) %>% as.numeric(), # start of sampling
         Period = as.factor(case_when(Week %in% 1:4 ~ 1,
                                      Week %in% 5:8 ~ 2,
                                      Week %in% 9:13 ~ 3))) %>% 
  # since each Point sampled only once per week, following code is okay
  # calculating relative abundance for species in each point count compared to that species
  # in all points that week
  group_by(Week, Spec_code) %>% 
  mutate(Spec_RelAbun = round(Spec_Abun / sum(Spec_Abun), 4),
         # mean across points that species present in (<= 32)
         Spec_HabSel = round(Spec_Abun / mean(Spec_Abun), 4), 
         # mean across all points
         Spec_HabSel2 = round(Spec_Abun / (sum(Spec_Abun) / 32), 4)) 

# HabSel2 makes more sense because it considers Points with 0 abundance, but it works
# only for the levels of All Birds and Guilds. For individual species, it might work
# for ones super abundant and widespread, but otherwise it degenerates into not such a
# meaningful measure. In such cases, using HabSel might be justified, as we already know
# the species is absent from many Points, then we can simply compare preferences among
# the Points it IS present in.


# Modifying for different levels of response:

# all species
birds_all <- birds2 %>% 
  group_by(Period, Week, Days, Observer, Point, CoD, Weather) %>% 
  summarise(All_Abun = sum(Spec_Abun)) %>% 
  group_by(Week) %>% 
  mutate(All_RelAbun = round(All_Abun / sum(All_Abun), 4),
         All_HabSel = round(All_Abun / mean(All_Abun), 4), 
         All_HabSel2 = round(All_Abun / (sum(All_Abun) / 32), 4)) %>% 
  ungroup()

# individual guilds
birds_guild <- birds2 %>% 
  group_by(Period, Week, Days, Observer, Point, CoD, Weather, GuildFeed) %>% 
  summarise(Guild_Abun = sum(Spec_Abun)) %>%
  group_by(Week) %>% 
  mutate(Guild_RelAbun = round(Guild_Abun / sum(Guild_Abun), 4),
         Guild_HabSel = round(Guild_Abun / mean(Guild_Abun), 4), 
         Guild_HabSel2 = round(Guild_Abun / (sum(Guild_Abun) / 32), 4)) %>% 
  ungroup()

# migrants and non-migrants
birds_mig <- birds2 %>% 
  group_by(Period, Week, Days, Observer, Point, CoD, Weather, Migration75) %>% 
  summarise(Migr_Abun = sum(Spec_Abun)) %>%
  group_by(Week) %>% 
  mutate(Migr_RelAbun = round(Migr_Abun / sum(Migr_Abun), 4),
         Migr_HabSel = round(Migr_Abun / mean(Migr_Abun), 4), 
         Migr_HabSel2 = round(Migr_Abun / (sum(Migr_Abun) / 32), 4)) %>% 
  ungroup()



# species detection information
birds_summary <- birds2 %>%
  group_by(Spec_code, Week) %>% 
  summarise(Tot_Points = n_distinct(Point),
            Tot_Det = sum(Spec_Abun)) %>% 
  ungroup()



# Habitat ---------------------------------------------------------------------------

# Habitat data from "data_HabVar_2021Feb26.xlsx"; Sheet "Point Descriptions"; 
# columns A:R (all columns); rows 1:33 (all rows)

# Removing unnecessary columns to leave only Latitude, Longitude, CCavg, TreeDens,
# UndDens, TreeRich, UndRich, TreeDiv, UndDiv, DOM, Moss

habvar <- readxl::read_xlsx("data/data_HabVar_2021Feb26.xlsx", "Point Descriptions") %>% 
  rename(TPD = TreePropDeci) %>% 
  # scaling variables where necessary
  mutate(CCavgscaled = scale(CCavg),
         CCavgsdscaled = scale(CCavgsd),
         TreeDensscaled = scale(TreeDens),
         TPDscaled = scale(TPD),
         UndDensscaled = scale(UndDens)) %>% 
  select(-c(CCavgvar, TreeTot, UndTot)) %>% 
  # removing cells in which I'd entered summaries
  filter(!is.na(Point))



# Merging bird and habitat data -----------------------------------------------------

# Merging bird and habitat variables into one tibble
m_all <- left_join(birds_all, habvar, by = "Point")
m_guild <- left_join(birds_guild, habvar, by = "Point")
m_mig <- left_join(birds_mig, habvar, by = "Point")


# Saving as RData to import later ---------------------------------------------------

# saving RData for use in analysis
rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("birds_all","birds_guild","birds_mig","habvar",
                    "m_all","m_guild","m_mig","birds_codes","birds_summary")), 
   pos = ".GlobalEnv")
save.image("data/01_dataimport.RData")

