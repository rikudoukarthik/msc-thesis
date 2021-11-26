library(tidyverse)
library(vegan)
library(data.table)

habdata <- read.delim("clipboard")
summary(forhabcalc)
head(habdata)

forhabcalc <- habdata %>% filter(Week == "13") %>% 
  select(-c(2,7:10))
head(forhabcalc)

# number and density ####
habcalcPart1 <- forhabcalc %>% select(-3) %>% filter(Layer != "H") %>%
  mutate(Number = as.numeric(Number)) %>%
           group_by(Point, Layer) %>%
           summarise(Abundance = sum(Number)) %>% 
  pivot_wider(id_cols = Point, names_from = Layer, values_from = Abundance) %>% 
  rename(TreeTot = C, UndTot = U) %>% 
  mutate(UndTot = replace_na(UndTot, 0)) %>% 
  mutate(TreeDens = (TreeTot/(pi*20^2))*100, UndDens = (UndTot/(pi*20^2))*100) # density x100
head(habcalc)

# richness and diversity ####
head(forhabcalc)
habrichcalc <- forhabcalc %>% select(-4) %>% filter(Layer != "H") %>%
  group_by(Point, Layer) %>%
  summarise(Richness = n()) %>% 
  pivot_wider(id_cols = Point, names_from = Layer, values_from = Richness) %>% 
  rename(TreeRich = C, UndRich = U) %>% 
  mutate(UndRich = replace_na(UndRich, 0))
head(habrichcalc)

# mutate(across(everything(), ~replace_na(.x,0))) %>% # or replace(is.na(.), 0) 

habdivcalc <- forhabcalc %>% select(-1) %>% filter(Layer != "H") %>% 
  complete(Point, Layer, fill = list(Number = '0')) %>%
  fill(Species) %>%
  pivot_wider(names_from = Species, values_from = Number,  
              values_fill = list(Number = '0'))

# earlier attempt which excludes some rows with 0 values... #
# habdivcalc <- forhabcalc %>% select(-1) %>% filter(Layer != "H") %>% 
#  group_by(Point, Layer, Species) %>% 
#  mutate(rn = 1:n()) %>% 
#  pivot_wider(names_from = Species, values_from = Number, 
#              values_fill = list(Number = "0")) %>% select(-3) 

habdivcalc_C <- habdivcalc %>% mutate(across(Lari_deci:Quer_rubr,.fns = as.numeric)) %>% 
  filter(Layer == "C") %>% ungroup %>% select(-2) %>% 
  column_to_rownames(var="Point") %>% as.matrix %>%  
  diversity(index = "simpson") %>% as.data.frame() %>% 
  rownames_to_column(var = "Point")

habdivcalc_U <- habdivcalc %>% mutate(across(Lari_deci:Quer_rubr,.fns = as.numeric)) %>%
  filter(Layer == "U") %>% ungroup %>% select(-2) %>% 
  column_to_rownames(var="Point") %>% as.matrix %>%  
  diversity(index = "simpson") %>% as.data.frame() %>% 
  rownames_to_column(var = "Point")

habcalcPart2 <- inner_join(habdivcalc_C, habdivcalc_U, by = "Point")
colnames(habcalcPart2) <- c("Point", "TreeDiv", "UndDiv")

# reprex 2 vegan::diversity() ####
reprex <- forhabcalc %>% select(-1) %>% filter(Layer != "H") %>% 
  group_by(Point, Layer, Species) %>% 
  mutate(rn = 1:n()) %>% 
  pivot_wider(names_from = Species, values_from = Number, 
              values_fill = list(Number = "0")) %>% select(-3)
dput(reprex[1:10,1:7])

reprex_C <- reprex %>% mutate(across(Lari_deci:Quer_rubr,.fns = as.numeric)) %>% 
  filter(Layer == "C") %>% ungroup %>% select(-2) %>% 
  column_to_rownames(var="Point") %>% as.matrix %>% 
  diversity(index = "simpson")
diversity(as.numeric(reprex_C[,1:14]))

# reprex 1 pivot_wider ####

rep_example <- forhabcalc %>% select(-1) %>% filter(Layer != "H") 
dput(rep_example[1:15])

# not wide
rep_example %>% group_by(Point, Layer) %>% 
  mutate(Number = as.numeric(Number)) %>% 
  distinct() %>% 
  mutate(rn = 1:n()) %>% 
  pivot_wider(id_cols = c(Point, Layer, rn), names_from = Species, values_from = Number)

# list-cols
rep_example %>% group_by(Point, Layer) %>% 
  mutate(Number = as.numeric(Number)) %>% 
  pivot_wider(id_cols = c(Point, Layer), names_from = Species, values_from = Number)

# ideal
r1 <- c("P03", "P03", "P06", "P06", "P07", "P07")
r2 <- rep(c("C", "U"), 3)
cbind(r1,r2)
r3 <- as.numeric(c(21, 0, 3, 0, 0, 0))
r4 <- as.numeric(c(17, 0, 28, 0, 3, 0))
r5 <- as.numeric(c(5, 0, 28, 0, 20, 0))
r6 <- as.numeric(c(0, 3, 0, 0, 1, 0))
r7 <- as.numeric(c(1, 0, 0, 0, 1, 0))
r8 <- as.numeric(c(0, 1, 0, 0, 0, 0))
rep_example_ideal <- cbind(r1,r2,r3,r4,r5,r6,r7,r8)
colnames(rep_example_ideal) <-  c("Point", "Layer", "Lari_deci", "Quer_rope", 
                              "Pinu_sylv", "Sorb_aucu", "Betu_pend", "Acer_pseu")
rep_example_ideal <- as.data.frame(rep_example_ideal)
dput(rep_example_ideal)

# Ronak Shah

rep_example %>% group_by(Point, Species) %>% mutate(row = row_number()) %>%
  pivot_wider(names_from = Species, values_from = Number,
              values_fill = list(Number = "0"))

# akrun

library(data.table)
rep_example %>% mutate(rn = rowid(Point, Species)) %>%
  pivot_wider(names_from = Species, values_from = Number, values_fill = list(Number = '0'))

rep_example %>% 
  complete(Point, Layer, fill = list(Number = '0')) %>%
  fill(Species) %>%  pivot_wider(names_from = Species,values_from = Number,
                                 values_fill = list(Number = '0'))

# herbs and dom ####

head(forhabcalc)

habcalcPart3 <- forhabcalc %>% select(-1) %>% 
  complete(Point, Layer, fill = list(Number = '0')) %>% fill(Species) %>% 
  filter(Layer == "H") %>% 
  mutate(Number = ifelse(Number=="DOM"|Number=="Y", Number, NA),
         Layer = ifelse(Layer=="H", Layer, NA)) %>% 
  pivot_wider(names_from = Number, values_from = Species) %>%
  select(-2,-4) %>% 
  mutate(across(everything(), ~replace_na(.x,0))) %>% 
  mutate(Moss3 = ifelse(Y=="Moss", 1, 0)) %>% 
  select(-3) %>% 
  mutate(Moss2 = ifelse(DOM=="Moss", 1, 0),
         Moss = Moss3 + Moss2) %>% 
  select(-3, -4) 

head(habcalcPart3)

# before finding out `complete` which retains those 0 rows...
# h1 <- c("P06","P07","P08", "P11", "P24", "P30", "P42")
# h2 <- as.character(rep(NA, 7))
# h3 <- as.numeric(rep(0,7))
# h0 <- cbind(h1, h2, h3)
# colnames(h0) <- c("Point", "DOM", "Moss")
# habherbs <- rbind(habherbs, h0) %>% arrange(Point)

# Final ####

habcalcFinal <- inner_join(habcalcPart1, habrichcalc, by = "Point") %>% 
  inner_join(habcalcPart2, by = "Point") %>% 
  inner_join(habcalcPart3, by = "Point")

write.csv(habcalcFinal, file = "habcalcFinal.csv", row.names = F)

# for ordination ####
# taking only Canopy right now
habord_C <- habdivcalc %>% filter(Layer == "C") %>% select(-2)

write.csv(habord_C, file = "habord_C.csv", row.names = F)
