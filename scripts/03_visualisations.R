library(tidyverse)
library(stargazer)
library(Hmisc)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(magrittr)
library(bbmle)
library(tictoc)

load("data/01_dataimport.RData")
load("data/02_analysis.RData")

source("scripts/functions.R")


# setting theme
theme_set(theme_classic())

# purples: #A204B4 #D91EFA
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072A2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072A2", "#D55E00", "#CC79A7")
custom_pal <- c("#d8b365", "#f5f5f5", "#5ab4ac")
Set1 <- brewer.pal(9, name = "Set1")
Set3 <- brewer.pal(12, name = "Set3")

Gauss_coeff <- 1.96


### Summaries ####

# summary statistics etc. of some results

sum_stat <- data.frame(
  
  UNIQUE.SPEC.RAW = birds_species_raw %>% filter(!is.na(Species)) %$% n_distinct(Spec_code),
  UNIQUE.SPEC.FILT = birds_species_filt %>% filter(!is.na(Species)) %$% n_distinct(Spec_code),
  
  FORAGING.RAW = birds_species_raw %>% filter(Foraging == 1) %>% nrow(),
  FORAGING.FILT = birds_species_filt %>% filter(Foraging == 1) %>% nrow(),
  FORAGING.FILT.FRUIT = birds_species_filt %>% filter(Fruits == 1) %>% nrow(),
  
  DET.TOT.RAW = birds_species_raw %>% nrow(),
  DET.AUD.RAW = birds_species_raw %>% filter(Heard == 1 & Seen == 0) %>% nrow(),
  DET.VIS.RAW = birds_species_raw %>% filter(Heard == 0 & Seen == 1) %>% nrow(),
  DET.BOTH.RAW = birds_species_raw %>% filter(Heard == 1 & Seen == 1) %>% nrow(),
  
  DET.TOT.FILT = birds_species_filt %>% nrow(),
  DET.AUD.FILT = birds_species_filt %>% filter(Heard == 1 & Seen == 0) %>% nrow(),
  DET.VIS.FILT = birds_species_filt %>% filter(Heard == 0 & Seen == 1) %>% nrow(),
  DET.BOTH.FILT = birds_species_filt %>% filter(Heard == 1 & Seen == 1) %>% nrow(),
  
  DET.IF.RAW = birds_species_raw %>% filter(Observer == "IF") %>% nrow(),
  DET.KT.RAW = birds_species_raw %>% filter(Observer == "KT") %>% nrow(),
  DET.IF.FILT = birds_species_filt %>% filter(Observer == "IF") %>% nrow(),
  DET.KT.FILT = birds_species_filt %>% filter(Observer == "KT") %>% nrow()
  
  ) %>% 
  mutate(
    
    DET.AUD.RAW.PC = (100*DET.AUD.RAW/DET.TOT.RAW) %>% round(1),
    DET.VIS.RAW.PC = (100*DET.VIS.RAW/DET.TOT.RAW) %>% round(1),
    DET.BOTH.RAW.PC = (100*DET.BOTH.RAW/DET.TOT.RAW) %>% round(1),
    
    DET.AUD.FILT.PC = (100*DET.AUD.FILT/DET.TOT.FILT) %>% round(1),
    DET.VIS.FILT.PC = (100*DET.VIS.FILT/DET.TOT.FILT) %>% round(1),
    DET.BOTH.FILT.PC = (100*DET.BOTH.FILT/DET.TOT.FILT) %>% round(1),
    
    DET.IF.RAW.PC = (100*DET.IF.RAW/DET.TOT.RAW) %>% round(1),
    DET.KT.RAW.PC = (100*DET.KT.RAW/DET.TOT.RAW) %>% round(1),
    DET.IF.FILT.PC = (100*DET.IF.FILT/DET.TOT.FILT) %>% round(1),
    DET.KT.FILT.PC = (100*DET.KT.FILT/DET.TOT.FILT) %>% round(1)
    
  )


### Tables  ####

### Table: model for all birds ###

allAIC <- AICctab(all1, all2, all3, all4, 
                  weights = T, base = T, logLik = T)
summary(all4)
# printing AICc table for import into Word in Methods section
stargazer(as.data.frame(allAIC), 
          summary = F, type = "html", out = "outputs/summary_all.html")


### Table: models for guilds ###

summary(inv3)
summary(omn3)
invAIC <- AICctab(inv1, inv2, inv3, 
                  weights = T, base = T, logLik = T)
omnAIC <- AICctab(omn1, omn2, omn3, 
                  weights = T, base = T, logLik = T)
stargazer(rbind(as.data.frame(invAIC), as.data.frame(omnAIC)), 
          summary = F, type = "html", out = "outputs/summary_guilds.html")


### Table of species detected within 30m  ###

birds_speclist <- birds_species_filt %>%
  # only considering those IDd to species, not slashes or spuhs
  filter(Species != is.na(Species)) %>%
  distinct(Spec_code) %>% 
  left_join(birds_codes, by = "Spec_code") %>% 
  arrange(Spec_code) %>% 
  rownames_to_column("S.No.")

stargazer(birds_speclist, type = "html", summary = F, out = "outputs/specieslist.html")


### Table of points and habitat variables  ###

stargazer(habvar %>% 
            dplyr::select(-c(scaleCC, scaleCH, scaleTPD, scaleUDens, logCH, logTDens)), 
          type = "html", summary = F, out = "outputs/pointhabvar.html")


### Table of species within 30m, with guild info and overall abundances ###

birds_speclist_overall <- birds_species_filt %>% 
  filter(!is.na(GuildFeed)) %>% # removing sampling events where no birds detected within 30m
  mutate(Spec_code = fct_collapse(
    
    Spec_code, 
    Acc_sp = c("Acc_nis", "Acc_sp"),
    Ant_sp = c("Ant_pra", "Ant_tri", "Ant_sp"),
    Cer_sp = c("Cer_bra", "Cer_fam", "Cer_sp"),
    Cor_sp = c("Cor_corax", "Cor_corone", "Cor_monedula", "Cor_sp"),
    Phy_sp = c("Phy_col", "Phy_tro", "Phy_sp"),
    Reg_sp = c("Reg_reg", "Reg_ign", "Reg_sp"),
    Den_sp = c("Den_maj", "Den_med", "Den_sp"),
    Fic_sp = c("Fic_sp"),
    Fal_sp = c("Fal_tin", "Fal_sp"),
    Mot_sp = c("Mot_alb", "Mot_sp"),
    Pas_sp = c("Pas_mon", "Pas_sp"),
    Poe_sp = c("Poe_mon", "Poe_pal", "Poe_sp"),
    Tur_sp = c("Tur_mer", "Tur_phi", "Tur_pil", "Tur_vis", "Tur_sp")
    
  )) %>% 
  group_by(GuildFeed, Spec_code) %>% 
  summarise(TotDet = n()) %>% 
  left_join(birds_codes %>% dplyr::select(-Genus, -Species)) %>% 
  arrange(GuildFeed, desc(TotDet)) %>% rownames_to_column("S.No.")

write.csv(birds_speclist_overall, file = "outputs/specieslist_overall.csv", row.names = F)


### Figure 1: all birds ####

tic.clearlog()

## Week & DOM main effects ##

# empty table to predict
pred_data1 <- data.frame(Week = 1:13) %>% 
  group_by(Week) %>% 
  # need to predict for all five DOM categories at each time-step
  reframe(DOM = unique(m_all$DOM)) %>% 
  # other variables at fixed values
  mutate(Observer = factor("KT", levels = levels(m_all$Observer)),
         logCH = mean(m_all$logCH),
         logTDens = mean(m_all$logTDens),
         # for predict.glmmTMB which needs empty column of random effects
         # see https://stackoverflow.com/a/72733566/13000254, https://github.com/glmmTMB/glmmTMB/issues/766
         Point = NA,
         Days = NA)

# predictions
tic("Bootstrapped prediction 1 for all-birds model")
prediction <- boot_conf_GLMM(all4,
                             new_data = pred_data1,
                             new_data_string = "pred_data1",
                             model_data_string = "m_all",
                             nsim = 1000)
save(prediction, file = "outputs/pred1.RData")
toc(quiet = TRUE, log = TRUE) 
# load("outputs/pred1.RData")

# calculating mean and SE from bootstrapped values
for (j in 1:nrow(pred_data1)) {
  
  pred_data1$PRED.LINK[j] <- median(na.omit(prediction[,j]))
  pred_data1$SE.LINK[j] <- sd(na.omit(prediction[,j]))
  
}

pred_data1 <- pred_data1 %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  mutate(CI.U = PRED + Gauss_coeff*SE,
         CI.L = PRED - Gauss_coeff*SE) %>% 
  dplyr::select(Week, DOM, PRED, CI.L, CI.U)

# plotting Week & DOM main effects
plot1 <- ggplot(pred_data1,
                aes(x = Week, y = PRED, col = DOM, fill = DOM)) +
  scale_fill_manual(values = cbbPalette) +
  scale_colour_manual(values = cbbPalette) +
  geom_point(size = 2, position = position_dodge(0.5)) +
  geom_line(linewidth = 1.5, position = position_dodge(0.5)) +
  geom_ribbon(aes(ymin = CI.L, ymax = CI.U), 
              alpha = 0.35, colour = NA, position = position_dodge(0.5)) +
  scale_y_continuous(breaks = seq(0, 40, 4)) +
  scale_x_continuous(breaks = 1:13) +
  coord_cartesian(ylim = c(2, 18)) +
  labs(y = "Bird detections per point count") +
  guides(col = guide_legend(title = "Ground vegetation",
                            title.position = "top"),
         fill = guide_legend(title = "Ground vegetation",
                             title.position = "top")) +
  theme(legend.position = c(0.4, 0.9), 
        legend.direction = "horizontal")


## logCH-logTDens interaction ##

temp1 <- data.frame(Week = c(1, 13)) %>% 
  group_by(Week) %>% 
  # need to predict for full range of values at each time-step
  reframe(logCH = seq(min(m_all$logCH), 
                      max(m_all$logCH), 
                      0.01)) %>% 
  group_by(Week, logCH) %>% 
  # for visualising interaction between CH and TDens
  reframe(logTDens = c(min(m_all$logTDens), mean(m_all$logTDens), max(m_all$logTDens))) %>% 
  # other variables at fixed values
  mutate(DOM = factor("Moss", levels = levels(m_all$DOM)), # moss has intermediate values
         Observer = factor("KT", levels = levels(m_all$Observer))) 

temp2 <- data.frame(Week = c(1, 13)) %>% 
  group_by(Week) %>% 
  # need to predict for full range of values at each time-step
  reframe(logTDens = seq(min(m_all$logTDens), 
                         max(m_all$logTDens), 
                         0.01)) %>% 
  group_by(Week, logTDens) %>% 
  # for visualising interaction between CH and TDens
  reframe(logCH = c(min(m_all$logCH), mean(m_all$logCH), max(m_all$logCH))) %>% 
  # other variables at fixed values
  mutate(DOM = factor("Moss", levels = levels(m_all$DOM)),
         Observer = factor("KT", levels = levels(m_all$Observer)))

pred_data2 <- bind_rows(temp1, temp2, .id = "ID") %>% 
  # for predict.glmmTMB which needs empty column of random effects
  mutate(Point = NA, Days = NA)


tic("Bootstrapped prediction 2 for all-birds model (CH and TDens)")
prediction <- boot_conf_GLMM(all4,
                             new_data = pred_data2,
                             new_data_string = "pred_data2",
                             model_data_string = "m_all",
                             nsim = 1000)
save(prediction, file = "outputs/pred2.RData")
toc(quiet = TRUE, log = TRUE) # 43 min
# load("outputs/pred2.RData")

# calculating mean and SE from bootstrapped values
for (j in 1:nrow(pred_data2)) {
  
  pred_data2$PRED.LINK[j] <- median(na.omit(prediction[,j]))
  pred_data2$SE.LINK[j] <- sd(na.omit(prediction[,j]))
  
}

pred_data2 <- pred_data2 %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  mutate(CI.U = PRED + Gauss_coeff*SE,
         CI.L = PRED - Gauss_coeff*SE) %>% 
  dplyr::select(ID, Week, logCH, logTDens, PRED, CI.L, CI.U) %>% 
  filter(Week == 13) %>% 
  group_by(logCH) %>% 
  # to plot three levels of one variable separately
  mutate(TYPE = case_when(ID == 1 & logTDens == min(m_all$logTDens) ~ min(m_all$logTDens),
                          ID == 1 & logTDens == mean(m_all$logTDens) ~ mean(m_all$logTDens),
                          ID == 1 & logTDens == max(m_all$logTDens) ~ max(m_all$logTDens),
                          ID == 2 & logCH == min(m_all$logCH) ~ min(m_all$logCH),
                          ID == 2 & logCH == mean(m_all$logCH) ~ mean(m_all$logCH),
                          ID == 2 & logCH == max(m_all$logCH) ~ max(m_all$logCH)),
         TYPE.LABEL = case_when(ID == 1 & logTDens == min(m_all$logTDens) ~ "min",
                                ID == 1 & logTDens == mean(m_all$logTDens) ~ "mean",
                                ID == 1 & logTDens == max(m_all$logTDens) ~ "max",
                                ID == 2 & logCH == min(m_all$logCH) ~ "min",
                                ID == 2 & logCH == mean(m_all$logCH) ~ "mean",
                                ID == 2 & logCH == max(m_all$logCH) ~ "max")) %>% 
  mutate(CH = exp(logCH),
         TDens = exp(logTDens),
         TYPE = exp(TYPE) %>% round(2) %>% factor(),
         TYPE.LABEL = factor(TYPE.LABEL, levels = c("min", "mean", "max"))) %>% 
  ungroup()


plot2 <- ggplot(pred_data2 %>% filter(ID == 1), 
                aes(x = CH, y = PRED, col = TYPE.LABEL, fill = TYPE.LABEL)) +
  scale_fill_manual(values = cbPalette[c(3, 1, 7)]) +
  scale_colour_manual(values = cbPalette[c(3, 1, 7)]) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = CI.L, ymax = CI.U), alpha = 0.45, colour = NA) +
  scale_y_continuous(breaks = seq(0, 40, 5)) +
  guides(col = guide_legend(title = expression(Tree~density~per~100~m^2), 
                            title.position = "top", reverse = T),
         fill = guide_legend(title = expression(Tree~density~per~100~m^2), 
                             title.position = "top", reverse = T)) +
  labs(y = "Bird detections per point count",
       x = "Canopy heterogeneity") +
  theme(legend.position = c(0.6, 0.9), 
        legend.direction = "horizontal")

plot3 <- ggplot(pred_data2 %>% filter(ID == 2), 
                aes(x = TDens, y = PRED, col = TYPE.LABEL, fill = TYPE.LABEL)) +
  scale_fill_manual(values = cbPalette[c(3, 1, 7)]) +
  scale_colour_manual(values = cbPalette[c(3, 1, 7)]) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = CI.L, ymax = CI.U), alpha = 0.45, colour = NA) +
  scale_y_continuous(breaks = seq(0, 40, 5)) +
  scale_x_continuous(breaks = seq(2, 42, 4)) +
  guides(col = guide_legend(title = "Canopy heterogeneity", 
                            title.position = "top", reverse = T),
         fill = guide_legend(title = "Canopy heterogeneity", 
                             title.position = "top", reverse = T)) +
  labs(y = "Bird detections per point count",
       x = expression(Tree~density~per~100~m^2)) +
  theme(legend.position = c(0.6, 0.9), 
        legend.direction = "horizontal")


## Figure 1 ##

fig_layout <- "
AAAAAA
AAAAAA
AAAAAA
CCCDDD
CCCDDD
"

fig1_allbirds <- plot1 + plot2 + 
  (plot3 + theme(axis.title.y = element_blank())) +
  plot_layout(design = fig_layout)  +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16, margin = margin(0, 7, 0, 0)),
        plot.margin = margin(10, 10, 10, 10),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 13, vjust = 3),
        axis.title.x = element_text(size = 13, vjust = -0.75)) 

ggsave("outputs/fig1_allbirds.png", fig1_allbirds, 
       width = 26, height = 25, units = "cm", dpi = 300)


tic.log()

### Figure 2: inv-feeders ####

summary(inv3)

m_guild_inv <- m_guild %>% filter(GuildFeed == "Invertebrate")

## Moss:Week interaction ##

# empty table to predict
pred_data3 <- data.frame(Week = 1:13) %>% 
  group_by(Week) %>% 
  reframe(Moss = unique(m_guild_inv$Moss)) %>% 
  # other variables at fixed values
  mutate(Observer = factor("KT", levels = levels(m_guild_inv$Observer)),
         CoD = 5, # somewhat intermediate, but also one of the higher-bird-activity periods
         DOM = factor("Moss", levels = levels(m_guild_inv$DOM)),
         Point = NA)

# predictions
tic("Bootstrapped prediction 1 for inv. model")
prediction <- boot_conf_GLMM(inv3,
                             new_data = pred_data3,
                             new_data_string = "pred_data3",
                             model_data_string = "m_guild_inv",
                             nsim = 1000)
save(prediction, file = "outputs/pred3.RData")
toc() # 
# load("outputs/pred3.RData")

# calculating mean and SE from bootstrapped values
for (j in 1:nrow(pred_data3)) {
  
  pred_data3$PRED.LINK[j] <- median(na.omit(prediction[,j]))
  pred_data3$SE.LINK[j] <- sd(na.omit(prediction[,j]))
  
}

pred_data3 <- pred_data3 %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  mutate(CI.U = PRED + Gauss_coeff*SE,
         CI.L = PRED - Gauss_coeff*SE) %>% 
  dplyr::select(Week, Moss, PRED, CI.L, CI.U)


plot_inv1 <- ggplot(pred_data3,
                    aes(x = Week, y = PRED, col = Moss, fill = Moss)) +
  scale_fill_manual(values = custom_pal[c(1, 3)], 
                    name = "Presence of\nmoss layer", labels = c("No", "Yes")) +
  scale_colour_manual(values = custom_pal[c(1, 3)], 
                      name = "Presence of\nmoss layer", labels = c("No", "Yes")) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = CI.L, ymax = CI.U), alpha = 0.35, colour = NA) +
  scale_y_continuous(breaks = seq(0, 40, 2)) +
  scale_x_continuous(breaks = 1:13) +
  # coord_cartesian(ylim = c(2, 18)) +
  labs(y = "Invertebrate-feeder detections per point count") +
  guides(col = guide_legend(title.position = "top"),
         fill = guide_legend(title.position = "top")) +
  theme(legend.position = c(0.7, 0.9), 
        legend.direction = "horizontal")


## DOM main effect ##

# empty table to predict
pred_data4 <- data.frame(DOM = unique(m_guild_inv$DOM)) %>% 
  group_by(DOM) %>% 
  # we need one, but just to investigate relationship between DOM-Moss and Moss
  # can later filter for Moss == 0
  reframe(Moss = unique(m_guild_inv$Moss)) %>% 
  # other variables at fixed values
  mutate(Observer = factor("KT", levels = levels(m_guild_inv$Observer)),
         CoD = 5, # somewhat intermediate, but also one of the higher-bird-activity periods
         Week = 13,
         Point = NA)

# predictions
tic("Bootstrapped prediction 2 for inv. model")
prediction <- boot_conf_GLMM(inv3,
                             new_data = pred_data4,
                             new_data_string = "pred_data4",
                             model_data_string = "m_guild_inv",
                             nsim = 1000)
save(prediction, file = "outputs/pred4.RData")
toc() # 
# load("outputs/pred4.RData")

# calculating mean and SE from bootstrapped values
for (j in 1:nrow(pred_data4)) {
  
  pred_data4$PRED.LINK[j] <- median(na.omit(prediction[,j]))
  pred_data4$SE.LINK[j] <- sd(na.omit(prediction[,j]))
  
}

pred_data4 <- pred_data4 %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  mutate(CI.U = PRED + Gauss_coeff*SE,
         CI.L = PRED - Gauss_coeff*SE) %>% 
  filter(Moss == 0) %>% 
  dplyr::select(DOM, PRED, CI.L, CI.U)


plot_inv2 <- ggplot(pred_data4,
                    aes(x = DOM, y = PRED)) +
  scale_colour_manual(values = custom_pal[3]) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                linewidth = 1, width = 0.6) +
  scale_y_continuous(breaks = seq(0, 40, 1), limits = c(1, 4)) +
  labs(y = "Invertebrate-feeder detections per point count",
       x = "Ground vegetation") 


## Figure 2 ##

fig2_invbirds <- plot_inv1 + 
  (plot_inv2 + theme(axis.title.y = element_blank())) +
  plot_layout(nrow = 1, widths = c(2, 1)) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16, margin = margin(0, 7, 0, 0)),
        plot.margin = margin(10, 10, 10, 10),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 13, vjust = 3),
        axis.title.x = element_text(size = 13, vjust = -0.75)) 

ggsave("outputs/fig2_invbirds.png", fig2_invbirds, 
       width = 30, height = 20, units = "cm", dpi = 300)

### ###


### Figure 3: omnivores ####

summary(omn3)

m_guild_omn <- m_guild %>% filter(GuildFeed == "Omnivore")


## logCH:Week interaction ##

# empty table to predict
pred_data5 <- data.frame(Week = c(1, 7, 13)) %>% 
  group_by(Week) %>% 
  reframe(logCH = seq(min(m_guild_omn$logCH), 
                      max(m_guild_omn$logCH), 
                      0.01)) %>% 
  mutate(Observer = factor("KT", levels = levels(m_guild_omn$Observer)),
         HabClass = factor("Road", levels = levels(m_guild_omn$HabClass)), # intermediate values
         logTDens = mean(m_guild_omn$logTDens),
         Point = NA)

# predictions
tic("Bootstrapped prediction 1 for omn. model")
prediction <- boot_conf_GLMM(omn3,
                             new_data = pred_data5,
                             new_data_string = "pred_data5",
                             model_data_string = "m_guild_omn",
                             nsim = 1000)
save(prediction, file = "outputs/pred5.RData")
toc() # 690 sec
# load("outputs/pred5.RData")

# calculating mean and SE from bootstrapped values
for (j in 1:nrow(pred_data5)) {
  
  pred_data5$PRED.LINK[j] <- median(na.omit(prediction[,j]))
  pred_data5$SE.LINK[j] <- sd(na.omit(prediction[,j]))
  
}

pred_data5 <- pred_data5 %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  mutate(CI.U = PRED + Gauss_coeff*SE,
         CI.L = PRED - Gauss_coeff*SE) %>% 
  dplyr::select(Week, logCH, PRED, CI.L, CI.U)


plot_omn1 <- ggplot(pred_data5,
                    aes(x = logCH, y = PRED, 
                        col = as.factor(Week), fill = as.factor(Week))) +
  scale_fill_manual(values = cbPalette[c(3, 1, 7)], name = "Week") +
  scale_colour_manual(values = cbPalette[c(3, 1, 7)], name = "Week") +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = CI.L, ymax = CI.U), alpha = 0.35, colour = NA) +
  scale_y_continuous(breaks = seq(0, 40, 4)) +
  labs(y = "Omnivore detections per point count",
       x = "Canopy heterogeneity") +
  guides(col = guide_legend(title.position = "top"),
         fill = guide_legend(title.position = "top")) +
  theme(legend.position = c(0.4, 0.9), 
        legend.direction = "horizontal")


## logTDens main effect ##

pred_data6 <- data.frame(logTDens = seq(min(m_guild_omn$logTDens), 
                                        max(m_guild_omn$logTDens), 
                                        0.01)) %>% 
  mutate(Observer = factor("KT", levels = levels(m_guild_omn$Observer)),
         Week = 13,
         HabClass = factor("Road", levels = levels(m_guild_omn$HabClass)), # intermediate values
         logCH = mean(m_guild_omn$logCH),
         Point = NA)

# predictions
tic("Bootstrapped prediction 2 for omn. model")
prediction <- boot_conf_GLMM(omn3,
                             new_data = pred_data6,
                             new_data_string = "pred_data6",
                             model_data_string = "m_guild_omn",
                             nsim = 1000)
save(prediction, file = "outputs/pred6.RData")
toc() # 730 sec
# load("outputs/pred6.RData")

# calculating mean and SE from bootstrapped values
for (j in 1:nrow(pred_data6)) {
  
  pred_data6$PRED.LINK[j] <- median(na.omit(prediction[,j]))
  pred_data6$SE.LINK[j] <- sd(na.omit(prediction[,j]))
  
}

pred_data6 <- pred_data6 %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  mutate(CI.U = PRED + Gauss_coeff*SE,
         CI.L = PRED - Gauss_coeff*SE) %>% 
  dplyr::select(logTDens, PRED, CI.L, CI.U)


plot_omn2 <- ggplot(pred_data6,
                    aes(x = logTDens, y = PRED)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = CI.L, ymax = CI.U), alpha = 0.35, colour = NA) +
  scale_y_continuous(breaks = seq(0, 40, 4), limits = c(4, 16)) +
  labs(y = "Omnivore detections per point count", 
       x = expression(Tree~density~per~100~m^2))


## HabClass main effect ##

pred_data7 <- data.frame(HabClass = unique(m_guild_omn$HabClass)) %>% 
  mutate(Observer = factor("KT", levels = levels(m_guild_omn$Observer)),
         Week = 13,
         logCH = mean(m_guild_omn$logCH),
         logTDens = mean(m_guild_omn$logTDens),
         Point = NA)

# predictions
tic("Bootstrapped prediction 3 for omn. model")
prediction <- boot_conf_GLMM(omn3,
                             new_data = pred_data7,
                             new_data_string = "pred_data7",
                             model_data_string = "m_guild_omn",
                             nsim = 1000)
save(prediction, file = "outputs/pred7.RData")
toc() # 880 sec
# load("outputs/pred7.RData")

# calculating mean and SE from bootstrapped values
for (j in 1:nrow(pred_data7)) {
  
  pred_data7$PRED.LINK[j] <- median(na.omit(prediction[,j]))
  pred_data7$SE.LINK[j] <- sd(na.omit(prediction[,j]))
  
}

pred_data7 <- pred_data7 %>% 
  mutate(PRED = exp(PRED.LINK),
         # to transform lower bound of SE (not CI.L! think "mean +- SE")
         SE.L = exp(PRED.LINK - SE.LINK)) %>% 
  mutate(SE = PRED - SE.L) %>% 
  mutate(CI.U = PRED + Gauss_coeff*SE,
         CI.L = PRED - Gauss_coeff*SE) %>% 
  dplyr::select(HabClass, PRED, CI.L, CI.U)


plot_omn3 <- ggplot(pred_data7,
                    aes(x = HabClass, y = PRED)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                linewidth = 1, width = 0.6) +
  scale_y_continuous(breaks = seq(0, 40, 4), limits = c(4, 16)) +
  labs(y = "Omnivore detections per point count", x = "Habitat class")


## Figure 3 ##

fig_layout <- "
AAA
BBC
"

fig3_omnbirds <- plot_omn1 + plot_omn2 +
     (plot_omn3 + theme(axis.title.y = element_blank())) +
  plot_layout(design = fig_layout) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16, margin = margin(0, 7, 0, 0)),
        plot.margin = margin(10, 10, 10, 10),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 13, vjust = 3),
        axis.title.x = element_text(size = 13, vjust = -0.75)) 

ggsave("outputs/fig3_omnbirds.png", fig3_omnbirds, 
       width = 28, height = 25, units = "cm", dpi = 300)

### ###


