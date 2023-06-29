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

### setting theme

theme_set(theme_classic())

# purples: #A204B4 #D91EFA
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072A2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072A2", "#D55E00", "#CC79A7")
# c("#d8b365", "#f5f5f5", "#5ab4ac")
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

## Week & DOM main effects ##

# empty table to predict
pred_data1 <- data.frame(Week = 1:13) %>% 
  group_by(Week) %>% 
  # need to predict for all five DOM categories at each time-step
  reframe(DOM = unique(m_all$DOM)) %>% 
  # other variables at fixed values
  mutate(Observer = factor("KT", levels = levels(m_all$Observer)),
         logCH = mean(m_all$logCH),
         logTDens = mean(m_all$logTDens))

# predictions
tic("Bootstrapped prediction 1 for all-birds model")
prediction <- boot_conf_GLMM(all4,
                             new_data = pred_data1,
                             new_data_string = "pred_data1",
                             nsim = 1000)
save(prediction, file = "outputs/pred1.RData")
toc() # 43 min
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
  mutate(DOM = factor("Bare", levels = levels(m_all$DOM)),
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
  mutate(DOM = factor("Bare", levels = levels(m_all$DOM)),
         Observer = factor("KT", levels = levels(m_all$Observer)))

pred_data2 <- bind_rows(temp1, temp2, .id = "ID")


tic("Bootstrapped predictions for all-birds model (CH and TDens)")
prediction <- boot_conf_GLMM(all4,
                             new_data = pred_data2,
                             new_data_string = "pred_data2",
                             nsim = 1000)
save(prediction, file = "outputs/pred2.RData")
toc() # 42 min
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

#

### Fig 2: inv-feeders ####


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


### Fig 3: omnivores ####


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


### secondary graphs ####

### Fig.3: average bird detections per point count over the weeks and points ###

fig3 <- 
  (ggplot(b_data_A1_all, aes(x=Week, y=All_Abun, group=Week)) +
     geom_point(col="#213B73", size=3.5, stat = "summary",
                fun.data = mean_se, fun.args = list(mult=2)) +
     stat_summary(col="#213B73", geom = "errorbar", size=1.5, width=0.2, 
                  fun.data = mean_se, fun.args = list(mult=2)) +
     scale_y_continuous(breaks = seq(0,40,4)) +
     scale_x_continuous(breaks = seq(1,13,1)) +
     coord_cartesian(ylim = c(0,16)) +
     labs(y = "Bird detections per point count")) /
  (ggplot(b_data_A1_all, aes(x=Point, y=All_Abun)) + 
     geom_boxplot(fill="#2E799E", alpha=0.6, size=0.9, outlier.alpha = 1) +
     scale_y_continuous(breaks = seq(0,40,4)) +
     labs(y = "Bird detections per point count") +
     theme(axis.text.x = element_text(angle = 45, hjust = 0.75, size = 8))) + # plots median
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) 

ggsave("Fig3.png", fig3, 
       width = 18, height = 24, units = "cm", dpi=300)



### Fig.9: correlations between canopy measures ###

# ( ggplot(habvar) + geom_point(aes(CCavg, TreeDens), size=2, colour="#655B1B") +
#   scale_y_continuous(breaks = seq(0,20,2)) + scale_x_continuous(breaks = seq(0,100,10)) +
#   labs(x = "Canopy cover", y = expression(Tree~density~per~100~m^2)) ) +
# ( ggplot(habvar) + geom_point(aes(CCavg, CCavgsd), size=2, colour="#655B1B") +
#   scale_y_continuous(breaks = seq(0,24,4)) + scale_x_continuous(breaks = seq(0,100,10)) +
#   labs(x = "Canopy cover", y = "Canopy heterogeneity") ) +
# ( ggplot(habvar) + geom_point(aes(CCavgsd, TreeDens), size=2, colour="#655B1B") +
#   scale_x_continuous(breaks = seq(0,24,4)) + scale_y_continuous(breaks = seq(0,20,2)) +
#   labs(x = "Canopy heterogeneity", y = expression(Tree~density~per~100~m^2)) ) +
# ( ggplot(habvar) + geom_point(aes(TreePropDeci, TreeDens), size=2, colour="#655B1B") +
#   scale_x_continuous(breaks = seq(0,100,10)) + scale_y_continuous(breaks = seq(0,20,2)) +
#   labs(x = "Proportion of deciduous trees", y = expression(Tree~density~per~100~m^2)) ) +
# plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) -> fig9

# with correlation lines
( ggscatter(habvar, x="CCavg", y="TreeDens", 
            cor.method="pearson", cor.coef=T,
            conf.int=T, add="reg.line", size=2, color="#655B1B") +
    scale_y_continuous(breaks = seq(0,20,2)) + scale_x_continuous(breaks = seq(0,100,10)) +
    labs(x = "Canopy cover", y = expression(Tree~density~per~100~m^2)) ) +
  ( ggscatter(habvar, x="CCavg", y="CCavgsd", 
              cor.method="pearson", cor.coef=T,
              conf.int=T, add="reg.line", size=2, color="#655B1B") +
      scale_y_continuous(breaks = seq(0,24,4)) + scale_x_continuous(breaks = seq(0,100,10)) +
      labs(x = "Canopy cover", y = "Canopy heterogeneity") ) +
  ( ggscatter(habvar, x="CCavgsd", y="TreeDens", 
              cor.method="pearson", cor.coef=T,
              conf.int=F, size=2, color="#655B1B") +
      scale_x_continuous(breaks = seq(0,24,4)) + scale_y_continuous(breaks = seq(0,20,2)) +
      labs(x = "Canopy heterogeneity", y = expression(Tree~density~per~100~m^2)) ) +
  ( ggscatter(habvar, x="TreePropDeci", y="TreeDens", 
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

h_data_DOMtrends <- inner_join(select(habvar, c(1,5,6,12,14,15)),
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

( ggplot(mdata, aes(Observer, All_Abun)) +
    geom_boxplot(fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
    labs(y = "Bird detections per point count") +
    scale_y_continuous(breaks = seq(0,40,4)) +
    coord_cartesian(ylim=c(0,28)) ) +
  ( ggplot(b_rawdata_all, aes(Observer, All_Abun)) +
      geom_boxplot(fill="#655B1B", alpha=0.6, size=0.7, outlier.alpha = 1) +
      labs(y = "Bird detections per point count") +
      scale_y_continuous(breaks = seq(0,60,8)) +
      coord_cartesian(ylim=c(0,56)) +
      theme(axis.title.y = element_blank()) ) +
  ( ggplot(mdata, aes(Observer, All_Abun, group=Observer)) +
      geom_bar(stat="summary", fun="sum", fill="#655B1B", alpha=0.6, size=0.7, col="black") +
      labs(y = "Overall bird detections") +
      coord_cartesian(ylim=c(0,4000)) +
      scale_y_continuous(breaks = seq(0,4000,500)) +
      geom_hline(yintercept = sum(mdata$All_Abun), linetype="dashed", size=1.2) ) +
  ( ggplot(b_rawdata_all, aes(Observer, All_Abun, group=Observer)) +
      geom_bar(stat="summary", fun="sum", fill="#655B1B", alpha=0.6, size=0.7, col="black") +
      labs(y = "Overall bird detections") +
      coord_cartesian(ylim=c(0,8000)) +
      scale_y_continuous(breaks = seq(0,9000,1000)) +
      geom_hline(yintercept = sum(b_rawdata_all$All_Abun), linetype="dashed", size=1.2) +
      theme(axis.title.y = element_blank()) ) +
  plot_layout(ncol=4) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)) -> fig11

ggsave("Fig11.png", fig11, 
       width = 20, height = 10, units = "cm", dpi=300)



( ggplot(mdata) +
    geom_col(aes(Week, All_Abun, fill=Observer)) +
    scale_fill_viridis_d(end=0.5, direction=-1) +
    scale_x_continuous(breaks = 1:13) +
    scale_y_continuous(breaks = seq(0,500,50)) +
    labs(y = "Number of detections") ) +
  ( ggplot(mdata) +
      geom_col(aes(Point, All_Abun, fill=Observer)) +
      scale_fill_viridis_d(end=0.5, direction=-1) +
      scale_y_continuous(breaks = seq(0,200,20)) +
      coord_cartesian(ylim = c(0,140)) +
      labs(y = "Number of detections") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) ) +
  ( ggplot(b_rawdata_all) +
      geom_col(aes(Week, All_Abun, fill=Observer)) +
      scale_fill_viridis_d(end=0.5, direction=-1) +
      scale_x_continuous(breaks = 1:13) +
      scale_y_continuous(breaks = seq(0,1000,100)) +
      labs(y = "Number of detections") ) +
  ( ggplot(b_rawdata_all) +
      geom_col(aes(Point, All_Abun, fill=Observer)) +
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


