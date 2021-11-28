library(tidyverse)
library(stargazer)
library(Hmisc)
library(patchwork)
library(RColorBrewer)
library(ggpubr)


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

stargazer(select(habvar, -c(16:20)), type="html", summary=F,
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
                        DOM = sample(unique(mdata$DOM), 1201, replace = T),
                        Observer = factor("KT", levels = levels(mdata$Observer)),
                        CCavgsd = mean(mdata$CCavgsd),
                        TreeDens = mean(mdata$TreeDens))
all.pred1$All_Abun <- predict(all.6, all.pred1, type="response", re.form=NA)
# no random effects because interested in population level patterns

# Ben Bolker's vignette for plotting confidence intervals:
#(https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions)
all.pred1mm <- model.matrix(terms(all.6), all.pred1)
all.pred1pvar <- diag(all.pred1mm %*% tcrossprod(vcov(all.6), all.pred1mm)) # fixef uncertainty only
# all.pred1tvar <- all.pred1pvar + VarCorr(all.6$Point[1]) + VarCorr(all.6$Days[1]) random uncertainty
all.pred1 <- data.frame(all.pred1,
                        lo = all.pred1$All_Abun - exp(cmult*sqrt(all.pred1pvar)),
                        hi = all.pred1$All_Abun + exp(cmult*sqrt(all.pred1pvar)))
# plotting
all.pred1plot <- ggplot(aes(x=Week, y=All_Abun, col=DOM, fill=DOM), data=all.pred1) +
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
                         DOM = sample(unique(mdata$DOM), 1201, replace = T),
                         Observer = factor("KT", levels = levels(mdata$Observer)),
                         CCavgsd = mean(mdata$CCavgsd),
                         TreeDens = mean(mdata$TreeDens))
all.pred2a$All_Abun <- predict(all.6, all.pred2a, type="response", re.form=NA)
all.pred2amm <- model.matrix(terms(all.6), all.pred2a)
all.pred2apvar <- diag(all.pred2amm %*% tcrossprod(vcov(all.6), all.pred2amm)) # fixef uncertainty only
all.pred2a <- data.frame(all.pred2a,
                         lo = all.pred2a$All_Abun - exp(cmult*sqrt(all.pred2apvar)),
                         hi = all.pred2a$All_Abun + exp(cmult*sqrt(all.pred2apvar)))
# Week 13
all.pred2b <- data.frame(Week = 13,
                         DOM = sample(unique(mdata$DOM), 1201, replace = T),
                         Observer = factor("KT", levels = levels(mdata$Observer)),
                         CCavgsd = mean(mdata$CCavgsd),
                         TreeDens = mean(mdata$TreeDens))
all.pred2b$All_Abun <- predict(all.6, all.pred2b, type="response", re.form=NA)
all.pred2bmm <- model.matrix(terms(all.6), all.pred2b)
all.pred2bpvar <- diag(all.pred2bmm %*% tcrossprod(vcov(all.6), all.pred2bmm)) # fixef uncertainty only
all.pred2b <- data.frame(all.pred2b,
                         lo = all.pred2b$All_Abun - exp(cmult*sqrt(all.pred2bpvar)),
                         hi = all.pred2b$All_Abun + exp(cmult*sqrt(all.pred2bpvar)))
# joining
all.pred2 <- rbind(all.pred2a, all.pred2b)
# plotting
all.pred2plot <- ggplot(data=all.pred2, mapping=aes(x=DOM, y=All_Abun, col=factor(Week))) +
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
                         DOM = factor("0", levels = levels(mdata$DOM)),
                         Observer = factor("KT", levels = levels(mdata$Observer)),
                         CCavgsd = seq(1, 20, 0.001),
                         TreeDens = mean(mdata$TreeDens))

all.pred3a$All_Abun <- predict(all.6, all.pred3a, type="response", re.form=NA)
all.pred3amm <- model.matrix(terms(all.6), all.pred3a)
all.pred3apvar <- diag(all.pred3amm %*% tcrossprod(vcov(all.6), all.pred3amm))
all.pred3a <- data.frame(all.pred3a,
                         lo = all.pred3a$All_Abun - exp(cmult*sqrt(all.pred3apvar)),
                         hi = all.pred3a$All_Abun + exp(cmult*sqrt(all.pred3apvar)))
# Week 13
all.pred3b <- data.frame(Week = 13,
                         DOM = factor("0", levels = levels(mdata$DOM)),
                         Observer = factor("KT", levels = levels(mdata$Observer)),
                         CCavgsd = seq(1, 20, 0.001),
                         TreeDens = mean(mdata$TreeDens))

all.pred3b$All_Abun <- predict(all.6, all.pred3b, type="response", re.form=NA)
all.pred3bmm <- model.matrix(terms(all.6), all.pred3b)
all.pred3bpvar <- diag(all.pred3bmm %*% tcrossprod(vcov(all.6), all.pred3bmm))
all.pred3b <- data.frame(all.pred3b,
                         lo = all.pred3b$All_Abun - exp(cmult*sqrt(all.pred3bpvar)),
                         hi = all.pred3b$All_Abun + exp(cmult*sqrt(all.pred3bpvar)))
# joining
all.pred3 <- rbind(all.pred3a, all.pred3b)
# plotting
all.pred3plot <- ggplot(data=all.pred3, mapping=aes(x=CCavgsd, y=All_Abun, 
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
                         DOM = factor("0", levels = levels(mdata$DOM)),
                         Observer = factor("KT", levels = levels(mdata$Observer)),
                         CCavgsd = mean(mdata$CCavgsd),
                         TreeDens = seq(1,15,0.001))

all.pred4a$All_Abun <- predict(all.6, all.pred4a, type="response", re.form=NA)
all.pred4amm <- model.matrix(terms(all.6), all.pred4a)
all.pred4apvar <- diag(all.pred4amm %*% tcrossprod(vcov(all.6), all.pred4amm))
all.pred4a <- data.frame(all.pred4a,
                         lo = all.pred4a$All_Abun - exp(cmult*sqrt(all.pred4apvar)),
                         hi = all.pred4a$All_Abun + exp(cmult*sqrt(all.pred4apvar)))
# Week 13
all.pred4b <- data.frame(Week = 13,
                         DOM = factor("0", levels = levels(mdata$DOM)),
                         Observer = factor("KT", levels = levels(mdata$Observer)),
                         CCavgsd = mean(mdata$CCavgsd),
                         TreeDens = seq(1,15,0.001))

all.pred4b$All_Abun <- predict(all.6, all.pred4b, type="response", re.form=NA)
all.pred4bmm <- model.matrix(terms(all.6), all.pred4b)
all.pred4bpvar <- diag(all.pred4bmm %*% tcrossprod(vcov(all.6), all.pred4bmm))
all.pred4b <- data.frame(all.pred4b,
                         lo = all.pred4b$All_Abun - exp(cmult*sqrt(all.pred4bpvar)),
                         hi = all.pred4b$All_Abun + exp(cmult*sqrt(all.pred4bpvar)))
# joining
all.pred4 <- rbind(all.pred4a, all.pred4b)
# plotting
all.pred4plot <- ggplot(data=all.pred4, mapping=aes(x=TreeDens, y=All_Abun, 
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


### Visualising results: secondary graphs ####

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


