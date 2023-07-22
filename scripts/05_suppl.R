# supplementary graphs ###

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


