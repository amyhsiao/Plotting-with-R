{
  dittoPlot(
    object = subset(liver_ss.ds, subset= Type %in% c("stromal cell", "hepatocyte")),
    plots = c('vlnplot', 'boxplot'),
    var = 'cs_signature',
    group.by = "health",
    split.by = "Type", 
    jitter.size = NA,
    boxplot.width=0.5,
    min=ymin,
    max=8200,
    colors = c(5,6))+
    ggtitle(paste("All samples [", sig_name, "]"))+
    stat_compare_means(comparisons = my_comparisons,
                       label.y = pvalue_y)+
    theme_bw()+
    theme(axis.text = element_text(size = 8),
          axis.text.x = element_text(angle = 0, size=14),
          axis.text.y = element_text(size=12),
          strip.text = element_text(size=12),
          legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          title = element_blank(),
          #axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("CS signature score")
}
