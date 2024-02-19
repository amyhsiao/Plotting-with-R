install.packages("ggpmisc")
library("ggpmisc")
library("ggplot2")
#Using Seurat V5

#formula and pval
{
  formula <- y~x
  Theme <- theme(axis.title.x = element_text(size=24),
                 axis.title.y = element_text(size=24),
                 axis.text.x = element_text(size=18),
                 axis.text.y = element_text(size=18),
                 plot.title = element_text(size=28, hjust=0.5)) 
  Pval <- stat_fit_glance(method = "lm", 
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "top",
                          aes(label = sprintf("italic(P)*\"-value = \"*%.3g", 
                                              after_stat(p.value))),
                          parse = TRUE)
}

#Turning single cell data into dataframe for ggplot2
DCN_df <- FetchData(object = liver_ss, vars = c("DCN")) #deleted layers = "counts"
CS <- liver_ss$cs_signature %>% data.frame() #cs_signature is a metadata column
colnames(CS) <- "CS"
head(CS)

DCN_df <- merge(DCN_df, CS, by = 'row.names', all = TRUE) 
head(DCN_df)

sp_DCN <- ggplot(data=DCN_df, aes(x=DCN, y=CS)) + geom_point(color='grey60')
sp_DCN 

#final plot
sp_DCN + stat_poly_line(formula=formula) + 
  ggtitle("") + theme_bw() + Theme + xlim(0.001, NA) + Pval +
  stat_poly_eq(formula = formula, label.x = "right", label.y = 0.995)
