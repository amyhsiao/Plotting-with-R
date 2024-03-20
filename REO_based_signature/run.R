# use REO_percentage to calculate the REO of all cell pairs in the reactome cs gene set
reactome_gs_cs <- geneIds(getGeneSets(species = "Homo sapiens",
                                      library =  c("C2"), 
                                      gene.sets = c("KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_CHONDROITIN_SULFATE")))

cs_genes <- unlist(reactome_gs_cs)
names(cs_genes) <- NULL
cs_genes


str_df <- data.frame(healthy = c(), Tumor = c())
hep_df <- data.frame(healthy = c(), Tumor = c())
for(i in 1:22){
  for(j in (i+1):22){
    if(j > 22){
      break
    }
    r <- REO_Percentage(subset(liver_HS, subset = Type %in% 'stromal cell'), c(cs_genes[i], cs_genes[j]), group_by = 'health')
    r_h <- REO_Percentage(subset(liver_HS, subset = Type %in% 'hepatocyte'), c(cs_genes[i], cs_genes[j]), group_by = 'health')
    
    str_df <- bind_rows(str_df, r)
    hep_df <- bind_rows(hep_df, r_h)

  }
}
str_df
hep_df
