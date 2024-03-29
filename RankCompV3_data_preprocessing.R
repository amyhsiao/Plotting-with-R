# Before using RankCompV3, download the original source code 
# and in RankCompV3.jl line 426 change "no change" to "nchange"
# use the modified  package instead of the original one

# subset(liver_HS, subset = Type %in% 'hepatocyte')
# meta of hep
{
cellnames_h <- colnames(subset(liver_HS, subset = Type %in% 'hepatocyte')) %>% data.frame()
rownames(cellnames_h) <- colnames(subset(liver_HS, subset = Type %in% 'hepatocyte'))
colnames(cellnames_h) <- "sample_name"
head(cellnames_h)

liver_hp <- as.character(levels(subset(liver_HS, subset = Type %in% 'hepatocyte')$health))[subset(liver_HS, subset = Type %in% 'hepatocyte')$health] %>%data.frame()
rownames(liver_hp) <- colnames(subset(liver_HS, subset = Type %in% 'hepatocyte'))
colnames(liver_hp) <- "group"
head(liver_hp)

p_h <- merge(cellnames_h, liver_hp, by = "row.names") %>% data.frame()
head(p_h)
p_h <- subset(p_h, select = -Row.names)
head(p_h)

write.table(p_h, "./RankComp/data/hep_meta.txt", row.names = F, quote = F)
}
#exp of hep
{
a <- subset(liver_HS, subset = Type %in% 'hepatocyte')[['RNA']]$counts[cs_genes,]

write.table(data.frame("gene_name"=rownames(a),a), 
            "./RankComp/data/hep_expr.txt",
            quote=F, row.names = F)
}

# meta of str
{
  cellnames_s <- colnames(subset(liver_HS, subset = Type %in% 'stromal cell')) %>% data.frame()
  rownames(cellnames_s) <- colnames(subset(liver_HS, subset = Type %in% 'stromal cell'))
  colnames(cellnames_s) <- "sample_name"
  head(cellnames_s)
  
  liver_sp <- as.character(levels(subset(liver_HS, subset = Type %in% 'stromal cell')$health))[subset(liver_HS, subset = Type %in% 'stromal cell')$health] %>%data.frame()
  rownames(liver_sp) <- colnames(subset(liver_HS, subset = Type %in% 'stromal cell'))
  colnames(liver_sp) <- "group"
  head(liver_sp)
  
  p_s <- merge(cellnames_s, liver_sp, by = "row.names") %>% data.frame()
  head(p_s)
  p_s <- subset(p_s, select = -Row.names)
  head(p_s)
  
  write.table(p_s, "./RankComp/data/str_meta.txt", row.names = F, quote = F)
}
#exp of str
{
  a <- subset(liver_HS, subset = Type %in% 'stromal cell')[['RNA']]$counts[cs_genes,]
  
  write.table(data.frame("gene_name"=rownames(a),a), 
              "./RankComp/data/str_expr.txt",
              quote=F, row.names = F)
}
#important: replace '.' or '-' with '_' in both of the files or there will be an error


###############for pg genes
pg_genes <- unlist(gs_pg)
names(pg_genes) <- NULL

#exp of hep
{
  a <- subset(liver_HS, subset = Type %in% 'hepatocyte')[['RNA']]$counts[pg_genes,]
  
  write.table(data.frame("gene_name"=rownames(a),a), 
              "./RankComp/data/hep_expr_pg.txt",
              quote=F, row.names = F)
}

#exp of str
{
  a <- subset(liver_HS, subset = Type %in% 'stromal cell')[['RNA']]$counts[pg_genes,]
  
  write.table(data.frame("gene_name"=rownames(a),a), 
              "./RankComp/data/str_expr_pg.txt",
              quote=F, row.names = F)
}

###############for hs genes
hs_genes <- unlist(reactome_gs_hs)
names(hs_genes) <- NULL

hs_genes

#exp of hep
{
  a <- subset(liver_HS, subset = Type %in% 'hepatocyte')[['RNA']]$counts[hs_genes,]
  
  write.table(data.frame("gene_name"=rownames(a),a), 
              "./RankComp/data/hep_expr_hs.txt",
              quote=F, row.names = F)
}

#exp of str
{
  a <- subset(liver_HS, subset = Type %in% 'stromal cell')[['RNA']]$counts[hs_genes,]
  
  write.table(data.frame("gene_name"=rownames(a),a), 
              "./RankComp/data/str_expr_hs.txt",
              quote=F, row.names = F)
}


