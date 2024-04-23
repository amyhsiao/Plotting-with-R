#SCP2154
#load umap coordinates
{
UMAP_coordinates <- read.table("./DATA/2022_SI_human_liver_allcells_cluster_file.txt", 
                               stringsAsFactors = F, header = T, row.names=1)
UMAP_coordinates <-  UMAP_coordinates[-1,]
colnames(UMAP_coordinates) <- c("UMAP_1", "UMAP_2")
UMAP_coordinates$UMAP_1 <- as.numeric(UMAP_coordinates$UMAP_1)
UMAP_coordinates$UMAP_2 <- as.numeric(UMAP_coordinates$UMAP_2)

str(UMAP_coordinates)

UMAP_coordinates_mat <- as(UMAP_coordinates, "matrix")
str(UMAP_coordinates_mat)
}

liver@reductions$UMAP <- CreateDimReducObject(embeddings = UMAP_coordinates_mat, 
                                                 key = "UMAP_", global = T, assay = "RNA")
