# seurat_object: seurat object that has been normalized
# features: a one dimentional array of two genes eg. c('gene A', 'gene B')

# calculates the percentage of gene A > gene B or gene A < gene B in each cell of each group

# example output
#                          healthy      Tumor
#B3GALT6 > B3GAT1        1.6829533  0.4494382
#B3GALT6 < B3GAT1        0.4343105  0.2247191

#modified from Percent_Expressing from scCustomize
#input features: array of two genes to compare
REO_Percentage <- function(
    seurat_object,
    features,
    threshold = 0,
    group_by = NULL,
    split_by = NULL,
    entire_object = FALSE,
    slot = deprecated(),
    layer = "data",
    assay = NULL
){
  #checks
  {
    # Check is slot is supplied
    if (lifecycle::is_present(slot)) {
      lifecycle::deprecate_warn(when = "2.0.0",
                                what = "Percent_Expressing(slot)",
                                with = "Percent_Expressing(layer)",
                                details = c("v" = "As of Seurat 5.0.0 the {.code slot} parameter is deprecated and replaced with {.code layer}.",
                                            "i" = "Please adjust code now to prepare for full deprecation.")
      )
      layer <- slot
    }
    
    # set assay (if null set to active assay)
    assay <- assay %||% DefaultAssay(object = seurat_object)
    
    # Check features exist in object
    features_list <- Feature_Present(data = seurat_object, features = features, print_msg = FALSE, case_check = TRUE, seurat_assay = 'SCT')[[1]]
    
    # Check group_by is in object
    if (!is.null(x = group_by) && group_by == "ident") {
      group_by <- NULL
    }
    
    if (!is.null(x = group_by)) {
      possible_groups <- colnames(x = seurat_object@meta.data)
      if (!group_by %in% possible_groups) {
        cli_abort("Grouping variable {.val {group_by}} was not found in Seurat Object.")
      }
    }
    
    # Check split_by is in object
    if (!is.null(x = split_by)) {
      possible_groups <- colnames(x = seurat_object@meta.data)
      if (!split_by %in% possible_groups) {
        cli_abort("Splitting variable {.val {split_by}} was not found in Seurat Object.")
      }
    }
    
  }
  # Pull Expression Info
  cells <- unlist(x = CellsByIdentities(object = seurat_object, idents = NULL))
  expression_info <- FetchData(object = seurat_object, vars = features_list, cells = cells, layer = 'data')
  
  # Add grouping variable
  if (isTRUE(x = entire_object)) {
    expression_info$id <- "All_Cells"
  } 
  else {
    expression_info$id <- if (is.null(x = group_by)) {
      Idents(object = seurat_object)[cells, drop = TRUE]
    } else {
      seurat_object[[group_by, drop = TRUE]][cells, drop = TRUE]
    }
  }
  
  if (!is.factor(x = expression_info$id)) {
    expression_info$id <- factor(x = expression_info$id)
  }
  id.levels <- levels(x = expression_info$id)
  expression_info$id <- as.vector(x = expression_info$id)
  
  # Calculate REO percentage
  
  #data.use: remove ident column of expression_info dataframe
  data.use <- expression_info[expression_info$id == 'healthy', 1:(ncol(x = expression_info) - 1), drop = FALSE]
  data.use_T <- expression_info[expression_info$id == 'Tumor', 1:(ncol(x = expression_info) - 1), drop = FALSE]
  AB <- sum(data.use[1] > data.use[2]) / nrow(data.use) * 100 # A > B in healthy
  BA <- sum(data.use[1] < data.use[2]) / nrow(data.use) * 100
  
  AB_T <- sum(data.use_T[1] > data.use_T[2]) / nrow(data.use_T) * 100 # A > B in tumor
  BA_T <- sum(data.use_T[1] < data.use_T[2]) / nrow(data.use_T) * 100
  
  row1 <- paste(colnames(data.use)[1],'>',colnames(data.use)[2])
  row2 <- paste(colnames(data.use)[1],'<',colnames(data.use)[2])
  
  #create returin df
  final_df <- data.frame('healthy' = c(AB,BA), 'Tumor' = c(AB_T,BA_T), row.names = c(row1, row2))
  return(final_df)
}
