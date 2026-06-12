deconv_spot <- function(sces,
                        spes,
                        sc_cell_id_column,
                        sc_celltype_column,
                        sc_gene_id_column,
                        filter.mt.genes = "^mt-",
                        st_spot_ID_column = "cell_ID",
                        st_gene_ID_column = "feat_ID",
                        project = "ddls",
                        threads = 8,
                        output_dir = NULL,
                        ...){
  # loading required packages
  library(qs)
  library(SpatialDDLS)
  library(SpatialExperiment)
  
  if(is.null(output_dir)){
    output_dir = "./"
  }

  # create spatialDDLS object
  message("create spatialDDLS object...")
  ddls <- createSpatialDDLSobject(
    sc.data = sces, 
    sc.cell.ID.column = sc_cell_id_column, 
    sc.cell.type.column = sc_celltype_column, 
    sc.gene.ID.column = sc_gene_id_column, 
    filter.mt.genes = filter.mt.genes,
    st.data = spes, 
    st.spot.ID.column = st_spot_ID_column, 
    st.gene.ID.column = st_gene_ID_column, 
    project = project
  )
  
  message("generate gene cell mix profile...")
  ddls <- genMixedCellProp(
    object = ddls, 
    cell.ID.column = sc_cell_id_column, 
    cell.type.column = sc_celltype_column, 
    num.sim.spots = 10000, 
    n.cells = 50
  )
  
  message("simulation mixed profile...")
  ddls <- simMixedProfiles(
    object = ddls, 
    type.data = "both", 
    threads = threads
  )
  message("training the prediction model...")
  ddls <- trainDeconvModel(
    object = ddls, 
    type.data.train = "mixed", 
    threads = threads
  )
  qsave(ddls, file = paste0(output_dir, "ddls.qs"))
  
  message("deconvolution spots...")
  ddls <- deconvSpatialDDLS(
    object = ddls, 
    index.st = 1, 
    k.spot = 10, 
    fast.pca = TRUE
  )
  
  return(ddls)
}
